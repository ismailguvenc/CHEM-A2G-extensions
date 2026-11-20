#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""E4: sample experiment on the laptop with the public V-USRP and the local CHEM build (no AERPAW access).

Two V-USRP B210 nodes, 'gt' (LW1, 10 m) and 'uav' (30 m), talk through the local CHEM process
(build_v2: a2g_doppler v2 + scm_mimo + CHEM-EXT-HOOK/MIMO). Control goes through CHEM's public
pychem client. The UAV flies a scripted radial leg (5 m/s, 30 m) through the extensions'
'trajectoryFiles' option. Cases:
  c0 stock CHEM (FREE_SPACE, no extension, 1 tone)
  a2g_<model>: a2g_doppler hook, PROPOSED | JAKES | TDL_A | TDL_D (1 tone)
  scm_copol / scm_xpol: scm_mimo hook, 2x2, two tones (one per TX port)
Run inside WSL with the 'vusrp' env (see run_e4.sh)."""
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', time.strftime('E4_%m%d_%H%M'))
os.makedirs(OUT, exist_ok=True)
CHEM_SRC = os.path.join(HERE, '..', '..', 'CHEM_source', 'CHEM')
sys.path.insert(0, os.path.join(CHEM_SRC, 'utils'))
import pychem  # noqa: E402

CHEM_BIN = os.path.expanduser('~/work/build_v2/bin/Release/CHEM')
LAT0, LON0 = 35.72750947, -78.69595819
GT = dict(lat=35.727489, lon=-78.695957, alt=10.0)
UAV0_ENU = (-320.0, 0.0, 30.0)


def ll(e, n):
    lat = LAT0 + n / 111319.5
    return lat, LON0 + e / (111319.5 * np.cos(np.radians(0.5 * (lat + LAT0))))


def log(*a):
    s = time.strftime('[%H:%M:%S] ') + ' '.join(str(x) for x in a)
    print(s, flush=True)
    open(os.path.join(OUT, 'timeline.txt'), 'a').write(s + '\n')


# ---- configs ----------------------------------------------------------------------------------
ulat, ulon = ll(UAV0_ENU[0], UAV0_ENU[1])
conf = f"""chem_address: 127.0.0.1
vusrps:
  - name: gt
    usrp_type: b210
    serial: VUSRP-GT
    num_channels: 2
    logical_node: 0
    location: {{lat: {GT['lat']}, lon: {GT['lon']}, alt: {GT['alt']}}}
  - name: uav
    usrp_type: b210
    serial: VUSRP-UAV
    num_channels: 2
    logical_node: 1
    location: {{lat: {ulat:.8f}, lon: {ulon:.8f}, alt: {UAV0_ENU[2]}}}
"""
os.makedirs(os.path.expanduser('~/.config/uhd'), exist_ok=True)
open(os.path.expanduser('~/.config/uhd/conf.yaml'), 'w').write(conf)
open(os.path.join(OUT, 'conf.yaml'), 'w').write(conf)
traj = os.path.join(OUT, 'uav_traj.csv')
with open(traj, 'w') as f:
    f.write('t,east,north,up,ve,vn,vu\n')
    for i in range(0, 1201):
        t = i * 0.1
        f.write(f'{t:.2f},{UAV0_ENU[0] + 5.0 * t:.3f},0,30,5,0,0\n')
cfg = json.load(open(os.path.join(CHEM_SRC, 'configs', 'config.json')))
cfg['logDirectory'] = os.path.join(OUT, 'chem.log')
cfg_path = os.path.join(OUT, 'chem_config.json')
json.dump(cfg, open(cfg_path, 'w'), indent=1)

env = dict(os.environ)
env['LD_LIBRARY_PATH'] = ':'.join([os.path.expanduser('~/deps/mavsdk_x/usr/lib'),
                                   os.path.expanduser('~/mm/root/envs/chembuild/lib')])
chem_log = open(os.path.join(OUT, 'chem_stdout.txt'), 'w')
chem = subprocess.Popen([CHEM_BIN, cfg_path], stdout=chem_log, stderr=subprocess.STDOUT, env=env)
log('CHEM started, pid', chem.pid)
c = pychem.ChemClient('127.0.0.1', 5000)
for i in range(60):
    try:
        v = c.get_version()
        if v:
            break
    except Exception:
        pass
    time.sleep(0.5)
log('CHEM version', c.get_version())
log('extensions', list((c.extension_list() or {}).get('extensions', {}).keys()))

PY = sys.executable
TOOLS = os.path.join(HERE, 'tools')
FREQ = 3400.0


def wait_link(src, dst, t=40):
    t0 = time.time()
    while time.time() - t0 < t:
        for L in (c.get_individual_channels() or {}).get(f'{FREQ:.6f}', []):
            if L.get('src') == src and L.get('dest') == dst:
                return True
        time.sleep(0.5)
    return False


def run_case(tag, tones, nch, ext=None, secs=25):
    log('--- case', tag)
    tx = subprocess.Popen([PY, os.path.join(TOOLS, 'tone_tx.py'), 'gt', str(secs + 25)] + [str(f) for f in tones],
                          stdout=open(os.path.join(OUT, tag + '_tx.txt'), 'w'), stderr=subprocess.STDOUT)
    time.sleep(4)
    rx = subprocess.Popen([PY, os.path.join(TOOLS, 'rx_dec.py'), 'uav', str(secs + 8), os.path.join(OUT, tag + '_rx.npz'),
                           str(nch)] + [str(f) for f in tones],
                          stdout=open(os.path.join(OUT, tag + '_rx.txt'), 'w'), stderr=subprocess.STDOUT)
    ok = wait_link('gt', 'uav')
    log('  link gt->uav', 'up' if ok else 'NOT FOUND')
    info = {}
    if ok:
        c.set_noise_model(FREQ, 'gt', 'uav', 'FIXED_SNR', 35.0)
        c.set_doppler(FREQ, 'gt', 'uav', False)
        if ext:
            name, params = ext
            r = c.extension(name, 'start', params)
            log('  ext', name, 'start:', r.get('status'), r.get('message', ''))
            time.sleep(2)
            info['ext_status'] = c.extension(name, 'status')
        info['links'] = [L for L in (c.get_individual_channels() or {}).get(f'{FREQ:.6f}', []) if L.get('src') == 'gt']
    rx.wait(timeout=secs + 60)
    if ext:
        c.extension(ext[0], 'stop')
    tx.wait(timeout=60)
    json.dump(info, open(os.path.join(OUT, tag + '_info.json'), 'w'), indent=1, default=str)
    log('  rx:', open(os.path.join(OUT, tag + '_rx.txt')).read().strip().splitlines()[-1:])
    time.sleep(20)  # CHEM node timeout before the next registration


base = dict(enabled=True, mode='hook', referenceOrigin=dict(lat=LAT0, lon=LON0, alt=0.0),
            links=[], nodes=dict(gt=dict(enu=[0.0, 0.0, 10.0])),
            trajectoryFiles=dict(uav=traj), trajectoryDelayMs=0, seed=1234,
            K_dB=8.5, shape='GAUSS', taps='ALL', paths=[dict(type='LOS'), dict(type='GROUND', reflCoeff=-0.3)])
ONLY = os.environ.get('E4_CASES', '').split()
_run = run_case


def run_case(tag, *a, **k):
    if ONLY and tag not in ONLY:
        return
    _run(tag, *a, **k)


try:
    run_case('c0_stock', [100e3], 1)
    for m in ('PROPOSED', 'JAKES', 'TDL_A', 'TDL_D'):
        run_case('a2g_' + m.lower(), [100e3], 1, ('a2g_doppler', dict(base, dopplerModel=m, tdlDelaySpreadNs=100)))
    for tag, xp in (('scm_copol', False), ('scm_xpol', True)):
        p = dict(base, scenario='A2G', doppler='PROPOSED',
                 groundArray=dict(elements=2, spacingWl=0.5, xpol=xp), airArray=dict(elements=2, spacingWl=0.5, xpol=xp))
        run_case(tag, [100e3, 230e3], 2, ('scm_mimo', p))
finally:
    chem.terminate()
    try:
        chem.wait(timeout=10)
    except Exception:
        chem.kill()
    log('CHEM stopped; results in', OUT)
