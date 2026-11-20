# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write the CIR timelines for the live-twin runs (chem_taps_timeline, same C++ models as the extensions).

T1 (SISO Doppler, tone sounder at 1 MS/s, b205 link): chem_today | proposed | jakes | tdl_a | tdl_d
    flight: GT 10 m, UAV 30 m, 5 m/s radial approach (as E1), 70 s, 1 kHz updates
T2 (2x2 LTE TM4 at 11.52 MS/s, b210 link): SCM A2G, Doppler PROPOSED co-pol 0.5 lambda, and
    cross-polarized arrays with PROPOSED | SCM | JAKES | TDL_D, plus identity (control)
    flight: UAV 30 m, 5 m/s radial outbound from 100 m, 130 s, 1 kHz updates
"""
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from chemlib import *  # noqa

OUT = os.path.join(os.path.dirname(__file__), '..', 'twin', 'timelines')
os.makedirs(OUT, exist_ok=True)
ENV = dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser('~/deps/mavsdk_x/usr/lib') + ':' +
           os.path.expanduser('~/mm/root/envs/chembuild/lib'))


def run(name, cfg):
    cfg['out'] = os.path.join(OUT, name + '.bin')
    cfg['infoOut'] = os.path.join(OUT, name + '_info.csv')
    jp = os.path.join(RUNDIR, name + '.json')
    json.dump(cfg, open(jp, 'w'), indent=1)
    r = subprocess.run([TIMELINE, jp], capture_output=True, text=True, env=ENV)
    if r.returncode:
        raise RuntimeError(r.stderr[-1500:])
    print(name, os.path.getsize(cfg['out']) // 1024, 'kB', flush=True)


gt = write_traj(traj_static([0.0, 0.0, 10.0]), os.path.join(RUNDIR, 'twin_gt.csv'))
SPEED = float(os.environ.get('SPEED', 5.0))
uav1 = write_traj(traj_line([-60.0 - 72.0 * SPEED, 0.0, 30.0], [SPEED, 0.0, 0.0], 72.0), os.path.join(RUNDIR, 'twin_uav1.csv'))
base1 = dict(model='a2g', fc=3.4e9, fs=1e6, rateHz=1000, durationS=70, nTx=1, nRx=1, maxTaps=4,
             trajTx=gt, trajRx=uav1)
T1 = {
    'T1_chem_today': dict(paths=[dict(type='LOS')], K_dB=300, dopplerMode='SHARED_LOS', updateMs=100.0,
                          phaseUpdate='START_HOLD'),
    'T1_proposed': dict(dopplerModel='PROPOSED'),
    'T1_jakes': dict(dopplerModel='JAKES'),
    'T1_tdl_a': dict(dopplerModel='TDL_A', tdlDelaySpreadNs=100),
    'T1_tdl_d': dict(dopplerModel='TDL_D', tdlDelaySpreadNs=100),
}
SUF = os.environ.get('SUFFIX', '')
for name, ext in T1.items():
    run(name + SUF, dict(base1, config=ext_a2g(**ext) if 'paths' not in ext else dict(ext, seed=1234)))

if os.environ.get('T1_ONLY') == '1':
    sys.exit(0)
uav2 = write_traj(traj_line([100.0, 20.0, 30.0], [5.0, 0.0, 0.0], 132.0), os.path.join(RUNDIR, 'twin_uav2.csv'))
base2 = dict(model='scm', fc=3.41e9, fs=11.52e6, rateHz=1000, durationS=130, nTx=2, nRx=2, maxTaps=8,
             trajTx=gt, trajRx=uav2)


def scm(dop, xpol, sp=0.5):
    return dict(scenario='A2G', doppler=dop, seed=2026, dropDistanceM=50.0,
                groundArray=dict(elements=2, spacingWl=sp, xpol=xpol), airArray=dict(elements=2, spacingWl=sp, xpol=xpol))


T2 = {'T2_prop_copol': scm('PROPOSED', False), 'T2_prop_xpol': scm('PROPOSED', True),
      'T2_scm_xpol': scm('SCM', True), 'T2_jakes_xpol': scm('JAKES', True), 'T2_tdld_xpol': scm('TDL_D', True)}
for name, c in T2.items():
    run(name, dict(base2, config=c))
# identity control (constant): write directly in the same format
K = 130 * 1000
hdr = dict(model='const', doppler='NONE', fc=3.41e9, fs=11.52e6, rateHz=1000, K=K, nRx=2, nTx=2, L=1)
with open(os.path.join(OUT, 'T2_identity.bin'), 'wb') as f:
    f.write((json.dumps(hdr) + '\n').encode())
    a = np.zeros((K, 2, 2, 1, 2), np.float32)
    a[:, 0, 0, 0, 0] = 1
    a[:, 1, 1, 0, 0] = 1
    f.write(a.tobytes())
print('done')
