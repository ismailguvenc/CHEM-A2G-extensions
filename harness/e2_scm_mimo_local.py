# SPDX-License-Identifier: AGPL-3.0-or-later
"""E2 (local CHEM): 2x2 SCM through the CHEM data path, with the Doppler options.

Two tones (TX port 0 at +20 kHz, TX port 1 at -30 kHz) go through chem::Intermediate with 2-channel
source and destination nodes (scm_harness). Per 1 ms block, the tone phasors on each RX port give
the narrowband 2x2 matrix H(t). Flight: GT 10 m, UAV 30 m, 5 m/s radial, 60 s.
Cases:
  chem_ones  : CHEM default real port matrix [[1,1],[1,1]] (twin default)
  chem_ident : CHEM real port matrix identity (set_mimo_matrix 'identity')
  scm_<dop>  : scm_mimo hook, scenario A2G, Doppler SCM | PROPOSED | JAKES | TDL_D, 0.5 lambda arrays
  scm_prop_4wl: PROPOSED with 4 lambda ground spacing
Also: timeline statistics (chem_taps_timeline) over a 300 s flight per option and spacing:
capacity at 10/20 dB SNR and condition number per 10 ms.
"""
import json
import os
import subprocess
import sys

import numpy as np
from scipy import signal as sps

sys.path.insert(0, os.path.dirname(__file__))
from chemlib import *  # noqa

OUT = os.path.join(RESULTS, 'e2_scm_mimo_local')
XPOL_ONLY = os.environ.get('XPOL_ONLY') == '1'
os.makedirs(OUT, exist_ok=True)
FS = 200e3
T = 60.0
F0, F1 = 20e3, -30e3
BLK = 200  # 1 ms
tr_gt = traj_static([0.0, 0.0, 10.0])
tr_uav = traj_line([-320.0, 0.0, 30.0], [5.0, 0.0, 0.0], T + 2.0)
n = np.arange(int(FS * T))
x0 = (np.exp(2j * np.pi * F0 * n / FS) / np.sqrt(2)).astype(np.complex64)
x1 = (np.exp(2j * np.pi * F1 * n / FS) / np.sqrt(2)).astype(np.complex64)
ENV = dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser('~/deps/mavsdk_x/usr/lib') + ':' +
           os.path.expanduser('~/mm/root/envs/chembuild/lib'))


def scm_ext(dop, spG=0.5, spA=0.5, xpol=False, **kw):
    e = dict(mode='hook', scenario='A2G', doppler=dop, updateMs=0.5, seed=kw.pop('seed', 77),
             groundArray=dict(elements=2, spacingWl=spG, xpol=xpol), airArray=dict(elements=2, spacingWl=spA, xpol=xpol))
    e.update(kw)
    return e


def run_scm(tag, chem, ext):
    files = [os.path.join(RUNDIR, f'{tag}_tx{i}.cf32') for i in range(2)]
    outs = [os.path.join(RUNDIR, f'{tag}_rx{i}.cf32') for i in range(2)]
    x0.tofile(files[0]); x1.tofile(files[1])
    cfg = {'fs': FS, 'fc': FC, 'frameSamples': 200, 'nTx': 2, 'nRx': 2, 'in': files, 'out': outs,
           'trajTx': write_traj(tr_gt, os.path.join(RUNDIR, tag + '_tx.csv')),
           'trajRx': write_traj(tr_uav, os.path.join(RUNDIR, tag + '_rx.csv')), 'chem': chem}
    if ext:
        cfg['ext'] = ext
    jp = os.path.join(RUNDIR, tag + '.json')
    json.dump(cfg, open(jp, 'w'), indent=1)
    r = subprocess.run([SCM_HARNESS, jp], capture_output=True, text=True, env=ENV)
    if r.returncode:
        raise RuntimeError(r.stdout[-1500:] + r.stderr[-1500:])
    y = [np.fromfile(o, np.complex64).astype(np.complex128) for o in outs]
    for p in files + outs:
        os.remove(p)
    Nb = len(y[0]) // BLK
    H = np.zeros((Nb, 2, 2), complex)
    for u in range(2):
        for s, f in enumerate((F0, F1)):
            z = y[u][:Nb * BLK] * np.exp(-2j * np.pi * f * n[:Nb * BLK] / FS)
            H[:, u, s] = z.reshape(Nb, BLK).mean(axis=1) * np.sqrt(2)
    return H, json.loads(r.stdout.strip().splitlines()[-1])


def mimo_stats(H):
    p = np.mean(np.abs(H) ** 2)
    Hn = H / np.sqrt(p)
    sv = np.linalg.svd(Hn, compute_uv=False)
    cond = 20 * np.log10(sv[:, 0] / np.maximum(sv[:, 1], 1e-12))
    cap = {}
    for snr in (10, 20):
        r = 10 ** (snr / 10)
        cap[snr] = np.log2(np.prod(1 + r / 2 * sv ** 2, axis=1))
    return sv, cond, cap


CASES = {
    'chem_ones': (dict(pl='NONE', noise=False, inputDbfs=-50.0, coeff=[[1, 1], [1, 1]]), None),
    'chem_ident': (dict(pl='NONE', noise=False, inputDbfs=-50.0, coeff=[[1, 0], [0, 1]]), None),
    'scm_SCM': (dict(pl='NONE', noise=False, inputDbfs=-50.0), scm_ext('SCM')),
    'scm_PROPOSED': (dict(pl='NONE', noise=False, inputDbfs=-50.0), scm_ext('PROPOSED')),
    'scm_JAKES': (dict(pl='NONE', noise=False, inputDbfs=-50.0), scm_ext('JAKES')),
    'scm_TDL_D': (dict(pl='NONE', noise=False, inputDbfs=-50.0), scm_ext('TDL_D')),
    'scm_PROPOSED_4wl': (dict(pl='NONE', noise=False, inputDbfs=-50.0), scm_ext('PROPOSED', spG=4.0)),
}
if XPOL_ONLY:
    CASES = {f'scm_{d}_xpol': (dict(pl='NONE', noise=False, inputDbfs=-50.0), scm_ext(d, xpol=True))
             for d in ('SCM', 'PROPOSED', 'JAKES', 'TDL_D')}
summ = {}
for tag, (chem, ext) in CASES.items():
    H, info = run_scm('e2_' + tag, chem, ext)
    sv, cond, cap = mimo_stats(H)
    # Doppler spectrum of H00 and of the sum over pairs (shows the Doppler option inside SCM)
    f, P = sps.welch(H[:, 0, 0], fs=1e3, window='hann', nperseg=8192, noverlap=4096,
                     return_onesided=False, detrend=False)
    o = np.argsort(f)
    np.savez_compressed(os.path.join(OUT, tag + '.npz'), H=H.astype(np.complex64), f=f[o], P=P[o])
    summ[tag] = dict(cond_med_db=float(np.median(cond)), cond_p90_db=float(np.percentile(cond, 90)),
                     cap10_med=float(np.median(cap[10])), cap20_med=float(np.median(cap[20])),
                     rank2_frac=float(np.mean(cond < 10.0)), wall_s=info['wall_s'],
                     buffers_ext=info['buffers_ext'])
    print(tag, summ[tag], flush=True)

# ---- timeline statistics (model only, many drops) ---------------------------------------------
TL = {}
traj_long = traj_line([-900.0, 50.0, 30.0], [5.0, 0.0, 0.0], 302.0)
ptx = write_traj(tr_gt, os.path.join(RUNDIR, 'tl_tx.csv'))
prx = write_traj(traj_long, os.path.join(RUNDIR, 'tl_rx.csv'))
for dop in ('SCM', 'PROPOSED', 'JAKES', 'TDL_D'):
    for spG in ((0.5, 4.0, 10.0) if not XPOL_ONLY else ('xpol',)):
        tag = f'{dop}_{spG:g}wl' if spG != 'xpol' else f'{dop}_xpol'
        ecfg = scm_ext(dop, spG=spG, dropDistanceM=20.0, seed=11) if spG != 'xpol' else scm_ext(dop, xpol=True, dropDistanceM=20.0, seed=11)
        cfg = dict(model='scm', config=ecfg, fc=FC, fs=11.52e6,
                   rateHz=200, durationS=300, nTx=2, nRx=2, maxTaps=8, trajTx=ptx, trajRx=prx,
                   out=os.path.join(RUNDIR, 'tl.bin'), infoOut=os.path.join(OUT, f'tl_{tag}.csv'))
        jp = os.path.join(RUNDIR, 'tl.json')
        json.dump(cfg, open(jp, 'w'))
        r = subprocess.run([TIMELINE, jp], capture_output=True, text=True, env=ENV)
        if r.returncode:
            raise RuntimeError(r.stderr[-1500:])
        d = np.genfromtxt(os.path.join(OUT, f'tl_{tag}.csv'), delimiter=',', names=True)
        TL[tag] = dict(cap10_med=float(np.median(d['cap10'])), cap20_med=float(np.median(d['cap20'])),
                       cond_med_db=float(np.median(d['cond_dB'])), rank2_frac=float(np.mean(d['cond_dB'] < 10)),
                       drops=int(d['drops'][-1]))
        print('TL', tag, TL[tag], flush=True)
save_json(dict(cases=summ, timeline=TL), os.path.join(OUT, 'e2_summary_xpol.json' if XPOL_ONLY else 'e2_summary.json'))
print('E2 done')
