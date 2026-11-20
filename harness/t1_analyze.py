# SPDX-License-Identifier: AGPL-3.0-or-later
"""Analyse T1 (live twin, deployed CHEM, V-USRP B205): Doppler models played by chem_tv_driver.
usage: t1_analyze.py <T1 results dir>"""
import json
import os
import sys

import numpy as np
from scipy import signal as sps

sys.path.insert(0, os.path.dirname(__file__))
from chemlib import lcr_afd, save_json  # noqa

D = sys.argv[1]
res = {}
for m in ('chem_today', 'proposed', 'jakes', 'tdl_a', 'tdl_d'):
    fn = os.path.join(D, f'T1_{m}_rx.npz')
    if not os.path.exists(fn):
        continue
    z = np.load(fn)
    y = z['y'][0].astype(complex)
    fs = float(z['fs'])
    drv = [json.loads(l) for l in open(os.path.join(D, f'{m}_driver.txt')) if l.startswith('{')]
    t0 = [x['t0_wall'] for x in drv if 't0_wall' in x][0]
    stats = [x for x in drv if 'pushed' in x][-1]
    i0 = int(round((t0 - float(z['wall0'])) * fs)) + int(2 * fs)  # skip 2 s
    seg = y[i0:i0 + int(60 * fs)]
    seg = seg / np.sqrt(np.mean(np.abs(seg) ** 2))
    hd = seg.reshape(-1, 5).mean(axis=1)  # 2 kHz
    fs2 = fs / 5
    f, P = sps.welch(hd, fs=fs2, window='hann', nperseg=16384, noverlap=8192, return_onesided=False, detrend=False)
    o = np.argsort(f)
    f, P = f[o], P[o]
    nu_c = np.sum(f * P) / np.sum(P)
    nu_rms = np.sqrt(np.sum((f - nu_c) ** 2 * P) / np.sum(P))
    env = np.abs(hd)
    lev = lcr_afd(env, fs2, [0.0, -10.0])
    nm = {}
    for lag_ms in (0.5, 1.0, 5.0, 10.0):
        L = int(round(lag_ms * 1e-3 * fs2)) or 1
        e = []
        for b in range(0, len(hd) - 2000, 2000):
            s2 = hd[b:b + 2000 + L]
            ph = np.angle(np.sum(s2[1:] * np.conj(s2[:-1])))
            s2 = s2 * np.exp(-1j * ph * np.arange(len(s2)))
            e.append(np.mean(np.abs(s2[L:] - s2[:-L]) ** 2) / np.mean(np.abs(s2) ** 2))
        nm[lag_ms] = 10 * np.log10(np.mean(e))
    floor = 10 * np.log10(np.sum(P[np.abs(f - nu_c) > 150]) / np.sum(P))
    res[m] = dict(nu_c=nu_c, nu_rms=nu_rms, lcr0=lev[0][0], afd0_ms=1e3 * lev[0][1], lcr10=lev[1][0],
                  afd10_ms=1e3 * lev[1][1], p_fade10=float(np.mean(env ** 2 < 0.1)), nmse_db=nm,
                  floor_beyond_150Hz_db=floor, driver=stats, rx_errs=int(z['errs']))
    np.savez_compressed(os.path.join(D, f'{m}_spec.npz'), f=f, P=P, env=env[::4], fs_env=fs2 / 4)
    print(m, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in res[m].items()}, flush=True)
save_json(res, os.path.join(D, 't1_summary.json'))
