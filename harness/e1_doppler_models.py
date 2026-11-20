# SPDX-License-Identifier: AGPL-3.0-or-later
"""E1 (local CHEM): the four Doppler models of a2g_doppler v2 through the CHEM data path.

Flight: GT (LW1) at 10 m; UAV at 30 m, 5 m/s, radial approach for 60 s (fD ~ +56 Hz), as Fig. 5 of
the Doppler report. A constant (DC) tone goes through chem::Intermediate (real Receiver, worker,
sc16 converters, capture Transmitter). Models:
  chem_today : CHEM core, one link Doppler shift updated every 100 ms, LOS only
  proposed   : a2g_doppler PROPOSED (LOS + ground, K = 8.5 dB, Gaussian diffuse)
  jakes      : a2g_doppler JAKES    (same taps, Rayleigh Jakes on every tap)
  tdl_a      : a2g_doppler TDL_A    (38.901 TDL-A, DS 100 ns, Jakes)
  tdl_d      : a2g_doppler TDL_D    (38.901 TDL-D, DS 100 ns, LOS line at 0.7 fm, K1 = 13.3 dB)
Outputs: Welch PSD of the channel gain, envelope CDF, LCR/AFD, channel-aging NMSE.
"""
import os
import sys

import numpy as np
from scipy import signal as sps

sys.path.insert(0, os.path.dirname(__file__))
from chemlib import *  # noqa

SPEED = float(os.environ.get('SPEED', 5.0))
OUT = os.path.join(RESULTS, 'e1_doppler_models' + ('' if SPEED == 5.0 else f'_{SPEED:g}ms'))
os.makedirs(OUT, exist_ok=True)
FS = 200e3
T = 60.0
tr_gt = traj_static([0.0, 0.0, 10.0])
tr_uav = traj_line([-60.0 - (T + 2.0) * SPEED, 0.0, 30.0] if SPEED != 5.0 else [-320.0, 0.0, 30.0], [SPEED, 0.0, 0.0], T + 2.0)
x = np.ones(int(FS * T), np.complex64)

MODELS = {
    'chem_today': (chem_today(pl='NONE', doppler=True), None),
    'proposed': (chem_today(pl='NONE', doppler=False), ext_a2g(dopplerModel='PROPOSED')),
    'jakes': (chem_today(pl='NONE', doppler=False), ext_a2g(dopplerModel='JAKES')),
    'tdl_a': (chem_today(pl='NONE', doppler=False), ext_a2g(dopplerModel='TDL_A', tdlDelaySpreadNs=100)),
    'tdl_d': (chem_today(pl='NONE', doppler=False), ext_a2g(dopplerModel='TDL_D', tdlDelaySpreadNs=100)),
}

res = {}
for name, (chem, ext) in MODELS.items():
    y, st, info = run_chem(x, FS, tr_gt, tr_uav, chem=chem, ext=ext, pipeline='intermediate',
                           frame=200, tag='e1_' + name, log_ms=100)
    h = y.astype(np.complex128)
    h = h / np.sqrt(np.mean(np.abs(h) ** 2))
    # decimate to 2 kHz for spectra and statistics
    hd = h.reshape(-1, 100).mean(axis=1)
    fs2 = FS / 100
    f, P = sps.welch(hd, fs=fs2, window='hann', nperseg=16384, noverlap=8192, return_onesided=False,
                     detrend=False)
    o = np.argsort(f)
    f, P = f[o], P[o]
    env = np.abs(hd)
    nu_c = np.sum(f * P) / np.sum(P)
    nu_rms = np.sqrt(np.sum((f - nu_c) ** 2 * P) / np.sum(P))
    lev = lcr_afd(env, fs2, [0.0, -10.0])
    # aging NMSE after ideal CFO removal at lag tau (LOS phase from the state log is not needed:
    # remove the mean Doppler of each 1 s block)
    nm = {}
    for lag_ms in (0.5, 1.0, 5.0, 10.0):
        L = int(round(lag_ms * 1e-3 * fs2)) or 1
        e = []
        for b in range(0, len(hd) - 2000, 2000):
            seg = hd[b:b + 2000 + L]
            ph = np.angle(np.sum(seg[1:] * np.conj(seg[:-1])))
            rot = np.exp(-1j * ph * np.arange(len(seg)))
            s2 = seg * rot
            e.append(np.mean(np.abs(s2[L:] - s2[:-L]) ** 2) / np.mean(np.abs(s2) ** 2))
        nm[lag_ms] = 10 * np.log10(np.mean(e))
    res[name] = dict(info=info, nu_c=nu_c, nu_rms=nu_rms, lcr0=lev[0][0], afd0_ms=1e3 * lev[0][1],
                     lcr10=lev[1][0], afd10_ms=1e3 * lev[1][1],
                     p_fade10=float(np.mean(env ** 2 < 0.1)), nmse_db=nm,
                     fD_geo=float(np.median(st['chemDopplerHz'])) if name == 'chem_today' else None)
    np.savez_compressed(os.path.join(OUT, name + '.npz'), f=f, P=P, env=env[::4], fs_env=fs2 / 4, hd=hd.astype(np.complex64))
    print(name, {k: v for k, v in res[name].items() if k != 'info'}, info['wall_s'], flush=True)

save_json(res, os.path.join(OUT, 'e1_summary.json'))
print('E1 done')
