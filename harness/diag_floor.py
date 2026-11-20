# SPDX-License-Identifier: AGPL-3.0-or-later
"""Find the source of the broadband floor of the a2g hook in E1."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from chemlib import *  # noqa
FS = 200e3; T = 10.0
tr_gt = traj_static([0.0, 0.0, 10.0]); tr_uav = traj_line([-320.0, 0.0, 30.0], [5.0, 0.0, 0.0], T + 2.0)
x = np.ones(int(FS * T), np.complex64)
def floor_db(y, f0):
    h = y.astype(complex); h /= np.sqrt(np.mean(abs(h) ** 2))
    n = np.arange(len(h)); z = h * np.exp(-2j * np.pi * f0 * n / FS)
    Z = np.fft.fft(z * np.hanning(len(z))); f = np.fft.fftfreq(len(z), 1 / FS); P = abs(Z) ** 2; P /= P.sum()
    m = abs(f) > 200
    top = np.argsort(P * m)[-4:]
    return round(10 * np.log10(P[m].sum()), 1), [(round(f[i], 1), round(10 * np.log10(P[i]), 1)) for i in top]
for tag, ext, pipe, frame in [
        ('prop_int_m50', ext_a2g(), 'intermediate', 200),
        ('prop_int_noLog', ext_a2g(paths=[dict(type='LOS')]), 'intermediate', 200),
        ('prop_int_f1000', ext_a2g(paths=[dict(type='LOS')]), 'intermediate', 1000),
        ('prop_int_upd5', ext_a2g(paths=[dict(type='LOS')], updateMs=5.0), 'intermediate', 200)]:
    y, st, info = run_chem(x, FS, tr_gt, tr_uav, chem=chem_today(pl='NONE', doppler=False), ext=ext,
                           pipeline=pipe, frame=frame, tag='diag_' + tag, log_ms=(1e6 if 'noLog' in tag else 10.0))
    print(tag, floor_db(y, 56.3), flush=True)
