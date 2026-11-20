# SPDX-License-Identifier: AGPL-3.0-or-later
"""Analyse E4 (laptop: public V-USRP + local CHEM). usage: e4_analyze.py <run dir with a2g cases> <run dir with scm cases>"""
import json
import os
import sys

import numpy as np
from scipy import signal as sps

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'harness'))
from chemlib import lcr_afd, save_json  # noqa

A, S = sys.argv[1], sys.argv[2]
OUT = os.path.join(os.path.dirname(__file__), 'results', 'e4_summary.json')
res = {}


def load(d, tag):
    z = np.load(os.path.join(d, tag + '_rx.npz'))
    return z['y'].astype(complex), float(z['fs']), int(z['errs'])


for tag in ('c0_stock', 'a2g_proposed', 'a2g_jakes', 'a2g_tdl_a', 'a2g_tdl_d'):
    y, fs, errs = load(A, tag)
    h = y[0, 0, int(8 * fs):int(32 * fs)]  # 24 s after the extension started
    p = np.abs(h) ** 2
    empty = float(np.mean(p < 0.01 * np.median(p)))
    hd = h[: len(h) // 5 * 5].reshape(-1, 5).mean(1)  # 2 kHz
    hd = hd / np.sqrt(np.mean(np.abs(hd) ** 2))
    f, P = sps.welch(hd, fs=fs / 5, window='hann', nperseg=8192, noverlap=4096, return_onesided=False, detrend=False)
    o = np.argsort(f)
    f, P = f[o], P[o]
    band = np.abs(f) < 150
    c = np.sum(f[band] * P[band]) / np.sum(P[band])
    rms = np.sqrt(np.sum((f[band] - c) ** 2 * P[band]) / np.sum(P[band]))
    lv = lcr_afd(np.abs(hd), fs / 5, [-10.0])
    e5 = []
    for b in range(0, len(hd) - 2010, 2000):
        s2 = hd[b:b + 2010]; ph = np.angle(np.sum(s2[1:] * np.conj(s2[:-1])))
        s2 = s2 * np.exp(-1j * ph * np.arange(len(s2))); e5.append(np.mean(np.abs(s2[10:] - s2[:-10]) ** 2))
    res[tag] = dict(centroid_hz=c, rms_hz_in_150=rms, lcr_m10=lv[0][0], p_fade10=float(np.mean(np.abs(hd) ** 2 < 0.1)),
                    nmse_5ms_db=10 * np.log10(np.mean(e5)), empty_frac=empty, rx_errors=errs,
                    floor_db=10 * np.log10(np.sum(P[~band]) / np.sum(P)))
    np.savez(os.path.join(os.path.dirname(OUT), tag + '_spec.npz'), f=f, P=P)
    print(tag, {k: round(v, 3) if isinstance(v, float) else v for k, v in res[tag].items()}, flush=True)

for tag in ('scm_copol', 'scm_xpol'):
    y, fs, errs = load(S, tag)
    H = np.transpose(y[:, :, int(8 * fs):int(32 * fs)], (2, 0, 1))  # [t, rx port, tx tone/port]
    H = H[: len(H) // 10 * 10].reshape(-1, 10, 2, 2).mean(1)  # 1 ms
    pw = np.sum(np.abs(H) ** 2, axis=(1, 2))
    keep = pw > 0.05 * np.median(pw)
    Hn = H[keep] / np.sqrt(np.mean(np.abs(H[keep]) ** 2))
    sv = np.linalg.svd(Hn, compute_uv=False)
    cond = 20 * np.log10(sv[:, 0] / np.maximum(sv[:, 1], 1e-12))
    cap = np.log2(np.prod(1 + 100 / 2 * sv ** 2, axis=1))
    res[tag] = dict(cond_med_db=float(np.median(cond)), cond_p90_db=float(np.percentile(cond, 90)),
                    rank2_frac=float(np.mean(cond < 10)), cap20_med=float(np.median(cap)),
                    empty_frac=float(1 - np.mean(keep)), rx_errors=errs,
                    cross_to_direct_db=float(10 * np.log10(np.mean(np.abs(Hn[:, 0, 1]) ** 2 + np.abs(Hn[:, 1, 0]) ** 2) /
                                                           np.mean(np.abs(Hn[:, 0, 0]) ** 2 + np.abs(Hn[:, 1, 1]) ** 2))))
    np.savez(os.path.join(os.path.dirname(OUT), tag + '_H.npz'), H=Hn.astype(np.complex64))
    print(tag, {k: round(v, 3) if isinstance(v, float) else v for k, v in res[tag].items()}, flush=True)
save_json(res, OUT)
