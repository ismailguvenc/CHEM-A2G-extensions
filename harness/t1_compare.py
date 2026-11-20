# SPDX-License-Identifier: AGPL-3.0-or-later
"""Compare the Doppler models on the live twin (T1) with the local CHEM (E1) at the same speed,
on 100 Hz channel samples (the twin applies one CIR per 10 ms V-USRP buffer).
usage: t1_compare.py <T1 dir> <E1 dir>"""
import json, os, sys
import numpy as np
from scipy import signal as sps
sys.path.insert(0, os.path.dirname(__file__))
from chemlib import lcr_afd, save_json  # noqa
T1, E1 = sys.argv[1], sys.argv[2]


def metrics(h100):
    h = h100 / np.sqrt(np.mean(np.abs(h100) ** 2))
    f, P = sps.welch(h, fs=100.0, window='hann', nperseg=1024, noverlap=512, return_onesided=False, detrend=False)
    o = np.argsort(f); f, P = f[o], P[o]
    c = np.sum(f * P) / np.sum(P); rms = np.sqrt(np.sum((f - c) ** 2 * P) / np.sum(P))
    env = np.abs(h); lv = lcr_afd(env, 100.0, [0.0, -10.0])
    out = dict(nu_c=c, nu_rms=rms, lcr0=lv[0][0], afd0_ms=1e3 * lv[0][1], lcr10=lv[1][0], p_fade10=float(np.mean(env ** 2 < 0.1)))
    for L in (1, 2):  # 10 and 20 ms
        e = []
        for b in range(0, len(h) - 100, 100):
            s = h[b:b + 100 + L]; ph = np.angle(np.sum(s[1:] * np.conj(s[:-1])))
            s = s * np.exp(-1j * ph * np.arange(len(s))); e.append(np.mean(np.abs(s[L:] - s[:-L]) ** 2))
        out[f'nmse_{10 * L}ms_db'] = 10 * np.log10(np.mean(e))
    return out, f, P


res = {}
for m in ('chem_today', 'proposed', 'jakes', 'tdl_a', 'tdl_d'):
    hd = np.load(os.path.join(E1, m + '.npz'))['hd'].astype(complex)          # 2 kHz
    loc, fl, Pl = metrics(hd[: len(hd) // 20 * 20].reshape(-1, 20).mean(1))
    z = np.load(os.path.join(T1, f'T1_{m}_rx.npz')); y = z['y'][0].astype(complex); fs = float(z['fs'])
    drv = [json.loads(l) for l in open(os.path.join(T1, f'{m}_driver.txt')) if l.startswith('{')]
    t0 = [x['t0_wall'] for x in drv if 't0_wall' in x][0]
    i0 = int(round((t0 - float(z['wall0'])) * fs)) + int(2 * fs)
    # align the 10 ms blocks with CHEM's CIR update instants (minimum within-block variance)
    best = None
    for o in range(0, 100, 2):
        sg = y[i0 + o:i0 + o + int(20 * fs)]
        B = sg[: len(sg) // 100 * 100].reshape(-1, 100)
        v = np.mean(np.abs(B[:, 10:90] - B[:, 10:90].mean(1, keepdims=True)) ** 2)
        if best is None or v < best[0]:
            best = (v, o)
    o = best[1]
    seg = y[i0 + o:i0 + o + int(60 * fs)]
    h100 = seg[: len(seg) // 100 * 100].reshape(-1, 100)[:, 10:90].mean(1)
    pw = np.abs(h100) ** 2
    drop = pw < 0.02 * np.median(pw)                     # empty buffers (dropped on the twin)
    tw_raw, _, _ = metrics(h100)
    idx = np.arange(len(h100))
    h100c = h100.copy()
    h100c[drop] = np.interp(idx[drop], idx[~drop], h100[~drop].real) + 1j * np.interp(idx[drop], idx[~drop], h100[~drop].imag)
    tw, ft, Pt = metrics(h100c)
    tw['block_offset'] = o
    tw['dropped_block_frac'] = float(np.mean(drop))
    tw['raw'] = tw_raw
    res[m] = dict(local=loc, twin=tw)
    np.savez(os.path.join(T1, m + '_cmp100.npz'), fl=fl, Pl=Pl, ft=ft, Pt=Pt)
    print(m, ' local:', {k: round(v, 2) for k, v in loc.items()}, '\n', ' ' * len(m), 'twin :', {k: round(v, 3) for k, v in tw.items() if k != 'raw'}, flush=True)
save_json(res, os.path.join(T1, 't1_compare100.json'))
