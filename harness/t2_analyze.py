# SPDX-License-Identifier: AGPL-3.0-or-later
"""Analyse T2 (live twin LTE 2x2 TM4): eNB console metrics (CQI, RI, MCS, rate, BLER), UE CSV, iPerf3.
usage: t2_analyze.py <T2 results dir>"""
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from chemlib import save_json  # noqa

D = sys.argv[1]
ANSI = re.compile(r'\x1b\[[0-9;]*m')


def num(t):
    t = t.strip().rstrip('%')
    m = {'k': 1e3, 'M': 1e6, 'G': 1e9}
    try:
        if t and t[-1] in m:
            return float(t[:-1]) * m[t[-1]]
        return float(t)
    except ValueError:
        return np.nan


def parse_enb(fn):
    rows = []
    hdr = None
    for line in open(fn, errors='ignore'):
        line = ANSI.sub('', line).replace('|', ' ')
        tk = line.split()
        if 'cqi' in tk and 'ri' in tk:
            hdr = tk
            continue
        if hdr is None or len(tk) < len(hdr) - 2:
            continue
        try:
            i = hdr.index('cqi')
            off = len(tk) - len(hdr)  # some columns (rat/pci) may be absent
            j = i + off if off < 0 else i
            cqi, ri, mcs, brate, ok, nok = (num(tk[j + k]) for k in range(6))
            pct = num(tk[j + 6])
        except (IndexError, ValueError):
            continue
        if np.isnan(cqi) or np.isnan(ri):
            continue
        rows.append((cqi, ri, mcs, brate, ok, nok, pct))
    return np.array(rows) if rows else np.zeros((0, 7))


summ = {}
for it in sorted(glob.glob(os.path.join(D, '*_iperf.json'))):
    S = os.path.basename(it)[:-len('_iperf.json')]
    case = re.match(r'T2_(.+)_s\d+_\d+$', S).group(1)
    try:
        ip = json.load(open(it))
        mbps = ip['end']['sum']['bits_per_second'] / 1e6
        loss = ip['end']['sum']['lost_percent']
        per_s = [x['sum']['bits_per_second'] / 1e6 for x in ip['intervals']]
    except Exception:
        mbps, loss, per_s = np.nan, np.nan, []
    e = parse_enb(os.path.join(D, S + '_enb.txt')) if os.path.exists(os.path.join(D, S + '_enb.txt')) else np.zeros((0, 7))
    act = e[e[:, 3] > 1e5] if len(e) else e  # reports with DL traffic
    r = dict(iperf_offered_mbps=mbps, loss_pct=loss, iperf_mbps=mbps * (1 - loss / 100) if loss == loss else np.nan, n_reports=int(len(act)))
    if len(act):
        r.update(cqi_med=float(np.median(act[:, 0])), ri2_frac=float(np.mean(act[:, 1] >= 2)),
                 mcs_med=float(np.median(act[:, 2])),
                 goodput_med_mbps=float(np.median(act[:, 3] * act[:, 4] / np.maximum(act[:, 4] + act[:, 5], 1)) / 1e6), dl_brate_med_mbps=float(np.median(act[:, 3]) / 1e6),
                 bler_med=float(np.nanmedian(act[:, 6])), ri_series=act[:, 1].tolist(), cqi_series=act[:, 0].tolist(),
                 brate_series=(act[:, 3] / 1e6).tolist())
    ue = os.path.join(D, S + '_ue_metrics.csv')
    if False and os.path.exists(ue):
        try:
            u = np.genfromtxt(ue, delimiter=';', names=True, invalid_raise=False)
            for k in ('dl_snr', 'dl_mcs', 'dl_bler'):
                if k in u.dtype.names:
                    v = u[k][np.isfinite(u[k])]
                    if len(v):
                        r['ue_' + k + '_med'] = float(np.median(v))
        except Exception as ex:
            r['ue_err'] = str(ex)
    r['iperf_series'] = per_s
    summ[case] = r
    print(case, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items() if not k.endswith('series')}, flush=True)
save_json(summ, os.path.join(D, 't2_summary.json'))
