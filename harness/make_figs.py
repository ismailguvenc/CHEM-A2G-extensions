# SPDX-License-Identifier: AGPL-3.0-or-later
"""Figures for the CHEM MIMO and Doppler models report (Nov 2026).
usage: make_figs.py <T1 dir or -> <T2 dir or ->"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from chemlib import PAL, RESULTS, mpl_setup  # noqa

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa
from matplotlib.patches import FancyBboxPatch  # noqa

mpl_setup()
FIG = os.path.join(os.path.dirname(__file__), '..', 'report', 'figures')
os.makedirs(FIG, exist_ok=True)
T1 = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != '-' else None
T2 = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != '-' else None
NAMES = {'chem_today': 'CHEM today', 'proposed': 'Proposed', 'jakes': 'Jakes (every tap)',
         'tdl_a': 'TDL-A', 'tdl_d': 'TDL-D'}


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(FIG, name + '.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)


# ---- F1 architecture ------------------------------------------------------------------------
def fig_arch():
    fig, ax = plt.subplots(figsize=(10, 4.6))
    ax.set_xlim(0, 100); ax.set_ylim(0, 50); ax.axis('off')

    def box(x, y, w, h, txt, c, fs=7.8, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=1.2', fc=c, ec='#444', lw=0.8))
        ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=fs, weight='bold' if bold else None)

    def arr(x0, y0, x1, y1):
        ax.annotate('', (x1, y1), (x0, y0), arrowprops=dict(arrowstyle='->', lw=0.9, color='#333'))
    box(1, 36, 22, 10, 'a2g model library\n(trajectories, geometry,\nSoS diffuse, Doppler presets)', '#dbe9f9')
    box(1, 21, 22, 10, 'scm model library\n(TR 25.996 drops, ULA / x-pol,\nDoppler options)', '#dbe9f9')
    box(30, 36, 20, 10, 'a2g_doppler (SISO)\nPROPOSED | JAKES |\nTDL_A | TDL_D', '#fde3cf', bold=False)
    box(30, 21, 20, 10, 'scm_mimo (MIMO)\nSCM | PROPOSED |\nJAKES | TDL_D', '#fde3cf')
    box(30, 5, 20, 10, 'chem_taps_timeline\n(1 kHz CIR per port pair)', '#e7f5ec')
    box(57, 28, 19, 18, 'local CHEM build\n(public 6b9e58f)\nCHEM-EXT-HOOK\nCHEM-EXT-MIMO (port sum)\n+ sionna client', '#f1f1f1')
    box(57, 5, 19, 12, 'chem_tv_driver.py\n(M1, CHG_CIR per\nport pair, FIXED_SNR)', '#e7f5ec')
    box(82, 28, 17, 18, 'offline harnesses\na2g_harness, scm_harness\n(Intermediate path)\nchem_sionna_probe', '#f1f1f1')
    box(82, 5, 17, 12, 'deployed CHEM v1.2.4\n+ V-USRP B205/B210\n(AERPAW X1546)', '#fff6c8')
    ax.text(66.5, 48.5, 'Sionna RT server (AERPAW-DT-SIONNA-EXTENSION, aerpaw scene) → CHEM sionna client',
            ha='center', fontsize=8, style='italic')
    for y in (41, 26):
        arr(23.6, y, 29.4, y)
    arr(50.6, 41, 56.4, 38); arr(50.6, 26, 56.4, 33); arr(40, 20.4, 40, 15.6)
    arr(50.6, 10, 56.4, 10); arr(76.6, 10, 81.4, 10); arr(76.6, 37, 81.4, 37)
    ax.text(40, 1.5, 'same C++ model in both paths', ha='center', fontsize=8, color='#2a7a4a')
    save(fig, 'fig01_architecture')


# ---- F2/F3 Doppler models: local CHEM (E1) and live twin (T1) ------------------------------------
def fig_doppler():
    e1 = os.path.join(RESULTS, 'e1_doppler_models')
    fig, axs = plt.subplots(1, 5, figsize=(13, 2.9), sharey=True)
    for ax, m in zip(axs, NAMES):
        d = np.load(os.path.join(e1, m + '.npz'))
        P = d['P'] / d['P'].max()
        ax.plot(d['f'], 10 * np.log10(P + 1e-12), color=PAL[0], lw=0.9, label='local CHEM (E1)')
        if False:
            t = np.load(os.path.join(T1, m + '_spec.npz'))
            Pt = t['P'] / t['P'].max()
            ax.plot(t['f'], 10 * np.log10(Pt + 1e-12), color=PAL[1], lw=0.9, alpha=0.85, label='live twin (T1)')
        ax.axvline(56.3, color='#999', ls=':', lw=0.8)
        ax.axvline(-56.7, color='#ccc', ls=':', lw=0.6)
        ax.set_xlim(-120, 120); ax.set_ylim(-60, 3)
        ax.set_title(NAMES[m], fontsize=9); ax.set_xlabel('Doppler (Hz)')
    axs[0].set_ylabel('PSD (dB, norm.)')
    axs[0].legend(fontsize=7, loc='lower left')
    save(fig, 'fig02_doppler_psd')

    s1 = json.load(open(os.path.join(e1, 'e1_summary.json')))
    st = {}
    fig, axs = plt.subplots(1, 3, figsize=(11, 2.9))
    x = np.arange(len(NAMES)); w = 0.38
    for k, (ax, key, lab, logy) in enumerate(zip(axs, ('nu_rms', 'lcr0', 'nmse05'),
                                              ('rms Doppler spread (Hz)', 'LCR at rms level (1/s)',
                                               'aging NMSE at 0.5 ms (dB)'), (False, False, False))):
        def val(s, m):
            if key == 'nmse05':
                return s[m]['nmse_db']['0.5'] if '0.5' in s[m]['nmse_db'] else s[m]['nmse_db'][0.5]
            return s[m][key]
        ax.bar(x - w / 2, [val(s1, m) for m in NAMES], w, color=PAL[0], label='local CHEM (E1)')
        if st:
            ax.bar(x + w / 2, [val(st, m) if m in st else np.nan for m in NAMES], w, color=PAL[1], label='live twin (T1)')
        ax.set_xticks(x); ax.set_xticklabels([NAMES[m].replace(' (every tap)', '') for m in NAMES], rotation=25, fontsize=8)
        ax.set_title(lab, fontsize=9)
    axs[0].legend(fontsize=7)
    save(fig, 'fig03_doppler_stats')


# ---- F4/F5 SCM MIMO (E2) ----------------------------------------------------------------------
def fig_mimo():
    e2 = os.path.join(RESULTS, 'e2_scm_mimo_local')
    cases = [('chem_ones', 'CHEM default [[1,1],[1,1]]', '#888', '-'), ('scm_PROPOSED', 'SCM A2G, co-pol 0.5λ', PAL[0], '-'),
             ('scm_PROPOSED_4wl', 'SCM A2G, co-pol 4λ (GT)', PAL[2], '-'),
             ('scm_PROPOSED_xpol', 'SCM A2G, x-pol, PROPOSED', PAL[1], '-'),
             ('scm_SCM_xpol', 'x-pol, SCM Doppler', PAL[3], '--'), ('scm_JAKES_xpol', 'x-pol, JAKES', PAL[4], '--'),
             ('scm_TDL_D_xpol', 'x-pol, TDL_D', PAL[5], '--'), ('chem_ident', 'CHEM identity', '#000', ':')]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
    for tag, lab, c, ls in cases:
        fn = os.path.join(e2, tag + '.npz')
        if not os.path.exists(fn):
            continue
        H = np.load(fn)['H'].astype(complex)
        H = H / np.sqrt(np.mean(np.abs(H) ** 2))
        sv = np.linalg.svd(H, compute_uv=False)
        cond = np.clip(20 * np.log10(sv[:, 0] / np.maximum(sv[:, 1], 1e-12)), 0, 60)
        cap = np.log2(np.prod(1 + 100 / 2 * sv ** 2, axis=1))
        for ax, v in zip(axs, (cond, cap)):
            vs = np.sort(v)
            ax.plot(vs, np.linspace(0, 1, len(vs)), color=c, ls=ls, lw=1.2, label=lab)
    axs[0].set_xlabel('condition number (dB, clipped at 60)'); axs[0].set_ylabel('CDF')
    axs[1].set_xlabel('capacity at 20 dB SNR (bit/s/Hz)')
    axs[0].axvline(10, color='#bbb', lw=0.7)
    axs[1].legend(fontsize=7, loc='upper left')
    save(fig, 'fig04_mimo_cdf')

    fig, axs = plt.subplots(1, 4, figsize=(12, 2.8), sharey=True)
    for ax, d in zip(axs, ('SCM', 'PROPOSED', 'JAKES', 'TDL_D')):
        z = np.load(os.path.join(e2, f'scm_{d}_xpol.npz'))
        P = z['P'] / z['P'].max()
        ax.plot(z['f'], 10 * np.log10(P + 1e-12), color=PAL[0], lw=0.9)
        ax.axvline(56.3, color='#999', ls=':', lw=0.8)
        ax.set_xlim(-120, 120); ax.set_ylim(-50, 3)
        ax.set_title(f'scm_mimo, Doppler {d}', fontsize=9); ax.set_xlabel('Doppler (Hz)')
    axs[0].set_ylabel('PSD of $h_{00}$ (dB)')
    save(fig, 'fig05_scm_doppler')


# ---- F6 LTE TM4 on the twin (T2) -----------------------------------------------------------------
def fig_lte():
    if not T2:
        return
    js = os.path.join(T2, 't2_summary.json')
    if not os.path.exists(js):
        return
    s = json.load(open(js))
    order = [c for c in ('default', 'identity', 'prop_copol', 'prop_xpol', 'scm_xpol', 'jakes_xpol', 'tdld_xpol') if c in s]
    lab = {'default': 'CHEM\ndef.', 'identity': 'ident.', 'prop_copol': 'co-pol\nPROP',
           'prop_xpol': 'x-pol\nPROP', 'scm_xpol': 'x-pol\nSCM', 'jakes_xpol': 'x-pol\nJAKES',
           'tdld_xpol': 'x-pol\nTDL-D'}
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.4))
    fig.subplots_adjust(wspace=0.35)
    x = np.arange(len(order))
    axs[0].bar(x, [s[c]['goodput_med_mbps'] for c in order], color=PAL[0])
    axs[0].axhline(30, color='#999', ls=':', lw=0.8)
    axs[0].set_ylabel('DL goodput, eNB MAC (Mbit/s, median)')
    r2 = [s[c].get('ri2_frac', np.nan) for c in order]
    axs[1].bar(x, r2, color=PAL[1]); axs[1].set_ylabel('fraction of eNB reports with RI = 2'); axs[1].set_ylim(0, 1.05)
    axs[2].bar(x - 0.2, [s[c].get('cqi_med', np.nan) for c in order], 0.4, color=PAL[2], label='CQI (median)')
    axs[2].bar(x + 0.2, [s[c].get('mcs_med', np.nan) for c in order], 0.4, color=PAL[3], label='DL MCS (median)')
    ax2 = axs[2].twinx(); ax2.plot(x, [s[c].get('bler_med', np.nan) for c in order], 'ko-', ms=3, lw=0.8)
    ax2.set_ylabel('DL BLER (%), median'); axs[2].legend(fontsize=7, loc='upper left')
    for ax in axs:
        ax.set_xticks(x); ax.set_xticklabels([lab[c] for c in order], fontsize=7)
    save(fig, 'fig06_lte_tm4')


def fig_twin_doppler():
    if not T1 or not os.path.exists(os.path.join(T1, 't1_compare100.json')):
        return
    c = json.load(open(os.path.join(T1, 't1_compare100.json')))
    fig, axs = plt.subplots(1, 5, figsize=(13, 2.9), sharey=True)
    for ax, m in zip(axs, NAMES):
        z = np.load(os.path.join(T1, m + '_cmp100.npz'))
        ax.plot(z['fl'], 10 * np.log10(z['Pl'] / z['Pl'].max() + 1e-12), color=PAL[0], lw=1.0, label='local CHEM')
        ax.plot(z['ft'], 10 * np.log10(z['Pt'] / z['Pt'].max() + 1e-12), color=PAL[1], lw=1.0, alpha=0.85, label='live twin')
        ax.axvline(22.6, color='#999', ls=':', lw=0.8)
        ax.set_xlim(-50, 50); ax.set_ylim(-40, 3); ax.set_title(NAMES[m], fontsize=9); ax.set_xlabel('Doppler (Hz)')
    axs[0].set_ylabel('PSD (dB, norm.)'); axs[0].legend(fontsize=7, loc='lower left')
    save(fig, 'fig07_twin_doppler_psd')
    fig, axs = plt.subplots(1, 3, figsize=(11, 2.8))
    x = np.arange(len(NAMES)); w = 0.38
    for ax, key, lab in zip(axs, ('nu_rms', 'lcr0', 'nmse_10ms_db'), ('rms Doppler spread (Hz)', 'LCR at rms level (1/s)', 'aging NMSE at 10 ms (dB)')):
        ax.bar(x - w / 2, [c[m]['local'][key] for m in NAMES], w, color=PAL[0], label='local CHEM')
        ax.bar(x + w / 2, [c[m]['twin']['raw'][key] if 'raw' in c[m]['twin'] else c[m]['twin'][key] for m in NAMES], w, color=PAL[1], label='live twin')
        ax.set_xticks(x); ax.set_xticklabels([NAMES[m].replace(' (every tap)', '') for m in NAMES], rotation=25, fontsize=8)
        ax.set_title(lab, fontsize=9)
    axs[0].legend(fontsize=7)
    save(fig, 'fig08_twin_doppler_stats')


def fig_sionna():
    e3 = os.path.join(RESULTS, 'e3_sionna_chem', 'e3_compare.npz')
    if not os.path.exists(e3):
        return
    z = np.load(e3)
    pr = np.genfromtxt(os.path.join(RESULTS, 'e3_sionna_chem', 'probe.csv'), delimiter=',', names=True)
    a = pr['active'][:len(z['t'])] > 0
    fig, ax = plt.subplots(figsize=(6.5, 3.0))
    ax.plot(z['t'], z['rt_tot_db'], color='#999', lw=1, label='offline RT: sum of path powers')
    ax.plot(z['t'], z['rt_los_db'], color=PAL[2], lw=1, ls='--', label='offline RT: LOS path')
    ax.plot(z['t'][a], z['chem_db'][a], 'o', color=PAL[1], ms=3.5, label='CIR energy applied by local CHEM (sionna ext.)')
    ax.axvspan(0, z['t'][a][0] - 1, color='#eee', zorder=0)
    ax.text(2, ax.get_ylim()[0] + 1, 'scene creation (unit CIR)', fontsize=7)
    ax.set_xlabel('flight time (s)'); ax.set_ylabel('channel gain (dB)'); ax.legend(fontsize=7)
    save(fig, 'fig09_sionna_chem')


def fig_vusrp():
    V = os.path.join(os.path.dirname(__file__), '..', 'vusrp_local', 'results')
    if not os.path.exists(os.path.join(V, 'e4_summary.json')):
        return
    e1 = os.path.join(RESULTS, 'e1_doppler_models')
    fig, axs = plt.subplots(1, 4, figsize=(12.5, 2.9), sharey=True)
    for ax, (m, tag) in zip(axs, (('proposed', 'a2g_proposed'), ('jakes', 'a2g_jakes'), ('tdl_a', 'a2g_tdl_a'), ('tdl_d', 'a2g_tdl_d'))):
        d = np.load(os.path.join(e1, m + '.npz'))
        ax.plot(d['f'], 10 * np.log10(d['P'] / d['P'].max() + 1e-12), color=PAL[0], lw=0.9, label='local CHEM, IQ files (E1)')
        z = np.load(os.path.join(V, tag + '_spec.npz'))
        ax.plot(z['f'], 10 * np.log10(z['P'] / z['P'].max() + 1e-12), color=PAL[2], lw=0.9, alpha=0.85, label='local CHEM + public V-USRP (E4)')
        ax.axvline(56.3, color='#999', ls=':', lw=0.8)
        ax.set_xlim(-120, 120); ax.set_ylim(-60, 3); ax.set_title(NAMES[m], fontsize=9); ax.set_xlabel('Doppler (Hz)')
    axs[0].set_ylabel('PSD (dB, norm.)'); axs[0].legend(fontsize=6.5, loc='lower left')
    save(fig, 'fig10_vusrp_doppler')


fig_arch()
fig_vusrp()
fig_twin_doppler()
fig_sionna()
fig_doppler()
fig_mimo()
fig_lte()
print('figures ->', os.path.abspath(FIG))
