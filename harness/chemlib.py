# SPDX-License-Identifier: AGPL-3.0-or-later
"""Shared code for the CHEM a2g_doppler test harness (runs in WSL, env 'chem').

- run_chem(): IQ -> a2g_harness (CHEM Channel / Intermediate + extension) -> IQ
- trajectory builders (ENU, GT at the origin)
- model presets: CHEM today, CHEM 2_RAY, CHEM 3GPP, Jakes TDL, TDL-D style,
  proposed A2G_MEAS (LOS + ground + scatterers, per-tap Doppler)
- estimators shared by the measured data and the emulated data
"""
import json
import os
import subprocess
import time

import numpy as np

C0 = 299792458.0
FC = 3.4e9
LAM = C0 / FC
HERE = os.path.dirname(os.path.abspath(__file__))
OUTROOT = os.path.abspath(os.path.join(HERE, '..'))
RESULTS = os.path.join(OUTROOT, 'results')
HARNESS = os.path.expanduser(os.environ.get('A2G_HARNESS', '~/work/build_v2/bin/Release/a2g_harness'))
SCM_HARNESS = os.path.expanduser(os.environ.get('SCM_HARNESS', '~/work/build_v2/bin/Release/scm_harness'))
TIMELINE = os.path.expanduser(os.environ.get('CHEM_TIMELINE', '~/work/build_v2/bin/Release/chem_taps_timeline'))
RUNDIR = os.path.expanduser('~/work/runs_v2')
os.makedirs(RUNDIR, exist_ok=True)
os.makedirs(RUNDIR, exist_ok=True)
LW1 = (35.727489, -78.695957)


def mpl_setup():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 9,
                         'axes.grid': True, 'grid.alpha': 0.3, 'figure.dpi': 150,
                         'savefig.bbox': 'tight', 'legend.fontsize': 7.5})
    return plt


PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948', '#555555']


# ---------------------------------------------------------------------------
# Trajectories (ENU metres; columns t, e, n, u, ve, vn, vu)
# ---------------------------------------------------------------------------
def traj_line(p0, v, T, dt=0.05, t0=-1.0):
    t = np.arange(t0, T + dt, dt)
    p = np.asarray(p0)[None, :] + t[:, None] * np.asarray(v)[None, :]
    return np.column_stack([t, p, np.repeat(np.asarray(v, float)[None, :], len(t), 0)])


def traj_static(p):
    return list(map(float, p))


def traj_waypoints(wps, v, amax, dt=0.01, hover=1.0, smooth_s=0.6):
    """Waypoint mission, trapezoidal speed per leg, Gaussian-smoothed velocity
    (same construction as report_doppler/sim/common.py make_traj)."""
    wps = np.asarray(wps, float)
    vel = [np.zeros((int(hover / dt), 3))]
    for i in range(len(wps) - 1):
        d = wps[i + 1] - wps[i]; L = np.linalg.norm(d)
        if L >= v * v / amax:
            ta, tc, vp = v / amax, (L - v * v / amax) / v, v
        else:
            vp = np.sqrt(L * amax); ta, tc = vp / amax, 0.0
        tt = np.arange(0, 2 * ta + tc, dt)
        s = np.where(tt < ta, amax * tt, np.where(tt < ta + tc, vp, np.maximum(vp - amax * (tt - ta - tc), 0)))
        vel.append(s[:, None] * (d / L)[None, :])
    vel.append(np.zeros((int(hover / dt), 3)))
    vel = np.concatenate(vel)
    sig = smooth_s / dt; m = int(4 * sig)
    k = np.exp(-0.5 * (np.arange(-m, m + 1) / sig) ** 2); k /= k.sum()
    vel = np.stack([np.convolve(np.pad(vel[:, j], m, mode='edge'), k, 'valid') for j in range(3)], 1)
    pos = wps[0] + np.concatenate([np.zeros((1, 3)), np.cumsum(0.5 * (vel[1:] + vel[:-1]) * dt, 0)])
    t = np.arange(len(vel)) * dt
    return np.column_stack([t, pos, vel])


def latlon_to_enu(lat, lon, lat0=LW1[0], lon0=LW1[1]):
    m = 111319.5
    e = (lon - lon0) * m * np.cos(np.deg2rad(0.5 * (lat + lat0)))
    n = (lat - lat0) * m
    return e, n


def write_traj(arr, path):
    if isinstance(arr, list):
        return arr
    np.savetxt(path, arr, delimiter=',', header='t,east,north,up,ve,vn,vu', comments='', fmt='%.6f')
    return path


def geometry(tr_tx, tr_rx, t):
    """positions/velocities of both ends at times t (linear interp of tables)."""
    def at(tr):
        if isinstance(tr, list):
            p = np.repeat(np.asarray(tr, float)[None, :], len(t), 0)
            return p, np.zeros_like(p)
        p = np.stack([np.interp(t, tr[:, 0], tr[:, 1 + j]) for j in range(3)], 1)
        v = np.stack([np.interp(t, tr[:, 0], tr[:, 4 + j]) for j in range(3)], 1)
        return p, v
    pa, va = at(tr_tx); pb, vb = at(tr_rx)
    r = pb - pa; d = np.linalg.norm(r, axis=1); u = r / d[:, None]
    fD = -np.sum((vb - va) * u, 1) / LAM
    return dict(d=d, fD=fD, pa=pa, pb=pb, va=va, vb=vb)


def path_doppler(tr_tx, tr_rx, t, point=None, ground=False, dt=1e-3):
    """Geometric Doppler of a scatterer / ground path (numerical derivative)."""
    def L(tt):
        g = geometry(tr_tx, tr_rx, tt)
        if ground:
            img = g['pb'].copy(); img[:, 2] *= -1
            return np.linalg.norm(img - g['pa'], axis=1)
        if point is not None:
            s = np.asarray(point)[None, :]
            return np.linalg.norm(g['pa'] - s, axis=1) + np.linalg.norm(s - g['pb'], axis=1)
        return g['d']
    return -(L(t + dt) - L(t - dt)) / (2 * dt) / LAM, L(t)


# ---------------------------------------------------------------------------
# Model presets
# ---------------------------------------------------------------------------
def chem_today(pl='FREE_SPACE', cir=None, doppler=True, **kw):
    c = dict(pl=pl, doppler=doppler, shadowStd=kw.get('shadowStd', 0.0), groundCoeff=kw.get('groundCoeff', -0.3),
             noise=False, channelUpdateMs=kw.get('channelUpdateMs', 100.0), inputDbfs=kw.get('inputDbfs', -50.0),
             scenario3gpp=kw.get('scenario3gpp', 'RMa'))
    if cir is not None:
        c['cir'] = [[float(np.real(v)), float(np.imag(v))] for v in cir]
    return c


def ext_a2g(paths=None, K_dB=8.5, shape='GAUSS', mode='hook', **kw):
    e = dict(mode=mode, updateMs=kw.pop('updateMs', 0.5), K_dB=K_dB, shape=shape,
             taps=kw.pop('taps', 'ALL'), seed=kw.pop('seed', 1234),
             paths=paths if paths is not None else [dict(type='LOS'), dict(type='GROUND', reflCoeff=-0.3)])
    e.update(kw)
    return e


# ---------------------------------------------------------------------------
# Running the harness
# ---------------------------------------------------------------------------
def run_chem(x, fs, tr_tx, tr_rx, chem=None, ext=None, pipeline='channel', frame=None,
             tag='run', fc=FC, log_ms=10.0, keep=False):
    tag = tag.replace(' ', '_')
    frame = frame or int(round(fs * 1e-3))
    xin = os.path.join(RUNDIR, tag + '_in.cf32')
    xout = os.path.join(RUNDIR, tag + '_out.cf32')
    st = os.path.join(RUNDIR, tag + '_state.csv')
    n = (len(x) // frame) * frame
    np.asarray(x[:n], np.complex64).tofile(xin)
    cfg = dict(**{'in': xin, 'out': xout, 'fs': fs, 'fc': fc, 'frameSamples': frame, 'pipeline': pipeline,
                  'trajTx': write_traj(tr_tx, os.path.join(RUNDIR, tag + '_tx.csv')),
                  'trajRx': write_traj(tr_rx, os.path.join(RUNDIR, tag + '_rx.csv')),
                  'chem': chem or chem_today(), 'stateOut': st, 'logStateEveryMs': log_ms})
    if ext is not None:
        cfg['ext'] = ext
    jp = os.path.join(RUNDIR, tag + '.json')
    json.dump(cfg, open(jp, 'w'), indent=1)
    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = os.path.expanduser('~/deps/mavsdk_x/usr/lib') + ':' + \
        os.path.expanduser('~/mm/root/envs/chembuild/lib') + ':' + env.get('LD_LIBRARY_PATH', '')
    t0 = time.time()
    r = subprocess.run([HARNESS, jp], capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError('harness failed: ' + r.stdout[-2000:] + r.stderr[-2000:])
    info = json.loads(r.stdout.strip().splitlines()[-1])
    info['py_wall_s'] = time.time() - t0
    y = np.fromfile(xout, np.complex64)
    state = np.genfromtxt(st, delimiter=',', names=True)
    if not keep:
        for p in (xin, xout):
            try:
                os.remove(p)
            except OSError:
                pass
    return y, state, info


# ---------------------------------------------------------------------------
# Signals and estimators
# ---------------------------------------------------------------------------
def zadoff_chu(N, u=1):
    n = np.arange(N)
    return np.exp(-1j * np.pi * u * n * (n + (N % 2)) / N)


def sounder_estimates(y, ref, P):
    """Circular cross-correlation of each period of y with one period ref:
    returns H[k, l] (period k, delay l), unit gain for an ideal channel."""
    K = len(y) // P
    Y = np.fft.fft(y[:K * P].reshape(K, P), axis=1)
    R = np.fft.fft(ref[:P])
    return np.fft.ifft(Y * np.conj(R)[None, :], axis=1) / np.sum(np.abs(ref[:P]) ** 2)


def awgn(x, snr_db, ref_power=None, rng=None):
    rng = rng or np.random.default_rng(0)
    p = np.mean(np.abs(x) ** 2) if ref_power is None else ref_power
    s = np.sqrt(p / 10 ** (snr_db / 10) / 2)
    return x + s * (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x)))


def spreading_function(H, fsnap, win=None):
    """Delay-Doppler spreading function |S(tau, nu)|^2 from snapshots H[k, l]."""
    K = H.shape[0]
    w = np.hanning(K) if win is None else win
    S = np.fft.fftshift(np.fft.fft(H * w[:, None], axis=0), axes=0)
    nu = np.fft.fftshift(np.fft.fftfreq(K, 1 / fsnap))
    return nu, np.abs(S) ** 2 / np.sum(w ** 2)


def doppler_moments(nu, P):
    p = P / max(P.sum(), 1e-30)
    m = np.sum(nu * p)
    return m, np.sqrt(max(np.sum((nu - m) ** 2 * p), 0))


def corr_abs(h, lag):
    a, b = h[lag:], h[:-lag]
    return np.abs(np.sum(a * np.conj(b))) / np.sqrt(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2))


def ks_distance(a, b):
    a = np.sort(np.asarray(a)[np.isfinite(a)]); b = np.sort(np.asarray(b)[np.isfinite(b)])
    if len(a) == 0 or len(b) == 0:
        return np.nan
    z = np.concatenate([a, b])
    return float(np.max(np.abs(np.searchsorted(a, z, 'right') / len(a) - np.searchsorted(b, z, 'right') / len(b))))


def k_moment(p):
    """Greenstein moment-method K (linear) from powers p."""
    m = np.mean(p); v = np.var(p)
    g = v / m ** 2
    if g >= 1:
        return 0.0
    s = np.sqrt(1 - g)
    return s / (1 - s)


def lcr_afd(env, fs, levels_db):
    """Level crossing rate (1/s, upward) and AFD (s) of an envelope, levels
    relative to the rms value."""
    r = env / np.sqrt(np.mean(env ** 2))
    T = len(r) / fs
    out = []
    for L in levels_db:
        th = 10 ** (L / 20)
        below = r < th
        up = np.sum(below[:-1] & ~below[1:])
        lcr = up / T
        out.append((lcr, (np.mean(below) / lcr) if lcr > 0 else np.nan))
    return np.array(out)


def save_json(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
