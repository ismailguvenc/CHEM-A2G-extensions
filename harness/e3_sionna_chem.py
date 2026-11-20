# SPDX-License-Identifier: AGPL-3.0-or-later
"""E3: Sionna RT with the local CHEM build (CHEM's own 'sionna' extension client + AERPAW Sionna server).

The UAV positions are those of the offline ray trace (sionna_rt/out/rt_paths_lw1_flight.json, LW1 ->
LW2 bearing, 30 m, 5 m/s). chem_sionna_probe moves the CHEM sandbox node along them; the sionna
extension posts the positions, the server ray-traces the aerpaw Lake Wheeler scene, and CHEM applies
the CIR. We compare the CIR energy that CHEM applies with the offline path power.
"""
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from chemlib import *  # noqa

OUT = os.path.join(RESULTS, 'e3_sionna_chem')
os.makedirs(OUT, exist_ok=True)
RT = os.path.join(os.path.dirname(__file__), '..', '..', 'sionna_rt', 'out', 'rt_paths_lw1_flight.json')
PROBE = os.path.expanduser('~/work/build_v2/bin/Release/chem_sionna_probe')
d = json.load(open(RT))
lat0, lon0 = d['origin']['lat'], d['origin']['lon']


def enu_to_ll(e, n):
    lat = lat0 + n / 111319.5
    lon = lon0 + e / (111319.5 * np.cos(np.radians(0.5 * (lat + lat0))))
    return lat, lon


rows = []
for s in d['samples'][::4]:
    e, n, u = s['rx_enu']
    lat, lon = enu_to_ll(e, n)
    g = np.array([p['gain_re'] + 1j * p['gain_im'] for p in s['paths']])
    los = [p for p in s['paths'] if p['type'] == 'LOS']
    rows.append((s['t'], lat, lon, u, 10 * np.log10(np.sum(np.abs(g) ** 2)),
                 10 * np.log10(abs(los[0]['gain_re'] + 1j * los[0]['gain_im']) ** 2) if los else np.nan,
                 10 * np.log10(abs(np.sum(g)) ** 2)))
pos = os.path.join(RUNDIR, 'e3_pos.csv')
with open(pos, 'w') as f:
    f.write('t,lat,lon,alt\n')
    for r in rows:
        f.write(f'{r[0]},{r[1]:.8f},{r[2]:.8f},{r[3]:.2f}\n')
txe = d['tx']['enu']
txlat, txlon = enu_to_ll(txe[0], txe[1])
cfg = dict(serverUrl='http://127.0.0.1:8000', sceneConfig='aerpaw', sceneOffset=dict(x=-3.7, y=3.4, z=0.0),
           referenceOrigin=dict(lat=lat0, lon=lon0, alt=112.0), tx=dict(lat=txlat, lon=txlon, alt=txe[2]),
           positions=pos, dwellMs=1500, updateRateMs=300, out=os.path.join(OUT, 'probe.csv'), fc=d['fc_hz'])
jp = os.path.join(RUNDIR, 'e3.json')
json.dump(cfg, open(jp, 'w'), indent=1)
env = dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser('~/deps/mavsdk_x/usr/lib') + ':' +
           os.path.expanduser('~/mm/root/envs/chembuild/lib'))
r = subprocess.run([PROBE, jp], capture_output=True, text=True, env=env)
open(os.path.join(OUT, 'probe_stderr.txt'), 'w').write(r.stderr)
if r.returncode:
    raise RuntimeError(r.stderr[-2000:])
pr = np.genfromtxt(cfg['out'], delimiter=',', names=True)
off = np.array([[x[0], x[4], x[5], x[6]] for x in rows])
m = min(len(pr), len(off))
diff_tot = pr['energy_db'][:m] - off[:m, 1]
res = dict(n=m, active_frac=float(np.mean(pr['active'][:m] > 0)),
           taps_med=float(np.median(pr['taps'][:m])),
           energy_vs_offline_total_db=dict(median=float(np.median(diff_tot)), p10=float(np.percentile(diff_tot, 10)),
                                           p90=float(np.percentile(diff_tot, 90))),
           energy_vs_offline_coherent_db=float(np.median(pr['energy_db'][:m] - off[:m, 3])),
           status=json.loads(r.stdout.strip().splitlines()[-1]))
np.savez(os.path.join(OUT, 'e3_compare.npz'), t=off[:m, 0], chem_db=pr['energy_db'][:m], rt_tot_db=off[:m, 1],
         rt_los_db=off[:m, 2], rt_coh_db=off[:m, 3], taps=pr['taps'][:m])
save_json(res, os.path.join(OUT, 'e3_summary.json'))
print(json.dumps(res, indent=1, default=str)[:1500])
