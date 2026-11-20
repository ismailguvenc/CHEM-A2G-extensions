# CHEM A2G Doppler and SCM MIMO extensions

This repository has two extensions for CHEM, the channel emulator of the AERPAW digital twin
([anilgurses/CHEM](https://github.com/anilgurses/CHEM)):

| Extension | Function |
|---|---|
| `extension/scm_mimo/` | 3GPP TR 25.996 spatial channel model (SCM) for MIMO links. Scenarios `SMA`, `UMA`, `UMI`, `UMI_LOS` and `A2G` (air-to-ground). Uniform linear arrays and cross-polarized pairs at the ground and air ends. Doppler options `SCM`, `PROPOSED`, `JAKES` and `TDL_D`. |
| `extension/a2g_doppler/` | Per-tap A2G Doppler for SISO links. Each multipath tap (LOS, ground reflection, point scatterers, fixed TDL taps) gets its own geometric Doppler from the node trajectories, plus an optional narrow diffuse Doppler spectrum. Doppler options `PROPOSED`, `JAKES`, `TDL_A` and `TDL_D` (3GPP TR 38.901). |

`scm_mimo` uses the geometry, trajectory and Doppler library of `a2g_doppler`. Therefore, if you enable
`scm_mimo`, you must also enable `a2g_doppler`.

Status: research code (November 2026). We tested it out of tree against the public CHEM snapshot
(commit `6b9e58f`), with unit tests, offline harnesses, a local V-USRP run, and a playback of CIR timelines into
the deployed CHEM of the AERPAW twin. The AERPAW deployment runs a newer CHEM version. Check the core
hook anchors before you merge into that version.

## Model summary

**a2g_doppler.** For path p, the tap phase is `psi_p(t) = -2 pi d_p(t) / lambda`, where `d_p(t)` is the path length
from the node trajectories (LOS, image-method ground reflection, or scatterer). The phase is continuous and in
double precision. The frequency in each packet (default 0.5 ms) uses the midpoint rule, so the phase error does
not accumulate. A diffuse part (sum of 32 sinusoids, Gaussian, exponential or Jakes shape) gives each tap a
Rician fading with K = 8.5 dB by default. The rms Doppler width comes from an altitude table measured at
3.51 GHz and 5 m/s (5.0 Hz at 30 m to 3.2 Hz at 110 m). The extension scales this width with the carrier
frequency and the speed. See `extension/a2g_doppler/README.md`.

**scm_mimo.** Each SCM drop (TR 25.996 steps 3 to 10) gives the path delays, powers and angle offsets. The
extension applies the angles relative to the LOS direction at each end, so the channel follows the UAV
trajectory. A new drop starts after `dropDistanceM` of travel, with a crossfade of `crossfadeS`. The sub-path
phase is the Doppler process and you select it with `doppler`:

- `SCM`: the TR 25.996 sub-path Doppler.
- `PROPOSED`: a LOS line at the geometric Doppler, plus the diffuse A2G spectrum of `a2g_doppler`.
- `JAKES`: the classical Jakes spectrum.
- `TDL_D`: the 3GPP TDL-D Doppler.

## Contents

| Folder | Content |
|---|---|
| `extension/a2g_doppler/` | Extension source, model library (pure C++17), 21 GoogleTest tests, offline harness `a2g_harness`, example configuration, and `tools/apply_core_hook.py` (adds the per-sample CIR hook to the CHEM core). |
| `extension/scm_mimo/` | Extension source, SCM model library, 10 GoogleTest tests, `scm_harness` (N IQ streams through CHEM), `chem_taps_timeline` (writes a CIR timeline for playback), `chem_sionna_probe`, example configuration, and `tools/apply_mimo_core.py` (adds the MIMO port summation to the CHEM core). |
| `extension/add_traj_files.py` | Adds the `trajectoryFiles` option to both extensions (already applied in this copy). |
| `build_scripts/` | `build_v2.sh` copies a clean CHEM clone, applies the core changes, adds both extensions, builds and runs the unit tests. `py.sh` runs a harness script. |
| `harness/` | Python experiments: E1 (Doppler options, SISO), E2 (SCM 2x2), E3 (Sionna RT with CHEM), timeline generation, analysis and figures. |
| `vusrp_local/` | E4: public V-USRP build and a sample over-the-emulator experiment with the local CHEM. No AERPAW access is necessary. |
| `results_summaries/` | JSON summaries of the experiments (the numbers in the report tables). |

All core changes are between `>>> CHEM-EXT` and `<<< CHEM-EXT` markers. The apply scripts are idempotent.

## Build (Linux or WSL2 Ubuntu)

Prerequisites:

- A clone of `github.com/anilgurses/CHEM` (tested at commit `6b9e58f`) in `~/work/CHEM_pristine`.
- micromamba at `~/mm` with the environments `chembuild` (CHEM build), `chem` (Python analysis) and `vusrp`
  (V-USRP build).
- MAVSDK at `~/deps/mavsdk_x/usr` (CHEM dependency).
- GCC 13 needs `-include cstdint` (the build script sets it).
- For E4: `github.com/anilgurses/VUSRP` (tested at commit `f87bce7`) and Boost 1.84 (not 1.87 or later).

Change the paths in the scripts if your setup is different. Then:

1. `bash build_scripts/build_v2.sh --fresh` builds CHEM with both extensions and runs the unit tests.
2. `bash build_scripts/py.sh e1_doppler_models.py` runs the Doppler experiment.
3. `bash build_scripts/py.sh e2_scm_mimo_local.py` runs the SCM MIMO experiment.
4. `bash vusrp_local/build_vusrp.sh`, then `bash vusrp_local/run_e4.sh`, runs the V-USRP experiment.

To add the extensions to another CHEM tree by hand:

1. Run `python3 extension/a2g_doppler/tools/apply_core_hook.py <CHEM dir>`.
2. Run `python3 extension/scm_mimo/tools/apply_mimo_core.py <CHEM dir>`.
3. Copy `extension/a2g_doppler` and `extension/scm_mimo` to `<CHEM dir>/extensions/`.
4. Configure with `-DCHEM_EXT_A2G_DOPPLER=ON -DCHEM_EXT_SCM_MIMO=ON`.

## Configuration

Each extension reads its block under `"extensions"` in the CHEM configuration. See
`extension/a2g_doppler/config/a2g_doppler.example.json` and `extension/scm_mimo/config/scm_mimo.example.json`.
An empty `links` list applies the model to all links. Fixed nodes can have fixed positions under `nodes`. Other
nodes use the positions and velocities that CHEM gets from the vehicles, or recorded trajectories through
`trajectoryFiles`.

## Authors and license

Ismail Guvenc, NC State University, with AI assistance (Claude). Developed for the
[AERPAW](https://aerpaw.org) digital twin.

AGPL-3.0-or-later, the same as CHEM. See `LICENSE`.
