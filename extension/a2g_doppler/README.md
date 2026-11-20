# a2g_doppler: per-tap A2G Doppler extension for CHEM

This extension adds a time-varying, per-tap Doppler channel model to CHEM.
It implements the model of the report *AERPAW A2G Doppler Model for CHEM*
(v1, Nov 2026). Each multipath tap has its own geometric Doppler. Each tap
can also have a narrow diffuse Doppler spectrum.

Today CHEM applies one Doppler shift per link. It applies that shift after
the static CIR, so all taps share it. CHEM updates the shift every 100 ms.
The extension replaces this for the links that it owns.

Status: tested out of tree against the public CHEM snapshot (commit
`6b9e58f`, 27 Jul 2026). The deployed twin runs CHEM v1.2.4, which is newer.
The hook anchors in `channel.h`/`channel.cpp` can differ in v1.2.4. Check
them before you merge (see *Core changes*).

## Model

For a link from node A to node B, with path p = 0 .. P-1:

```
y[n] = sum_p  H_p[k] * exp(j psi_p[n]) * x[n - D_p[k]]
psi_p(t) = -2 pi kappa d_p(t) / lambda      (continuous phase, double precision)
H_p[k]   = a_p ( sqrt(K_p/(K_p+1)) + sqrt(1/(K_p+1)) g_p(t_k) )
```

* `d_p(t)` is the path length. It is computed from the node trajectories:
  * LOS: `|r_A - r_B|`.
  * GROUND: specular reflection on the plane `z = groundZ`, by the image
    method.
  * SCATTERER: a point at a fixed ENU position: `|r_A - s| + |s - r_B|`.
  * FIXED: a TDL tap with a constant excess delay and a constant Doppler
    (`dopplerHz`, or `f_m cos(aoaDeg)`).
* The tap frequency in packet k uses the *midpoint* (range-difference) rule:
  `f_k = -(d(t_k+Delta) - d(t_k)) / (lambda Delta)`. With this rule the
  phase is exact at every packet boundary, and the error does not
  accumulate. The `START_HOLD` rule (what CHEM does now) is also available.
* The packet interval is `Delta = updateMs` (default 0.5 ms). Inside the
  packet, the phase advances per sample.
* `g_p` is an independent, unit-power diffuse process for each path. It is a
  sum of 32 sinusoids with stratified quantile frequencies, so the rms width
  is exact. Available shapes:
  * `GAUSS` (default),
  * `EXP` (two-sided exponential),
  * `JAKES` (for comparison).

  Each sinusoid has its own phase accumulator, so the width can change during
  a flight without a jump.
* The rms width `sigma_f` comes from the altitude table of the letter
  (5.0, 4.5, 4.0, 3.1, 3.2 Hz at 30 .. 110 m, measured at 3.51 GHz and
  5 m/s). It is scaled by `fc/3.51 GHz`. It is also scaled by `speed/5 m/s`
  when `widthScaling` is `carrier_speed`.
* The LOS default is `K = 8.5 dB`. With `"taps": "ALL"`, the other paths get
  a diffuse part with the same K. With `"LOS_ONLY"`, the other paths are
  specular.
* Path amplitudes are relative to the LOS:
  * GROUND: `Gamma d_0/d_g`.
  * SCATTERER: `10^(gainDb/20) d_0/d_s`.
  * FIXED: `10^(gainDb/20)`.

  The LOS amplitude is 1, so CHEM's path-loss model stays the LOS loss.
  With `normalizeTotal` the sum of the path powers is 1.
* Delays are integer samples (as the CHEM CIR uses). The option
  `fractionalDelay` uses a 16-tap windowed sinc and adds a bulk latency of
  7 samples.
* Optional oscillator impairment (report eq. 12):
  * a drift of 125 Hz per 30 min,
  * plus a random walk,
  * plus 2 Hz rms jitter held for 0.1 s.

## Files

```
extensions/a2g_doppler/
  CMakeLists.txt            own targets and options (added by one marked line in the top CMakeLists)
  include/a2g_doppler.h     chem::A2GDopplerExtension (ChannelExtension)
  src/a2g_doppler.cpp       config parser, live poll loop, hook / taps modes, EXT commands
  model/include/a2g/model.h model library (pure C++17, no CHEM dependency)
  model/src/model.cpp
  tests/test_model.cpp      16 model tests (GoogleTest)
  tests/test_hook.cpp       5 tests through chem::Channel, Intermediate, node map, ExtensionRegistry context
  tools/a2g_harness.cpp     offline harness: IQ file -> CHEM Channel or Intermediate -> IQ file
  tools/apply_core_hook.py  applies the marked core edits (idempotent)
  config/a2g_doppler.example.json
patch/
  chem_a2g_doppler_full.patch   git diff vs the snapshot: core edits + extension
  chem_core_only.patch          only the marked core edits
```

## Core changes (all between `>>> CHEM-EXT` and `<<< CHEM-EXT` markers)

| File | Change | Needed for |
|---|---|---|
| `CMakeLists.txt` | `option(CHEM_EXT_A2G_DOPPLER ON)` + `add_subdirectory(extensions/a2g_doppler)` (6 lines) | build |
| `src/chem/channel/coordinator.cpp` | `#include "a2g_doppler.h"` + one `registerExtension(...)` under `#ifdef` (CHEM's documented extension registration) | registration |
| `configs/config.json` | default `"a2g_doppler"` block, `enabled: false` | config |
| `include/chem/channel/channel.h`, `src/chem/channel/channel.cpp` | optional time-varying CIR hook: `setTimeVaryingCIR(fn)`, `clearTimeVaryingCIR()`, `hasTimeVaryingCIR()`. In `processChannel()`, the hook replaces the static CIR stage when it is set (about 25 lines) | `mode: "hook"` only |

The extension also builds and runs **without** the channel hook. Build with
`apply_core_hook.py --no-hook`; the extension then uses `mode: "taps"`. In
that mode it pushes the CIR with `Intermediate::updateCIR()` every
`tapUpdateMs`, and CHEM holds the taps for each buffer (zero-order hold).
The hook is still recommended. In taps mode, each tap is held constant for
one CHEM buffer and one update interval (zero-order hold). This adds
spectral images at multiples of the update rate. In the report (E6), a
-150 Hz tap with 1 ms buffers has images at -11.6 dB. With updates of 5 ms
or more, the tap Doppler aliases. The taps mode is accurate only if
|f_D| x (update interval) << 1, so the buffers must also be short.

## Build and test

```bash
python3 extensions/a2g_doppler/tools/apply_core_hook.py .     # or --no-hook
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCHEM_ENABLE_UNIT_TESTING=ON   # CHEM_EXT_A2G_DOPPLER=ON by default
make -j$(nproc)
ctest --output-on-failure          # a2g_model, a2g_hook + CHEM's own 3 suites
```

Turn the extension off with `-DCHEM_EXT_A2G_DOPPLER=OFF`. CHEM is then
unchanged in behaviour. The hook members still exist but they are never set.

## Configure

In `~/.config/CHEM/config.json` (see `config/a2g_doppler.example.json`):

```json
"extensions": {
  "a2g_doppler": {
    "enabled": true, "mode": "hook", "updateMs": 0.5,
    "phaseUpdate": "MIDPOINT", "dopplerMode": "PER_TAP",
    "K_dB": 8.5, "shape": "GAUSS", "taps": "ALL", "widthScaling": "carrier_speed",
    "paths": [ {"type": "LOS"}, {"type": "GROUND", "reflCoeff": -0.3},
               {"type": "SCATTERER", "point": [260, -140, 15], "gainDb": -6} ],
    "referenceOrigin": {"lat": 35.72750947, "lon": -78.69595819, "alt": 0.0},
    "links": [ {"src": "<node id>", "dest": "<node id>"} ],
    "nodes": { "LW1": {"lat": 35.727489, "lon": -78.695957, "alt": 10.0} },
    "trajectoryDelayMs": 200, "seed": 1234,
    "oscImpairment": {"enable": false}
  }
}
```

* An empty `links` list means all links.
* `nodes` overrides the position of fixed nodes. Other nodes use CHEM's node
  map: the vehicle position and velocity.
* Positions are converted to ENU with the same flat-earth formula as
  `Coordinator::configChannels()`.
* The trajectory is evaluated `trajectoryDelayMs` late. The model then
  interpolates between position reports (cubic Hermite) and does not
  extrapolate.

At run time (pychem / coordinator `EXT` command):

```json
{"CMD": "EXT", "extension": "a2g_doppler", "action": "configure", "K_dB": 6.0, "shape": "GAUSS"}
{"CMD": "EXT", "extension": "a2g_doppler", "action": "status"}
{"CMD": "EXT", "extension": "a2g_doppler", "action": "stop"}
```

The `status` action returns the following per link: the Doppler, the
excess delay, the amplitude and `sigma_f` of each path.

## Things to know

* CHEM compiles with `-ffast-math`. The model library builds with
  `-fno-fast-math -fcx-limited-range`, and it uses finite sentinels
  (K >= 300 dB means specular, K <= -300 dB means Rayleigh). Do not use inf
  or NaN in the config. Without `-fcx-limited-range`, GCC calls `__muldc3`
  for every complex multiply, and the hook is several times slower.
* Cost: see the report (E5) for the throughput per path at 23.04 MS/s.
* CHEM finding: a node without RF characteristics gives NaN samples from
  `processChannel()`. Its `source_power_dbfs` is NaN by default, and
  `-ffast-math` removes the `std::isfinite()` test. Twin nodes have
  characteristics. The harness sets them.
* The extension turns off the CHEM link Doppler (`setDopplerEnabled(false)`)
  on the links that it owns, so the shift is not applied twice. CHEM path
  loss, antenna gains and noise stay as they are.
* Thread safety: the hook state has its own mutex. Hook installation uses
  `std::atomic_store` on a `shared_ptr`.
* Time base: in live mode the stream clock re-syncs to `steady_clock` when
  the two drift apart by more than 50 ms. The midpoint phase is anchored to
  the geometry, so a re-sync does not break the phase model.

## Offline harness

`a2g_harness run.json` pushes an IQ file (cf32) through
`chem::Channel::processChannel()` (`"pipeline": "channel"`). With
`"pipeline": "intermediate"` it uses the full `chem::Intermediate` worker
path: a real `Receiver`, UHD sc16 converters, and a capturing `Transmitter`
in place of UDP. Every 100 ms of stream time, the harness does what
`Coordinator::configChannels()` does from ENU trajectories. The Python
drivers are in `../harness/` (`chemlib.run_chem()`).
