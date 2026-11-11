# Recorded-data attitude benchmark

This example compares the original covariance/square-root UKF model and three
robust-estimator stages against Madgwick, Mahony, and online VQF on BROAD and
Sassari. Both six-axis
(gyro + accelerometer) and nine-axis (+ magnetometer) configurations are run.
See [the comparison and methodology](../../docs/attitude-benchmark.md).

| Filter | Model and enabled features |
|---|---|
| `ukf`, `srukf` | Original joint quaternion/gyro-bias model with fixed-noise raw direction updates |
| `robust-rest` | Quaternion SR-UKF, separate rest-only bias SR-UKF, raw acceleration and isolated magnetic heading |
| `robust-acc` | Same architecture plus acceleration filtering in the gyro-propagated frame |
| `robust` | Same filtered architecture plus acceleration norm/innovation gating and magnetic norm/dip rejection with reacquisition |

All robust stages use two-component tangent-plane gravity observations, constrain
the gravity-update gain to the world-horizontal tilt subspace so it cannot claim
absolute heading information, and transport attitude covariance after a correction.
The covariance update accounts for that constrained gain. The first stage changes the
estimator architecture as well as adding rest-based bias learning; its difference
from the baseline is not a one-feature ablation. The robust stages estimate bias
at rest, while VQF also has motion-based bias estimation. None uses optical ground
truth in its sensor updates. Runtime parameters and feature flags can be inspected
with `attitude_benchmark --describe-filter robust` (also accepts the other stages).

From the repository root, using Python 3.12 and a C++17 compiler:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r examples/attitude_benchmark/requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target attitude_benchmark --parallel 2
build/examples/attitude_benchmark/robust_geometry_test
python -m unittest discover -s examples/attitude_benchmark -p 'test_*.py' -v
python examples/attitude_benchmark/run.py --smoke
python examples/attitude_benchmark/run.py
```

The first run downloads datasets into `.benchmark-cache` (~215 MB for the full
suite). Every file is checked against its SHA-256 and size in
`datasets.lock.json`, including on cache hits. Do not commit dataset binaries.
Network access is needed for initial dependency and data downloads; a verified
populated dataset cache supports offline replay.

The geometry executable is built as a dependency of `attitude_benchmark`. It
checks tangent observations, quaternion sign handling, and covariance transport;
it is run explicitly so the repository's existing unittest-only CTest workflow
does not require building this example.

Full runs use all 39 BROAD trials and all 18 Sassari speed/sensor combinations,
both modes and all eight algorithms/variants: 912 evaluations. Smoke runs use the
complete BROAD 01 and Sassari slow/XS1 recordings (32 evaluations); they are
integration checks, not representative accuracy rankings. Select algorithms with
`--algorithms ukf srukf robust-rest robust-acc robust vqf`, datasets with
`--datasets broad`, or an axis configuration with `--axes 6`. Use
`--cache`, `--output`, and `--executable` to change paths.

New-estimator development uses BROAD 01, 08, 10, 24, 28 and Sassari slow (all six
sensors): 11 recordings. Validation uses the other 34 BROAD trials plus Sassari
medium/fast (12 recordings). All sensors sharing a physical Sassari motion stay
in the same split. These are development/validation labels for the new
implementation; earlier baseline results were already observed on every trial.

```sh
python examples/attitude_benchmark/run.py --split development --output benchmark-results/development
python examples/attitude_benchmark/run.py --split validation --output benchmark-results/validation
python examples/attitude_benchmark/run.py --datasets broad --trials '10_*' \
  --algorithms robust --axes 6 --diagnostics --output benchmark-results/translation-diagnostic
```

Quoted `--trials` patterns match either the recording name or its qualified name,
such as `'sassari/slow_*'`. Explicit selections are marked custom in the report;
unmatched patterns fail. `--smoke` cannot be combined with a split or selectors.

The runner writes `benchmark-results/comparison.md`, `results.json`, and
`per_trial.csv`. JSON contains source and executable hashes, package versions,
the executable's robust-stage parameter descriptions, exact recording selections,
per-trial RMSE/P95, sample counts, alignment windows, and diagnostic timing.
It also includes summaries by split/condition and robust-stage accuracy ratios
against the original UKF and VQF when those references were selected.
Failed filters, malformed datasets, and checksum failures stop the run instead
of disappearing from the average. Neither sources nor the executable may change
during a run.

`--diagnostics` writes robust-stage traces under `OUTPUT/diagnostics/`, including
gyro bias (rad/s), rest/rejection decisions, innovations (rad), covariance traces,
acceleration magnitude and magnetic norm/dip. These writes occur inside the timed
loop and are marked `cpp_update_loop_with_diagnostics_io` in per-trial results.
Run without diagnostics for update-loop timing. C++ and Python reference timings
include different overheads and are not a cross-language speed comparison.

The C++ executable is independently usable:

```sh
build/examples/attitude_benchmark/attitude_benchmark \
  --input sensors.csv --output attitude.csv --filter robust --axes 9 \
  --diagnostics states.csv
```

Input header: `time,gx,gy,gz,ax,ay,az,mx,my,mz`. Time is seconds; gyro is rad/s;
acceleration is specific force in m/s², including gravity; magnetic values may
use any consistent scale. Supply all ten columns even in six-axis mode (zeros
are acceptable for the unused magnetometer). Samples must have finite values
and strictly increasing timestamps. Baselines and intermediate robust stages
require nonzero used direction vectors. The final `robust` stage can reject zero
acceleration and propagate through that invalid gravity observation. Nine-axis
mode still requires a nonzero magnetometer vector for every sample.

Output header: `time,qw,qx,qy,qz`. Quaternions rotate body to a right-handed
z-up world; nine-axis initialization puts horizontal magnetic north on +X.
Six-axis world yaw is arbitrary. Initialization uses current sensor data;
the final robust stage defers gravity initialization if startup acceleration is
implausible and initializes when a usable observation arrives. Ordinarily the
first output is sensor-only initialization;
each subsequent output incorporates the sample at that timestamp. Timing JSON
is written to stderr; input and attitude-output CSV I/O is outside the timed
update loop, while optional diagnostic writes are inside it.

Data attribution and upstream licenses are recorded in `datasets.lock.json`.
BROAD is CC BY 4.0; the Sassari repository provides an MIT license. Reference
implementations remain installed dependencies rather than vendored source:
AHRS 0.4.0 (MIT) and VQF 2.0.1 (MIT).
