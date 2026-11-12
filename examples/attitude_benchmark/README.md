# Recorded-data attitude benchmark

This example estimates attitude from recorded gyroscope, accelerometer and
magnetometer data, and scores the result against optical ground truth. One
estimator model is run on both filter cores, `UKF::Core` (`ukf`) and
`UKF::SquareRootCore` (`srukf`), and compared with Madgwick, Mahony and online
VQF on the BROAD and Sassari datasets. Both six-axis (gyroscope +
accelerometer) and nine-axis (+ magnetometer) configurations are run. See
[the comparison and methodology](../../docs/attitude-benchmark.md).

## Where to look

| File | Contents |
|---|---|
| [`attitude_estimator.h`](attitude_estimator.h) | The estimator: state and measurement vectors, process and measurement models, tuning, and the update step. Start here. |
| [`sensor_checks.h`](sensor_checks.h) | Rest detection and magnetic disturbance detection. They are not part of the UKF; they decide which measurements it is given. |
| [`imu_csv.h`](imu_csv.h) | Reading sensor CSV files and writing attitude CSV files. |
| [`estimator.cpp`](estimator.cpp) | The command-line program. |
| `run.py`, `datasets.py`, `metrics.py`, `references.py` | The benchmark harness: pinned downloads, scoring and the reference algorithms. |

## How the estimator works

The filter is a small strapdown inertial navigator with three state vector
fields:

| State field | Type | Process model |
|---|---|---|
| `Attitude` | Quaternion, body to world | Turned by the gyroscope reading, less the bias |
| `Velocity` | Vector, world frame | Changed by the accelerometer reading, rotated into the world frame, plus gravity |
| `GyroBias` | Vector, body frame | Constant; process noise lets it wander |

Nothing measures the velocity, but a sensor moved by hand never gets far. So
the filter is given a pseudo-measurement saying that the velocity is zero, with
a generous uncertainty. That is enough to find which way is down: if the
attitude estimate is tilted, some of gravity is mistaken for acceleration, the
velocity estimate runs away, and the measurement of zero pulls both the
velocity and the tilt back. Real accelerations do little harm, because they
average out within a few seconds while the effect of a tilt error keeps
growing.

| Measurement field | Set when | Measurement model |
|---|---|---|
| `ZeroVelocity` | Every step | The `Velocity` state |
| `RestGyroscope` | The sensor is at rest | The `GyroBias` state |

The measurement vector is a `UKF::DynamicMeasurementVector`, because the
gyroscope only measures its own bias while the sensor is at rest. The bias is
also learned during motion, through its effect on the velocity.

None of that can observe heading. With a magnetometer, a second filter with a
single state (`HeadingOffset`) estimates the rotation about the vertical which
points the attitude at magnetic north. Keeping it separate means that a
disturbed magnetic field can spoil the heading but never the tilt. Readings
whose field strength or dip angle have strayed from the Earth's are ignored.

Because heading is unobservable, the attitude filter's uncertainty about it
would grow without limit, and the unscented transform cannot represent
attitude errors approaching half a turn. The estimator therefore caps that
one uncertainty (`limit_heading_sigma`).

`ukf` and `srukf` share all of this; the few lines where the two cores differ
(covariances versus root covariances) are collected in six small helper
functions. The tuning is in `Attitude_Tuning`, and can be printed with
`attitude_benchmark --describe-filter srukf`.

## Build and replay

From the repository root, using Python 3.12 and a C++17 compiler:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r examples/attitude_benchmark/requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target attitude_benchmark --parallel 2
python -m unittest discover -s examples/attitude_benchmark -p 'test_*.py' -v
python examples/attitude_benchmark/run.py --smoke
python examples/attitude_benchmark/run.py
```

The first run downloads datasets into `.benchmark-cache` (~215 MB for the full
suite). Every file is checked against its SHA-256 and size in
`datasets.lock.json`, including on cache hits. Do not commit dataset binaries.
Network access is needed for initial dependency and data downloads; a verified
populated dataset cache supports offline replay.

Full runs use all 39 BROAD trials and all 18 Sassari speed/sensor combinations,
both modes and all five algorithms: 570 evaluations. Smoke runs use the
complete BROAD 01 and Sassari slow/XS1 recordings (20 evaluations); they are
integration checks, not representative accuracy rankings. Select algorithms with
`--algorithms ukf srukf vqf`, datasets with `--datasets broad`, or an axis
configuration with `--axes 6`. Use `--cache`, `--output`, and `--executable` to
change paths.

The estimator is tuned on a fixed development subset: BROAD 01, 08, 10, 24, 28
and Sassari slow (all six sensors), 11 recordings. Validation uses the other 34
BROAD trials plus Sassari medium/fast (12 recordings). All sensors sharing a
physical Sassari motion stay in the same split.

```sh
python examples/attitude_benchmark/run.py --split development --output benchmark-results/development
python examples/attitude_benchmark/run.py --split validation --output benchmark-results/validation
python examples/attitude_benchmark/run.py --datasets broad --trials '10_*' \
  --algorithms srukf --axes 6 --diagnostics --output benchmark-results/translation-diagnostic
```

Quoted `--trials` patterns match either the recording name or its qualified name,
such as `'sassari/slow_*'`. Explicit selections are marked custom in the report;
unmatched patterns fail. `--smoke` cannot be combined with a split or selectors.

The runner writes `benchmark-results/comparison.md`, `results.json`, and
`per_trial.csv`. JSON contains source and executable hashes, package versions,
the executable's tuning, exact recording selections, per-trial RMSE/P95, sample
counts, alignment windows, and diagnostic timing. It also includes summaries by
split and condition, and accuracy ratios of the UKF estimators against VQF when
VQF was selected. Failed filters, malformed datasets, and checksum failures stop
the run instead of disappearing from the average. Neither sources nor the
executable may change during a run.

`--diagnostics` writes one trace per UKF run under `OUTPUT/diagnostics/`: the
gyroscope bias (rad/s) and velocity (m/s) estimates, whether the sensor was
judged to be at rest, the heading offset (rad), and whether the magnetometer
reading was rejected. These writes occur inside the timed loop and are marked
`cpp_update_loop_with_diagnostics_io` in per-trial results. Run without
diagnostics for update-loop timing. C++ and Python reference timings include
different overheads and are not a cross-language speed comparison.

## Using the executable directly

```sh
build/examples/attitude_benchmark/attitude_benchmark \
  --input sensors.csv --output attitude.csv --filter srukf --axes 9 \
  --diagnostics states.csv
```

Input header: `time,gx,gy,gz,ax,ay,az,mx,my,mz`. Time is seconds; gyroscope is
rad/s; acceleration is specific force in m/s², including gravity; magnetic
values may use any consistent scale. Supply all ten columns even in six-axis
mode (zeros are acceptable for the unused magnetometer). Samples must have
finite values and strictly increasing timestamps. The first accelerometer
reading must not be zero, because it sets the initial tilt. A magnetometer
reading with no horizontal component is skipped.

Output header: `time,qw,qx,qy,qz`. Quaternions rotate body to a right-handed
z-up world; in nine-axis mode horizontal magnetic north is along +X, and in
six-axis mode the heading is arbitrary. The first output comes from the first
sample alone; each subsequent output incorporates the sample at that timestamp.
Timing JSON is written to stderr; input and attitude-output CSV I/O is outside
the timed update loop, while optional diagnostic writes are inside it.

Data attribution and upstream licenses are recorded in `datasets.lock.json`.
BROAD is CC BY 4.0; the Sassari repository provides an MIT license. Reference
implementations remain installed dependencies rather than vendored source:
AHRS 0.4.0 (MIT) and VQF 2.0.1 (MIT).
