# EuRoC MAV attitude benchmark

This example adds a model of how a multirotor flies to the attitude estimator
of the [attitude benchmark](../attitude_benchmark), and scores it on the EuRoC
recordings of an AscTec Firefly hexrotor. The results are compared with the
unmodified `ukf` and `srukf` estimators of that benchmark, Madgwick, Mahony and
online VQF. Two controls separate the effect of the vehicle model from the
effect of merely allowing for this IMU's large gyroscope bias:

- `mav-srukf-gravity` is the same estimator with its rotor model switched off.
- `vqf-mav` is VQF with its gyro-bias clipping limit raised to 10°/s and its
  initial bias standard deviation raised to 5°/s; other VQF settings remain at
  their online defaults.

EuRoC publishes gyroscope and accelerometer samples at nominally 200 Hz, but no
magnetometer or rotor-command stream. This benchmark therefore runs **six-axis
estimation only**. It does not synthesize magnetic readings, use camera images,
or feed reference position, velocity, orientation or bias into any estimator.
A physically appropriate process model does not guarantee an accuracy advantage;
the reports retain all results without a must-win test.

## Where to look

| File | Contents |
|---|---|
| [`mav_estimator.h`](mav_estimator.h) | The estimator: state and measurement vectors, process and measurement models, tuning, and the update step. Start here. |
| [`rotor_vibration.h`](rotor_vibration.h) | Detecting from accelerometer vibration whether the rotors are turning. Not part of the UKF. |
| [`estimator.cpp`](estimator.cpp) | The command-line program. |
| [`calibrate_mount.py`](calibrate_mount.py) | Estimates the direction of rotor thrust in the IMU frame from the development recordings. |
| `run.py`, `datasets.py` | The benchmark harness. Metrics, reference algorithms, CSV reading and the rest detector are shared with the attitude benchmark. |

## Vehicle model

The filter has the same three state vector fields as the hand-held estimator,
except that the velocity is expressed in the IMU frame:

| State field | Type | Process model |
|---|---|---|
| `Attitude` | Quaternion, IMU to world | Turned by the gyroscope reading, less the bias |
| `Velocity` | Vector, IMU frame | Changed by the specific force plus gravity |
| `GyroBias` | Vector, IMU frame | Constant; process noise lets it wander |

A multirotor's rotors can only push along one axis, the thrust axis. The only
other force on the airframe is rotor drag, which opposes the velocity across
the rotor plane and is proportional to it. So in powered flight

```text
lateral accelerometer reading = -drag rate * lateral velocity
```

which turns the accelerometer into a velocity sensor, and through the velocity
pins down the tilt. In the process model the specific force along the thrust
axis is whatever the accelerometer measures, and the lateral specific force is
the rotor drag which follows from the velocity state. The lateral accelerometer
reading is used as a measurement instead:

| Measurement field | Set when | Measurement model |
|---|---|---|
| `RotorDrag` | Powered flight | Minus the drag rate times the lateral `Velocity` |
| `ZeroVelocity` | Otherwise | The `Velocity` state |
| `RestGyroscope` | At rest, rotors off | The `GyroBias` state |

The model follows [Leishman et al., *Quadrotors and Accelerometers*, IEEE
Control Systems Magazine, 2014](https://doi.org/10.1109/MCS.2013.2287362). It
represents the Firefly's aggregate translation rather than six individually
commanded motors. EuRoC does not publish motor commands, so the measured axial
specific force stands in for the thrust. Velocity is an internal state, not a
reference input or an independently validated navigation output.

The rotor model applies while the rotors are turning, detected from the
vibration they cause in the accelerometer signal, and the thrust is between
0.5 g and 2 g. Otherwise (while the vehicle is carried, or standing on the
ground) the filter behaves exactly like the hand-held estimator: the complete
accelerometer reading is integrated, and the velocity is given a
pseudo-measurement of zero. Vibration is an approximate operating-condition
detector, not an externally supplied flight mask. The drag rate and thrust
direction are fixed across sequences. Wind, aggressive maneuvers and unmodeled
forces can invalidate the approximation.

Velocity is expressed in the IMU frame because rotor drag acts in that frame:
the rotor model then never involves the heading, which a six-axis IMU cannot
observe. For the same reason the filter's uncertainty about its heading would
grow without limit, and the unscented transform cannot represent attitude
errors approaching half a turn, so that one uncertainty is capped.

The IMU is not mounted square to the airframe. The thrust direction,
`[cos(20°), 0, -sin(20°)]` in the IMU frame, is an empirical calibration from
the sensor measurements in the two development recordings
(`calibrate_mount.py`), not a published mechanical calibration or a fit to
reference orientation.

The tuning is in `MAV_Tuning`, and can be printed with
`euroc_benchmark --describe-filter mav-srukf`. The gyroscope bias starts with
a 5°/s standard deviation, and the rest detector accepts average rates up to
10°/s as a possible bias; this IMU's bias exceeds 4°/s.

## Build and replay

From the repository root, using Python 3.12 and a C++17 compiler:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r examples/attitude_benchmark/requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target attitude_benchmark euroc_benchmark --parallel 2
python -m unittest discover -s examples/attitude_benchmark -p 'test_*.py' -v
python -m unittest discover -s examples/euroc_benchmark -p 'test_*.py' -v
python examples/euroc_benchmark/run.py --smoke
python examples/euroc_benchmark/run.py
```

`--smoke` replays the complete MH_01_easy recording. It checks integration and
does not provide an independent optical-attitude ranking. A full run evaluates
all 11 recordings with eight algorithms/variants: **88 sequence evaluations**.
Inputs are verified against the pinned sizes and SHA-256 hashes in
`datasets.lock.json`, including cached files. Only sensor/reference/calibration
files are needed; camera images are not estimator inputs. The default cache is
`.benchmark-cache/euroc` and generated reports go to `benchmark-results/euroc`.

The development set is fixed to MH_01_easy and V1_02_medium. The other nine
complete sequences form the validation set; do not tune on their results.
Selectors apply before loading, so development runs need not fetch validation
recordings.

```sh
python examples/euroc_benchmark/run.py --split development \
  --output benchmark-results/euroc-development
python examples/euroc_benchmark/run.py --split validation \
  --output benchmark-results/euroc-validation
python examples/euroc_benchmark/run.py --sequences 'V1_02_*' \
  --algorithms srukf mav-srukf vqf vqf-mav --diagnostics \
  --output benchmark-results/euroc-diagnostic
```

Quoted `--sequences` patterns match complete recording names. Unmatched patterns
fail. Custom selections are marked as such and cannot be confused with the full
suite. `--smoke` cannot be combined with a split or sequence selectors. Use
`--cache`, `--output`, `--executable` or `--mav-executable` to change paths.

The native MAV executable accepts the same sensor-only CSV format as the
attitude benchmark:

```sh
build/examples/euroc_benchmark/euroc_benchmark \
  --input sensors.csv --output attitude.csv --filter mav-srukf --axes 6 \
  --diagnostics states.csv
build/examples/euroc_benchmark/euroc_benchmark --describe-filter mav-srukf
```

Input columns are `time,gx,gy,gz,ax,ay,az,mx,my,mz`: seconds, rad/s, m/s² and
three unused zero magnetic columns. Output is `time,qw,qx,qy,qz`, scalar-first
IMU-to-world orientation with world +Z up. World yaw is unobservable from the
IMU alone. The first accelerometer reading must not be zero, because it sets
the initial tilt. The optional diagnostics trace holds the gyroscope bias and
velocity estimates, whether the vehicle was judged to be at rest, whether the
rotor model applied, and the rotor vibration level; writing it is included in
that run's elapsed time.

## Reference and scoring limits

The primary attitude comparison uses the five Vicon recordings other than
V1_01_easy, whose canonical orientation has a documented quality issue. Four of
those recordings are in the validation split. Machine Hall reference orientation
is reconstructed with IMU data and has no independent optical orientation
observation. Results for all 11 sequences remain available as a descriptive
comparison, with environment-specific summaries and sequence-level warnings.
All published state references include offline IMU fusion; Vicon supplies an
additional optical orientation observation, not an entirely IMU-independent
processed trajectory.
The primary vehicle-model target is **inclination RMSE** on that Vicon subset.
Heading and total errors remain reported, but six-axis yaw is unobservable and
their improvements can come from bias handling rather than aircraft dynamics.

The loader preserves integer nanoseconds before subtracting the timestamp
origin. EuRoC's source timestamps have sub-microsecond clock quantization. VQF
requires a fixed update interval, so its adapter uses nominal 200 Hz only after
checking every timestamp and interval is within 1 µs of that grid. No IMU values
are interpolated or dropped. All other estimators use the measured intervals.
The maximum clock deviation is recorded in each sequence's provenance. Actual
sampling gaps or variable-rate inputs are rejected by this VQF protocol.

Every filter processes all IMU samples. Scores include valid ground-truth
samples after a five-second burn-in and the initial heading-alignment window.
There is no reference-derived flight selection. Evaluation applies one constant
yaw offset over the first second of valid ground truth, with no roll/pitch fit
or full-trajectory alignment. Heading scores measure subsequent drift. Summary
tables average per-sequence RMSEs with equal sequence weight; a longer recording
does not dominate the aggregate.

Reports contain `comparison.md`, `results.json` and `per_sequence.csv`. JSON
includes configurations, source and binary hashes, package versions, exact
selection, reference limitations, full/primary/validation summaries, paired
MAV-versus-VQF differences, and per-sequence RMSE/P95 errors. Sources and
executables must remain unchanged throughout a replay. Failures are reported
instead of removed from averages. C++ and Python timing scopes differ and are
not a cross-language runtime ranking.

The [CI workflow](../../.github/workflows/euroc-benchmark.yml) runs the smoke
comparison on relevant pull requests, the full suite on relevant `master`
pushes, and either scope through manual dispatch. It publishes comparison and
provenance artifacts and never asserts that the vehicle model must beat VQF.
