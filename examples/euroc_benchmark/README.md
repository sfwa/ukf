# EuRoC MAV attitude benchmark

This example compares a vehicle-model SR-UKF for EuRoC's AscTec Firefly hexrotor
with the original UKF/SR-UKF, the three robust SR-UKF variants, Madgwick, Mahony
and online VQF from the [previous attitude benchmark](../attitude_benchmark).
The additional `vqf-mav` control changes VQF's gyro-bias clipping limit to 10°/s
and its initial bias standard deviation to 5°/s, widening its permitted bias
range and matching the MAV estimator's initial prior. It tests whether an
apparent advantage comes from
accommodating sensor bias rather than from vehicle dynamics; other VQF settings
remain at their online defaults.
The `mav-srukf-gravity` control uses the same gyro propagation, initial bias
prior, rest observations and attitude correction conventions as the MAV variant,
substituting filtered gravity observations for the drag/velocity model. Its measurement units and
noise differ accordingly; it is an estimator-model comparison, not a single
process-coefficient toggle.

EuRoC publishes gyro and accelerometer samples at nominally 200 Hz, but no
magnetometer or rotor-command stream. This benchmark therefore runs **six-axis
estimation only**. It does not synthesize magnetic readings, use camera images,
or feed reference position, velocity, orientation or bias into any estimator.
A physically appropriate process model does not guarantee an accuracy advantage;
the reports retain all results without a must-win test.

## Vehicle process model

The native estimator maintains a joint quaternion, body-frame velocity and
three-axis gyro-bias state, with a nine-dimensional error covariance. Let `q`
rotate the IMU frame to the world, `v` be velocity expressed in the IMU frame,
`b_g` be gyro bias, `s` the nominal rotor-thrust axis in that frame, and `P` an
orthonormal basis perpendicular to `s`. The multirotor approximation is

```text
omega = measured_gyro - b_g
q_dot = 0.5 * q * [0, omega]
v_dot = q^-1 * [0, 0, -g] - omega × v + u * s - d * P * P^T * v
b_g_dot = gyro-bias random walk
u = s^T * measured_acceleration
P^T * measured_acceleration = -d * P^T * v + observation noise
```

The process includes gravity, rotating-frame transport, axial collective force
and lateral rotor drag. It represents the Firefly's aggregate multirotor
translation rather than six individually commanded motors. EuRoC does not
publish motor commands, so axial specific force supplies `u`; its uncertainty
is absorbed into process noise. Lateral accelerometer components constrain
velocity and, through the process model, attitude. Velocity is an internal
state, not a reference input or an independently validated navigation output.

A causal accelerometer-vibration detector and plausible collective force enable
the rotor-drag model. Outside that condition, gravity observations are used.
This is an approximate operating-condition detector, not an externally supplied
flight mask. Model coefficients and the IMU-frame thrust direction are fixed
across sequences and reported by `--describe-filter`. The model assumes modest
drag around a nominal rotor operating point; wind, aggressive maneuvers and
unmodeled forces can invalidate that approximation.

Both MAV stages use the same gyro-bias random walk, 5°/s initial bias standard
deviation, and direct gyro observations during detected rest. Motion updates may
also correct gyro bias through the joint covariance in both stages. The
observations driving that correction differ: rotor drag in the MAV model and
filtered gravity in the control. The rest detector permits filtered gyro rates
up to 10°/s; this is a rest-classification threshold, not a hard clip on the
native bias state. No zero-velocity observation is inferred from quiet IMU data.

Corrections suppress direct world-yaw gain and transport the entire covariance,
including attitude/velocity/bias cross terms, into the corrected quaternion
chart. The compiled geometry test checks the constrained covariance against
error samples and the chart transport against finite differences.

## Build and replay

From the repository root, using Python 3.12 and a C++17 compiler:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r examples/attitude_benchmark/requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target attitude_benchmark euroc_benchmark --parallel 2
build/examples/attitude_benchmark/robust_geometry_test
build/examples/euroc_benchmark/mav_geometry_test
python -m unittest discover -s examples/attitude_benchmark -p 'test_*.py' -v
python -m unittest discover -s examples/euroc_benchmark -p 'test_*.py' -v
python examples/euroc_benchmark/run.py --smoke
python examples/euroc_benchmark/run.py
```

`--smoke` replays the complete MH_01_easy recording. It checks integration and
does not provide an independent optical-attitude ranking. A full run evaluates
all 11 recordings with eleven algorithms/variants: **121 sequence evaluations**.
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
  --algorithms robust mav-srukf vqf vqf-mav --diagnostics \
  --output benchmark-results/euroc-diagnostic
```

Quoted `--sequences` patterns match complete recording names. Unmatched patterns
fail. Custom selections are marked as such and cannot be confused with the full
suite. `--smoke` cannot be combined with a split or sequence selectors. Use
`--cache`, `--output`, `--executable` or `--mav-executable` to change paths.

The native MAV executable accepts the same sensor-only CSV format as the
previous benchmark:

```sh
build/examples/euroc_benchmark/euroc_benchmark \
  --input sensors.csv --output attitude.csv --filter mav-srukf --axes 6 \
  --diagnostics states.csv
build/examples/euroc_benchmark/euroc_benchmark --describe-filter mav-srukf
```

Input columns are `time,gx,gy,gz,ax,ay,az,mx,my,mz`: seconds, rad/s, m/s² and
three unused zero magnetic columns. Output is `time,qw,qx,qy,qz`, scalar-first
body-to-world orientation with world +Z up. World yaw is unobservable from the
IMU alone. Diagnostics are optional and are written only for supported native
estimators; their I/O is included in the diagnostic run's elapsed time.

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
