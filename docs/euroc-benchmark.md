# EuRoC MAV attitude benchmark

The [EuRoC example](../examples/euroc_benchmark) compares an aircraft-specific
SR-UKF with the eight estimators/variants in the BROAD/Sassari benchmark. It also
includes a matched gravity-observation control and an adapted VQF configuration
so that changing gyro-bias assumptions is not mistaken for a process-model gain.
The goal is to measure whether rotor dynamics improve attitude estimation; the
benchmark does not require the new estimator to win.

## Aircraft and available measurements

[EuRoC](https://ethz-asl.github.io/datasets/euroc-mav/) was recorded using an
AscTec Firefly **hexrotor**, carrying a stereo visual-inertial sensor with an
ADIS16448 IMU. The published IMU stream contains angular velocity and specific
force at nominally 200 Hz. It does not contain magnetometer samples, rotor
speeds, commanded thrust, or control torques. This is therefore a **six-axis**
attitude comparison. Camera images are not estimator inputs.

A torque-driven rigid-body simulation cannot be evaluated from these inputs
without inventing actuator commands. Instead, the MAV estimator combines
measured angular velocity with a multirotor rotor-drag model. Collective
specific force is inferred from the measured acceleration along the thrust
axis, and transverse specific force supplies a drag observation. The model is
based on [Leishman et al., *Quadrotors and Accelerometers*, IEEE Control Systems
Magazine, 2014](https://doi.org/10.1109/MCS.2013.2287362), whose rotor-drag argument
also applies to multirotors with more than four rotors.

The IMU coordinate frame is not a conventional thrust-Z aircraft frame. The
nominal thrust direction is approximately `[cos(20°), 0, -sin(20°)]` in the
published IMU frame. The 20° angle is an empirical calibration from the **sensor
measurements in the two development recordings**, not a published mechanical
calibration or a fit to reference orientation. Estimated directions are 19.22°
on MH_01_easy and 20.15° on V1_02_medium. The same fixed mounting assumption is
used for every evaluation sequence.

## Process and observation model

The SR-UKF jointly estimates a quaternion, three body-velocity components and
three gyro-bias components, with a nine-dimensional error covariance. Each
sigma point propagates with its own bias-corrected gyro. This allows the motion
model to learn observable gyro-bias components during flight, including when
there is no useful initial rest interval. Detected rest adds a direct gyro-bias
observation; stable IMU measurements never impose a zero-velocity observation.

Let `q` rotate IMU/body coordinates into a z-up world, `v_B` denote body-frame
velocity relative to still air, `t` denote the unit thrust direction, and
`P = I - t*tᵀ` project onto the rotor plane. With gravity `g_W = [0,0,-9.80665]`,
angular rate `omega = gyro - gyro_bias`, collective specific force
`u = tᵀ*a_measured`, and drag
rate `d`, the multirotor model is

```text
q_next = q * Exp_quaternion(omega * dt)
dv_B/dt = R(q)ᵀ*g_W - omega × v_B + u*t - d*P*v_B
Eᵀ*a_measured = -d*Eᵀ*v_B + observation_error
```

`E` is an orthonormal basis for the rotor plane. Axial acceleration is a process
input; lateral acceleration is an observation. This introduces real coupling
between attitude and velocity through gravity and body rotation. It does not
reinterpret the complete accelerometer vector as a gravity direction during
powered flight.

The fixed nominal drag rate is 0.25/s. The default propagates all three body
velocity components; the optional planar approximation is not used in the
recorded comparison. Outside the rotor model's applicability gate, velocity
prediction is held while its uncertainty grows. This is an attitude estimator,
not a validated position/velocity navigation solution.

The unknown absolute world yaw remains a gauge. Measurement corrections cannot
directly rotate the quaternion about world Z. The constrained gain uses the
full covariance expression `P - K*Cᵀ - C*Kᵀ + K*S*Kᵀ`; after a quaternion
correction the complete covariance, including attitude/velocity cross terms,
is transported to the new local Gibbs chart.

Rotor operation is inferred causally from high-frequency accelerometer energy.
The model falls back to gravity filtering when rotor operation is not detected.
This matters because EuRoC recordings include handling, support on the ground,
spin-up and flight. Vibration indicates rotor operation, **not proof of free
flight**. Motor-on ground contact, wind, varying rotor speed, sensor lever-arm
acceleration and imperfect mounting calibration remain model errors. Axial
velocity is weakly observable; the recorded configuration specifies how it is
handled. No reference position, velocity or sensor-bias column is supplied to an
estimator.

## Reference configuration

`ukf`, `srukf`, `robust-rest`, `robust-acc`, `robust`, `madgwick`, `mahony`, and
`vqf` retain the configurations of the
[previous benchmark](attitude-benchmark.md). `mav-srukf-gravity` is a matched
control for the new estimator's bias/noise assumptions. `vqf-mav` uses online
VQF with a 10°/s bias limit and 5°/s initial bias standard deviation; its other
parameters retain their defaults. The native estimator also starts with a
5°/s bias prior and permits a 10°/s filtered gyro magnitude in rest detection;
that detection threshold is not a hard bias-state clip. These settings address gyro biases exceeding
4°/s in the development data, compared with the original 2°/s VQF limit. The
original default VQF result is also retained.

Both new SR-UKF variants permit motion-bias updates. Before rotor operation,
they use a short gravity-observation average that accounts for rotation during
the averaging interval. The control uses gyro-frame filtered gravity during
rotor operation; the MAV variant uses rotor drag and body-velocity dynamics.
Their observation-noise settings reflect different models. Filtered observation
noise is approximate; covariance traces are health diagnostics, not proof of
statistically calibrated uncertainty.

The executable's `--describe-filter` output and every report record the exact
fixed numerical settings. Sensor-only development calibration is distinguished
from evaluation. No reference orientation initializes a filter, no reference
bias is subtracted, and no whole-trajectory alignment corrects an estimate.

## Data and reference quality

The manifest pins 55 small sensor/calibration/reference files by repository
commit, byte count and SHA-256: approximately 85 MB for all eleven recordings.
Dataset binaries remain outside Git. Independent copies were checked for every
IMU and reference CSV. One mirror's V1_01_easy reference was a corrected
OpenVINS trajectory; the manifest deliberately retrieves and labels the
**canonical original** instead, verified against the original archive.

The supplied `state_groundtruth_estimate0` trajectories are offline batch
estimates. Their author-provided sensor metadata makes the distinction explicit:

| Group | Reference observations | Interpretation |
|---|---|---|
| Machine Hall | Leica positions and IMU | No independent optical orientation observation; attitude results are descriptive |
| Vicon rooms | Vicon poses and IMU | External orientation observations constrain the reference, but it still shares IMU data |
| V1_01_easy | Canonical Vicon/IMU estimate | [Known orientation-quality issue](https://docs.openvins.com/gs-datasets.html); retained separately |

The primary model comparison is **inclination RMSE on the five Vicon sequences
excluding V1_01_easy**. Heading and total errors are retained; six-axis heading
is drift after initial yaw alignment. All eleven sequences, environment groups,
and the known-problem sequence remain visible in the report. These are not
fully independent ground-truth attitude measurements or a visual-odometry
leaderboard.

Integer nanosecond timestamps are preserved before subtracting a local origin.
Reference quaternions are interpolated along the shortest arc only within
reference support and across gaps no larger than 20 ms. Missing reference
samples do not remove IMU samples. The recorded IMU clock has 256 ns
quantization. VQF receives a nominal 200 Hz clock only after verifying that
all timestamps remain within 1 microsecond of that grid; sensor samples are
neither interpolated nor resampled. Other estimators use measured intervals.

## Evaluation and reproduction

Development recordings are **MH_01_easy** and **V1_02_medium**. The other nine
recordings are reserved for validation of the new model; four belong to the
primary Vicon subset. Parameters are fixed before validation replay. Complete
recordings stay together, with no sample-level random split.

Every estimator processes every IMU sample. Scoring uses valid reference
orientations after a five-second burn-in and a one-second initial yaw-alignment
window beginning at the first valid reference. Alignment affects scoring only.
There is no reference-derived flight mask or trajectory-wide fit. Means give
each recording equal weight. Paired differences, ratios and sequence wins are
reported without claiming statistical significance from a small sample.

## Recorded outcome

The implementation completed all **121 evaluations** again on 2026-11-12 after
the CMake integration repair and a fix that eagerly evaluates the midpoint
velocity, avoiding a dangling Eigen expression. Model parameters were unchanged.
On the four primary Vicon validation sequences, mean inclination RMSE is
**1.130°** for `mav-srukf`, compared with **1.615°** for adapted `vqf-mav`,
**1.631°** for the joint-state gravity control, and **5.483°** for default VQF.
This is a 30.0% reduction versus adapted VQF and 30.7% versus the gravity control.
The MAV model still has lower inclination error on all four reserved sequences.
That comparison supports a process-model benefit beyond changing the bias prior
or allowable bias range; it is an observed effect on a small benchmark, not a
claim of statistical significance or universal superiority.

Across all eleven recordings, the MAV model's mean inclination/total RMSE is
1.393°/8.858°, versus adapted VQF's 1.946°/61.127°. The large heading differences
also reflect startup and motion-bias handling, so they should not be attributed
entirely to rotor dynamics. MH_01 and MH_02 retain substantial heading offsets;
adapted VQF has lower inclination error on MH_02 and the known-problem V1_01
reference, and some per-sequence heading or total errors favor the gravity
control. All results and both VQF configurations
are retained in the recorded tables. No model parameter changed after this
validation replay.

See the [example README](../examples/euroc_benchmark/README.md) for build and
replay commands and [the recorded comparison](euroc/comparison.md) for results.
The [CI workflow](../.github/workflows/euroc-benchmark.yml) runs complete-recording
smoke replays on relevant PRs and full replays on master or manual dispatch.
Failures, invalid output and missing coverage fail the run; an estimator losing
an accuracy comparison does not. Timing is diagnostic and is not a fair
cross-language speed ranking.

Please cite [Burri et al., *The EuRoC micro aerial vehicle datasets*,
2016](https://doi.org/10.1177/0278364915620033). The ETH dataset record specifies
[In Copyright – Non-Commercial Use Permitted](https://doi.org/10.3929/ethz-b-000690084);
the example's source-code license does not relicense the dataset.
