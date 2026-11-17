# EuRoC MAV attitude benchmark

The [EuRoC example](../examples/euroc_benchmark) compares an aircraft-specific
SR-UKF with the five estimators in the BROAD/Sassari benchmark. It also
includes a matched control without the rotor model and an adapted VQF
configuration, so that changing gyro-bias assumptions is not mistaken for a
process-model gain. The goal is to measure whether rotor dynamics improve
attitude estimation; the benchmark does not require the new estimator to win.

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
calibration or a fit to reference orientation. Estimated directions are 19.05°
on MH_01_easy and 20.11° on V1_02_medium, rounded to the nearest 5°
([`calibrate_mount.py`](../examples/euroc_benchmark/calibrate_mount.py)). The
same fixed mounting assumption is used for every evaluation sequence.

## Process and observation model

[`mav_estimator.h`](../examples/euroc_benchmark/mav_estimator.h) is the
[attitude benchmark's estimator](attitude-benchmark.md#ukf-estimator) with a
rotor model added. It uses `UKF::SquareRootCore` with the library's own
prediction and correction steps. The filter jointly estimates a quaternion,
three velocity components and three gyro-bias components, with a
nine-dimensional error covariance. Each sigma point propagates with its own
bias-corrected gyro. This allows the motion model to learn observable gyro-bias
components during flight, including when there is no useful initial rest
interval. Detected rest adds a direct gyro-bias observation.

Let `q` rotate IMU coordinates into a z-up world, `v` denote velocity expressed
in the IMU frame, `t` denote the unit thrust direction, and `P = I - t*tᵀ`
project onto the rotor plane. With gravity `g_W = [0,0,-9.80665]`, angular rate
`omega = gyro - gyro_bias`, measured specific force `a`, and drag rate `d`, the
model in powered flight is

```text
q_dot = 0.5 * q * [0, omega]
v_dot = R(q)ᵀ*g_W - omega × v + (tᵀ*a)*t - d*P*v
Eᵀ*a  = -d*Eᵀ*v + observation_error
```

`E` is an orthonormal basis for the rotor plane. Axial acceleration is a process
input; lateral acceleration is an observation. This introduces real coupling
between attitude and velocity through gravity and body rotation. It does not
reinterpret the complete accelerometer vector as a gravity direction during
powered flight.

When the rotor model does not apply, the complete accelerometer reading is a
process input and the velocity receives a zero pseudo-measurement, exactly as
in the hand-held estimator:

```text
v_dot = R(q)ᵀ*g_W - omega × v + a
0     = v + observation_error
```

Velocity is expressed in the IMU frame because rotor drag acts in that frame.
The observation then depends on velocity alone, and the unobservable heading
never enters the rotor model. The fixed nominal drag rate is 0.25/s. This is an
attitude estimator, not a validated position/velocity navigation solution.

The unknown absolute world yaw remains a gauge, and the filter's uncertainty
about it grows without limit. Once that uncertainty approaches a radian, the
sigma points spread towards half a turn, where the attitude error
representation fails. The estimator therefore scales back the heading row of
the root covariance whenever its standard deviation exceeds 0.3 rad. This
matters here: when a recording starts with the vehicle being carried, the gyro
bias (over 4°/s on this IMU) is unknown until the first rest, and the heading
uncertainty reaches the limit within seconds.

Rotor operation is inferred causally from high-frequency accelerometer energy:
the RMS of the signal above 20 Hz, averaged over 0.5 s, switches the rotor
state on above 0.8 m/s² and off below 0.5 m/s². The rotor model applies while
the rotors are on and the axial specific force is between 0.5 g and 2 g.
This matters because EuRoC recordings include handling, support on the ground,
spin-up and flight. Vibration indicates rotor operation, **not proof of free
flight**. Motor-on ground contact, wind, varying rotor speed, sensor lever-arm
acceleration and imperfect mounting calibration remain model errors. Axial
velocity is weakly observable. No reference position, velocity or sensor-bias
column is supplied to an estimator.

## Reference configuration

`ukf`, `srukf`, `madgwick`, `mahony`, and `vqf` retain the configurations of the
[attitude benchmark](attitude-benchmark.md). `mav-srukf-gravity` is the MAV
estimator with its rotor model switched off, and so a matched control for its
bias and noise assumptions. `vqf-mav` uses online VQF with a 10°/s bias limit
and 5°/s initial bias standard deviation; its other parameters retain their
defaults. The MAV estimator also starts with a 5°/s bias prior and permits a
10°/s average gyro magnitude in rest detection; that detection threshold is not
a hard bias-state clip. These settings address gyro biases exceeding 4°/s in
the development data, compared with the original 2°/s VQF limit. The original
default VQF result is also retained. The unmodified `ukf` and `srukf` keep
their 0.5°/s bias prior and 2°/s rest limit, so they never detect rest on this
IMU and their heading drifts accordingly.

The MAV estimator's gyro, accelerometer and zero-velocity noise densities are
those of the attitude benchmark. Its bias random walk is 0.01 (°/s)/√s. The
rest gyroscope noise is 0.09 (°/s)·√s, looser than for a hand-held sensor: the
bias measured with the rotors stopped shifts when they start. In powered
flight the velocity process noise is 0.2 m/s/√s and the rotor-drag observation
noise 0.067 m/s²·√s. The rest detector sees the IMU signals low-pass filtered
at 2 Hz, requires an accelerometer magnitude between 0.7 g and 1.3 g, and never
reports rest while the rotors are on.

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
recordings are reserved for validation of the model; four belong to the
primary Vicon subset. Complete recordings stay together, with no sample-level
random split.

Every estimator processes every IMU sample. Scoring uses valid reference
orientations after a five-second burn-in and a one-second initial yaw-alignment
window beginning at the first valid reference. Alignment affects scoring only.
There is no reference-derived flight mask or trajectory-wide fit. Means give
each recording equal weight. Paired differences, ratios and sequence wins are
reported without claiming statistical significance from a small sample.

## Recorded outcome

The 2026-11-12 run completed all **88 evaluations**. On the four primary Vicon
validation sequences, mean inclination RMSE is **1.172°** for `mav-srukf`,
compared with **1.615°** for adapted `vqf-mav`, **1.543°** for the control
without the rotor model, 1.511° for the unmodified `srukf`, and **5.483°** for
default VQF. This is 27% lower than adapted VQF and 24% lower than the control.
The MAV model has lower inclination error than adapted VQF on all four reserved
sequences, and than the control on three of them; on V2_01_easy the control is
better (1.132° versus 1.255°). That comparison supports a process-model benefit
beyond changing the bias prior or allowable bias range; it is an observed effect
on a small benchmark, not a claim of statistical significance or universal
superiority.

Across all eleven recordings, the MAV model's mean inclination/heading/total
RMSE is 1.228°/3.461°/3.710°, versus the control's 1.801°/8.116°/8.433° and
adapted VQF's 1.946°/60.936°/61.127°. Heading here is drift from an
unobservable initial value, and the large differences reflect startup and
motion-bias handling as much as rotor dynamics. Adapted VQF has lower
inclination error on the known-problem V1_01 reference (2.847° versus 3.025°).
The unmodified `ukf` and `srukf` reach 1.717° inclination with their hand-held
tuning, but their heading drifts by tens of degrees because they never measure
this IMU's bias at rest. All results and both VQF configurations are retained
in the recorded tables.

Two structural choices were made after validation results had been seen, so
the validation recordings are not untouched by this estimator's development.
Velocity was first expressed in the world frame; heading drift on validation
flights prompted a return to the IMU frame. The heading uncertainty limit was
added after a synthetic test exposed the failure which it prevents, and was
checked on the development recordings before the validation replay. No noise
parameter was tuned on a validation recording.

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
