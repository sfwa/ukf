# Attitude benchmark: BROAD and Sassari

The [benchmark example](../examples/attitude_benchmark) compares two baseline
estimators and three robust variants built with this library against three
published reference algorithms.
It measures attitude accuracy, not the general superiority of a Kalman-filter
implementation. Motion models, bias estimation, disturbance handling, and tuning
all affect the result.

## Estimators

| Label | Implementation | Configuration |
|---|---|---|
| `ukf` | `UKF::Core` | Quaternion + three gyro-bias states, double precision |
| `srukf` | `UKF::SquareRootCore` | Identical state, model, initialization and noise to `ukf` |
| `robust-rest` | Two `UKF::SquareRootCore` instances | Separate attitude and rest-only gyro-bias estimation; raw gravity direction; isolated magnetic heading |
| `robust-acc` | Same architecture | Adds causal acceleration filtering in a gyro-propagated frame |
| `robust` | Same architecture | Adds magnetic disturbance rejection and sensor-validity guards |
| `madgwick` | AHRS 0.4.0 | Fixed gain 0.05 in both modes |
| `mahony` | AHRS 0.4.0 | Fixed proportional gain 1.0, integral gain 0.3 |
| `vqf` | VQF 2.0.1 | Full online defaults, including bias estimation and magnetic disturbance rejection |

All estimators are causal. No optical reference, proprietary sensor orientation,
or future sample enters their updates. Propagation uses the mean of the previous
and current gyro readings over the interval ending at the current sample. VQF's
offline forward/backward algorithm is not used. The adapters preserve reference
implementation behavior, including AHRS's zero-gyro early return.

The baseline `ukf` and `srukf` examples propagate a body-to-world quaternion using the exponential of
the bias-corrected angular increment. Bias follows a random walk. Measurement
updates compare normalized acceleration with world +Z and, in nine-axis mode,
normalized magnetometer readings with a fixed world magnetic direction. That
direction, including dip, is initialized from the first sensor sample. These
are deliberately simple fixed-noise estimators: there is no acceleration or
magnetic disturbance rejection. They are not the more elaborate dual-filter
model in `examples/ahrs`.

Both baselines initialize attitude from the first accelerometer/magnetometer sample
and gyro bias to zero. Their initial attitude and bias standard deviations are
0.1 rad and 0.01 rad/s. Unit-vector measurement standard deviations are 0.1.
Gyro noise density is 0.02 rad/√s and bias random-walk density is
0.0002 (rad/s)/√s. The first-order discrete process covariance is
`diag(density²) * dt`; the SR-UKF receives its square root. Unscented parameters
are `alpha²=1`, `beta=2`, and `kappa=0`. These settings are fixed across datasets,
axes and trials; no trial-specific optimization or test-set tuning is performed.
The references use their corresponding sensor-only initialization, including
VQF's native initialization procedure.

### Robust estimator

[`robust_estimator.cpp`](../examples/attitude_benchmark/robust_estimator.cpp)
uses this library for a quaternion-only attitude SR-UKF and a separate three-axis
bias SR-UKF. Separating the bias state prevents translational acceleration or
magnetic disturbances from being interpreted as gyro bias. The attitude filter
propagates using the estimated bias; the bias filter receives direct gyro
observations only during detected rest. This is a different estimator model
from the baseline, not a change to the generic UKF implementation.

Rest detection requires stable accelerometer and gyro readings for 1.5 s.
Deviations are measured against causal 0.5 s low-pass references. Gyro deviation
and filtered absolute gyro limits are 2°/s; acceleration deviation is limited
to 0.5 m/s². Applying the absolute limit to the filtered gyro avoids rejecting
rest whenever measurement noise crosses the bias limit. Slow constant rotation
below the detection limits can still be ambiguous with gyro bias.

Bias process noise is 0.01 (°/s)/√s. The rest observation variance is chosen
from a target converged posterior standard deviation of 0.03°/s:
`R = P_target² / Q_step + P_target`, where `P_target = (0.03°/s)²` and
`Q_step` is the bias process variance over the observation interval. The
observation interval is nominally 0.05 s. This specifies the intended smoothing
behavior across sampling rates rather than treating 0.03°/s as the raw sensor
measurement noise.

`robust-acc` and `robust` rotate unnormalized acceleration into an auxiliary
frame propagated by the bias-corrected gyro. A causal second-order Butterworth
filter with a 3 s time constant estimates gravity there. It uses the exact
zero-order-hold continuous-state transition for each measured `dt`; it does not
assume a constant sampling rate. Filtering in this frame avoids lagging actual
body rotations. This design is inspired by the
[VQF paper](https://doi.org/10.1016/j.inffus.2022.10.014); the implementation is
independent and does not call VQF.

Gravity updates use two tangent-plane components at about 20 Hz. The effective
measurement noise is inflated for overlapping filtered observations. This is an
approximation to correlated measurement noise, not a whitening transformation.
The gravity gain is constrained to world-horizontal attitude corrections, so
gravity cannot directly correct yaw through covariance correlations. With this
constrained gain, the posterior covariance uses
`P - K*Cᵀ - C*Kᵀ + K*S*Kᵀ`, where `C` is the state/measurement cross covariance
and `S` is the innovation covariance. After applying the quaternion correction,
the covariance is transported into the new local attitude chart. Tests check
this transport against finite differences and the constrained update against
independent error samples.
The attitude covariance is conditional on the separately estimated bias;
diagnostic covariance values are filter-health information, not a claim of
fully calibrated joint attitude/bias confidence.

The final stage skips gravity and rest-bias updates when acceleration norm is
outside 0.7–1.3 g, and skips gravity corrections above 60° innovation. All
nonzero physical acceleration samples still enter the linear low-pass filter:
discarding large samples before filtering can bias cancellation over a motion
cycle. Zero acceleration holds the previous filter input while gyro propagation
continues. Invalid startup acceleration postpones gravity initialization until
a plausible reading arrives. These guards cannot distinguish every sustained
translation from gravity.

Magnetometer measurements update a separate scalar world-yaw correction:
`q_9D = q_yaw(delta) * q_6D`. They cannot change inclination, the acceleration
filter, or gyro bias. The final stage checks magnetic field norm and dip against
an accepted reference, requires an undisturbed dwell before resuming correction,
and can accept a persistent new field after sufficient rotation. This preserves
field magnitude for disturbance detection instead of normalizing it away.
A field change preserving norm and dip can remain undetectable. Rejection can
persist until the accepted field returns or a new reference passes the rotation
test; this implementation does not add VQF's prolonged-disturbance fallback.
Bias learning also requires rest, so recordings without rest and changing gyro
bias remain a limitation compared with VQF's additional motion-bias estimator.

The complete fixed configuration is available from
`attitude_benchmark --describe-filter robust` and captured in every report's
`cpp_parameters`. The first robust stage changes bias and heading architecture
as well as adding rest detection, so its difference from the baseline is not a
single-feature rest ablation. The three robust stages are replayable cumulative
variants.

## Datasets and coordinates

| Dataset | Full coverage | Input/reference handling |
|---|---|---|
| [BROAD](https://github.com/dlaidig/broad) | All 39 trials, including the 16 disturbance trials | Author-synchronized and aligned HDF5 data; use each file's sampling-rate attribute (~285.714 Hz), not a rounded 287 Hz |
| [Sassari](https://github.com/marcocaruso/mimu_optical_dataset_caruso_sassari) | Slow, medium and fast × XS1, XS2, AP1, AP2, SH1, SH2 = 18 recordings | v5 synchronized rows at 100 Hz; sensor gyro/acc/mag columns and optical `Qs` quaternion |

Downloads are pinned by commit/release, byte count, and SHA-256 in
[`datasets.lock.json`](../examples/attitude_benchmark/datasets.lock.json).
BROAD data is CC BY 4.0; Sassari supplies an MIT repository license. Please cite
the [BROAD paper](https://doi.org/10.3390/data6070072) and
[Sassari paper](https://doi.org/10.1109/JSEN.2020.3024806) when reusing the data.

All comparisons use scalar-first body-to-world quaternions and a z-up world.
BROAD supplies aligned ENU optical ground truth. Sassari uses the supplied `Qs`
unchanged: the authors nominally aligned local sensor/optical axes by mounting
them on a common board. Residual mounting/calibration errors are retained;
we do not fit a mounting transform to the evaluated trajectories. Sassari's
first sensor column mixes native timestamps and counters, so the documented
100 Hz synchronized sample clock is used instead. Proprietary IMU orientation
columns are ignored. Mahony and VQF's native magnetic-north conventions are
converted to the examples' +X convention.

## Evaluation protocol

Filters process every sensor sample, including rest and optical gaps. Scores
include only the authors' movement intervals with valid optical references,
after the first five seconds and after the heading-alignment window. Missing
optical references never remove sensor samples from filter propagation.

Optical world yaw differs from magnetic north and six-axis heading has an
arbitrary gauge. For **both** modes, compute one fixed yaw offset over the first
one-second window beginning at the first valid optical sample, then hold it
constant for the remainder of the recording. No roll/pitch alignment or
whole-trajectory fit is performed. Some BROAD trials have long leading optical
gaps (trial 37 first has a reference at about 31 s); the actual alignment start
is reported per trial. Drift before that window cannot be assessed by this
heading metric. The alignment is evaluation-only and never initializes or
corrects the estimator.

The relative quaternion is `q_est * conjugate(q_reference)`. Its z-axis twist
gives heading error, the remaining swing gives inclination error, and its
shortest rotation angle gives total error. Quaternion sign does not matter.
Reports contain RMSE and absolute 95th-percentile errors in degrees for each
recording. Heading is undefined at an exact 180° tilt; its scored count is
reported separately, while inclination/total still include that sample.
Raw unaligned heading and total errors are retained in JSON/CSV for diagnosis;
they are not generally magnetic-heading accuracy measurements.

Summary tables take the arithmetic mean of **per-trial RMSEs**, giving each
recording equal weight. Sassari's six sensors share the same physical motion;
these are 18 sensor recordings, not 18 independent motion experiments. Dataset
scores are kept separate. Six-axis heading measures drift after alignment;
absolute heading is unobservable without an additional reference. Nine-axis
heading scores also describe stability relative to initial alignment rather
than absolute magnetic-north calibration. This protocol is explicit and
reproducible but is not a reproduction of either publication's leaderboard.

New-estimator iteration uses a fixed development subset: BROAD 01, 08, 10, 24
and 28, plus Sassari slow at all six sensors. The remaining 34 BROAD and 12
Sassari recordings form the validation subset. All sensors from the same
Sassari motion stay together. Baseline results on all recordings were already
known when this split was established; validation means reserved from iteration
of the new estimator. Parameters are frozen before its first validation run.
Use `--split development` or `--split validation` to replay these subsets.

## Recorded comparison

The complete measured tables and provenance are in
[the recorded comparison](benchmarks/attitude-comparison.md). Values are results
of this particular model and fixed parameter configuration. UKF and SR-UKF share
the model and settings, and their nine-axis aggregates match at the displayed
precision. Their six-axis BROAD results differ on several recordings; these
differences are retained rather than requiring equality under large heading
drift. The synthetic equivalence test covers benign motion. The robust variants
test the effect of a disturbance-aware estimator model while retaining these
baselines for comparison.

To refresh the checked-in snapshot after running the two datasets separately:

```sh
python examples/attitude_benchmark/run.py --datasets broad --output benchmark-results/broad
python examples/attitude_benchmark/run.py --datasets sassari --output benchmark-results/sassari
python docs/benchmarks/update.py benchmark-results/broad/results.json benchmark-results/sassari/results.json
```

The snapshot updater rejects smoke, incomplete, or differently configured runs.

## Reproduction and CI

Follow the [example commands](../examples/attitude_benchmark/README.md) to build
and replay locally. `--smoke` runs complete BROAD 01 and Sassari slow/XS1
recordings, not truncated excerpts. It checks the integration and is not a
representative accuracy benchmark. Full replay comprises 57 recordings ×
2 modes × 8 algorithms/variants = **912 evaluations**.

The [GitHub Actions workflow](../.github/workflows/attitude-benchmark.yml) runs
smoke comparisons on relevant pull requests and full comparisons on relevant
pushes to `master`. Manual dispatch supports either scope. Separate BROAD and
Sassari jobs build the C++ example in Release mode, run synthetic/data/geometry
tests, cache verified downloads, and publish Markdown, JSON, and CSV artifacts
plus a job-summary table. Dataset files and reports generated locally are
ignored by Git; the recorded documentation snapshot is updated deliberately.

CI fails on invalid input, non-finite/non-unit filter output, failed synthetic
checks, missing scored reference data, or checksum mismatch. It does not assert
that one algorithm must win or impose an arbitrary error threshold. Review
accuracy changes through the per-trial artifacts, especially disturbance trials.
The geometry tests cover signs, coordinate frames, optical gaps and yaw drift;
executable tests cover tilted rest, multi-axis rotation, bias convergence,
UKF/SR equivalence and malformed input. Robust-estimator tests also exercise
prefix causality, irregular sampling, rest discrimination, acceleration
filtering, magnetic isolation/recovery and covariance geometry.

Timing is diagnostic only: the C++ measurement normally excludes CSV I/O;
`--diagnostics` includes state-trace writes in the measured loop. The Python
reference-call measurement includes wrapper overhead. It is not a fair
cross-language speed ranking, and shared CI CPU timing is not a performance
regression gate. For embedded runtime comparisons, use equivalent native
implementations on the same target and measure their update loops separately.
