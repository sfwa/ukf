# Attitude benchmark: BROAD and Sassari

The [benchmark example](../examples/attitude_benchmark) runs one attitude
estimator, built with this library, on both filter cores and compares it with
three published reference algorithms. It measures attitude accuracy, not the
general superiority of a Kalman-filter implementation. Motion models, bias
estimation, disturbance handling, and tuning all affect the result.

## Estimators

| Label | Implementation | Configuration |
|---|---|---|
| `ukf` | `UKF::Core` | Quaternion, velocity and gyro-bias states; separate one-state magnetic heading filter; double precision |
| `srukf` | `UKF::SquareRootCore` | Identical model, initialization and noise to `ukf` |
| `madgwick` | AHRS 0.4.0 | Fixed gain 0.05 in both modes |
| `mahony` | AHRS 0.4.0 | Fixed proportional gain 1.0, integral gain 0.3 |
| `vqf` | VQF 2.0.1 | Full online defaults, including bias estimation and magnetic disturbance rejection |

All estimators are causal. No optical reference, proprietary sensor orientation,
or future sample enters their updates. Propagation uses the mean of the previous
and current gyro readings over the interval ending at the current sample. VQF's
offline forward/backward algorithm is not used. The adapters preserve reference
implementation behavior, including AHRS's zero-gyro early return.

### UKF estimator

[`attitude_estimator.h`](../examples/attitude_benchmark/attitude_estimator.h)
is written once and instantiated for both cores. It uses the library as the
other examples do: a state vector, a dynamic measurement vector, a process
model integrated with `UKF::IntegratorRK4`, per-field measurement models, and
the cores' own `a_priori_step`, `innovation_step` and `a_posteriori_step`.

The attitude filter is a small strapdown inertial navigator. Its state is a
body-to-world quaternion, the velocity in the world frame, and the gyroscope
bias. The gyroscope reading, less the bias, turns the attitude; the
accelerometer reading, rotated into the world frame with gravity added, changes
the velocity. Both are process model inputs rather than measurements.

Tilt becomes observable through a pseudo-measurement that the velocity is zero.
A hand-held sensor never moves far, so its velocity stays small; a tilt error
makes gravity leak into the velocity estimate, where the pseudo-measurement
finds it. Translational accelerations average out within a few seconds, and are
integrated rather than mistaken for gravity. The gyroscope bias is learned from
the same velocity errors during motion. No acceleration low-pass filter,
auxiliary frame or measurement gating is involved; the response is set by the
process and measurement noise alone.

Rest detection requires stable accelerometer and gyro readings for 1.5 s.
Deviations are measured against 0.5 s running averages. Gyro deviation and
average gyro limits are 2°/s; acceleration deviation is limited to 0.5 m/s².
While the sensor is at rest, the gyroscope reading is added to the measurement
vector as a direct measurement of the bias. Slow constant rotation below the
detection limits can still be ambiguous with gyro bias.

Magnetometer measurements update a separate filter whose single state is a
world-yaw correction: `q_9D = q_yaw(offset) * q_6D`. They cannot change
inclination, velocity, or gyro bias. The magnetic field is rotated into the
world frame by the six-axis attitude; its horizontal direction is the
measurement. Field strength and dip angle are compared with an accepted
reference, and readings are ignored while either has strayed by more than 10%
or 10°, until the field has matched again for 0.5 s. A persistent new field is
adopted as the reference after 20 s of rotation above 20°/s. A field change
preserving strength and dip remains undetectable.

A six-axis estimator cannot observe heading, so the attitude filter's
uncertainty about its own heading would grow without limit. Once that
uncertainty approaches a radian the sigma points spread towards half a turn,
where the attitude error representation fails. The estimator therefore scales
back the heading row of the covariance whenever its standard deviation exceeds
0.3 rad. The limit changes no result on these recordings; a synthetic test
covers the case it exists for.

Noise is specified as densities and scaled to each measured time step, so the
tuning is independent of sample rate. The initial tilt, velocity and bias
standard deviations are 0.05 rad, 0.1 m/s and 0.5°/s. Gyro noise density is
0.003 rad/√s, combined in quadrature with 0.002 rad/√s per rad/s of rotation
rate; accelerometer noise is 0.05 m/s/√s; bias random walk is 0.001 (°/s)/√s. The zero-velocity
pseudo-measurement noise is 0.09 m/s·√s, and the rest gyroscope noise
0.012 (°/s)·√s. The heading offset has a random walk of 0.17°/√s and the
magnetic heading a noise of 3°·√s. The unscented transform uses the library's
default parameters. These settings are fixed across datasets, axes and trials.
The complete configuration is available from
`attitude_benchmark --describe-filter srukf` and captured in every report's
`cpp_parameters`.

The first accelerometer sample sets the initial tilt, and the first usable
magnetometer sample the initial heading. The references use their corresponding
sensor-only initialization, including VQF's native initialization procedure.

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

The UKF estimator's model and tuning were chosen on a fixed development subset:
BROAD 01, 08, 10, 24 and 28, plus Sassari slow at all six sensors. The remaining
34 BROAD and 12 Sassari recordings form the validation subset. All sensors from
the same Sassari motion stay together. The tuning was frozen before the
estimator's first validation run and has not been changed since. The magnetic
disturbance thresholds were carried over unchanged from the previous estimator.
The heading uncertainty limit was added after that run.
Reference-algorithm results on all recordings were already known. Use
`--split development` or `--split validation` to replay these subsets.

## Recorded comparison

The complete measured tables and provenance are in
[the recorded comparison](benchmarks/attitude-comparison.md). Values are results
of this particular model and fixed parameter configuration. `ukf` and `srukf`
share the model and settings, and agree at the displayed precision on every
aggregate.

The 2026-11-12 run completed all 570 evaluations. `ukf` and `srukf` differ by
less than 10⁻¹¹ degrees on every recording. Over the full datasets the UKF
estimator's total RMSE is 6.404° versus VQF's 6.448° on BROAD six-axis and
2.472° versus 2.593° on BROAD nine-axis. On Sassari it is 4.423° versus 4.323°
six-axis and 4.199° versus 4.235° nine-axis. On the validation recordings
alone, which played no part in tuning, the totals are 6.374° versus 6.277°
(BROAD six-axis), 2.486° versus 2.589° (BROAD nine-axis), 4.514° versus 4.343°
(Sassari six-axis) and 4.511° versus 4.580° (Sassari nine-axis). These results
demonstrate near parity with VQF, not superiority.

The differences are not uniform. The UKF estimator's inclination error is
higher than VQF's on BROAD (0.936° versus 0.857°) and lower on Sassari (1.861°
versus 2.095°), most clearly in the fast recordings (2.294° versus 2.820°). Its
six-axis heading drifts more on Sassari (3.862° versus 3.637°). On disturbed
nine-axis BROAD recordings its total error is 2.689° versus VQF's 3.185°; on
undisturbed ones it is 2.322° versus 2.181°. Madgwick and Mahony, run with
fixed gains, are well behind both.

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
2 modes × 5 algorithms = **570 evaluations**.

The [GitHub Actions workflow](../.github/workflows/attitude-benchmark.yml) runs
smoke comparisons on relevant pull requests and full comparisons on relevant
pushes to `master`. Manual dispatch supports either scope. Separate BROAD and
Sassari jobs build the C++ example in Release mode, run the synthetic and data
tests, cache verified downloads, and publish Markdown, JSON, and CSV artifacts
plus a job-summary table. Dataset files and reports generated locally are
ignored by Git; the recorded documentation snapshot is updated deliberately.

CI fails on invalid input, non-finite/non-unit filter output, failed synthetic
checks, missing scored reference data, or checksum mismatch. It does not assert
that one algorithm must win or impose an arbitrary error threshold. Review
accuracy changes through the per-trial artifacts, especially disturbance trials.
The metric tests cover signs, coordinate frames, optical gaps and yaw drift.
Executable tests on synthetic recordings cover tilted rest, multi-axis rotation,
bias convergence, agreement between the two cores, malformed input, prefix
causality, irregular sampling, rest discrimination, translation, accelerometer
dropout, magnetic isolation and recovery, and unobservable heading.

Timing is diagnostic only: the C++ measurement normally excludes CSV I/O;
`--diagnostics` includes state-trace writes in the measured loop. The Python
reference-call measurement includes wrapper overhead. It is not a fair
cross-language speed ranking, and shared CI CPU timing is not a performance
regression gate. For embedded runtime comparisons, use equivalent native
implementations on the same target and measure their update loops separately.
