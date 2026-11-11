# Attitude benchmark: BROAD and Sassari

The [benchmark example](../examples/attitude_benchmark) compares two concrete
estimators built with this library against three published reference algorithms.
It measures attitude accuracy, not the general superiority of a Kalman-filter
implementation. Motion models, bias estimation, disturbance handling, and tuning
all affect the result.

## Estimators

| Label | Implementation | Configuration |
|---|---|---|
| `ukf` | `UKF::Core` | Quaternion + three gyro-bias states, double precision |
| `srukf` | `UKF::SquareRootCore` | Identical state, model, initialization and noise to `ukf` |
| `madgwick` | AHRS 0.4.0 | Fixed gain 0.05 in both modes |
| `mahony` | AHRS 0.4.0 | Fixed proportional gain 1.0, integral gain 0.3 |
| `vqf` | VQF 2.0.1 | Full online defaults, including bias estimation and magnetic disturbance rejection |

All estimators are causal. No optical reference, proprietary sensor orientation,
or future sample enters their updates. Propagation uses the mean of the previous
and current gyro readings over the interval ending at the current sample. VQF's
offline forward/backward algorithm is not used. The adapters preserve reference
implementation behavior, including AHRS's zero-gyro early return.

The UKF examples propagate a body-to-world quaternion using the exponential of
the bias-corrected angular increment. Bias follows a random walk. Measurement
updates compare normalized acceleration with world +Z and, in nine-axis mode,
normalized magnetometer readings with a fixed world magnetic direction. That
direction, including dip, is initialized from the first sensor sample. These
are deliberately simple fixed-noise estimators: there is no acceleration or
magnetic disturbance rejection. They are not the more elaborate dual-filter
model in `examples/ahrs`.

Both UKFs initialize attitude from the first accelerometer/magnetometer sample
and gyro bias to zero. Their initial attitude and bias standard deviations are
0.1 rad and 0.01 rad/s. Unit-vector measurement standard deviations are 0.1.
Gyro noise density is 0.02 rad/√s and bias random-walk density is
0.0002 (rad/s)/√s. The first-order discrete process covariance is
`diag(density²) * dt`; the SR-UKF receives its square root. Unscented parameters
are `alpha²=1`, `beta=2`, and `kappa=0`. These settings are fixed across datasets,
axes and trials; no trial-specific optimization or test-set tuning is performed.
The references use their corresponding sensor-only initialization, including
VQF's native initialization procedure.

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

## Recorded comparison

The complete measured tables and provenance are in
[the recorded comparison](benchmarks/attitude-comparison.md). Values are results
of this particular model and fixed parameter configuration. UKF and SR-UKF share
the model and settings, and their nine-axis aggregates match at the displayed
precision. Their six-axis BROAD results differ on several recordings; these
differences are retained rather than requiring equality under large heading
drift. The synthetic equivalence test covers benign motion. A stronger
disturbance-aware attitude model may improve the UKF example without changing
the filtering library.

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
UKF/SR equivalence and malformed input.

Timing is diagnostic only: the C++ measurement excludes CSV I/O, while the Python
reference-call measurement includes wrapper overhead. It is not a fair
cross-language speed ranking, and shared CI CPU timing is not a performance
regression gate. For embedded runtime comparisons, use equivalent native
implementations on the same target and measure their update loops separately.
