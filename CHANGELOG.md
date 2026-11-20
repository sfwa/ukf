# Changelog

User-facing changes to UKF, newest first. Each release has a section named
for its version; [RELEASING.md](RELEASING.md) says what a version promises and
how numerical changes are reported. Changes on master that are not yet
released are under **Unreleased**.

Sections within a release: **Added**, **Changed**, **Deprecated**,
**Removed**, **Fixed**, and **Numerical changes** (anything that can move a
filter's output, with the configurations affected, the expected effect, the
evidence, and whether to replay or retune).

## [Unreleased]

Nothing yet.

## [1.0.0-rc.1] - 2026-11-19

The first candidate for `1.0.0`, which will be the first numbered release.
UKF had no releases before it, so this section summarises what has changed
since the library's long-standing master of 2016 to 2025, for users of that
code. With no earlier release to compare against,
[RELEASING.md](RELEASING.md#the-baseline-for-the-first-release) records the
replay and code-size figures that later releases are compared with.

### Added

- `UKF/Version.h`, with `UKF_VERSION_MAJOR`, `_MINOR`, `_PATCH`,
  `_PRERELEASE`, `_STRING` and `UKF_VERSION`. The CMake package carries the
  same version, and `find_package(ukf 1.0)` accepts any later 1.x
  release (a pre-release such as `1.0.0-rc.1` comes before `1.0.0`) ([#99]).
- Install rules and a CMake package: `find_package(ukf)` gives `ukf::ukf`,
  `ukf::single` and `ukf::double` ([#96]). The install includes UKF's
  `LICENSE` under `share/licenses/ukf` and, with the downloaded Eigen, its
  `COPYING.*` licence files under `share/licenses/ukf/eigen3` ([#109]).
- CMake options `UKF_EIGEN_VERSION` and `UKF_EIGEN_URL` to download another
  Eigen (3.3 or newer, unpatched), and `UKF_USE_SYSTEM_EIGEN` now accepts an
  installed Eigen 3.3 to 5 ([#96]). `UKF_EIGEN_UNROLLING_LIMIT` sets
  `EIGEN_UNROLLING_LIMIT` for every target that links UKF ([#94]).
- Precision targets `ukf::single` and `ukf::double`, and opt-in project tools
  (`UKF_BUILD_TESTS`, `UKF_BUILD_BENCHMARKS`, `UKF_BUILD_EXAMPLES`) so a
  parent project gets only the header targets ([#73], [#88]).
- Quaternion measurement fields, with three-dimensional innovation and noise,
  and `InnovationVector` / `InnovationVectorType` ([#98]). *Experimental in
  1.0: see [RELEASING.md](RELEASING.md#experimental-interfaces).*
- `get_sigma_points()` on all three filters ([#97]).
- Dense measurement noise covariance and root factors
  (`UKF::DenseCovarianceMatrix`) ([#84]).
- Innovation gating: `step_with_gate` and `a_posteriori_step(gate, ...)`
  ([#85]). `step_with_gate_result` and `a_posteriori_step_with_gate_result`
  take the same arguments and return a `UKF::GateResult` (`Applied`,
  `Rejected`, `Empty` or `Failed`), which tells a rejected outlier from an
  empty measurement and from a numerical failure ([#104]).
- Measurement-only updates for independent observations made at the same
  time: `measurement_step(z, inputs...)` draws fresh sigma points around the
  current estimate and corrects with one observation, and
  `prepare_measurement_step()` does the first half for a split or gated
  correction ([#103]).
- `UKF::Angle`, a scalar state or measurement field with period 2π, for
  headings and bearings that cross ±π ([#105]).
- Fixed-lag and offline unscented Rauch–Tung–Striebel smoothing
  (`UKF/Smoother.h`), and `a_priori_step_with_cross_covariance` on `Core`
  and `SquareRootCore` to capture what it needs from a prediction ([#106]).
  *Experimental in 1.0: see
  [RELEASING.md](RELEASING.md#experimental-interfaces).*
- Scalar states and scalar measurements in all three filters ([#80]).
- A direct (discrete) process model as an alternative to an integrated
  derivative ([#63], [#75]).
- `is_field_set` on dynamic measurement vectors ([#61]).
- Single-precision unit tests (`unittest_single`) alongside double ([#94]).
- Runnable tracking and attitude examples, a usage guide ([#95]), attitude and
  EuRoC estimator benchmarks ([#87], [#89], [#91]) and embedded timing and
  code-size benchmarks for STM32H735 and ESP32-S31 ([#93]).

### Changed

- C++17 is required ([#74]).
- Filter steps report failure by returning `false` and leave the filter
  unchanged; the library no longer throws, and builds with exceptions
  disabled ([#90]). Code that caught exceptions must check return values.
- The standard UKF treats process noise as an additive covariance, as
  documented ([#77]).
- A successful correction with a nonempty measurement uses up its prediction
  and innovation. A second `innovation_step` or `a_posteriori_step` returns
  `false` until a new prediction, or `prepare_measurement_step()`, draws
  sigma points again; a rejected gate leaves the prediction available
  ([#101]). Code that corrected twice from one prediction should call
  `measurement_step` for each observation instead.
- A split-step caller's edit of the public innovation covariance, or its
  root, is checked before it is used. `Core` refuses a covariance that is not
  symmetric and positive semidefinite to within roundoff ([#101]); the
  square-root filters refuse a root that is not finite, and every filter
  refuses one that is not the size of the innovation ([#109]). A refused
  posterior returns `false`, or `GateResult::Failed`, and changes nothing.
- CI moved from Travis to GitHub Actions and now covers GCC 13 and 15, Clang
  18 and 22, macOS, Eigen 3.2.8, 3.4.0 and 5.0.0 ([#64], [#65], [#96]), and
  a build with the Ninja generator ([#109]).

### Fixed

- A square-root filter whose public innovation root had been set to NaN
  reported a successful correction, and `GateResult::Applied`, while leaving
  the estimate unchanged. Both square-root filters now refuse such a root
  ([#109]).
- `unittest`, `unittest_single` and `benchmark` build with the Ninja
  generator; the GoogleTest and Google Benchmark libraries they link are now
  declared as products of their external projects ([#109]).
- `Core` solves its gain in unit-variance coordinates, so a correction no
  longer depends on the units of each sensor, and a measurement channel with
  exactly zero variance cannot turn into a false independent observation by
  rounding. A field vector measured opposite to its prediction, exactly or
  nearly, gives a consistent residual and noise Jacobian, where an angled,
  reversed field gave a zero residual ([#101]).
- The `ahrs` and `sfwa_ukf` example C APIs: sensor innovations are read by
  key, missing accelerometer data is handled, the selected dynamics model
  and the controls reach every sigma point, the custom dynamics callback is
  honoured, and the X8 model's pitch inertia term is corrected ([#101]).
- The hand-held estimator of the attitude benchmarks keeps its velocity in
  the body frame, so that its zero-velocity pseudo-measurement no longer
  turns the heading ([#108]).
- GCC 15 and Clang 22 builds ([#96]).
- Field-vector rotations are taken between unit vectors, with the Jacobian
  scaled by the measured magnitude; a stored covariance is the repaired one;
  an exactly zero innovation covariance makes no correction; quaternions are
  renormalised after updates and sign-invariant ([#94]).
- Antiparallel field vectors use a deterministic perpendicular axis ([#81]);
  field-vector MRP rotations and their Jacobian are normalised ([#82]).
- Empty dynamic measurements skip the correction ([#83]); failed square-root
  downdates and exact observations are handled ([#78]); the covariance used
  for prediction is the one validated ([#76]); parameter noise keeps
  covariance correlations ([#79]).
- Linker errors from non-inline definitions ([#62]); GCC builds ([#40]);
  return type inference ([#43]); use as a CMake subdirectory ([#45]); the
  Eigen download URL ([#52]); unnecessary `template` keywords ([#56]).

### Numerical changes

Against the 2016 to 2025 master, estimates change in every configuration:

- **Additive process noise in `UKF::Core`** ([#77]): the predicted covariance
  of the standard filter now includes the process noise as documented. Users
  who tuned `process_noise_covariance` against the old behaviour should
  retune.
- **Field vectors and quaternions** ([#81], [#82], [#94]): filters with
  `UKF::FieldVector` measurements or quaternion states give different, and
  correct, innovations and corrections, independent of field magnitude.
  Replay recorded data before relying on old tuning.
- **Square-root filter speed-ups** ([#92], [#94]): results change only by
  rounding.

Against master of 17 November 2026 (83995d0, where this changelog was
first written), for anyone who has tuned against master since then:

- **Innovation solve in `UKF::Core`** ([#101]): the standard filter with
  sensors of very different units or variances in one measurement, a
  redundant measurement, or a channel of exactly zero variance. Gains,
  estimates and covariances change, toward the values the equations give:
  with unit prior and noise, a second channel rescaled by 10^10 recovers
  its gain of 0.5. Well-scaled measurements move by rounding only: the
  EuRoC, BROAD and Sassari replays reproduce the earlier figures to within
  4e-11 degrees. Replay recorded data if the sensors in one measurement
  differ in scale by many orders of magnitude.
- **Reversed field vectors** ([#101]): a `UKF::FieldVector` measured at or very
  near a half turn from its prediction. The residual and its noise
  Jacobian change there; other innovations do not.
- **Dense measurement noise without a field vector** ([#102]): the same
  values, computed without an identity transform.
- **The invalid innovation root fix** ([#109]): no filter given finite
  values of the right size changes; the unit tests compare a filter that
  went through a refused edit with one that did not.
- **The benchmark estimator, not the library** ([#108]): the hand-held
  rows of the attitude and EuRoC benchmark documents change (mean heading
  error on BROAD from 6.218 to 6.206 degrees without a magnetometer and
  from 2.244 to 2.156 with one; on Sassari from 3.862 to 3.502 and from
  3.689 to 3.511). The EuRoC MAV estimator's rows do not change. A project
  that copied `examples/attitude_benchmark/attitude_estimator.h` should
  take the change and replay.

### Contributors

Thanks to [cwecht](https://github.com/cwecht), [DieHertz](https://github.com/DieHertz),
[aronhoff](https://github.com/aronhoff), [tomcreutz](https://github.com/tomcreutz),
[johnschlag](https://github.com/johnschlag) and [remipch](https://github.com/remipch),
and to Erik Andreasson and Niclas Bäckman ([#73]), for their contributions.

[#40]: https://github.com/sfwa/ukf/pull/40
[#43]: https://github.com/sfwa/ukf/pull/43
[#45]: https://github.com/sfwa/ukf/pull/45
[#52]: https://github.com/sfwa/ukf/pull/52
[#56]: https://github.com/sfwa/ukf/pull/56
[#61]: https://github.com/sfwa/ukf/pull/61
[#62]: https://github.com/sfwa/ukf/pull/62
[#63]: https://github.com/sfwa/ukf/pull/63
[#64]: https://github.com/sfwa/ukf/pull/64
[#65]: https://github.com/sfwa/ukf/pull/65
[#73]: https://github.com/sfwa/ukf/pull/73
[#74]: https://github.com/sfwa/ukf/pull/74
[#75]: https://github.com/sfwa/ukf/pull/75
[#76]: https://github.com/sfwa/ukf/pull/76
[#77]: https://github.com/sfwa/ukf/pull/77
[#78]: https://github.com/sfwa/ukf/pull/78
[#79]: https://github.com/sfwa/ukf/pull/79
[#80]: https://github.com/sfwa/ukf/pull/80
[#81]: https://github.com/sfwa/ukf/pull/81
[#82]: https://github.com/sfwa/ukf/pull/82
[#83]: https://github.com/sfwa/ukf/pull/83
[#84]: https://github.com/sfwa/ukf/pull/84
[#85]: https://github.com/sfwa/ukf/pull/85
[#87]: https://github.com/sfwa/ukf/pull/87
[#88]: https://github.com/sfwa/ukf/pull/88
[#89]: https://github.com/sfwa/ukf/pull/89
[#90]: https://github.com/sfwa/ukf/pull/90
[#91]: https://github.com/sfwa/ukf/pull/91
[#92]: https://github.com/sfwa/ukf/pull/92
[#93]: https://github.com/sfwa/ukf/pull/93
[#94]: https://github.com/sfwa/ukf/pull/94
[#95]: https://github.com/sfwa/ukf/pull/95
[#96]: https://github.com/sfwa/ukf/pull/96
[#97]: https://github.com/sfwa/ukf/pull/97
[#98]: https://github.com/sfwa/ukf/pull/98
[#99]: https://github.com/sfwa/ukf/pull/99
[#101]: https://github.com/sfwa/ukf/pull/101
[#102]: https://github.com/sfwa/ukf/pull/102
[#103]: https://github.com/sfwa/ukf/pull/103
[#104]: https://github.com/sfwa/ukf/pull/104
[#105]: https://github.com/sfwa/ukf/pull/105
[#106]: https://github.com/sfwa/ukf/pull/106
[#108]: https://github.com/sfwa/ukf/pull/108
[#109]: https://github.com/sfwa/ukf/pull/109
