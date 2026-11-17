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

The first numbered release, planned as `1.0.0` after the candidate
`1.0.0-rc.1`. UKF had no releases before it, so this section summarises what
has changed since the library's long-standing master of 2016 to 2025, for
users of that code.

### Added

- `UKF/Version.h`, with `UKF_VERSION_MAJOR`, `_MINOR`, `_PATCH`,
  `_PRERELEASE`, `_STRING` and `UKF_VERSION`. The CMake package carries the
  same version, and `find_package(ukf 1.0)` accepts any later 1.x
  release (a pre-release such as `1.0.0-rc.1` comes before `1.0.0`).
- Install rules and a CMake package: `find_package(ukf)` gives `ukf::ukf`,
  `ukf::single` and `ukf::double` ([#96]).
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
  ([#85]).
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
- CI moved from Travis to GitHub Actions and now covers GCC 13 and 15, Clang
  18 and 22, macOS, and Eigen 3.2.8, 3.4.0 and 5.0.0 ([#64], [#65], [#96]).

### Fixed

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
