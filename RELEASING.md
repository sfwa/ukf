# Releases and compatibility

UKF follows [Semantic Versioning](https://semver.org). This file says what a
version number promises, what is and is not covered by that promise, and how a
release is made. [CONTRIBUTING.md](CONTRIBUTING.md) says what a contributor
needs to do, which is much less. User-facing changes are listed in
[CHANGELOG.md](CHANGELOG.md).

## The compatibility contract

From `1.0.0`, the maintainers support the interfaces and behaviour listed
here. UKF is header-only, so the promise is **source compatibility**: code that
compiles and behaves as documented against one release of a major version
compiles and behaves as documented against every later release of that major
version, in a supported configuration. There is no promise of binary
compatibility between releases, and no promise of bitwise-identical numerical
results (see [Numerical changes](#numerical-changes)).

### Covered

- **Filters.** `UKF::Core`, `UKF::SquareRootCore` and
  `UKF::SquareRootParameterEstimationCore`: their template parameters and
  defaults, their public member functions (`step`, `step_with_gate`,
  `a_priori_step`, `innovation_step`, `a_posteriori_step`,
  `get_sigma_points`, `prepare_measurement_step`, `measurement_step`) and
  their public fields (`state`, `covariance` or `root_covariance`,
  `process_noise_covariance` or `process_noise_root_covariance`,
  `measurement_covariance` or `measurement_root_covariance`, `innovation`,
  `innovation_covariance` or `innovation_root_covariance`), and the member type
  `InnovationVectorType`.
- **State and measurement types.** `UKF::StateVector`,
  `UKF::FixedMeasurementVector`, `UKF::DynamicMeasurementVector`, `UKF::Field`,
  `UKF::Vector`, `UKF::Quaternion`, `UKF::FieldVector`,
  `UKF::DenseCovarianceMatrix` and the members documented in their headers and
  in [docs/usage-guide.md](docs/usage-guide.md), including
  `CovarianceMatrix`, `covariance_size()`, `InnovationVector`, `set_field`,
  `get_field` and `is_field_set`.
- **Specialisation points.** The templates an application specialises:
  `StateVector::process_model`, `StateVector::derivative`,
  `MeasurementVector::expected_measurement`, and the sigma-point and MRP parameters in `UKF::Parameters` (`AlphaSquared`,
  `Beta`, `Kappa`, `MRP_A`, `MRP_F`). Their names, signatures and meanings
  stay fixed within a major version.
- **Integrators.** `UKF::NotUsedIntegrator`, `UKF::IntegratorEuler`,
  `UKF::IntegratorHeun` and `UKF::IntegratorRK4`.
- **Gates.** The gate concept accepted by `step_with_gate` and
  `a_posteriori_step(gate, ...)`, as [docs/innovation-gating.md](docs/innovation-gating.md)
  describes it. The outcome-returning forms `step_with_gate_result` and
  `a_posteriori_step_with_gate_result` expose `UKF::GateResult`.
- **Conventions.** What every documented quantity means: covariance versus
  lower root factor (`P = S Sᵀ`), the additive process-noise model, the
  three-dimensional uncertainty of a quaternion and its left-perturbation MRP
  coordinates, the innovation of a quaternion or field-vector measurement,
  declared-field order for noise storage, and the sigma-point weights derived
  from `AlphaSquared`, `Beta` and `Kappa`. Default values of every parameter.
- **Error handling.** No exceptions: each step returns `true` on success, and
  on failure returns `false` and leaves the public state, covariance,
  innovation and innovation covariance as they were. A dependent step fails
  until the next successful prediction. `step()` is all or nothing.
- **Allocation.** With fixed-size and bounded measurement storage, filter
  updates make no Eigen heap allocation, and they pass `EIGEN_NO_MALLOC`
  checks with every supported Eigen.
- **Build interface.** The CMake targets `ukf::ukf` (and `ukf`), `ukf::single`
  and `ukf::double`; the installed package `find_package(ukf)`; the options
  `UKF_USE_SYSTEM_EIGEN`, `UKF_EIGEN_VERSION`, `UKF_EIGEN_URL`,
  `UKF_EIGEN_UNROLLING_LIMIT`, `UKF_WARNINGS_AS_ERRORS`, `UKF_BUILD_TESTS`,
  `UKF_BUILD_BENCHMARKS`, `UKF_BUILD_EXAMPLES` and `UKF_INSTALL` and their
  defaults; the precision macros `UKF_SINGLE_PRECISION` and
  `UKF_DOUBLE_PRECISION` and the type `real_t`; and the version macros in
  [`UKF/Version.h`](include/UKF/Version.h).
- **Supported configurations**, below.

### Not covered

- Anything in `namespace UKF::Detail`, and protected or private members of any
  class. Use `get_sigma_points()` rather than a derived class's view of
  `sigma_points`.
- Header layout beyond the documented `#include "UKF/..."` lines: which header
  happens to declare what, and what each one includes transitively.
- Numerical results to the last bit, timing, code size and stack use. They are
  measured and reported (see [Numerical changes](#numerical-changes) and
  [Release gates](#release-gates)) but not promised.
- The examples, benchmarks, benchmark harnesses, Python packages and
  datasets under `examples/`, `benchmark/` and `docs/`. They show how to use
  the library and are tested, but may change in any release.
- Configurations outside the supported list. They may work, and fixes are
  welcome, but a release may break them.

### Experimental interfaces

An interface marked experimental is in a release so that users can try it, but
it may change in a minor release. Experimental interfaces are listed here and
say so in their header; nothing else is experimental. An interface leaves this
list in a minor release, after which the contract covers it.

| Interface | Since | Why it is experimental |
| --- | --- | --- |
| Quaternion measurement fields (`UKF::Quaternion` in a measurement vector), and `InnovationVector` / `InnovationVectorType` where they differ from the measurement type | 1.0.0 | Added just before 1.0; the innovation and noise coordinates should see downstream use before they are fixed. |
| Smoothing (`UKF/Smoother.h`: `UKF::SmoothingEstimate`, `UKF::unscented_rts_step`, `UKF::FixedLagSmoother`) and `a_priori_step_with_cross_covariance` on `Core` and `SquareRootCore` | 1.0.0 | Added just before 1.0; the snapshot and cross-covariance interface should see downstream use before it is fixed. |

### Supported configurations

What CI tests is what is supported. A newer compiler or Eigen is not supported
until CI tests it, and a release notes when one is added.

| | Supported |
| --- | --- |
| Language | C++17 |
| Precision | `float` and `double` (`UKF_SINGLE_PRECISION`, `UKF_DOUBLE_PRECISION`) |
| Compilers | GCC 13 and 15, Clang 18 and 22 on Linux; Apple Clang on the current macOS runner |
| Eigen | 3.2.8 with this repository's patches (the default download), 3.4.0 and 5.0.0, downloaded or installed |
| CMake | 3.10 or newer to configure; 3.14 or newer to install the package version file as architecture-independent |
| Build types | Debug and Release |

MSVC and embedded cross-compilers are not in CI. The maintainers measure
Cortex-M7 (STM32H735) and ESP32-S31 builds for timing and size, and report what
they find, but these are not supported configurations. Raising a minimum in
this table (dropping GCC 13, Clang 18, Eigen 3.2.8 or C++17) is a major
change; adding a newer one is a minor change.

## What each part of the version means

| Release | Policy | Example |
| --- | --- | --- |
| **Patch: `1.2.3`** | Compatible corrections and implementation improvements | Fix a numerical error while keeping the documented model and interface; make an update faster |
| **Minor: `1.3.0`** | Compatible features and deprecations | Add a measurement type, an accessor or an optional algorithm; support a newer compiler or Eigen; deprecate an interface that has a replacement |
| **Major: `2.0.0`** | Changes that need downstream adaptation | Change a noise convention, failure handling or an existing default; remove a deprecated interface; raise a supported minimum |

The maintainers classify every change and choose each version.

### Numerical changes

A program that still compiles does not show that its estimator behaves as it
did. Every change that can move a filter's output (estimates, covariance,
innovation, or which steps succeed) gets its own entry under **Numerical
changes** in the changelog and the release notes, whatever kind of release it
is in. The entry says:

- which filters and configurations are affected (filter, precision, Eigen,
  field types);
- what the effect is expected to be, and its size;
- the regression evidence: tests, replays or benchmarks, with links;
- whether downstream users should replay recorded data or retune.

A numerical bug fix can be a patch even when estimates change: the incorrect
mathematics is not kept to preserve old output. A deliberate change to
documented model semantics or to an existing default is a major change. The
process-noise correction ([#77](https://github.com/sfwa/ukf/pull/77)) is an
example of the first; the move from exceptions to failure returns
([#90](https://github.com/sfwa/ukf/pull/90)) of the second.

### Deprecation

1. A replacement arrives in a minor release, with migration notes in the
   changelog.
2. The old interface is marked `[[deprecated("use ...")]]` in the same or a
   later minor release, and keeps working.
3. It is removed only in the next major release, and no sooner than six
   months after it was first marked deprecated.

### Support

The latest minor release of the current major version gets fixes. A previous
major version gets fixes only for a maintenance window announced in the
release notes of the major release that replaces it; without such an
announcement there is none.

## Cadence and branches

- A patch release follows when useful fixes are ready. A minor release is due
  every six to eight weeks while development is active, and is skipped when
  there is nothing to release. A major release follows a compatibility change
  that justifies it, not the calendar.
- `master` holds the next release. Contributions arrive on short-lived
  branches. Between releases, `UKF_VERSION_PRERELEASE` on master is `"dev"`
  and the version is the next one planned, so a build of master is never
  mistaken for a release.
- When a fix must ship without what is on master, a maintenance branch
  `release/X.Y` is cut from the tag `vX.Y.Z` and the fix is cherry-picked to it.
- Incompatible work waits on a branch until the maintainers open the next
  major version; it does not merge into master before then.
- A GitHub milestone named for each planned version shows what is meant to be
  in it.

## Release gates

A release candidate or release is made from a commit on which all of these
have passed:

- the [CMake workflow](.github/workflows/cmake.yml), every configuration in
  Debug and Release: the unit tests in both precisions (`unittest` and
  `unittest_single`), with `EIGEN_NO_MALLOC`; the runnable examples; the
  consumer fixtures as a subdirectory, with Eigen unrolling disabled, and from
  the installed package;
- for changes to filter mathematics since the last release, the estimator
  replays in the [attitude](.github/workflows/attitude-benchmark.yml) and
  [EuRoC](.github/workflows/euroc-benchmark.yml) benchmark workflows, compared
  with the last release;
- for performance-sensitive changes, the embedded timing or code-size figures
  ([examples/stm32h735_benchmark](examples/stm32h735_benchmark),
  [examples/embedded_code_size](examples/embedded_code_size)) compared with
  the last release, in the release notes;
- the downstream projects the maintainers know of (such as sfwa/models) built
  and tested against the candidate.

## Making a release

1. Check the gates above on the commit to release.
2. In one pull request:
   - set the version in [`include/UKF/Version.h`](include/UKF/Version.h):
     the numbers, `UKF_VERSION_PRERELEASE` (`"rc.1"` for a candidate, `""`
     for a release) and `UKF_VERSION_STRING` to match. CMake reads the version
     from this header and refuses to configure if the string and the numbers
     disagree, so the package version cannot drift from it;
   - rename `## [Unreleased]` in [CHANGELOG.md](CHANGELOG.md) to
     `## [X.Y.Z] - YYYY-MM-DD` (or `X.Y.Z-rc.N`) and start a new empty
     `## [Unreleased]`;
   - check that the README describes this release.
3. After it merges, tag that commit with an annotated tag and push the tag:

   ```sh
   git tag -a vX.Y.Z -m "UKF X.Y.Z"
   git push origin vX.Y.Z
   ```

   The tag push runs the CMake workflow again, including a check that the tag
   names the version in the header and that the changelog has its section.
4. Create a GitHub Release from the tag. Mark a candidate as a pre-release.
   The notes are curated, not generated: start from the changelog section and
   add known limitations, migration guidance, the numerical-change entries and
   the contributors to credit.
5. Close the milestone. On master, set the version to the next planned
   version with `UKF_VERSION_PRERELEASE "dev"`.

A published tag is never moved, deleted and re-created, or force-pushed. A
mistake in a release is fixed by the next release. The documentation for a
release is the repository as of its tag.

### Release notes template

```markdown
## UKF X.Y.Z

One paragraph: what this release is for and who should upgrade.

### Changes
(From CHANGELOG.md.)

### Numerical changes
(Each entry: affected configurations, expected effect, evidence, replay or retune advice. "None" if none.)

### Migration
(What a user must change, if anything.)

### Known limitations

### Contributors
Thanks to @... for ...
```

## Using a release

Pin a release by its tag. As a subdirectory or with FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(ukf GIT_REPOSITORY https://github.com/sfwa/ukf.git GIT_TAG v1.0.0)
FetchContent_MakeAvailable(ukf)
target_link_libraries(my_model PRIVATE ukf::double)
```

From an installed copy, ask for the major version you were written against;
the package accepts any later release with the same major version:

```cmake
find_package(ukf 1.0 REQUIRED CONFIG)
```

As in SemVer, a pre-release comes before its release: `1.0.0-rc.1` does not
meet a request for `1.0` or `1.0.0`, even with `EXACT`. To try a release
candidate, or master, find the package without a version and check
`ukf_VERSION_PRERELEASE`.

In code, `UKF_VERSION` is `MAJOR * 10000 + MINOR * 100 + PATCH`, and
`UKF_VERSION_STRING` includes any pre-release:

```cpp
#include "UKF/Version.h"
#if UKF_VERSION < 10100
/* Before 1.1.0. */
#endif
```

## The path to 1.0.0

1.0.0 has a fixed scope: the interfaces above, as they are on master. New
features go to 1.1.0.

1. **Interface review.** Check every covered interface above against its
   documentation, and decide which interfaces are experimental (the table
   above is the proposal).
2. **Installation and packaging**, done in [#96](https://github.com/sfwa/ukf/pull/96):
   install rules, the `ukf` CMake package and its version file, and the
   installed-consumer test.
3. **`v1.0.0-rc.1`**, then two to four weeks for downstream testing. Release
   blockers found in that time are fixed and a further candidate is made.
4. **`v1.0.0`**, from the last candidate with no blockers.
