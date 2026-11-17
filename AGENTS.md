# AGENTS.md

Guidance for coding agents, and for people, working in this repository. The
[README](README.md) and the [usage guide](docs/usage-guide.md) are for users of
the library, [CONTRIBUTING.md](CONTRIBUTING.md) says what a pull request needs,
and [RELEASING.md](RELEASING.md) holds the compatibility contract and the
release process. Where this file and those differ, they win.

## What this is

UKF is a header-only C++17 library of unscented Kalman filters for small
processors: the standard filter (`UKF::Core`), the square-root filter
(`UKF::SquareRootCore`) and the square-root parameter estimation filter
(`UKF::SquareRootParameterEstimationCore`), with state vectors, fixed and
dynamic measurement vectors, and quaternion and field-vector types. It depends
only on Eigen.

It has two siblings. [NMPC](https://github.com/sfwa/nmpc) shares UKF's build
options, Eigen selection and CI matrix, and one program may use both.
[models](https://github.com/sfwa/models) builds the vehicles and controllers on
the two libraries.

## Layout

| Path | Contents |
| --- | --- |
| `include/UKF/` | The library. `Version.h` is the only place the version is written; `Config.h` selects `real_t`. |
| `test/` | GoogleTest tests, built twice: `unittest` in double precision and `unittest_single` in single precision. |
| `test/cmake/` | Consumer fixtures: as a subdirectory, and installed (`installed/`). |
| `benchmark/` | Google Benchmark micro-benchmarks. |
| `examples/` | Runnable `tracking` and `attitude` examples; the `ahrs` and `sfwa_ukf` libraries with C interfaces; the attitude and EuRoC estimator benchmarks (each with its own workflow); embedded timing and code-size benchmarks for STM32H735 and ESP32-S31, built with those boards' toolchains, not in CI. |
| `docs/` | The usage guide, feature notes and benchmark results. |
| `cmake/` | `ReadVersion.cmake` (identical in NMPC), the package configuration, and `EigenGccWarnings.h`. |
| `patch/` | Patches applied to the bundled Eigen 3.2.8. |

## Build and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target unittest unittest_single ukf_examples --parallel 2
ctest --test-dir build --output-on-failure
```

The default build downloads Eigen 3.2.8 and patches it, and downloads
GoogleTest and Google Benchmark from GitHub. The
[README](README.md#build-and-run) lists the options, including another Eigen
(`UKF_USE_SYSTEM_EIGEN`, `UKF_EIGEN_VERSION`) and the unrolling limit. Project
targets build with `-Wall -Wextra -Werror` and without exceptions, and tests,
benchmarks and examples define `EIGEN_NO_MALLOC`.

If you change CMake files, also run the consumer fixtures, as CI does:

```sh
cmake -S test/cmake -B build-consumer && cmake --build build-consumer && ctest --test-dir build-consumer
```

The CI workflow, `.github/workflows/cmake.yml`, has the installed-package
check.

## Formatting and lint

Both are standard tools with stock settings, checked in CI by the `lint` job.

- **clang-format 18**, Google style, unmodified (`.clang-format`). Format
  everything before committing:

  ```sh
  git ls-files '*.h' '*.cpp' '*.c' | xargs clang-format-18 -i
  ```

  Other clang-format versions format a little differently, so use 18. Where
  include order matters (a vendor header that must come first), wrap the
  block in `// clang-format off` and `// clang-format on`, with a comment
  saying why, as the embedded examples do.

- **clang-tidy 18**, with the `bugprone-*`, `performance-*` and
  `clang-analyzer-*` groups as errors (`.clang-tidy` lists the four checks
  that are off, and why). It runs over the tests, benchmarks and examples,
  and through them over the headers:

  ```sh
  cmake -S . -B build-tidy -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
  cmake --build build-tidy --target eigen3 googletest googlebenchmark
  run-clang-tidy-18 -p build-tidy -quiet "$PWD/(test|benchmark|examples)/"
  ```

  Fix a finding rather than silencing it. If a `NOLINT` is unavoidable, name
  the check and give the reason on the same line.

## Rules of the code

- **No exceptions or heap in a filter step.** A step that fails returns
  `false` and leaves the filter unchanged. Project targets build with
  `-fno-exceptions`, and `EIGEN_NO_MALLOC` catches allocations in Debug.
- **Both precisions.** `real_t` is `float` or `double`, chosen by
  `UKF_SINGLE_PRECISION` or `UKF_DOUBLE_PRECISION` (the `ukf::single` and
  `ukf::double` targets). Write literals as `real_t(0.5)` and use `std::`
  maths functions, so that single precision stays in single precision. A
  change to filter behaviour is tested in both `unittest` and
  `unittest_single`.
- **Eigen 3.2.8 to 5.0.** The bundled Eigen 3.2.8 is the oldest supported. It
  has no `Eigen::Index` (use `std::ptrdiff_t`).
- **CMake 3.10** is the minimum. Keep new CMake within it.
- **Keep step with NMPC.** The build options, the Eigen selection, the CI
  matrix, the lint setup, `cmake/ReadVersion.cmake` and `EigenGccWarnings.h`
  match NMPC's. A change to one of them belongs in both repositories.
- **Comments** explain why. Write them as full sentences.

## Versions, compatibility and the changelog

- Read [RELEASING.md](RELEASING.md) before changing anything under
  `include/UKF/`. Covered interfaces cannot change incompatibly outside a
  major release. `UKF::Detail` is not covered. Quaternion measurement fields
  are experimental and may change in a minor release.
- Say in the pull request if a change can move a filter's estimates,
  covariance or innovation, even when it is a fix.
- Fill in the pull request template's changelog line and compatibility
  guess. The maintainers write `CHANGELOG.md` and set versions. Do not edit
  `Version.h` or a released section of the changelog unless you are making a
  release.
- Never create, move or delete tags or GitHub Releases. Release tags are
  immutable and made only by the maintainers, following RELEASING.md.

## Before you finish

1. Formatted with clang-format 18, and clang-tidy 18 reports nothing.
2. Tests and examples build with warnings as errors, and `ctest` passes in
   both precisions.
3. Consumer fixtures, if CMake files changed.
4. Documentation for anything a user can see: the README, the usage guide or
   the header comment.
