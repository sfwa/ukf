# UKF

A C++17 Unscented Kalman Filter library for state and parameter estimation,
using [Eigen](https://eigen.tuxfamily.org) for linear algebra. Fixed-size and
bounded measurement storage support embedded applications without Eigen heap
allocation during filter updates. The library does not use C++ exceptions and
builds with them disabled: every filter step returns whether it succeeded.

| Filter | Purpose | State uncertainty |
| --- | --- | --- |
| `UKF::Core` | State estimation [[1]](#ref1), including quaternion states [[2]](#ref2) | Covariance `P` |
| `UKF::SquareRootCore` | Square-root state estimation [[3]](#ref3) | Lower factor `S`, where `P = S Sᵀ` |
| `UKF::SquareRootParameterEstimationCore` | Parameter estimation without a process model [[3]](#ref3) | Lower factor `S` |

All three support scalar, vector, [angle](docs/usage-guide.md#heading-and-bearing) and quaternion measurements, fixed or dynamic measurement
sets, optional [dense measurement noise](docs/dense-measurement-covariance.md),
and [innovation gating](docs/innovation-gating.md). Each can be compiled in
single or double precision.

## Releases

[![Latest release](https://img.shields.io/github/v/release/sfwa/ukf?include_prereleases&sort=semver)](https://github.com/sfwa/ukf/releases)

UKF uses [Semantic Versioning](https://semver.org). The first release will be
`1.0.0`, after a release candidate `1.0.0-rc.1`; until then master is the
only version. Use the [latest release](https://github.com/sfwa/ukf/releases/latest)
and [pin it by its tag](RELEASING.md#using-a-release).

- [RELEASING.md](RELEASING.md): the compatibility contract (which interfaces,
  conventions, build options and configurations a release supports, and
  which are experimental), what patch, minor and major releases may change,
  how numerical changes are reported, and how releases are made.
- [CHANGELOG.md](CHANGELOG.md): what has changed in each release.
- [CONTRIBUTING.md](CONTRIBUTING.md): what a pull request needs.

## Build and run

Use a C++17 compiler, CMake and `patch`. The commands below assume CMake 3.20
or newer; the project files retain a CMake 3.10 minimum. The default build
fetches pinned Eigen 3.2.8 (with the repository's allocation and C++17
patches), Google Test and Google Benchmark sources as their targets are built;
[another Eigen](#choosing-an-eigen-version) can be used instead. CMake 4 is
supported. Visual Studio 2015 does not provide the required C++17 support; CI
checks GCC (13 and 15) and Clang (18 and 22) on Linux and Clang on macOS, in
Debug and Release, with Eigen 3.2.8, 3.4.0 and 5.0.0.

```sh
git clone https://github.com/sfwa/ukf.git
cd ukf
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target unittest unittest_single ukf_examples --parallel 2
ctest --test-dir build --output-on-failure
./build/examples/tracking/ukf_tracking_example
./build/examples/attitude/ukf_attitude_example
```

The commands above use a single-configuration generator such as Unix Makefiles
or Ninja. For a multi-configuration generator, also pass `--config Debug` to
`cmake --build` and `-C Debug` to `ctest`; executable locations include the
configuration directory. Always use a separate build directory.

`unittest` runs the tests in double precision and `unittest_single` runs the
same tests in single precision, with tolerances for each. The `ukf_examples`
target builds both runnable examples and the larger `ahrs` and `sfwa_ukf` C++
libraries. Tests, benchmarks and examples enable `EIGEN_NO_MALLOC`; use Debug
to retain Eigen's allocation assertions. Build the micro-benchmark separately
if needed:

```sh
cmake --build build --target benchmark --parallel 2
./build/benchmark/benchmark
```

## Using UKF in a CMake project

```cmake
add_subdirectory(path/to/ukf)
add_executable(my_model src/my_model.cpp)
target_link_libraries(my_model PRIVATE ukf::double)
```

`ukf::double` selects `UKF_DOUBLE_PRECISION`; use `ukf::single` for
`UKF_SINGLE_PRECISION`. Both targets propagate the UKF/Eigen include paths,
the C++17 requirement and any Eigen download dependency. If your application
already defines its precision, link `ukf::ukf` (also available as `ukf`)
instead: this target does not define either precision macro. Choose exactly one
precision consistently across interacting translation units; adding the
opposite macro does not override a precision target.

To use an installed copy instead, build and install UKF, then find its package:

```sh
cmake -S path/to/ukf -B build-ukf -DCMAKE_INSTALL_PREFIX=/opt/ukf
cmake --build build-ukf
cmake --install build-ukf
```

```cmake
find_package(ukf 1.0 REQUIRED CONFIG) # With CMAKE_PREFIX_PATH including /opt/ukf.
target_link_libraries(my_model PRIVATE ukf::double)
```

The package accepts a request for any earlier version with the same major
version, and sets `ukf_VERSION` and `ukf_VERSION_PRERELEASE`. A pre-release,
such as master's `1.0.0-dev` or `1.0.0-rc.1`, comes before `1.0.0`, so it does
not meet a request for `1.0`: find it with `find_package(ukf REQUIRED CONFIG)`. In code,
`UKF/Version.h` defines `UKF_VERSION_STRING` and `UKF_VERSION` (10000 times
the major version, plus 100 times the minor, plus the patch). To build
against a release rather than master, check out or fetch its tag, such as
`v1.0.0`; [RELEASING.md](RELEASING.md#using-a-release) has a FetchContent
example.

The package provides the same `ukf::ukf`, `ukf::single` and `ukf::double`
targets, carrying the unrolling limit the installed copy was configured with.
With the downloaded Eigen, its headers (patched, for 3.2.8) are installed
under `include/ukf/eigen3`; with `UKF_USE_SYSTEM_EIGEN`, the package finds
Eigen again when it is loaded. Install rules are generated when UKF is built
on its own, or with `UKF_INSTALL=ON` as a subdirectory.

Without CMake, add `include/` and the Eigen headers to the include path,
compile as C++17, and define exactly one of `UKF_SINGLE_PRECISION` or
`UKF_DOUBLE_PRECISION` for every translation unit, before including UKF
headers.

| CMake option | Default | Purpose |
| --- | --- | --- |
| `UKF_USE_SYSTEM_EIGEN` | `OFF` | Use an installed Eigen 3.3 or newer (`Eigen3::Eigen`) instead of downloading one |
| `UKF_EIGEN_VERSION` | `3.2.8` | Eigen release to download when not using an installed one: 3.2.8 (patched) or 3.3 and newer (unpatched) |
| `UKF_EIGEN_URL` | empty (GitLab release) | Eigen source archive, as a URL or local file, to download instead; set `UKF_EIGEN_VERSION` to its version |
| `UKF_EIGEN_UNROLLING_LIMIT` | empty (Eigen's 100) | Set `EIGEN_UNROLLING_LIMIT` for every target that links `ukf`; see below |
| `UKF_WARNINGS_AS_ERRORS` | `ON` | Treat warnings in UKF's own targets as errors |
| `UKF_BUILD_TESTS` | `ON` standalone, `OFF` as a subdirectory | Configure the unit tests and `check` target |
| `UKF_BUILD_BENCHMARKS` | `ON` standalone, `OFF` as a subdirectory | Configure the benchmark target |
| `UKF_BUILD_EXAMPLES` | `ON` standalone, `OFF` as a subdirectory | Configure the examples and benchmark executables |
| `UKF_INSTALL` | `ON` standalone, `OFF` as a subdirectory | Generate install rules and the `ukf` CMake package |

### Choosing an Eigen version

By default UKF downloads Eigen 3.2.8 and applies the repository's patches,
which keep its LLT free of heap allocation and let it compile as C++17; the
embedded benchmarks and code-size figures use it. To use another Eigen:

- Download another release with `-DUKF_EIGEN_VERSION=3.4.0` (or `5.0.0`, or
  any release from 3.3). Releases other than 3.2.8 are used unpatched. To
  build offline or from a mirror, also give `-DUKF_EIGEN_URL=` with the
  archive's URL or path.
- Use an installed Eigen 3.3 or newer, including Eigen 5, with
  `-DUKF_USE_SYSTEM_EIGEN=ON`. Choose between installations with
  `-DEigen3_DIR=` (the directory holding `Eigen3Config.cmake`) or
  `CMAKE_PREFIX_PATH`.

An installed Eigen package is used as installed; the bundled patches are not
applied to it. Square-root rank updates use UKF-owned bounded workspace, so
`EIGEN_NO_MALLOC` checks also apply with an unpatched Eigen. Warnings are
enabled privately on UKF's own targets and are not imposed on downstream
consumers. CI builds a [downstream consumer fixture](test/cmake/) that checks
both precision targets, manual precision, C++17 propagation, the unrolling
option and isolation from parent project settings.

### Single precision

Single precision is what most embedded targets run, and the full test suite
runs in both precisions. In single precision keep the sigma-point parameter
`AlphaSquared` well above 1e-6 (the tests use 1e-2): the outer weights,
1 / (2 alpha² (n + kappa)), multiply the rounding of every sigma point, so
with 1e-6 a mean is good to only a few tenths of a percent of its value. See
[sigma-point scaling](docs/usage-guide.md#process-models-and-inputs).

### Eigen unrolling and flash size

Eigen unrolls fixed-size expressions up to `EIGEN_UNROLLING_LIMIT` (100 by
default). A flash-limited application can lower it with
`-DUKF_EIGEN_UNROLLING_LIMIT=N`; `0` disables unrolling. Without CMake, define
`EIGEN_UNROLLING_LIMIT` on the command line. Either way it must be the same in
every translation unit of the program that includes Eigen, so set it for the
whole build rather than in a source file.

It is a trade between flash and time. For the float attitude estimator
(`examples/attitude_benchmark`) on a Cortex-M7 at `-O2` with GCC 13.2, `0`
saved 12.9 KB of text (63 KB rather than 76 KB for the whole program) and took
15 % more instructions an update for the UKF and 13 % more for the SR-UKF,
counted under QEMU; cycles on hardware may differ. Results change only by
roundoff. For comparison, `-Os` on the STM32H735 made updates 2.2 to 2.4 times
as long. To measure an application's own estimator,
`examples/embedded_code_size/measure.py --eigen-unrolling-limit N` and the
STM32H735 benchmark's `UKF_EIGEN_UNROLLING_LIMIT` cache variable apply the same
setting.

## Quick start

This complete scalar-position example uses a direct, discrete-time model. State
keys and sensor keys are deliberately different: `expected_measurement` defines
their relationship.

```cpp
#define UKF_DOUBLE_PRECISION
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"

enum StateKey { Position, Velocity };
enum MeasurementKey { PositionSensor = 100 };
using State = UKF::StateVector<
    UKF::Field<Position, real_t>, UKF::Field<Velocity, real_t>>;
using Measurement = UKF::FixedMeasurementVector<
    UKF::Field<PositionSensor, real_t>>;

namespace UKF {
template <> template <>
State State::process_model<>(real_t dt) const {
    State next = *this;
    next.set_field<Position>(get_field<Position>() + dt * get_field<Velocity>());
    return next;
}
template <> template <>
real_t Measurement::expected_measurement<State, PositionSensor>(const State& state) {
    return state.get_field<Position>();
}
}

int main() {
    UKF::Core<State, Measurement> filter; // Default: direct process_model.
    filter.state << 0, 1;
    filter.covariance.setIdentity();
    // Discrete process covariance for this one-second step.
    filter.process_noise_covariance = State::CovarianceMatrix::Identity() * real_t(0.01);
    filter.measurement_covariance.set_field<PositionSensor>(real_t(0.25)); // Variance.
    Measurement measurement;
    measurement.set_field<PositionSensor>(real_t(0.9));
    if (!filter.step(real_t(1), measurement)) {
        return 1; // The filter is unchanged.
    }
    std::cout << filter.state.get_field<Position>() << '\n'; // Estimated position.
}
```

Include Eigen geometry before UKF headers that use quaternions, as above. The
[usage guide](docs/usage-guide.md) explains continuous-time models, covariance
initialisation, error handling, prediction, confidence, missing measurements
and quaternion conventions.

## Runnable examples

| Example | Demonstrates |
| --- | --- |
| [Tracking](examples/tracking/) | Direct constant-velocity model; distinct state/sensor keys; dynamic correlated observations; calibrated synthetic confidence; outlier gating; raw versus estimated error; forecasting on a copy |
| [Attitude](examples/attitude/) | Quaternion state and angular velocity; user-defined kinematics and measurement mapping; directional observations and gyroscope measurements; square-root uncertainty |
| [Attitude benchmark](examples/attitude_benchmark/) | A complete gyroscope/accelerometer/magnetometer attitude estimator written once for both state cores, compared with Madgwick, Mahony and VQF on recorded data ([results](docs/attitude-benchmark.md)) |
| [EuRoC MAV benchmark](examples/euroc_benchmark/) | That estimator with a rotor-drag model of a hexrotor ([results](docs/euroc-benchmark.md)) |
| [ESP32-S31](examples/esp32s31_benchmark/) and [STM32H735](examples/stm32h735_benchmark/) timing | The attitude estimator on microcontrollers, with update-time percentiles against a 2 ms budget, code size, and float or double precision on the STM32H735 |
| [AHRS](examples/ahrs/) | A larger attitude filter and parameter estimator with C and Python wrappers |
| [SFWA UKF](examples/sfwa_ukf/) | A larger navigation filter with C and Python wrappers |

The tracking and attitude examples are deterministic and return a failure
status if their consistency checks fail; CTest runs them. Start with them, or
with the attitude benchmark's estimator, before adapting the larger wrappers.
CI builds the C++ targets; the Python wrappers are not part of the automated
tests.

Filtering estimates the current state from observations so far; the library does
not currently implement an RTS or fixed-lag smoother that revises past estimates
using later observations.

## References

<a name="ref1">[1]</a> "A New Extension of the Kalman Filter to Nonlinear Systems",
S. J. Julier and J. K. Uhlmann,
https://www.cs.unc.edu/~welch/kalman/media/pdf/Julier1997_SPIE_KF.pdf

<a name="ref2">[2]</a> "Unscented Filtering for Spacecraft Attitude Estimation",
John L. Crassidis and F. Landis Markley,
http://www.acsu.buffalo.edu/~johnc/uf_att.pdf

<a name="ref3">[3]</a> "The Square-Root Unscented Kalman Filter for State and Parameter-Estimation",
Rudolph van der Merwe and Eric A. Wan,
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf
