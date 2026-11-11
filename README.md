# UKF

A C++17 Unscented Kalman Filter library for state and parameter estimation,
using [Eigen](https://eigen.tuxfamily.org) for linear algebra. Fixed-size and
bounded measurement storage support embedded applications without Eigen heap
allocation during filter updates.

| Filter | Purpose | State uncertainty |
| --- | --- | --- |
| `UKF::Core` | State estimation [[1]](#ref1), including quaternion states [[2]](#ref2) | Covariance `P` |
| `UKF::SquareRootCore` | Square-root state estimation [[3]](#ref3) | Lower factor `S`, where `P = S Sᵀ` |
| `UKF::SquareRootParameterEstimationCore` | Parameter estimation without a process model [[3]](#ref3) | Lower factor `S` |

All three support scalar/vector measurements, fixed or dynamic measurement
sets, optional [dense measurement noise](docs/dense-measurement-covariance.md),
and [innovation gating](docs/innovation-gating.md).

## Build and run

Use a C++17 compiler, CMake and `patch`. The commands below assume CMake 3.20
or newer; the project files retain a CMake 3.10 minimum. The supplied build
fetches pinned Eigen, Google Test and Google Benchmark sources as their targets
are built. CMake 4 is supported. Visual Studio 2015 does not provide the required
C++17 support; the current CI checks GCC and Clang on Linux and Clang on macOS.

```sh
git clone https://github.com/sfwa/ukf.git
cd ukf
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target unittest ukf_examples --parallel 2
ctest --test-dir build --output-on-failure
./build/examples/tracking/ukf_tracking_example
./build/examples/attitude/ukf_attitude_example
```

The commands above use a single-configuration generator such as Unix Makefiles
or Ninja. For a multi-configuration generator, also pass `--config Debug` to
`cmake --build` and `-C Debug` to `ctest`; executable locations include the
configuration directory. Always use a separate build directory.

The `ukf_examples` target builds both runnable examples and the larger `ahrs`
and `sfwa_ukf` C++ libraries. The unit tests and examples enable
`EIGEN_NO_MALLOC`; use Debug to retain Eigen's allocation assertions. Build the
benchmark separately if needed:

```sh
cmake --build build --target benchmark --parallel 2
./build/benchmark/benchmark
```

For an application, add `include/` and the selected Eigen headers to its include
path, compile as C++17, and define exactly one of `UKF_SINGLE_PRECISION` or
`UKF_DOUBLE_PRECISION` consistently across its translation units. Set the macro
before including UKF headers, or supply it as a compiler definition. The example
CMake targets set it themselves. The current build does not provide an installed
CMake package or a system-Eigen selection option.

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
    filter.step(real_t(1), measurement);
    std::cout << filter.state.get_field<Position>() << '\n'; // Estimated position.
}
```

Include Eigen geometry before UKF headers that use quaternions, as above. The
[usage guide](docs/usage-guide.md) explains continuous-time models, covariance
initialisation, prediction, confidence, missing measurements and quaternion
conventions.

## Runnable examples

| Example | Demonstrates |
| --- | --- |
| [Tracking](examples/tracking/) | Direct constant-velocity model; distinct state/sensor keys; dynamic correlated observations; calibrated synthetic confidence; outlier gating; raw versus estimated error; forecasting on a copy |
| [Attitude](examples/attitude/) | Quaternion state and angular velocity; user-defined kinematics and measurement mapping; directional observations and gyroscope measurements; square-root uncertainty |
| [AHRS](examples/ahrs/) | A larger attitude filter and parameter estimator with C and Python wrappers |
| [SFWA UKF](examples/sfwa_ukf/) | A larger navigation filter with C and Python wrappers |

The two small examples are deterministic and return a failure status if their
consistency checks fail. Start with them before adapting the larger wrappers.
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
