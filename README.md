# UKF

Unscented Kalman filter library. Several different UKF implementations are
provided:

  * Standard Unscented Kalman Filter for state estimation, as originally
    described in [[1]](#ref1), with extensions for quaternions as described
    in [[2]](#ref2).
  * Square-root Unscented Kalman Filter for state estimation, implemented as
    described in [[3]](#ref3).
  * Optimised form of square-root Unscented Kalman filter for parameter
    estimation, implemented as described in [[3]](#ref3).

This library makes use of the [Eigen](https://eigen.tuxfamily.org) library
for linear algebra routines and matrix and vector operations. Heavy use is
made of C++11 and C++14 features in an attempt to avoid any dynamic memory
allocations and maximise opportunities for compile-time optimisations.

A primary goal of the library is to provide efficient UKF implementations for
use on embedded systems, so there is a strong focus on having minimal
dependencies and avoiding non-deterministic operations.

The filter can be compiled using either single or double precision by
choosing one of the following preprocessor definitions:

  * UKF_SINGLE_PRECISION
  * UKF_DOUBLE_PRECISION

## Usage

The library contains a number of class templates which are to be specialised
for the particular application. There are three main classes which need to be
specialised to make up an implementation:

  * State vector
  * Measurement vector
  * Core

The following examples are all derived from the [unit tests](test/), so have
a look at them for more detail.

### State vector

The state vector and measurement vector are made up of a number of fields,
each of which contains a key and a type. The reason for this is to allow
the filter to handle quaternions in the state vector transparently, so
something like the following is allowed:

```C++
enum MyFields {
    LatLon,
    Altitude,
    Velocity,
    Attitude
};

using MyStateVector = UKF::StateVector<
    UKF::Field<LatLon, UKF::Vector<2>>,
    UKF::Field<Altitude, real_t>,
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>
>;
```

Internally, these fields are all stored together as one contiguous Eigen
column vector, and all key lookups are done at compile time.

UKF scaling parameters can be adjusted in the following way:

```C++
template <> constexpr real_t UKF::Parameters::AlphaSquared<MyStateVector> = 1e-6;
```

For a description of what the scaling parameters do, see [[2]](#ref2) or read
the comments in the [code](include/UKF/StateVector.h). The weights and sigma
point spread are derived from `AlphaSquared`, `Kappa` and `Beta` (via
`ScaledDimension`, which is n + lambda), so adjust those rather than `Lambda`.

### Measurement vector

The measurement vector can be specialised in a similar way, but with the
choice of a fixed or dynamic measurement vector:

```C++
enum MyFields {
    StaticPressure,
    DynamicPressure,
    Accelerometer,
    Gyroscope
};

using MyMeasurementVector = UKF::FixedMeasurementVector<
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Gyroscope, UKF::Vector<3>>,
    UKF::Field<StaticPressure, real_t>,
    UKF::Field<DynamicPressure, real_t>
>;
```

or:

```C++
enum MyFields {
    StaticPressure,
    DynamicPressure,
    Accelerometer,
    Gyroscope
};

using MyMeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Gyroscope, UKF::Vector<3>>,
    UKF::Field<StaticPressure, real_t>,
    UKF::Field<DynamicPressure, real_t>
>;
```

For the fixed measurement vector, all measurements have to be provided every
filter iteration. The dynamic measurement vector allows for a filter where
some measurements are not available at every iteration, and so should only be
fused when they are available.

There is a small performance overhead for the dynamic measurement vector,
but it does not do any dynamic memory allocation.

### Core

The Core class contains the filter state and step function. It can be
specialised as follows:

```C++
using MyCore = UKF::Core<
    MyStateVector,
    MyMeasurementVector,
    UKF::IntegratorRK4
>;
```

Here, the user-specialised state vector and measurement vector classes are
provided as template parameters, along with the integrator to use for the
process model. For a list of available integrators, see
[Integrator.h](include/UKF/Integrator.h). If the integrator is omitted (or is
`UKF::NotUsedIntegrator`), both `Core` and `SquareRootCore` call a direct
`process_model(delta, input...)` specialisation instead of integrating a
derivative.

The square-root state estimation filter can be used instead like this:

```C++
using MyCore = UKF::SquareRootCore<
    MyStateVector,
    MyMeasurementVector,
    UKF::IntegratorRK4
>;
```

Specialisations of the process model (for state estimation filters) and
measurement model must also be provided.

### Process model

The process model is implemented as an ODE, with a user-provided function to
calculate the derivative. The ODE is then solved using the integrator method
specified in the Core class specialisation.

Here is an example of a process model for a simple state vector:

```C++
using ProcessModelTestStateVector = UKF::StateVector<
    UKF::Field<Position, UKF::Vector<3>>,
    UKF::Field<Velocity, UKF::Vector<3>>
>;

template <> template <>
ProcessModelTestStateVector ProcessModelTestStateVector::derivative<>() const {
    ProcessModelTestStateVector temp;
    /* Position derivative. */
    temp.set_field<Position>(get_field<Velocity>());

    /* Velocity derivative. */
    temp.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));

    return temp;
}
```

Also, the process model can take an arbitrary number of user-specified
inputs, like this:

```C++
template <> template <>
ProcessModelTestStateVector ProcessModelTestStateVector::derivative<UKF::Vector<3>>(
        const UKF::Vector<3>& acceleration) const {
    ProcessModelTestStateVector temp;
    /* Position derivative. */
    temp.set_field<Position>(get_field<Velocity>());

    /* Velocity derivative. */
    temp.set_field<Velocity>(acceleration);

    return temp;
}
```

If the Core class specialisation has no integrator, the process model is
instead a function which returns the next state directly. Leave its template
argument list empty; the input types are deduced from the parameters:

```C++
template <> template <>
ProcessModelTestStateVector ProcessModelTestStateVector::process_model<>(
        real_t delta, const UKF::Vector<3>& acceleration) const {
    ProcessModelTestStateVector next = *this;
    next.set_field<Position>(get_field<Position>() + delta * get_field<Velocity>());
    next.set_field<Velocity>(get_field<Velocity>() + delta * acceleration);

    return next;
}
```

### Measurement model

The measurement model is specified per field, in order to allow the expected
measurement vector to be constructed for the dynamic measurement vector where
not all measurements may be available each iteration. Each measurement model
function takes a state vector as an input.

The state vector type is provided to the measurement model specialisation as
a template parameter; this allows a measurement vector class to be shared
across multiple state vectors, with difference process model defined for
each.

Here is an example of a measurement model for a simple measurement vector and
state vector:

```C++
using MyMeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<StaticPressure, real_t>,
    UKF::Field<DynamicPressure, real_t>,
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Gyroscope, UKF::Vector<3>>
>;

using MyStateVector = UKF::StateVector<
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<AngularVelocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>,
    UKF::Field<Altitude, real_t>
>;

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -9.8);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope>(const MyStateVector& state) {
    return state.get_field<AngularVelocity>();
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, StaticPressure>(const MyStateVector& state) {
    return 101.3 - 1.2*(state.get_field<Altitude>() / 100.0);
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, DynamicPressure>(const MyStateVector& state) {
    return 0.5 * 1.225 * state.get_field<Velocity>().squaredNorm();
}
```

As with the process model, the measurement model can take an arbitrary number
of user-specified inputs:

```C++
template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer, UKF::Vector<3>>(const MyStateVector& state, const UKF::Vector<3>& input) {
    return state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -9.8) + input;
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope, UKF::Vector<3>>(const MyStateVector& state, const UKF::Vector<3>& input) {
    return state.get_field<AngularVelocity>();
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, StaticPressure, UKF::Vector<3>>(const MyStateVector& state, const UKF::Vector<3>& input) {
    return 101.3 - 1.2*(state.get_field<Altitude>() / 100.0);
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, DynamicPressure, UKF::Vector<3>>(const MyStateVector& state, const UKF::Vector<3>& input) {
    return 0.5 * 1.225 * state.get_field<Velocity>().squaredNorm();
}
```

### Initialisation

The filter state, covariance, process noise covariance and measurement noise
covariance should be initialised to appropriate values, e.g.:

```C++
MyCore test_filter;
test_filter.state.set_field<Position>(UKF::Vector<3>(0, 0, 0));
test_filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
test_filter.state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
test_filter.state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));
test_filter.covariance = MyStateVector::CovarianceMatrix::Zero();
test_filter.covariance.diagonal() << 10000, 10000, 10000, 100, 100, 100, 1, 1, 5, 10, 10, 10;
test_filter.process_noise_covariance = MyStateVector::CovarianceMatrix::Identity()*1e-5;
test_filter.measurement_covariance << 10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;

```

Or, for the SR-UKF:
```C++
MyCore test_filter;
test_filter.state.set_field<Position>(UKF::Vector<3>(0, 0, 0));
test_filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
test_filter.state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
test_filter.state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));
test_filter.root_covariance = MyStateVector::CovarianceMatrix::Zero();
test_filter.root_covariance.diagonal() << 100, 100, 100, 10, 10, 10, 1, 1, 2.2, 3.2, 3.2, 3.2;
test_filter.process_noise_root_covariance = MyStateVector::CovarianceMatrix::Identity()*3e-2;
test_filter.measurement_root_covariance << 10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;
test_filter.measurement_root_covariance = test_filter.measurement_root_covariance.array().sqrt();
```

Currently, only a diagonal measurement noise covariance matrix is supported.

Both state estimation filters use additive process noise:
`x_next = process_model(x, dt, inputs) + noise`. The process noise covariance
is added after propagation, and the measurement update includes that noise.
Supply the discrete-time covariance for each step; continuous-time noise
densities must be discretised by the application. The square-root filter
takes a factor whose product with its transpose is this same covariance.

Square-root correction steps support positive semidefinite posterior
covariances, including exact observations.

### Error handling

The library does not use C++ exceptions and can be built with them disabled.
Each of `step()`, `a_priori_step()`, `innovation_step()` and
`a_posteriori_step()` returns `true` on success. On failure it returns `false`
and leaves the public state, (root) covariance, innovation and innovation
(root) covariance exactly as they were. A step fails when a covariance (or a
covariance downdate) is indefinite beyond roundoff, or when a result is not
finite. With a small `AlphaSquared`, the central sigma point weight is negative,
and a strongly nonlinear model can make the predicted covariance indefinite;
both cores report this the same way. A step that depends on a failed (or
missing) earlier step also returns `false` until the next successful
`a_priori_step()`. `step()` is all-or-nothing: if its correction fails after a
successful prediction, the prediction is rolled back too. To keep the
prediction when a measurement cannot be applied, call the split steps below.

### Iteration

The general steps for carrying out a filter iteration are something like:

```C++
MyMeasurementVector meas;
meas.set_field<Accelerometer>(UKF::Vector<3>(0.0, 0.0, 9.8));
meas.set_field<Gyroscope>(UKF::Vector<3>(0.0, 0.0, 0.0));
meas.set_field<StaticPressure>(101300.0);
meas.set_field<DynamicPressure>(101300.0);

if (!test_filter.step(0.01, meas)) {
    /* The filter is unchanged; reinitialise or skip this measurement. */
}
```

Or, if it's necessary to do things with the internal filter state (e.g.
filter health monitoring), then iteration can be split up into three steps:

```C++
MyMeasurementVector meas;
meas.set_field<Accelerometer>(UKF::Vector<3>(0.0, 0.0, 9.8));
meas.set_field<Gyroscope>(UKF::Vector<3>(0.0, 0.0, 0.0));
meas.set_field<StaticPressure>(101300.0);
meas.set_field<DynamicPressure>(101300.0);

test_filter.a_priori_step(0.01);
/*
At this point, the state and covariance (or root_covariance) variables
reflect the a priori state and (root_)covariance.
*/
test_filter.innovation_step(meas);
/* Innovation and innovation_(root_)covariance variables are now set. */
test_filter.a_posteriori_step();
/* State and (root_)covariance variables are set to the a priori values.
```

## Building and CMake integration

The public headers require C++17. The commands below use CMake 3.20 or newer;
the project files retain a CMake 3.10 minimum.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target unittest benchmark ahrs sfwa_ukf --parallel 2
ctest --test-dir build --output-on-failure
./build/benchmark/benchmark
```

These commands assume a single-configuration generator. For multi-configuration
generators, add `--config Debug` when building and `-C Debug` when running CTest;
executables are placed in the configuration subdirectory.

The default build downloads Eigen 3.2.8 and applies the allocation and C++17
compatibility patches. To use an installed Eigen 3.3 or newer CMake package,
configure with `-DUKF_USE_SYSTEM_EIGEN=ON`. That package supplies
`Eigen3::Eigen` and is used as installed; the bundled patches are not applied
to it. Tests, benchmarks and examples enable `EIGEN_NO_MALLOC`; Debug retains
Eigen's allocation assertions.
Square-root rank updates use UKF-owned bounded workspace, so these checks also
apply when using an unpatched installed Eigen package.

To include UKF in another CMake project:

```cmake
add_subdirectory(path/to/ukf)
add_executable(my_model src/my_model.cpp)
target_link_libraries(my_model PRIVATE ukf::double)
```

`ukf::double` selects `UKF_DOUBLE_PRECISION`; use `ukf::single` for
`UKF_SINGLE_PRECISION`. Both targets propagate the UKF/Eigen include paths,
C++17 requirement and any Eigen download dependency. If your application already
defines its precision, link `ukf::ukf` (also available as `ukf`) instead: this
target does not define either precision macro. Choose exactly one precision
consistently across interacting translation units; adding the opposite macro
does not override a precision target. An installed `find_package(ukf)` package
is not provided.

| CMake option | Default | Purpose |
| --- | --- | --- |
| `UKF_USE_SYSTEM_EIGEN` | `OFF` | Use an installed Eigen package instead of patched Eigen 3.2.8 |
| `UKF_WARNINGS_AS_ERRORS` | `ON` | Treat warnings in UKF's own targets as errors |
| `UKF_BUILD_TESTS` | `ON` standalone, `OFF` as a subdirectory | Configure the unit tests and `check` target |
| `UKF_BUILD_BENCHMARKS` | `ON` standalone, `OFF` as a subdirectory | Configure the benchmark target |
| `UKF_BUILD_EXAMPLES` | `ON` standalone, `OFF` as a subdirectory | Configure the example libraries and benchmark executables |

Warnings are enabled privately on UKF's own targets, with flags selected for the
compiler. They are not imposed on downstream consumers or third-party
dependencies. Set `UKF_WARNINGS_AS_ERRORS=OFF` to keep warnings without making
them fatal.

CI builds the project and its [downstream consumer fixture](test/cmake/) in
Debug and Release. The fixture checks both precision targets, manual precision,
C++17 propagation, parallel dependency ordering and isolation from parent
project settings.

## Examples

Some examples are provided [here](examples/).

The [recorded-data attitude benchmark](examples/attitude_benchmark/) is the
best place to start: [`attitude_estimator.h`](examples/attitude_benchmark/attitude_estimator.h)
is a complete attitude estimator for a gyroscope, accelerometer and optional
magnetometer, written once for both `UKF::Core` and `UKF::SquareRootCore`. It
is compared with Madgwick, Mahony and VQF on the BROAD and Sassari recordings
in six- and nine-axis modes. See
[methodology and measured results](docs/attitude-benchmark.md)
for local reproduction and CI artifacts.

The [EuRoC MAV benchmark](examples/euroc_benchmark/) extends that estimator
with a rotor-drag model of a Firefly hexrotor, and compares six-axis attitude
estimates with the same references. See its
[model, evaluation protocol and results](docs/euroc-benchmark.md).

The [interactive explainer](explainer/) is an article that explains the
unscented Kalman filter with figures you can change, from a single noisy
sensor up to that EuRoC estimator. The estimator runs in the browser, on the
recorded flights. See [its README](explainer/README.md) for how to build and
publish it.

## References

<a name="ref1">[1]</a> "A New Extension of the Kalman Filter to Nonlinear Systems:,
S. J. Julier and J. K. Uhlmann,
https://www.cs.unc.edu/~welch/kalman/media/pdf/Julier1997_SPIE_KF.pdf

<a name="ref2">[2]</a> "Unscented Filtering for Spacecraft Attitude Estimation", John L.
Crassidis and F. Landis Markley, http://www.acsu.buffalo.edu/~johnc/uf_att.pdf

<a name="ref3">[3]</a> "The Square-Root Unscented Kalman Filter for State and Parameter-Estimation",
Rudolph van der Merwe and Eric A. Wan,
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf
