# UKF

Unscented Kalman filter library. Several different UKF implementations are
provided:

  * Standard Unscented Kalman Filter for state estimation, as originally
    described in [^1], with extensions for quaternions as described in [^2].
  * Square-root Unscented Kalman Filter for state estimation, implemented as
    described in [^3].
  * Optimised form of square-root Unscented Kalman filter for parameter
    estimation, implemented as described in [^3].

This library makes use of the Eigen library for linear algebra routines and
matrix and vector operations. Heavy use is made of C++11 and C++14 features
in an attempt to avoid any dynamic memory allocations and maximise
opportunities for compile-time optimisations.

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

For a description of what the scaling parameters do, see [^2] or read the
comments in the [code](include/UKF/StateVector.h).

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
[Integrator.h](include/UKF/Integrator.h).

The square-root state estimation filter can be used instead like this:

```C++
using MyCore = UKF::SquareRootCore<
    MyStateVector,
    MyMeasurementVector,
    UKF::IntegratorRK4
>;
```

Specialisations of the process noise covariance (or process noise root
covariance for the SR-UKF), measurement noise covariance (or measurement
noise root covariance for the SR-UKF), process model (for state estimation
filters) and measurement model must also be provided.

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

### Process noise covariance

The process noise covariance function should return a covariance matrix,
which can be dependent upon the time step size if desired:

```C++
template <>
ProcessModelTestStateVector::CovarianceMatrix ProcessModelTestStateVector::process_noise_covariance(real_t delta) {
    ProcessModelTestStateVector::CovarianceMatrix temp;
    temp << 0.1*delta*delta,               0,               0,         0,         0,         0,
                          0, 0.1*delta*delta,               0,         0,         0,         0,
                          0,               0, 0.1*delta*delta,         0,         0,         0,
                          0,               0,               0, 0.1*delta,         0,         0,
                          0,               0,               0,         0, 0.1*delta,         0,
                          0,               0,               0,         0,         0, 0.1*delta;
    return temp;
}
```

For the SR-UKF, the process noise root covariance should be provided instead,
which could be something like:

```C++
template <>
ProcessModelTestStateVector::CovarianceMatrix ProcessModelTestStateVector::process_noise_root_covariance(real_t delta) {
    ProcessModelTestStateVector::CovarianceMatrix temp;
    real_t a = std::sqrt(0.1*delta*delta);
    real_t b = std::sqrt(0.1*delta);
    temp << a, 0, 0, 0, 0, 0,
            0, a, 0, 0, 0, 0,
            0, 0, a, 0, 0, 0,
            0, 0, 0, b, 0, 0,
            0, 0, 0, 0, b, 0,
            0, 0, 0, 0, 0, b;
    return temp;
}
```

### Measurement model

The measurement model is specified per field, in order to allow the expected
measurement vector to be constructed for the dynamic measurement vector,
where not all measurements may be available each iteration. Each measurement
model function takes a state vector as an input.

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

```
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

### Measurement noise covariance

The measurement noise covariance is a static member variable of the
specialised measurement vector class, and represents the diagonal of the
measurement noise covariance:

```C++
template <>
MyMeasurementVector::CovarianceVector MyMeasurementVector::measurement_covariance(
    (MyMeasurementVector::CovarianceVector() <<
    	10*10,
    	10*10,
        0.5*0.5 * UKF::Vector<3>::Ones(),
        0.004*0.004 * UKF::Vector<3>::Ones()).finished());
```

For the SR-UKF, the measurement noise root covariance should be provided
instead, which could be something like:

```C++
template <>
MyMeasurementVector::CovarianceVector MyMeasurementVector::measurement_root_covariance(
    (MyMeasurementVector::CovarianceVector() <<
    	10,
    	10,
        0.5 * UKF::Vector<3>::Ones(),
        0.004 * UKF::Vector<3>::Ones()).finished());
```

### Initialisation

The filter state and covariance should be initialised to appropriate values,
e.g.:

```C++
MyCore test_filter;
test_filter.state.set_field<Position>(UKF::Vector<3>(0, 0, 0));
test_filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
test_filter.state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
test_filter.state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));
test_filter.covariance = MyStateVector::CovarianceMatrix::Zero();
test_filter.covariance.diagonal() << 10000, 10000, 10000, 100, 100, 100, 1, 1, 5, 10, 10, 10;
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
```

### Iteration

The general steps for carrying out a filter iteration are something like:

```C++
MyMeasurementVector meas;
meas.set_field<Accelerometer>(UKF::Vector<3>(0.0, 0.0, 9.8));
meas.set_field<Gyroscope>(UKF::Vector<3>(0.0, 0.0, 0.0));
meas.set_field<StaticPressure>(101300.0);
meas.set_field<DynamicPressure>(101300.0);

test_filter.step(0.01, meas);
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

## Examples

Some examples are provided [here](examples/).

[^1]: "A New Extension of the Kalman Filter to Nonlinear Systems:,
S. J. Julier and J. K. Uhlmann,
https://www.cs.unc.edu/~welch/kalman/media/pdf/Julier1997_SPIE_KF.pdf
[^2]: "Unscented Filtering for Spacecraft Attitude Estimation", John L.
Crassidis and F. Landis Markley, http://www.acsu.buffalo.edu/~johnc/uf_att.pdf
[^3]: "The Square-Root Unscented Kalman Filter for State and Parameter-Estimation",
Rudolph van der Merwe and Eric A. Wan,
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf