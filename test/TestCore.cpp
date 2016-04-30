#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Types.h"
#include "Integrator.h"
#include "StateVector.h"
#include "MeasurementVector.h"
#include "Core.h"
#include "comparisons.h"

/* Set up state vector class. */
enum MyStateFields {
    Position,
    Velocity,
    Attitude,
    AngularVelocity
};

using MyStateVector = UKF::StateVector<
    UKF::Field<Position, UKF::Vector<3>>,
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>,
    UKF::Field<AngularVelocity, UKF::Vector<3>>
>;

/*
State vector process model. One version takes body frame kinematic
acceleration and angular acceleration as inputs, the other doesn't (assumes
zero accelerations).
*/
template <> template <>
MyStateVector MyStateVector::derivative<UKF::Vector<3>, UKF::Vector<3>>(
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) const {
    MyStateVector temp;
    
    /* Position derivative. */
    temp.set_field<Position>(get_field<Velocity>());

    /* Velocity derivative. */
    temp.set_field<Velocity>(get_field<Attitude>().conjugate() * acceleration);

    /* Attitude derivative. */
    UKF::Quaternion temp_q;
    temp_q.vec() = get_field<AngularVelocity>();
    temp_q.w() = 0;

    /* Angular velocity derivative. */
    temp.set_field<AngularVelocity>(angular_acceleration);

    return temp;
}

template <> template <>
MyStateVector MyStateVector::derivative<>() const {
    return derivative(UKF::Vector<3>(0, 0, 0), UKF::Vector<3>(0, 0, 0));
}

/*
State vector process noise covariance. These are just completely arbitrary
for this test.
*/
template <>
MyStateVector::CovarianceMatrix MyStateVector::process_noise_covariance(real_t dt) {
    MyStateVector::CovarianceMatrix temp;
    temp << 0.1*dt*dt,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                    0, 0.1*dt*dt,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                    0,         0, 0.1*dt*dt,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                    0,         0,         0,    0.1*dt,         0,         0,         0,         0,         0,         0,         0,         0,
                    0,         0,         0,         0,    0.1*dt,         0,         0,         0,         0,         0,         0,         0,
                    0,         0,         0,         0,         0,    0.1*dt,         0,         0,         0,         0,         0,         0,
                    0,         0,         0,         0,         0,         0, 0.1*dt*dt,         0,         0,         0,         0,         0,
                    0,         0,         0,         0,         0,         0,         0, 0.1*dt*dt,         0,         0,         0,         0,
                    0,         0,         0,         0,         0,         0,         0,         0, 0.1*dt*dt,         0,         0,         0,
                    0,         0,         0,         0,         0,         0,         0,         0,         0,    0.1*dt,         0,         0,
                    0,         0,         0,         0,         0,         0,         0,         0,         0,         0,    0.1*dt,         0,
                    0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,    0.1*dt;
    return temp;
}

/* Set up measurement vector class. */
enum MyMeasurementFields {
    GPS_Position,
    GPS_Velocity,
    Accelerometer,
    Magnetometer,
    Gyroscope
};

using MyMeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<GPS_Position, UKF::Vector<3>>,
    UKF::Field<GPS_Velocity, UKF::Vector<3>>,
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Magnetometer, UKF::Vector<3>>,
    UKF::Field<Gyroscope, UKF::Vector<3>>
>;

using MyCore = UKF::Core<
    MyStateVector,
    MyMeasurementVector,
    UKF::IntegratorRK4
>;

/*
Define measurement model to be used in tests. NOTE: These are just for
testing, don't expect them to make any physical sense whatsoever.
*/
template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, GPS_Position>(const MyStateVector& state) {
    return state.get_field<Position>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, GPS_Velocity>(const MyStateVector& state) {
    return state.get_field<Velocity>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -9.8);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * UKF::Vector<3>(1, 0, 0);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope>(const MyStateVector& state) {
    return state.get_field<AngularVelocity>();
}

/*
These versions of the predicted measurement functions take kinematic
acceleration and angular acceleration as inputs. Note that in reality, the
inputs would probably be a control vector and the accelerations would be
calculated using the state vector and a dynamics model.
*/
template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, GPS_Position>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Position>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, GPS_Velocity>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Velocity>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer, UKF::Vector<3>>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -9.8) + acceleration;
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer, UKF::Vector<3>>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Attitude>() * UKF::Vector<3>(1, 0, 0);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope, UKF::Vector<3>>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<AngularVelocity>();
}

/* Set the measurement covariance vector. */
template <>
MyMeasurementVector::CovarianceVector MyMeasurementVector::measurement_covariance =
    MyMeasurementVector::CovarianceVector();

TEST(CoreTest, Initialisation) {
    
}
