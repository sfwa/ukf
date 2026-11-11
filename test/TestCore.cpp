#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
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
namespace UKF {
namespace Parameters {
template <> constexpr real_t AlphaSquared<MyStateVector> = 1e-6;
}

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
    temp.set_field<Attitude>(temp_q);

    /* Angular velocity derivative. */
    temp.set_field<AngularVelocity>(angular_acceleration);

    return temp;
}

template <> template <>
MyStateVector MyStateVector::derivative<>() const {
    return derivative(UKF::Vector<3>(0, 0, 0), UKF::Vector<3>(0, 0, 0));
}

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
    UKF::Field<Magnetometer, UKF::FieldVector>,
    UKF::Field<Gyroscope, UKF::Vector<3>>
>;

using MyCore = UKF::Core<
    MyStateVector,
    MyMeasurementVector,
    UKF::IntegratorRK4
>;

namespace UKF {
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
UKF::FieldVector MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * UKF::FieldVector(1, 0, 0);
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
        [[maybe_unused]] const UKF::Vector<3>& acceleration, [[maybe_unused]] const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Position>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, GPS_Velocity>(const MyStateVector& state,
        [[maybe_unused]] const UKF::Vector<3>& acceleration, [[maybe_unused]] const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Velocity>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer, UKF::Vector<3>>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, [[maybe_unused]] const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -9.8) + acceleration;
}

template <> template <>
UKF::FieldVector MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer, UKF::Vector<3>>(const MyStateVector& state,
        [[maybe_unused]] const UKF::Vector<3>& acceleration, [[maybe_unused]] const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Attitude>() * UKF::FieldVector(1, 0, 0);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope, UKF::Vector<3>>(const MyStateVector& state,
        [[maybe_unused]] const UKF::Vector<3>& acceleration, [[maybe_unused]] const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<AngularVelocity>();
}

}

MyCore create_initialised_test_filter() {
    MyCore test_filter;
    test_filter.state.set_field<Position>(UKF::Vector<3>(0, 0, 0));
    test_filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
    test_filter.state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_filter.state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));
    test_filter.covariance = MyStateVector::CovarianceMatrix::Zero();
    test_filter.covariance.diagonal() << 10000, 10000, 10000, 100, 100, 100, 1, 1, 5, 10, 10, 10;
    test_filter.measurement_covariance << 10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;

    real_t a, b;
    real_t dt = 0.01;
    a = std::sqrt(0.1*dt*dt);
    b = std::sqrt(0.1*dt);
    test_filter.process_noise_covariance << a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                            0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b;

    return test_filter;
}

TEST(CoreTest, Initialisation) {
    MyCore test_filter = create_initialised_test_filter();
    ASSERT_EQ(13, test_filter.state.size());
}

/*
All these tests check that the estimated state matches the 'actual' state to
within 2-sigma.
*/
TEST(CoreTest, APrioriStep) {
    MyCore test_filter = create_initialised_test_filter();

    test_filter.a_priori_step(0.01);

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

TEST(CoreTest, APrioriStepWithInputs) {
    MyCore test_filter = create_initialised_test_filter();

    test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

TEST(CoreTest, InnovationStep) {
    MyCore test_filter = create_initialised_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01);
    test_filter.innovation_step(m);

    /*
    With the field vector, we expect the determinant to be approximately zero,
    so allow for it to be slightly negative due to numerical precision.
    */
    EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

TEST(CoreTest, InnovationStepPartialMeasurement) {
    MyCore test_filter = create_initialised_test_filter();
    MyMeasurementVector m;

    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01);
    test_filter.innovation_step(m);

    EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

TEST(CoreTest, InnovationStepWithInputs) {
    MyCore test_filter = create_initialised_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -15));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));
    test_filter.innovation_step(m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

TEST(CoreTest, InnovationStepPartialMeasurementWithInputs) {
    MyCore test_filter = create_initialised_test_filter();
    MyMeasurementVector m;

    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -15));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));
    test_filter.innovation_step(m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GE(test_filter.innovation_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

TEST(CoreTest, APosterioriStep) {
    MyCore test_filter = create_initialised_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01);
    test_filter.innovation_step(m);
    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

TEST(CoreTest, FullStep) {
    MyCore test_filter = create_initialised_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -14.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.step(0.01, m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GT(test_filter.covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(0).norm())*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(3).norm())*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(6).norm())*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        std::sqrt(test_filter.covariance.diagonal().segment<3>(9).norm())*2);
}

using ZeroPriorState = UKF::StateVector<UKF::Field<110, UKF::Vector<2>>>;
using ZeroPriorMeasurement = UKF::FixedMeasurementVector<UKF::Field<111, UKF::Vector<2>>>;
namespace UKF {
template <> template <>
ZeroPriorState ZeroPriorState::process_model<>(real_t) const { return *this; }
}

TEST(CoreTest, ZeroPriorWithPositiveProcessNoise) {
    UKF::Core<ZeroPriorState, ZeroPriorMeasurement> filter;
    filter.state << 2, -3;
    filter.covariance.setZero();
    filter.process_noise_covariance << 2, 0.5, 0.5, 1;
    filter.a_priori_step(real_t(1));
    EXPECT_LT((filter.state - UKF::Vector<2>(2, -3)).norm(), real_t(1e-5));
    EXPECT_LT((filter.covariance - filter.process_noise_covariance).norm(), real_t(1e-5));
}

using AdditiveNoiseState = UKF::StateVector<UKF::Field<150, UKF::Vector<2>>>;
using AdditiveNoiseMeasurement = UKF::FixedMeasurementVector<UKF::Field<151, UKF::Vector<2>>>;
namespace UKF {
template <> template <>
AdditiveNoiseState AdditiveNoiseState::process_model<>(real_t dt) const {
    UKF::Matrix<2, 2> transition;
    transition << 1, dt, 0, 2;
    return transition * *this;
}
template <> template <>
UKF::Vector<2> AdditiveNoiseMeasurement::expected_measurement<AdditiveNoiseState, 151>(const AdditiveNoiseState& state) {
    UKF::Matrix<2, 2> observation;
    observation << 1, 0.5, -0.3, 1;
    return observation * state;
}
}
struct AdditiveNoiseIntegrator {
    template <typename S> static S integrate(real_t dt, const S& state) { return state.process_model(dt); }
};

TEST(AdditiveNoiseTest, MatchesLinearKalmanFilterAndSquareRootFilter) {
    UKF::Core<AdditiveNoiseState, AdditiveNoiseMeasurement> filter;
    UKF::SquareRootCore<AdditiveNoiseState, AdditiveNoiseMeasurement, AdditiveNoiseIntegrator> sr;
    filter.state << 2, -1;
    filter.covariance << 2, 0.4, 0.4, 1;
    filter.process_noise_covariance << 0.5, 0.1, 0.1, 0.3;
    filter.measurement_covariance << 0.7, 0.9;
    sr.state = filter.state;
    sr.root_covariance = filter.covariance.llt().matrixL();
    sr.process_noise_root_covariance = filter.process_noise_covariance.llt().matrixL();
    sr.measurement_root_covariance = filter.measurement_covariance.array().sqrt();
    UKF::Vector<2> mean = filter.state;
    UKF::Matrix<2, 2> covariance = filter.covariance, F, H, R;
    F << 1, 0.5, 0, 2;
    H << 1, 0.5, -0.3, 1;
    R = filter.measurement_covariance.asDiagonal();
    AdditiveNoiseMeasurement z;
    z << 1, -2;
    for (int iteration = 0; iteration < 3; ++iteration) {
        mean = (F * mean).eval();
        covariance = (F * covariance * F.transpose() + filter.process_noise_covariance).eval();
        filter.a_priori_step(real_t(0.5));
        sr.a_priori_step(real_t(0.5));
        EXPECT_LT((filter.state - mean).norm(), real_t(1e-4));
        EXPECT_LT((filter.covariance - covariance).norm(), real_t(1e-4));
        const UKF::Matrix<2, 2> innovation_covariance = H * covariance * H.transpose() + R;
        const UKF::Matrix<2, 2> gain = covariance * H.transpose() * innovation_covariance.inverse();
        mean += gain * (z - H * mean);
        covariance = (covariance - gain * innovation_covariance * gain.transpose()).eval();
        filter.innovation_step(z);
        sr.innovation_step(z);
        EXPECT_LT((filter.innovation_covariance - innovation_covariance).norm(), real_t(1e-4));
        filter.a_posteriori_step();
        sr.a_posteriori_step();
        EXPECT_LT((filter.state - mean).norm(), real_t(1e-4));
        EXPECT_LT((filter.covariance - covariance).norm(), real_t(1e-4));
        EXPECT_LT((sr.state - mean).norm(), real_t(1e-4));
        EXPECT_LT((sr.root_covariance * sr.root_covariance.transpose() - covariance).norm(), real_t(1e-4));
    }
}

TEST(AdditiveNoiseTest, HandlesSemidefinitePriorAndRejectsIndefinitePrior) {
    UKF::Core<AdditiveNoiseState, AdditiveNoiseMeasurement> filter;
    filter.state.setZero();
    filter.covariance << 0, 0, 0, 1;
    filter.process_noise_covariance.setZero();
    filter.a_priori_step(real_t(0));
    UKF::Matrix<2, 2> expected;
    expected << 0, 0, 0, 4;
    EXPECT_LT((filter.covariance - expected).norm(), real_t(1e-5));
    filter.covariance << 1, 2, 2, 1;
    filter.state << 1, 2;
    const AdditiveNoiseState original_state = filter.state;
    const UKF::Matrix<2, 2> original_covariance = filter.covariance;
    EXPECT_FALSE(filter.a_priori_step(real_t(0)));
    EXPECT_EQ(real_t(0), (filter.state - original_state).norm());
    EXPECT_EQ(real_t(0), (filter.covariance - original_covariance).norm());

    /* Steps that depend on the failed prediction do nothing. */
    AdditiveNoiseMeasurement z;
    z << 1, -2;
    EXPECT_FALSE(filter.innovation_step(z));
    EXPECT_FALSE(filter.a_posteriori_step());
    EXPECT_FALSE(filter.step(real_t(0), z));
    EXPECT_EQ(real_t(0), (filter.state - original_state).norm());
    EXPECT_EQ(real_t(0), (filter.covariance - original_covariance).norm());
}

TEST(AdditiveNoiseTest, StepsRequireTheirPrerequisites) {
    UKF::Core<AdditiveNoiseState, AdditiveNoiseMeasurement> filter;
    filter.state << 2, -1;
    filter.covariance << 2, 0.4, 0.4, 1;
    filter.process_noise_covariance.setZero();
    filter.measurement_covariance << 0.7, 0.9;
    AdditiveNoiseMeasurement z;
    z << 1, -2;

    /* A fresh filter has no sigma points to correct with. */
    EXPECT_FALSE(filter.innovation_step(z));
    EXPECT_FALSE(filter.a_posteriori_step());
    EXPECT_TRUE(filter.step(real_t(0.5), z));
}

TEST(AdditiveNoiseTest, RejectsNonFiniteStateWithoutChangingIt) {
    UKF::Core<AdditiveNoiseState, AdditiveNoiseMeasurement> filter;
    filter.state << std::numeric_limits<real_t>::quiet_NaN(), 0;
    filter.covariance.setIdentity();
    filter.process_noise_covariance.setIdentity();
    const UKF::Matrix<2, 2> original_covariance = filter.covariance;
    EXPECT_FALSE(filter.a_priori_step(real_t(0.5)));
    EXPECT_TRUE(std::isnan(filter.state(0)));
    EXPECT_EQ(real_t(0), (filter.covariance - original_covariance).norm());
}

TEST(AdditiveNoiseTest, FactorsCorrelatedSemidefiniteCovariancesWithoutAllocation) {
    /* Columns are deliberately dependent, with different zero-pivot orders. */
    for (int i = 0; i < 20; ++i) {
        UKF::Matrix<3, 3> basis;
        basis << real_t(i)*real_t(0.1), 0, 0,
                 real_t(0.3), real_t(i % 3)*real_t(0.2), 0,
                 real_t(-0.7), real_t(0.6), 0;
        const UKF::Matrix<3, 3> covariance = basis * basis.transpose();
        UKF::Matrix<3, 3> root;
        ASSERT_TRUE(UKF::Detail::covariance_root(covariance, root));
        EXPECT_LT((root * root.transpose() - covariance).norm(), real_t(1e-5));
        EXPECT_EQ(real_t(0), root.template triangularView<Eigen::StrictlyUpper>().toDenseMatrix().norm());
    }
}

TEST(AdditiveNoiseTest, CovarianceRootReportsInvalidInputsWithoutChangingOutput) {
    UKF::Matrix<2, 2> root = UKF::Matrix<2, 2>::Constant(real_t(7));
    UKF::Matrix<2, 2> indefinite;
    indefinite << 1, 2, 2, 1;
    EXPECT_FALSE(UKF::Detail::covariance_root(indefinite, root));
    UKF::Matrix<2, 2> non_finite = UKF::Matrix<2, 2>::Identity();
    non_finite(1, 0) = std::numeric_limits<real_t>::infinity();
    EXPECT_FALSE(UKF::Detail::covariance_root(non_finite, root));
    EXPECT_EQ(real_t(0), (root.array() - real_t(7)).abs().maxCoeff());
}

/*
A small-alpha scalar filter whose quadratic process model makes the unscented
covariance indefinite (negative kappa). Both cores must reject the prediction
without changing the filter; the square-root core previously ignored the failed
downdate and kept a partially downdated factor.
*/
using IndefiniteUTState = UKF::StateVector<UKF::Field<170, real_t>>;
using IndefiniteUTMeasurement = UKF::FixedMeasurementVector<UKF::Field<171, real_t>>;
namespace UKF {
namespace Parameters {
template <> constexpr real_t Kappa<IndefiniteUTState> = -0.5;
}
template <> template <>
IndefiniteUTState IndefiniteUTState::process_model<>(real_t) const {
    IndefiniteUTState next;
    next << (*this)(0) * (*this)(0);
    return next;
}
template <> template <>
real_t IndefiniteUTMeasurement::expected_measurement<IndefiniteUTState, 171>(const IndefiniteUTState& state) {
    return state(0);
}
}

TEST(AdditiveNoiseTest, IndefiniteUnscentedPredictionFailsInBothCores) {
    ASSERT_LT(UKF::Parameters::Sigma_WC0<IndefiniteUTState>, real_t(0));
    UKF::Core<IndefiniteUTState, IndefiniteUTMeasurement> standard;
    standard.state << 0;
    standard.covariance << 1;
    standard.process_noise_covariance << real_t(1e-4);
    EXPECT_FALSE(standard.a_priori_step(real_t(0)));
    EXPECT_EQ(real_t(0), standard.state(0));
    EXPECT_EQ(real_t(1), standard.covariance(0, 0));

    UKF::SquareRootCore<IndefiniteUTState, IndefiniteUTMeasurement> root;
    root.state << 0;
    root.root_covariance << 1;
    root.process_noise_root_covariance << real_t(1e-2);
    EXPECT_FALSE(root.a_priori_step(real_t(0)));
    EXPECT_EQ(real_t(0), root.state(0));
    EXPECT_EQ(real_t(1), root.root_covariance(0, 0));
}

/*
With small alpha the central and outer mean weights are huge and of opposite
sign. The mean must still reproduce a component that is identical in every
sigma point exactly, even at ECEF-like magnitudes.
*/
TEST(CoreTest, SigmaPointMeanIsExactForConstantComponents) {
    ASSERT_LT(UKF::Parameters::Sigma_WM0<MyStateVector>, real_t(-1e5));
    MyStateVector state;
    state.set_field<Position>(UKF::Vector<3>(6378137.123, -4.5e6, 1.0/3.0));
    state.set_field<Velocity>(UKF::Vector<3>(0.1, -7.25, 1e-3));
    state.set_field<Attitude>(UKF::Quaternion(0.5, -0.5, 0.5, 0.5));
    state.set_field<AngularVelocity>(UKF::Vector<3>(0.3, 0.2, 0.1));
    MyStateVector::SigmaPointDistribution sigma_points;
    for (std::size_t i = 0; i < MyStateVector::num_sigma(); ++i) {
        sigma_points.col(i) = state;
    }
    const MyStateVector mean = MyStateVector::calculate_sigma_point_mean(sigma_points);
    const UKF::Vector<3> mean_position = mean.get_field<Position>();
    const UKF::Vector<3> mean_velocity = mean.get_field<Velocity>();
    const UKF::Vector<3> mean_angular_velocity = mean.get_field<AngularVelocity>();
    EXPECT_EQ(UKF::Vector<3>(6378137.123, -4.5e6, 1.0/3.0), mean_position);
    EXPECT_EQ(UKF::Vector<3>(0.1, -7.25, 1e-3), mean_velocity);
    EXPECT_EQ(UKF::Vector<3>(0.3, 0.2, 0.1), mean_angular_velocity);
    EXPECT_LT((state.get_field<Attitude>().coeffs() - mean.get_field<Attitude>().coeffs()).norm(),
        real_t(4) * std::numeric_limits<real_t>::epsilon());
}

TEST(CoreTest, ScaledDimensionAvoidsCancellation) {
    constexpr real_t expected = real_t(1e-6) * real_t(MyStateVector::covariance_size() + 3);
    EXPECT_EQ(expected, UKF::Parameters::ScaledDimension<MyStateVector>);
    EXPECT_NEAR(real_t(1), UKF::Parameters::Sigma_WM0<MyStateVector> +
        real_t(2 * MyStateVector::covariance_size()) * UKF::Parameters::Sigma_WMI<MyStateVector>,
        real_t(1e-6));
}
