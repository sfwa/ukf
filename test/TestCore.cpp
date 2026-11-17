#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include <type_traits>
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
/*
A small alpha, as much of the literature uses. Single precision cannot carry
alpha^2 = 1e-6 (see AlphaSquared), so its build uses 1e-2, whose central
weight is still negative.
*/
template <> constexpr real_t AlphaSquared<MyStateVector> =
    sizeof(real_t) == sizeof(double) ? real_t(1e-6) : real_t(1e-2);
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

    ASSERT_TRUE(test_filter.a_priori_step(0.01));

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

    ASSERT_TRUE(test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0)));

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

    ASSERT_TRUE(test_filter.a_priori_step(0.01));
    ASSERT_TRUE(test_filter.innovation_step(m));

    /*
    With the field vector, the innovation covariance is singular, so allow
    for an eigenvalue that is slightly negative due to numerical precision.
    */
    EXPECT_POSITIVE_SEMIDEFINITE(test_filter.innovation_covariance);

    ASSERT_TRUE(test_filter.a_posteriori_step());

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

    ASSERT_TRUE(test_filter.a_priori_step(0.01));
    ASSERT_TRUE(test_filter.innovation_step(m));

    EXPECT_POSITIVE_SEMIDEFINITE(test_filter.innovation_covariance);

    ASSERT_TRUE(test_filter.a_posteriori_step());

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

    ASSERT_TRUE(test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0)));
    ASSERT_TRUE(test_filter.innovation_step(m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0)));

    EXPECT_POSITIVE_SEMIDEFINITE(test_filter.innovation_covariance);

    ASSERT_TRUE(test_filter.a_posteriori_step());

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

    ASSERT_TRUE(test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0)));
    ASSERT_TRUE(test_filter.innovation_step(m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0)));

    EXPECT_POSITIVE_SEMIDEFINITE(test_filter.innovation_covariance);

    ASSERT_TRUE(test_filter.a_posteriori_step());

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

    ASSERT_TRUE(test_filter.a_priori_step(0.01));
    ASSERT_TRUE(test_filter.innovation_step(m));
    ASSERT_TRUE(test_filter.a_posteriori_step());

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

    ASSERT_TRUE(test_filter.step(0.01, m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0)));

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
        EXPECT_GE(root.diagonal().minCoeff(), real_t(0));
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
    ASSERT_LT(UKF::Parameters::Sigma_WM0<MyStateVector>, precision_tolerance(-1e5, -50));
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
    constexpr real_t expected = UKF::Parameters::AlphaSquared<MyStateVector> *
        real_t(MyStateVector::covariance_size() + 3);
    EXPECT_EQ(expected, UKF::Parameters::ScaledDimension<MyStateVector>);
    /*
    The weights sum to one to within the rounding of the central weight,
    which is about -7.5e5 in double precision here. The means are
    accumulated as offsets, so they do not depend on that sum.
    */
    EXPECT_NEAR(real_t(1), UKF::Parameters::Sigma_WM0<MyStateVector> +
        real_t(2 * MyStateVector::covariance_size()) * UKF::Parameters::Sigma_WMI<MyStateVector>,
        real_t(4) * std::numeric_limits<real_t>::epsilon() * std::abs(UKF::Parameters::Sigma_WM0<MyStateVector>));
}

/*
An identity process with a quadratic measurement and negative kappa: the
prediction is valid, but the unscented innovation covariance is indefinite.
*/
using IndefiniteInnovationState = UKF::StateVector<UKF::Field<172, real_t>>;
using IndefiniteInnovationMeasurement = UKF::FixedMeasurementVector<UKF::Field<173, real_t>>;
namespace UKF {
namespace Parameters {
template <> constexpr real_t Kappa<IndefiniteInnovationState> = -0.5;
}
template <> template <>
IndefiniteInnovationState IndefiniteInnovationState::process_model<>(real_t) const {
    return *this;
}
template <> template <>
real_t IndefiniteInnovationMeasurement::expected_measurement<IndefiniteInnovationState, 173>(
        const IndefiniteInnovationState& state) {
    return state(0) * state(0);
}
}

TEST(AdditiveNoiseTest, IndefiniteInnovationCovarianceFailsWithoutChangingFilter) {
    UKF::Core<IndefiniteInnovationState, IndefiniteInnovationMeasurement> filter;
    filter.state << 0;
    filter.covariance << 1;
    filter.process_noise_covariance << 0;
    filter.measurement_covariance << real_t(1e-4);
    IndefiniteInnovationMeasurement z;
    z << real_t(0.5);

    ASSERT_TRUE(filter.a_priori_step(real_t(0)));
    const real_t predicted_state = filter.state(0);
    const real_t predicted_covariance = filter.covariance(0, 0);
    EXPECT_NEAR(real_t(1), predicted_covariance, real_t(16) * std::numeric_limits<real_t>::epsilon());
    EXPECT_FALSE(filter.innovation_step(z));
    EXPECT_FALSE(filter.a_posteriori_step());
    EXPECT_EQ(real_t(0), filter.innovation(0));
    EXPECT_EQ(real_t(0), filter.innovation_covariance(0, 0));
    EXPECT_EQ(predicted_state, filter.state(0));
    EXPECT_EQ(predicted_covariance, filter.covariance(0, 0));

    /* As one step, the prediction is rolled back as well. */
    filter.state << 0;
    filter.covariance << 1;
    EXPECT_FALSE(filter.step(real_t(0), z));
    EXPECT_EQ(real_t(0), filter.state(0));
    EXPECT_EQ(real_t(1), filter.covariance(0, 0));
}

/* step() rolls back a successful prediction when its correction fails. */
TEST(AdditiveNoiseTest, StepIsAllOrNothing) {
    UKF::Core<AdditiveNoiseState, AdditiveNoiseMeasurement> filter;
    filter.state << 2, -1;
    filter.covariance << 2, 0.4, 0.4, 1;
    filter.process_noise_covariance << 0.5, 0.1, 0.1, 0.3;
    filter.measurement_covariance << 0.7, 0.9;
    AdditiveNoiseMeasurement z;
    z << 1, -2;
    ASSERT_TRUE(filter.step(real_t(0.5), z));

    const AdditiveNoiseState state = filter.state;
    const UKF::Matrix<2, 2> covariance = filter.covariance;
    const AdditiveNoiseMeasurement innovation = filter.innovation;
    const UKF::Matrix<2, 2> innovation_covariance = filter.innovation_covariance;
    filter.measurement_covariance(1) = std::numeric_limits<real_t>::quiet_NaN();
    EXPECT_FALSE(filter.step(real_t(0.5), z));
    EXPECT_EQ(real_t(0), (filter.state - state).norm());
    EXPECT_EQ(real_t(0), (filter.covariance - covariance).norm());
    EXPECT_EQ(real_t(0), (filter.innovation - innovation).norm());
    EXPECT_EQ(real_t(0), (filter.innovation_covariance - innovation_covariance).norm());

    /* The rolled-back prediction cannot be corrected without predicting again. */
    filter.measurement_covariance << 0.7, 0.9;
    EXPECT_FALSE(filter.innovation_step(z));
    EXPECT_TRUE(filter.step(real_t(0.5), z));
}

/*
An exact observation of correlated states leaves a posterior which is zero
only to within roundoff, with off-diagonal entries of either sign. The
correction is accepted against the prior's scale, and the covariance stored
must be the repaired one, or the next prediction, which judges it against
its own scale, fails.
*/
using ExactObservationState = UKF::StateVector<UKF::Field<218, real_t>, UKF::Field<219, real_t>>;
using ExactObservation = UKF::FixedMeasurementVector<UKF::Field<220, real_t>, UKF::Field<221, real_t>>;
namespace UKF {
template <> template <>
ExactObservationState ExactObservationState::process_model<>(real_t) const {
    return *this;
}
template <> template <>
real_t ExactObservation::expected_measurement<ExactObservationState, 220>(const ExactObservationState& state) {
    return state(0);
}
template <> template <>
real_t ExactObservation::expected_measurement<ExactObservationState, 221>(const ExactObservationState& state) {
    return state(1);
}
}

TEST(CoreTest, ExactCorrelatedObservationKeepsPredicting) {
    UKF::Core<ExactObservationState, ExactObservation> filter;
    filter.state << 1, 2;
    filter.covariance << 0.25, 0.05, 0.05, 1.01;
    filter.process_noise_covariance.setZero();
    filter.measurement_covariance << 0, 0;
    ExactObservation z;
    z << 1.1, 1.9;
    for (int i = 0; i < 4; ++i) {
        ASSERT_TRUE(filter.step(real_t(0), z)) << i;
        EXPECT_NEAR(real_t(1.1), filter.state(0), precision_tolerance(1e-12, 1e-6));
        EXPECT_NEAR(real_t(1.9), filter.state(1), precision_tolerance(1e-12, 1e-6));
        EXPECT_LT(filter.covariance.cwiseAbs().maxCoeff(), precision_tolerance(1e-12, 1e-5));
        EXPECT_POSITIVE_SEMIDEFINITE(filter.covariance);
    }
    EXPECT_TRUE(filter.a_priori_step(real_t(0)));
}

/*
A process model which keeps w >= 0, as some attitude integrators do, flips
the sigma points that cross 180 degrees. The filters must estimate what they
estimate with a model whose quaternions are continuous.
*/
using ContinuousTurnState = UKF::StateVector<UKF::Field<212, UKF::Quaternion>, UKF::Field<213, UKF::Vector<3>>>;
using CanonicalTurnState = UKF::StateVector<UKF::Field<214, UKF::Quaternion>, UKF::Field<215, UKF::Vector<3>>>;
using TurnMeasurement = UKF::FixedMeasurementVector<UKF::Field<216, UKF::Vector<3>>>;

template <typename S, int Attitude, int Rate>
S turn(const S& state, real_t dt, bool canonical) {
    const UKF::Vector<3> rate = state.template get_field<Rate>();
    UKF::Quaternion q = UKF::Quaternion(Eigen::AngleAxis<real_t>(rate.norm() * dt,
        rate.norm() > 0 ? UKF::Vector<3>(rate.normalized()) : UKF::Vector<3>::UnitZ())) *
        state.template get_field<Attitude>();
    if (canonical && q.w() < 0) {
        q.coeffs() = -q.coeffs();
    }
    S next = state;
    next.template set_field<Attitude>(q);
    return next;
}

namespace UKF {
template <> template <>
ContinuousTurnState ContinuousTurnState::process_model<>(real_t dt) const {
    return turn<ContinuousTurnState, 212, 213>(*this, dt, false);
}
template <> template <>
CanonicalTurnState CanonicalTurnState::process_model<>(real_t dt) const {
    return turn<CanonicalTurnState, 214, 215>(*this, dt, true);
}
template <> template <>
UKF::Vector<3> TurnMeasurement::expected_measurement<ContinuousTurnState, 216>(const ContinuousTurnState& state) {
    return state.get_field<212>() * UKF::Vector<3>(1, 0, 0);
}
template <> template <>
UKF::Vector<3> TurnMeasurement::expected_measurement<CanonicalTurnState, 216>(const CanonicalTurnState& state) {
    return state.get_field<214>() * UKF::Vector<3>(1, 0, 0);
}
}

template <template <typename...> class Filter>
void check_canonical_process_model() {
    Filter<ContinuousTurnState, TurnMeasurement, UKF::NotUsedIntegrator> continuous;
    Filter<CanonicalTurnState, TurnMeasurement, UKF::NotUsedIntegrator> canonical;
    const UKF::Quaternion start(Eigen::AngleAxis<real_t>(real_t(2.9), UKF::Vector<3>::UnitZ()));
    continuous.state.template set_field<212>(start);
    continuous.state.template set_field<213>(UKF::Vector<3>(0, 0, 0.5));
    canonical.state << continuous.state;
    if constexpr (std::is_same_v<Filter<ContinuousTurnState, TurnMeasurement, UKF::NotUsedIntegrator>,
            UKF::Core<ContinuousTurnState, TurnMeasurement, UKF::NotUsedIntegrator>>) {
        continuous.covariance = (UKF::Vector<6>() << 0.09, 0.09, 0.09, 0.01, 0.01, 0.04).finished().asDiagonal();
        continuous.process_noise_covariance = UKF::Matrix<6, 6>::Identity() * real_t(1e-4);
        continuous.measurement_covariance << 0.01, 0.01, 0.01;
        canonical.covariance = continuous.covariance;
        canonical.process_noise_covariance = continuous.process_noise_covariance;
        canonical.measurement_covariance = continuous.measurement_covariance;
    } else {
        continuous.root_covariance = (UKF::Vector<6>() << 0.3, 0.3, 0.3, 0.1, 0.1, 0.2).finished().asDiagonal();
        continuous.process_noise_root_covariance = UKF::Matrix<6, 6>::Identity() * real_t(1e-2);
        continuous.measurement_root_covariance << 0.1, 0.1, 0.1;
        canonical.root_covariance = continuous.root_covariance;
        canonical.process_noise_root_covariance = continuous.process_noise_root_covariance;
        canonical.measurement_root_covariance = continuous.measurement_root_covariance;
    }

    /* The truth turns through 180 degrees at 0.4 rad/s. */
    for (int i = 0; i < 20; ++i) {
        const UKF::Quaternion truth(Eigen::AngleAxis<real_t>(real_t(3.0) + real_t(0.04) * real_t(i + 1),
            UKF::Vector<3>::UnitZ()));
        TurnMeasurement z;
        z << truth * UKF::Vector<3>(1, 0, 0);
        ASSERT_TRUE(continuous.step(real_t(0.1), z)) << i;
        ASSERT_TRUE(canonical.step(real_t(0.1), z)) << i;
        const UKF::Quaternion a = continuous.state.template get_field<212>();
        const UKF::Quaternion b = canonical.state.template get_field<214>();
        EXPECT_LT(real_t(2) * (a * b.conjugate()).vec().norm(), precision_tolerance(1e-12, 1e-5)) << i;
        EXPECT_LT((continuous.state.template get_field<213>() - canonical.state.template get_field<215>()).norm(),
            precision_tolerance(1e-12, 1e-5)) << i;
    }
}

TEST(CoreTest, CanonicalQuaternionProcessModel) {
    check_canonical_process_model<UKF::Core>();
}

TEST(CoreTest, CanonicalQuaternionProcessModelSquareRoot) {
    check_canonical_process_model<UKF::SquareRootCore>();
}

/*
The sigma points of the last a priori step are readable through a const
reference: the central point is the predicted state, and the outer points
reproduce the predicted covariance, process noise included. An a posteriori
step leaves them as the points the measurement update used.
*/
template <typename Filter>
void check_sigma_points(Filter& filter, const UKF::Matrix<2, 2>& covariance, bool after_update) {
    using Points = typename AdditiveNoiseState::SigmaPointDistribution;
    static_assert(std::is_same<decltype(filter.get_sigma_points()), const Points&>::value,
        "Sigma points must be exposed by const reference");
    const Points& points = filter.get_sigma_points();
    const UKF::Vector<2> centre = points.col(0);
    if (!after_update) {
        EXPECT_LT((centre - filter.state).norm(), real_t(1e-5));
    }
    UKF::Matrix<2, 2> spread = UKF::Matrix<2, 2>::Zero();
    for (std::size_t i = 1; i < AdditiveNoiseState::num_sigma(); ++i) {
        spread += (points.col(i) - centre) * (points.col(i) - centre).transpose();
    }
    spread /= real_t(2) * UKF::Parameters::ScaledDimension<AdditiveNoiseState>;
    EXPECT_LT((spread - covariance).norm(), real_t(1e-4));
}

TEST(SigmaPointAccessTest, AllCoresExposeTheirLastSigmaPoints) {
    UKF::Core<AdditiveNoiseState, AdditiveNoiseMeasurement> filter;
    UKF::SquareRootCore<AdditiveNoiseState, AdditiveNoiseMeasurement, AdditiveNoiseIntegrator> sr;
    UKF::SquareRootParameterEstimationCore<AdditiveNoiseState, AdditiveNoiseMeasurement> parameter;
    EXPECT_TRUE(filter.get_sigma_points().isZero());
    EXPECT_TRUE(sr.get_sigma_points().isZero());
    EXPECT_TRUE(parameter.get_sigma_points().isZero());

    filter.state << 2, -1;
    filter.covariance << 2, 0.4, 0.4, 1;
    filter.process_noise_covariance << 0.5, 0.1, 0.1, 0.3;
    filter.measurement_covariance << 0.7, 0.9;
    sr.state = filter.state;
    sr.root_covariance = filter.covariance.llt().matrixL();
    sr.process_noise_root_covariance = filter.process_noise_covariance.llt().matrixL();
    sr.measurement_root_covariance = filter.measurement_covariance.array().sqrt();
    parameter.state = filter.state;
    parameter.root_covariance = sr.root_covariance;
    parameter.process_noise_root_covariance = sr.process_noise_root_covariance;
    parameter.measurement_root_covariance = sr.measurement_root_covariance;
    AdditiveNoiseMeasurement z;
    z << 1, -2;

    ASSERT_TRUE(filter.a_priori_step(real_t(0.5)));
    ASSERT_TRUE(sr.a_priori_step(real_t(0.5)));
    ASSERT_TRUE(parameter.a_priori_step());
    const UKF::Matrix<2, 2> predicted = filter.covariance;
    const UKF::Matrix<2, 2> sr_predicted = sr.root_covariance * sr.root_covariance.transpose();
    const UKF::Matrix<2, 2> parameter_predicted = parameter.root_covariance * parameter.root_covariance.transpose();
    check_sigma_points(filter, predicted, false);
    check_sigma_points(sr, sr_predicted, false);
    check_sigma_points(parameter, parameter_predicted, false);

    const auto before = filter.get_sigma_points();
    ASSERT_TRUE(filter.innovation_step(z) && filter.a_posteriori_step());
    ASSERT_TRUE(sr.innovation_step(z) && sr.a_posteriori_step());
    ASSERT_TRUE(parameter.innovation_step(z) && parameter.a_posteriori_step());
    EXPECT_TRUE(before == filter.get_sigma_points());
    check_sigma_points(filter, predicted, true);
    check_sigma_points(sr, sr_predicted, true);
    check_sigma_points(parameter, parameter_predicted, true);
}
