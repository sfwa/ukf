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

/*
Set up state vector class. The order of these is changed to prevent
linker collisions with the ones in TestCore.cpp.
*/
enum MyStateFields {
    Attitude,
    AngularVelocity,
    Position,
    Velocity
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

/*
Set up measurement vector class. The order of these is changed to prevent
linker collisions with the ones in TestCore.cpp.
*/
enum MyMeasurementFields {
    Accelerometer,
    Magnetometer,
    Gyroscope,
    GPS_Position,
    GPS_Velocity
};

using MyMeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<GPS_Position, UKF::Vector<3>>,
    UKF::Field<GPS_Velocity, UKF::Vector<3>>,
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Magnetometer, UKF::FieldVector>,
    UKF::Field<Gyroscope, UKF::Vector<3>>
>;

using MyCore = UKF::SquareRootCore<
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

/*
Initialise covariances as square roots to be comparable with TestCore.cpp.
*/
MyCore create_initialised_sr_test_filter() {
    MyCore test_filter;
    test_filter.state.set_field<Position>(UKF::Vector<3>(0, 0, 0));
    test_filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
    test_filter.state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_filter.state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));
    test_filter.root_covariance = MyStateVector::CovarianceMatrix::Zero();
    test_filter.root_covariance.diagonal() <<
        10000, 10000, 10000, 100, 100, 100, 1, 1, 5, 10, 10, 10;
    test_filter.root_covariance = test_filter.root_covariance.llt().matrixU();
    test_filter.measurement_root_covariance << 
        10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;
    test_filter.measurement_root_covariance = test_filter.measurement_root_covariance.array().sqrt();

    real_t a, b;
    real_t dt = 0.01;
    a = std::sqrt(0.1*dt*dt);
    b = std::sqrt(0.1*dt);
    test_filter.process_noise_root_covariance << a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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
    test_filter.process_noise_root_covariance = test_filter.process_noise_root_covariance.llt().matrixU();

    return test_filter;
}

using ParameterNoiseState = UKF::StateVector<UKF::Field<120, UKF::Vector<2>>>;
using ParameterNoiseMeasurement = UKF::FixedMeasurementVector<UKF::Field<121, UKF::Vector<2>>>;
using ParameterNoiseFilter = UKF::SquareRootParameterEstimationCore<ParameterNoiseState, ParameterNoiseMeasurement>;

TEST(ParameterNoiseTest, CorrelatedPriorAndDiagonalNoise) {
    ParameterNoiseFilter filter;
    filter.state << 2, -3;
    filter.root_covariance << 1, 0, 1, 1;
    filter.process_noise_root_covariance.setIdentity();
    for (int i = 1; i <= 3; ++i) {
        filter.a_priori_step();
        UKF::Matrix<2, 2> expected;
        expected << 1+i, 1, 1, 2+i;
        EXPECT_LT((filter.root_covariance * filter.root_covariance.transpose() - expected).norm(), real_t(1e-5));
        EXPECT_LT((filter.state - UKF::Vector<2>(2, -3)).norm(), real_t(1e-5));
    }
}

TEST(ParameterNoiseTest, CorrelatedNoiseAndZeroNoise) {
    ParameterNoiseFilter filter;
    filter.state.setZero();
    filter.root_covariance << -1, 0, 0.5, -2;
    filter.process_noise_root_covariance << 0.2, 0, 0.3, 0.4;
    UKF::Matrix<2, 2> expected = filter.root_covariance * filter.root_covariance.transpose() +
        filter.process_noise_root_covariance * filter.process_noise_root_covariance.transpose();
    filter.a_priori_step();
    EXPECT_LT((filter.root_covariance * filter.root_covariance.transpose() - expected).norm(), real_t(1e-5));
    filter.process_noise_root_covariance.setZero();
    filter.a_priori_step();
    EXPECT_LT((filter.root_covariance * filter.root_covariance.transpose() - expected).norm(), real_t(1e-5));
}

TEST(SquareRootCoreTest, Initialisation) {
    MyCore test_filter = create_initialised_sr_test_filter();
    ASSERT_EQ(13, test_filter.state.size());
}

/*
All these tests check that the estimated state matches the 'actual' state to
within 2-sigma.
*/
TEST(SquareRootCoreTest, APrioriStep) {
    MyCore test_filter = create_initialised_sr_test_filter();

    test_filter.a_priori_step(0.01);

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

TEST(SquareRootCoreTest, APrioriStepWithInputs) {
    MyCore test_filter = create_initialised_sr_test_filter();

    test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

TEST(SquareRootCoreTest, InnovationStep) {
    MyCore test_filter = create_initialised_sr_test_filter();
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
    EXPECT_GE(test_filter.innovation_root_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

TEST(SquareRootCoreTest, InnovationStepPartialMeasurement) {
    MyCore test_filter = create_initialised_sr_test_filter();
    MyMeasurementVector m;

    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01);
    test_filter.innovation_step(m);

    EXPECT_GE(test_filter.innovation_root_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

TEST(SquareRootCoreTest, InnovationStepWithInputs) {
    MyCore test_filter = create_initialised_sr_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -15));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));
    test_filter.innovation_step(m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GE(test_filter.innovation_root_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

TEST(SquareRootCoreTest, InnovationStepPartialMeasurementWithInputs) {
    MyCore test_filter = create_initialised_sr_test_filter();
    MyMeasurementVector m;

    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -15));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));
    test_filter.innovation_step(m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GE(test_filter.innovation_root_covariance.determinant(), -std::numeric_limits<real_t>::epsilon());

    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

TEST(SquareRootCoreTest, APosterioriStep) {
    MyCore test_filter = create_initialised_sr_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01);
    test_filter.innovation_step(m);
    test_filter.a_posteriori_step();

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

TEST(SquareRootCoreTest, FullStep) {
    MyCore test_filter = create_initialised_sr_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -14.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.step(0.01, m, UKF::Vector<3>(0, 0, -5), UKF::Vector<3>(1, 0, 0));

    EXPECT_GT(test_filter.root_covariance.determinant(), std::numeric_limits<real_t>::epsilon());
    EXPECT_LT((UKF::Vector<3>(100, 10, -50) - test_filter.state.get_field<Position>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(0).norm()*2);
    EXPECT_LT((UKF::Vector<3>(20, 0, 0) - test_filter.state.get_field<Velocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(3).norm()*2);
    EXPECT_LT(2*std::acos(std::abs(UKF::Quaternion(0.7071, 0, 0, -0.7071).dot(test_filter.state.get_field<Attitude>()))),
        test_filter.root_covariance.diagonal().segment<3>(6).norm()*2);
    EXPECT_LT((UKF::Vector<3>(0.5, 0, 0) - test_filter.state.get_field<AngularVelocity>()).norm(),
        test_filter.root_covariance.diagonal().segment<3>(9).norm()*2);
}

using DowndateState = UKF::StateVector<UKF::Field<180, UKF::Vector<2>>>;
using DowndateMeasurement = UKF::FixedMeasurementVector<UKF::Field<181, UKF::Vector<2>>>;
struct DowndateIdentityIntegrator {
    template <typename S> static S integrate(real_t, const S& state) { return state; }
};
namespace UKF {
template <> template <>
UKF::Vector<2> DowndateMeasurement::expected_measurement<DowndateState, 181>(const DowndateState& state) { return state; }
}
using DowndateFilter = UKF::SquareRootCore<DowndateState, DowndateMeasurement, DowndateIdentityIntegrator>;
using DowndateParameterFilter = UKF::SquareRootParameterEstimationCore<DowndateState, DowndateMeasurement>;

template <typename Filter>
void check_semidefinite_downdate(Filter filter, bool exact_both) {
    DowndateMeasurement z;
    z << 2, -3;
    filter.innovation_step(z);
    filter.a_posteriori_step();
    UKF::Matrix<2, 2> expected;
    expected << 0, 0, 0, exact_both ? real_t(0) : real_t(0.5);
    EXPECT_LT((filter.root_covariance * filter.root_covariance.transpose() - expected).norm(), real_t(1e-5));
    EXPECT_NEAR(real_t(2), filter.state(0), real_t(1e-5));
    EXPECT_NEAR(exact_both ? real_t(-3) : real_t(-1.5), filter.state(1), real_t(1e-5));
}

template <typename Filter>
void check_invalid_downdate(Filter filter) {
    DowndateMeasurement z;
    z << 2, -3;
    filter.innovation_step(z);
    filter.innovation_root_covariance *= real_t(0.1);
    const DowndateState original_state = filter.state;
    const UKF::Matrix<2, 2> original_root = filter.root_covariance;
    EXPECT_FALSE(filter.a_posteriori_step());
    EXPECT_EQ(real_t(0), (filter.state - original_state).norm());
    EXPECT_EQ(real_t(0), (filter.root_covariance - original_root).norm());
}

TEST(SquareRootDowndateTest, ExactAndPartiallyExactObservations) {
    for (bool exact_both : {false, true}) {
        DowndateFilter filter;
        filter.state.setZero();
        filter.root_covariance.setIdentity();
        filter.process_noise_root_covariance.setZero();
        filter.measurement_root_covariance << 0, exact_both ? 0 : 1;
        filter.a_priori_step(real_t(0));
        check_semidefinite_downdate(filter, exact_both);
        DowndateParameterFilter parameter;
        parameter.state.setZero();
        parameter.root_covariance.setIdentity();
        parameter.process_noise_root_covariance.setZero();
        parameter.measurement_root_covariance = filter.measurement_root_covariance;
        parameter.a_priori_step();
        check_semidefinite_downdate(parameter, exact_both);
    }
}

TEST(SquareRootDowndateTest, FailureDoesNotPartiallyApplyCorrection) {
    DowndateFilter filter;
    filter.state.setZero();
    filter.root_covariance.setIdentity();
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance.setOnes();
    filter.a_priori_step(real_t(0));
    check_invalid_downdate(filter);
    DowndateParameterFilter parameter;
    parameter.state.setZero();
    parameter.root_covariance.setIdentity();
    parameter.process_noise_root_covariance.setZero();
    parameter.measurement_root_covariance.setOnes();
    parameter.a_priori_step();
    check_invalid_downdate(parameter);
}

template <typename Scalar, int UpLo>
void check_bounded_rank_updates() {
    using Matrix = Eigen::Matrix<Scalar, 3, 3>;
    using Vector = Eigen::Matrix<Scalar, 3, 1>;
    Matrix lower;
    lower << -2, 0, 0, Scalar(0.5), 3, 0, Scalar(-0.3), Scalar(0.2), Scalar(-1.5);
    Vector update;
    update << Scalar(0.25), Scalar(-0.5), Scalar(0.125);
    for (Scalar weight : {Scalar(1.75), Scalar(-0.75)}) {
        const Matrix expected = lower * lower.transpose() + weight * update * update.transpose();
        Matrix factor = lower;
        if constexpr (UpLo == Eigen::Upper) {
            factor.transposeInPlace();
        }
        // Sentinel values outside the stored triangle must remain untouched.
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                if ((UpLo == Eigen::Lower && row < column) ||
                        (UpLo == Eigen::Upper && row > column)) {
                    factor(row, column) = Scalar(17 + row + column);
                }
            }
        }
        ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(factor, update, weight));
        const Matrix triangle = factor.template triangularView<UpLo>();
        Matrix actual;
        if constexpr (UpLo == Eigen::Lower) {
            actual = triangle * triangle.transpose();
        } else {
            actual = triangle.transpose() * triangle;
        }
        EXPECT_LT((actual - expected).norm(),
            Scalar(128) * std::numeric_limits<Scalar>::epsilon() * expected.norm());
        for (int row = 0; row < 3; ++row) {
            EXPECT_GT(factor(row, row), Scalar(0));
            for (int column = 0; column < 3; ++column) {
                if ((UpLo == Eigen::Lower && row < column) ||
                        (UpLo == Eigen::Upper && row > column)) {
                    EXPECT_EQ(Scalar(17 + row + column), factor(row, column));
                }
            }
        }
    }

    // Runtime sizes retain their compile-time storage bound, including zero.
    using BoundedMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, 3, 3>;
    using BoundedVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, 0, 3, 1>;
    for (int size = 0; size <= 3; ++size) {
        BoundedMatrix factor = BoundedMatrix::Zero(size, size);
        BoundedVector vector = update.head(size);
        if (size > 1) {
            vector(0) = Scalar(0); // Exercise a leading zero pivot.
        }
        ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(factor, vector, Scalar(0)));
        EXPECT_EQ(Scalar(0), factor.squaredNorm());
        const BoundedMatrix expected = Scalar(2) * vector * vector.transpose();
        ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(factor, vector, Scalar(2)));
        BoundedMatrix actual;
        if constexpr (UpLo == Eigen::Lower) {
            actual = factor * factor.transpose();
        } else {
            actual = factor.transpose() * factor;
        }
        EXPECT_LE((actual - expected).norm(), Scalar(128) *
            std::numeric_limits<Scalar>::epsilon() * std::max(Scalar(1), expected.norm()));
        const BoundedMatrix unchanged = factor;
        ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(factor, vector, Scalar(0)));
        EXPECT_EQ(Scalar(0), (factor - unchanged).norm());
    }

    // Fixed scalar roots are also used by both square-root filter variants.
    Eigen::Matrix<Scalar, 1, 1> scalar_factor, scalar_vector;
    scalar_factor << Scalar(-2);
    scalar_vector << Scalar(1);
    ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(scalar_factor, scalar_vector, Scalar(-1)));
    EXPECT_LT(std::abs(scalar_factor.squaredNorm() - Scalar(3)),
        Scalar(16) * std::numeric_limits<Scalar>::epsilon());
    EXPECT_GT(scalar_factor(0), Scalar(0));
    scalar_factor << Scalar(0);
    ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(scalar_factor, scalar_vector, Scalar(2)));
    EXPECT_LT(std::abs(scalar_factor.squaredNorm() - Scalar(2)),
        Scalar(16) * std::numeric_limits<Scalar>::epsilon());
    scalar_factor << Scalar(1);
    EXPECT_EQ(0, UKF::Detail::rank_update<UpLo>(scalar_factor, scalar_vector, Scalar(-1)));

    // A zero pivot with a zero entry is a valid (identity) downdate step, so
    // a semidefinite factor must be downdated through it, not rejected.
    Matrix semidefinite;
    semidefinite << 0, 0, 0, 0, 2, 0, 0, Scalar(0.5), 1;
    if constexpr (UpLo == Eigen::Upper) {
        semidefinite.transposeInPlace();
    }
    Vector downdate;
    downdate << 0, 1, Scalar(0.25);
    Matrix expected_semidefinite;
    if constexpr (UpLo == Eigen::Lower) {
        expected_semidefinite = semidefinite * semidefinite.transpose();
    } else {
        expected_semidefinite = semidefinite.transpose() * semidefinite;
    }
    expected_semidefinite -= downdate * downdate.transpose();
    ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(semidefinite, downdate, Scalar(-1)));
    Matrix actual_semidefinite;
    if constexpr (UpLo == Eigen::Lower) {
        actual_semidefinite = semidefinite * semidefinite.transpose();
    } else {
        actual_semidefinite = semidefinite.transpose() * semidefinite;
    }
    EXPECT_LT((actual_semidefinite - expected_semidefinite).norm(),
        Scalar(128) * std::numeric_limits<Scalar>::epsilon() * expected_semidefinite.norm());

    Matrix invalid = Matrix::Identity();
    update << Scalar(0.1), Scalar(2), Scalar(0);
    EXPECT_EQ(1, UKF::Detail::rank_update<UpLo>(invalid, update, Scalar(-1)));

    // The input may reference the factor: workspace must be copied first.
    Matrix aliased = lower;
    if constexpr (UpLo == Eigen::Upper) {
        aliased.transposeInPlace();
    }
    const Vector aliased_vector = aliased.col(0);
    const Matrix expected_alias = lower * lower.transpose() +
        Scalar(0.5) * aliased_vector * aliased_vector.transpose();
    ASSERT_EQ(-1, UKF::Detail::rank_update<UpLo>(aliased, aliased.col(0), Scalar(0.5)));
    Matrix actual_alias;
    if constexpr (UpLo == Eigen::Lower) {
        actual_alias = aliased * aliased.transpose();
    } else {
        actual_alias = aliased.transpose() * aliased;
    }
    EXPECT_LT((actual_alias - expected_alias).norm(), Scalar(128) *
        std::numeric_limits<Scalar>::epsilon() * expected_alias.norm());
}

TEST(SquareRootRankUpdateTest, FixedAndBoundedDoublePrecision) {
    check_bounded_rank_updates<double, Eigen::Lower>();
    check_bounded_rank_updates<double, Eigen::Upper>();
}

TEST(SquareRootRankUpdateTest, FixedAndBoundedSinglePrecision) {
    check_bounded_rank_updates<float, Eigen::Lower>();
    check_bounded_rank_updates<float, Eigen::Upper>();
}

/* The square-root cores accept direct process models, like Core. */
using DirectRootState = UKF::StateVector<UKF::Field<190, UKF::Vector<2>>>;
using DirectRootMeasurement = UKF::DynamicMeasurementVector<
    UKF::Field<191, UKF::Vector<2>>,
    UKF::Field<192, real_t>>;
namespace UKF {
template <> template <>
DirectRootState DirectRootState::process_model<>(real_t dt) const {
    DirectRootState next;
    next << (*this)(0) + dt * (*this)(1), (*this)(1);
    return next;
}
template <> template <>
UKF::Vector<2> DirectRootMeasurement::expected_measurement<DirectRootState, 191>(const DirectRootState& state) {
    return state;
}
template <> template <>
real_t DirectRootMeasurement::expected_measurement<DirectRootState, 192>(const DirectRootState& state) {
    return state(0);
}
}

TEST(SquareRootDirectModelTest, PredictsWithDirectProcessModel) {
    UKF::SquareRootCore<DirectRootState, DirectRootMeasurement> filter;
    filter.state << 1, 2;
    filter.root_covariance.setIdentity();
    filter.process_noise_root_covariance.setZero();
    ASSERT_TRUE(filter.a_priori_step(real_t(0.5)));
    EXPECT_LT((filter.state - UKF::Vector<2>(2, 2)).norm(), real_t(1e-9));
    UKF::Matrix<2, 2> transition, expected;
    transition << 1, 0.5, 0, 1;
    expected = transition * transition.transpose();
    EXPECT_LT((filter.root_covariance * filter.root_covariance.transpose() - expected).norm(), real_t(1e-9));
}

TEST(SquareRootDirectModelTest, EmptyMeasurementPredictsOnly) {
    UKF::SquareRootCore<DirectRootState, DirectRootMeasurement> filter;
    filter.state << 1, 2;
    filter.root_covariance.setIdentity();
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance << 1, 1, 1;
    ASSERT_TRUE(filter.step(real_t(0.5), DirectRootMeasurement()));
    EXPECT_EQ(0, filter.innovation.size());
    EXPECT_LT((filter.state - UKF::Vector<2>(2, 2)).norm(), real_t(1e-9));

    DirectRootMeasurement z;
    z.set_field<192>(real_t(2.5));
    const UKF::Matrix<2, 2> prior = filter.root_covariance * filter.root_covariance.transpose();
    ASSERT_TRUE(filter.step(real_t(0), z));
    const UKF::Matrix<2, 2> posterior = filter.root_covariance * filter.root_covariance.transpose();
    EXPECT_LT(posterior(0, 0), prior(0, 0));

    UKF::SquareRootParameterEstimationCore<DirectRootState, DirectRootMeasurement> parameter;
    parameter.state << 1, 2;
    parameter.root_covariance.setIdentity();
    parameter.process_noise_root_covariance.setZero();
    parameter.measurement_root_covariance << 1, 1, 1;
    ASSERT_TRUE(parameter.step(DirectRootMeasurement()));
    EXPECT_EQ(UKF::Vector<2>(1, 2), UKF::Vector<2>(parameter.state));
}

TEST(SquareRootDirectModelTest, StepsRequireTheirPrerequisites) {
    UKF::SquareRootCore<DirectRootState, DirectRootMeasurement> filter;
    filter.state << 1, 2;
    filter.root_covariance.setIdentity();
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance << 1, 1, 1;
    DirectRootMeasurement z;
    z.set_field<192>(real_t(2.5));
    EXPECT_FALSE(filter.innovation_step(z));
    EXPECT_FALSE(filter.a_posteriori_step());
    EXPECT_EQ(UKF::Vector<2>(1, 2), UKF::Vector<2>(filter.state));
}

/*
After an exact observation the posterior factor has a zero row. The next
correction must downdate through that zero pivot rather than rebuilding the
factor or failing, and must keep the observed component exact.
*/
TEST(SquareRootDowndateTest, SemidefinitePosteriorKeepsCorrecting) {
    DowndateFilter filter;
    filter.state.setZero();
    filter.root_covariance.setIdentity();
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance << 0, 1;
    DowndateMeasurement z;
    z << 2, -3;
    for (int i = 0; i < 3; ++i) {
        ASSERT_TRUE(filter.step(real_t(0), z));
        const UKF::Matrix<2, 2> covariance = filter.root_covariance * filter.root_covariance.transpose();
        EXPECT_NEAR(real_t(2), filter.state(0), real_t(1e-9));
        EXPECT_LT(std::abs(covariance(0, 0)), real_t(1e-12));
        EXPECT_NEAR(real_t(1) / real_t(i + 2), covariance(1, 1), real_t(1e-9));
    }
}
