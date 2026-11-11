#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
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
UKF::FieldVector MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer, UKF::Vector<3>>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) {
    return state.get_field<Attitude>() * UKF::FieldVector(1, 0, 0);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope, UKF::Vector<3>>(const MyStateVector& state,
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) {
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
    EXPECT_THROW(filter.a_posteriori_step(), std::domain_error);
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
