#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "UKF/Types.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "comparisons.h"

enum MyFields {
    StaticPressure,
    DynamicPressure,
    Accelerometer,
    Gyroscope
};

using MyMeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<StaticPressure, real_t>,
    UKF::Field<DynamicPressure, real_t>,
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Gyroscope, UKF::Vector<3>>
>;

TEST(DynamicMeasurementVectorTest, Instantiation) {
    MyMeasurementVector test_measurement;

    EXPECT_EQ(8, MyMeasurementVector::MaxRowsAtCompileTime);
    EXPECT_EQ(8, test_measurement.max_size());
}

TEST(DynamicMeasurementVectorTest, Assignment) {
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.get_field<Gyroscope>());

    test_measurement.set_field<DynamicPressure>(4);

    EXPECT_EQ(4, test_measurement.size());
    EXPECT_EQ(4, test_measurement.get_field<DynamicPressure>());

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));

    EXPECT_EQ(7, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(5, 6, 7), test_measurement.get_field<Accelerometer>());

    test_measurement.set_field<StaticPressure>(8);

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(8, test_measurement.get_field<StaticPressure>());

    UKF::Vector<8> expected;
    expected << 1, 2, 3, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(DynamicMeasurementVectorTest, Reassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement);
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.get_field<Gyroscope>());

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(4, 5, 6));

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement);
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.get_field<Gyroscope>());
}

TEST(DynamicMeasurementVectorTest, MultipleReassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    test_measurement.set_field<DynamicPressure>(4);
    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    test_measurement.set_field<StaticPressure>(8);

    EXPECT_EQ(8, test_measurement.size());
    UKF::Vector<8> expected;
    expected << 1, 2, 3, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(4, 5, 6));

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.get_field<Gyroscope>());
    expected << 4, 5, 6, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(7, 8, 9));

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(7, 8, 9), test_measurement.get_field<Accelerometer>());
    expected << 4, 5, 6, 4, 7, 8, 9, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<DynamicPressure>(1);

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(1, test_measurement.get_field<DynamicPressure>());
    expected << 4, 5, 6, 1, 7, 8, 9, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<StaticPressure>(3);

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(3, test_measurement.get_field<StaticPressure>());
    expected << 4, 5, 6, 1, 7, 8, 9, 3;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(DynamicMeasurementVectorTest, CopyConstructor) {
    MyMeasurementVector test_measurement_1, test_measurement_2;

    test_measurement_1.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    test_measurement_1.set_field<DynamicPressure>(4);
    test_measurement_1.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    test_measurement_1.set_field<StaticPressure>(8);

    test_measurement_2 = test_measurement_1;

    EXPECT_EQ(test_measurement_1.size(), test_measurement_2.size());
    EXPECT_VECTOR_EQ(test_measurement_1, test_measurement_2);
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement_2.get_field<Gyroscope>());
    EXPECT_EQ(4, test_measurement_2.get_field<DynamicPressure>());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(5, 6, 7), test_measurement_2.get_field<Accelerometer>());
    EXPECT_EQ(8, test_measurement_2.get_field<StaticPressure>());
}

TEST(DynamicMeasurementVectorTest, Arithmetic) {
    
}

enum MyStateVectorFields {
    AngularVelocity,
    Altitude,
    Velocity,
    Attitude
};

using MyStateVector = UKF::StateVector<
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<AngularVelocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>,
    UKF::Field<Altitude, real_t>
>;

/*
Define measurement model to be used in tests. NOTE: These are just for
testing, don't expect them to make any physical sense whatsoever.
*/
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

/*
These versions of the predicted measurement functions have non-state inputs.
This could be used to add predicted kinematic acceleration by feeding control
inputs into a dynamics model, for example.
*/
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

TEST(DynamicMeasurementVectorTest, SigmaPointGeneration) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<StaticPressure>(0);
    test_measurement.set_field<DynamicPressure>(0);

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> target_sigma_points(
        test_measurement.size(), MyStateVector::num_sigma());

    target_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,
                            0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,
                         -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,
                            1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                            0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                            0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                         89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.257,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.343,
                        8.575, 20.954, 25.371, 29.788,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575, 12.121,  7.704,  3.287,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575;
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points;
    measurement_sigma_points = test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    EXPECT_VECTOR_EQ(target_sigma_points.col(0),  measurement_sigma_points.col(0));
    EXPECT_VECTOR_EQ(target_sigma_points.col(1),  measurement_sigma_points.col(1));
    EXPECT_VECTOR_EQ(target_sigma_points.col(2),  measurement_sigma_points.col(2));
    EXPECT_VECTOR_EQ(target_sigma_points.col(3),  measurement_sigma_points.col(3));
    EXPECT_VECTOR_EQ(target_sigma_points.col(4),  measurement_sigma_points.col(4));
    EXPECT_VECTOR_EQ(target_sigma_points.col(5),  measurement_sigma_points.col(5));
    EXPECT_VECTOR_EQ(target_sigma_points.col(6),  measurement_sigma_points.col(6));
    EXPECT_VECTOR_EQ(target_sigma_points.col(7),  measurement_sigma_points.col(7));
    EXPECT_VECTOR_EQ(target_sigma_points.col(8),  measurement_sigma_points.col(8));
    EXPECT_VECTOR_EQ(target_sigma_points.col(9),  measurement_sigma_points.col(9));
    EXPECT_VECTOR_EQ(target_sigma_points.col(10), measurement_sigma_points.col(10));
    EXPECT_VECTOR_EQ(target_sigma_points.col(11), measurement_sigma_points.col(11));
    EXPECT_VECTOR_EQ(target_sigma_points.col(12), measurement_sigma_points.col(12));
    EXPECT_VECTOR_EQ(target_sigma_points.col(13), measurement_sigma_points.col(13));
    EXPECT_VECTOR_EQ(target_sigma_points.col(14), measurement_sigma_points.col(14));
    EXPECT_VECTOR_EQ(target_sigma_points.col(15), measurement_sigma_points.col(15));
    EXPECT_VECTOR_EQ(target_sigma_points.col(16), measurement_sigma_points.col(16));
    EXPECT_VECTOR_EQ(target_sigma_points.col(17), measurement_sigma_points.col(17));
    EXPECT_VECTOR_EQ(target_sigma_points.col(18), measurement_sigma_points.col(18));
    EXPECT_VECTOR_EQ(target_sigma_points.col(19), measurement_sigma_points.col(19));
    EXPECT_VECTOR_EQ(target_sigma_points.col(20), measurement_sigma_points.col(20));
}

TEST(DynamicMeasurementVectorTest, SigmaPointGenerationWithInput) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<StaticPressure>(0);
    test_measurement.set_field<DynamicPressure>(0);

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> target_sigma_points(
        test_measurement.size(), MyStateVector::num_sigma());

    target_sigma_points <<  1,      1,      1,      1,      1,      1,      1,      1, -7.313,      1,      1,      1,      1,      1,      1,      1,      1,      1,  9.313,      1,      1,
                            2,      2,      2,      2,      2,      2,      2, 10.313,      2,      2,      2,      2,      2,      2,      2,      2,      2, -6.313,      2,      2,      2,
                         -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,  8.188,  8.188,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,  8.188,  8.188,   -6.8,   -6.8,
                            1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                            0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                            0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                         89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.257,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.343,
                        8.575, 20.954, 25.371, 29.788,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575, 12.121,  7.704,  3.287,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575;
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points;
    measurement_sigma_points = test_measurement.calculate_sigma_point_distribution<MyStateVector>(
        sigma_points, UKF::Vector<3>(1, 2, 3));

    EXPECT_VECTOR_EQ(target_sigma_points.col(0),  measurement_sigma_points.col(0));
    EXPECT_VECTOR_EQ(target_sigma_points.col(1),  measurement_sigma_points.col(1));
    EXPECT_VECTOR_EQ(target_sigma_points.col(2),  measurement_sigma_points.col(2));
    EXPECT_VECTOR_EQ(target_sigma_points.col(3),  measurement_sigma_points.col(3));
    EXPECT_VECTOR_EQ(target_sigma_points.col(4),  measurement_sigma_points.col(4));
    EXPECT_VECTOR_EQ(target_sigma_points.col(5),  measurement_sigma_points.col(5));
    EXPECT_VECTOR_EQ(target_sigma_points.col(6),  measurement_sigma_points.col(6));
    EXPECT_VECTOR_EQ(target_sigma_points.col(7),  measurement_sigma_points.col(7));
    EXPECT_VECTOR_EQ(target_sigma_points.col(8),  measurement_sigma_points.col(8));
    EXPECT_VECTOR_EQ(target_sigma_points.col(9),  measurement_sigma_points.col(9));
    EXPECT_VECTOR_EQ(target_sigma_points.col(10), measurement_sigma_points.col(10));
    EXPECT_VECTOR_EQ(target_sigma_points.col(11), measurement_sigma_points.col(11));
    EXPECT_VECTOR_EQ(target_sigma_points.col(12), measurement_sigma_points.col(12));
    EXPECT_VECTOR_EQ(target_sigma_points.col(13), measurement_sigma_points.col(13));
    EXPECT_VECTOR_EQ(target_sigma_points.col(14), measurement_sigma_points.col(14));
    EXPECT_VECTOR_EQ(target_sigma_points.col(15), measurement_sigma_points.col(15));
    EXPECT_VECTOR_EQ(target_sigma_points.col(16), measurement_sigma_points.col(16));
    EXPECT_VECTOR_EQ(target_sigma_points.col(17), measurement_sigma_points.col(17));
    EXPECT_VECTOR_EQ(target_sigma_points.col(18), measurement_sigma_points.col(18));
    EXPECT_VECTOR_EQ(target_sigma_points.col(19), measurement_sigma_points.col(19));
    EXPECT_VECTOR_EQ(target_sigma_points.col(20), measurement_sigma_points.col(20));
}

TEST(DynamicMeasurementVectorTest, PartialSigmaPointGeneration) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> target_sigma_points(
        test_measurement.size(), MyStateVector::num_sigma());

    target_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,
                            0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,
                         -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8;
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points;
    measurement_sigma_points = test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    EXPECT_VECTOR_EQ(target_sigma_points.col(0),  measurement_sigma_points.col(0));
    EXPECT_VECTOR_EQ(target_sigma_points.col(1),  measurement_sigma_points.col(1));
    EXPECT_VECTOR_EQ(target_sigma_points.col(2),  measurement_sigma_points.col(2));
    EXPECT_VECTOR_EQ(target_sigma_points.col(3),  measurement_sigma_points.col(3));
    EXPECT_VECTOR_EQ(target_sigma_points.col(4),  measurement_sigma_points.col(4));
    EXPECT_VECTOR_EQ(target_sigma_points.col(5),  measurement_sigma_points.col(5));
    EXPECT_VECTOR_EQ(target_sigma_points.col(6),  measurement_sigma_points.col(6));
    EXPECT_VECTOR_EQ(target_sigma_points.col(7),  measurement_sigma_points.col(7));
    EXPECT_VECTOR_EQ(target_sigma_points.col(8),  measurement_sigma_points.col(8));
    EXPECT_VECTOR_EQ(target_sigma_points.col(9),  measurement_sigma_points.col(9));
    EXPECT_VECTOR_EQ(target_sigma_points.col(10), measurement_sigma_points.col(10));
    EXPECT_VECTOR_EQ(target_sigma_points.col(11), measurement_sigma_points.col(11));
    EXPECT_VECTOR_EQ(target_sigma_points.col(12), measurement_sigma_points.col(12));
    EXPECT_VECTOR_EQ(target_sigma_points.col(13), measurement_sigma_points.col(13));
    EXPECT_VECTOR_EQ(target_sigma_points.col(14), measurement_sigma_points.col(14));
    EXPECT_VECTOR_EQ(target_sigma_points.col(15), measurement_sigma_points.col(15));
    EXPECT_VECTOR_EQ(target_sigma_points.col(16), measurement_sigma_points.col(16));
    EXPECT_VECTOR_EQ(target_sigma_points.col(17), measurement_sigma_points.col(17));
    EXPECT_VECTOR_EQ(target_sigma_points.col(18), measurement_sigma_points.col(18));
    EXPECT_VECTOR_EQ(target_sigma_points.col(19), measurement_sigma_points.col(19));
    EXPECT_VECTOR_EQ(target_sigma_points.col(20), measurement_sigma_points.col(20));
}

TEST(DynamicMeasurementVectorTest, SigmaPointMean) {
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points(
        MyMeasurementVector::max_size(), MyStateVector::num_sigma());
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<StaticPressure>(0);
    test_measurement.set_field<DynamicPressure>(0);

    measurement_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,
                                 0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,
                              -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,
                                 1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                                 0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                                 0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                              89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.257,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.343,
                             8.575, 20.954, 25.371, 29.788,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575, 12.121,  7.704,  3.287,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575;

    MyMeasurementVector expected_mean;

    expected_mean.set_field<Accelerometer>(UKF::Vector<3>(0.0, 0.0, -7.494));
    expected_mean.set_field<Gyroscope>(UKF::Vector<3>(1.0, 0.0, 0.0));
    expected_mean.set_field<StaticPressure>(89.3);
    expected_mean.set_field<DynamicPressure>(10.4125);

    EXPECT_VECTOR_EQ(expected_mean, test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points));
}

TEST(DynamicMeasurementVectorTest, PartialSigmaPointMean) {
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points(
        3, MyStateVector::num_sigma());
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));

    measurement_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,
                                 0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,
                              -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8;

    MyMeasurementVector expected_mean;

    expected_mean.set_field<Accelerometer>(UKF::Vector<3>(0.0, 0.0, -7.494));

    EXPECT_VECTOR_EQ(expected_mean, test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points));
}

TEST(DynamicMeasurementVectorTest, SigmaPointDeltas) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<StaticPressure>(0);
    test_measurement.set_field<DynamicPressure>(0);

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> sigma_point_deltas(mean_measurement.size(), MyStateVector::num_sigma());
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> target_sigma_point_deltas(mean_measurement.size(), MyStateVector::num_sigma());

    target_sigma_point_deltas <<       0,       0,       0,       0,       0,       0,       0,       0,  -8.314,       0,       0,      0,       0,       0,       0,       0,       0,       0,   8.314,       0,       0,
                                       0,       0,       0,       0,       0,       0,       0,   8.314,       0,       0,       0,      0,       0,       0,       0,       0,       0,  -8.314,       0,       0,       0,
                                  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  12.682,  12.682,  -2.306,  -2.306, -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  12.682,  12.682,  -2.306,  -2.306,
                                       0,       0,       0,       0,   3.606,       0,       0,       0,       0,       0,       0,      0,       0,       0,  -3.606,       0,       0,       0,       0,       0,       0,
                                       0,       0,       0,       0,       0,   3.606,       0,       0,       0,       0,       0,      0,       0,       0,       0,  -3.606,       0,       0,       0,       0,       0,
                                       0,       0,       0,       0,       0,       0,   3.606,       0,       0,       0,       0,      0,       0,       0,       0,       0,  -3.606,       0,       0,       0,       0,
                                       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -0.043,      0,       0,       0,       0,       0,       0,       0,       0,       0,   0.043,
                                 -1.8375,  10.542,  14.959,  19.376, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, 1.7082, -2.7086,  -7.126, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375;
    sigma_point_deltas = mean_measurement.calculate_sigma_point_deltas<MyStateVector>(measurement_sigma_points);

    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(0),  sigma_point_deltas.col(0));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(1),  sigma_point_deltas.col(1));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(2),  sigma_point_deltas.col(2));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(3),  sigma_point_deltas.col(3));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(4),  sigma_point_deltas.col(4));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(5),  sigma_point_deltas.col(5));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(6),  sigma_point_deltas.col(6));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(7),  sigma_point_deltas.col(7));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(8),  sigma_point_deltas.col(8));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(9),  sigma_point_deltas.col(9));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(10), sigma_point_deltas.col(10));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(11), sigma_point_deltas.col(11));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(12), sigma_point_deltas.col(12));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(13), sigma_point_deltas.col(13));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(14), sigma_point_deltas.col(14));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(15), sigma_point_deltas.col(15));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(16), sigma_point_deltas.col(16));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(17), sigma_point_deltas.col(17));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(18), sigma_point_deltas.col(18));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(19), sigma_point_deltas.col(19));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(20), sigma_point_deltas.col(20));
}

TEST(DynamicMeasurementVectorTest, PartialSigmaPointDeltas) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> sigma_point_deltas(3, MyStateVector::num_sigma());
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> target_sigma_point_deltas(3, MyStateVector::num_sigma());

    target_sigma_point_deltas <<       0,       0,       0,       0,       0,       0,       0,       0,  -8.314,       0,       0,      0,       0,       0,       0,       0,       0,       0,   8.314,       0,       0,
                                       0,       0,       0,       0,       0,       0,       0,   8.314,       0,       0,       0,      0,       0,       0,       0,       0,       0,  -8.314,       0,       0,       0,
                                  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  12.682,  12.682,  -2.306,  -2.306, -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  12.682,  12.682,  -2.306,  -2.306;
    sigma_point_deltas = mean_measurement.calculate_sigma_point_deltas<MyStateVector>(measurement_sigma_points);

    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(0),  sigma_point_deltas.col(0));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(1),  sigma_point_deltas.col(1));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(2),  sigma_point_deltas.col(2));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(3),  sigma_point_deltas.col(3));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(4),  sigma_point_deltas.col(4));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(5),  sigma_point_deltas.col(5));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(6),  sigma_point_deltas.col(6));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(7),  sigma_point_deltas.col(7));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(8),  sigma_point_deltas.col(8));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(9),  sigma_point_deltas.col(9));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(10), sigma_point_deltas.col(10));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(11), sigma_point_deltas.col(11));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(12), sigma_point_deltas.col(12));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(13), sigma_point_deltas.col(13));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(14), sigma_point_deltas.col(14));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(15), sigma_point_deltas.col(15));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(16), sigma_point_deltas.col(16));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(17), sigma_point_deltas.col(17));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(18), sigma_point_deltas.col(18));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(19), sigma_point_deltas.col(19));
    EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(20), sigma_point_deltas.col(20));
}

TEST(DynamicMeasurementVectorTest, SigmaPointCovariance) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<StaticPressure>(0);
    test_measurement.set_field<DynamicPressure>(0);

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> sigma_point_deltas =
        mean_measurement.calculate_sigma_point_deltas<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::CovarianceMatrix calculated_covariance =
        mean_measurement.calculate_sigma_point_covariance<MyStateVector>(sigma_point_deltas);

    MyMeasurementVector::CovarianceMatrix expected_covariance(mean_measurement.size(), mean_measurement.size());

    expected_covariance << 5.31709,       0,       0,       0,       0,       0,       0,       0,
                                 0, 5.31709,       0,       0,       0,       0,       0,       0,
                                 0,       0, 29.2440,       0,       0,       0,       0,  -4.237,
                                 0,       0,       0,       1,       0,       0,       0,       0,
                                 0,       0,       0,       0,       1,       0,       0,       0,
                                 0,       0,       0,       0,       0,       1,       0,       0,
                                 0,       0,       0,       0,       0,       0, 1.44e-4,       0,
                                 0,       0,  -4.237,       0,       0,       0,       0,  32.263;

    EXPECT_VECTOR_EQ(expected_covariance.col(0),  calculated_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_covariance.col(1),  calculated_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_covariance.col(2),  calculated_covariance.col(2));
    EXPECT_VECTOR_EQ(expected_covariance.col(3),  calculated_covariance.col(3));
    EXPECT_VECTOR_EQ(expected_covariance.col(4),  calculated_covariance.col(4));
    EXPECT_VECTOR_EQ(expected_covariance.col(5),  calculated_covariance.col(5));
    EXPECT_VECTOR_EQ(expected_covariance.col(6),  calculated_covariance.col(6));
    EXPECT_VECTOR_EQ(expected_covariance.col(7),  calculated_covariance.col(7));
}

TEST(DynamicMeasurementVectorTest, PartialSigmaPointCovariance) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 0));

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> sigma_point_deltas =
        mean_measurement.calculate_sigma_point_deltas<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::CovarianceMatrix calculated_covariance =
        mean_measurement.calculate_sigma_point_covariance<MyStateVector>(sigma_point_deltas);

    MyMeasurementVector::CovarianceMatrix expected_covariance(mean_measurement.size(), mean_measurement.size());

    expected_covariance << 5.31709,       0,       0,
                                 0, 5.31709,       0,
                                 0,       0, 29.2440;

    EXPECT_VECTOR_EQ(expected_covariance.col(0),  calculated_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_covariance.col(1),  calculated_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_covariance.col(2),  calculated_covariance.col(2));
}

template <>
MyMeasurementVector::CovarianceVector MyMeasurementVector::measurement_covariance = MyMeasurementVector::CovarianceVector();

TEST(DynamicMeasurementVectorTest, MeasurementCovariance) {
    MyMeasurementVector test_measurement;

    MyMeasurementVector::measurement_covariance.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    MyMeasurementVector::measurement_covariance.set_field<DynamicPressure>(4);
    MyMeasurementVector::measurement_covariance.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    MyMeasurementVector::measurement_covariance.set_field<StaticPressure>(8);

    MyMeasurementVector::CovarianceMatrix expected_measurement_covariance = MyMeasurementVector::CovarianceMatrix::Zero(8, 8);
    expected_measurement_covariance.diagonal() << 5, 6, 7, 1, 2, 3, 8, 4;

    MyMeasurementVector::CovarianceMatrix measurement_covariance = test_measurement.calculate_measurement_covariance();

    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(0),  measurement_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(1),  measurement_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(2),  measurement_covariance.col(2));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(3),  measurement_covariance.col(3));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(4),  measurement_covariance.col(4));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(5),  measurement_covariance.col(5));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(6),  measurement_covariance.col(6));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(7),  measurement_covariance.col(7));
}

TEST(DynamicMeasurementVectorTest, PartialMeasurementCovariance) {
    MyMeasurementVector test_measurement;

    MyMeasurementVector::measurement_covariance.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));

    MyMeasurementVector::CovarianceMatrix expected_measurement_covariance = MyMeasurementVector::CovarianceMatrix::Zero(3, 3);
    expected_measurement_covariance.diagonal() << 1, 2, 3;

    MyMeasurementVector::CovarianceMatrix measurement_covariance = test_measurement.calculate_measurement_covariance();

    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(0),  measurement_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(1),  measurement_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(2),  measurement_covariance.col(2));
}
