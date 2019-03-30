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
    Gyroscope,
    Magnetometer
};

using MyMeasurementVector = UKF::FixedMeasurementVector<
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Gyroscope, UKF::Vector<3>>,
    UKF::Field<Magnetometer, UKF::FieldVector>,
    UKF::Field<StaticPressure, real_t>,
    UKF::Field<DynamicPressure, real_t>
>;

TEST(FixedMeasurementVectorTest, Instantiation) {
    MyMeasurementVector test_measurement;

    EXPECT_EQ(11, MyMeasurementVector::MaxRowsAtCompileTime);
    EXPECT_EQ(11, test_measurement.size());
}

TEST(FixedMeasurementVectorTest, Assignment) {
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    test_measurement.set_field<DynamicPressure>(4);
    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    test_measurement.set_field<StaticPressure>(8);
    test_measurement.set_field<Magnetometer>(UKF::FieldVector(9, 10, 11));

    EXPECT_EQ(11, test_measurement.size());

    EXPECT_EQ(8, test_measurement.get_field<StaticPressure>());
    EXPECT_EQ(4, test_measurement.get_field<DynamicPressure>());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.get_field<Gyroscope>());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(5, 6, 7), test_measurement.get_field<Accelerometer>());
    EXPECT_VECTOR_EQ(UKF::FieldVector(9, 10, 11), test_measurement.get_field<Magnetometer>());

    UKF::Vector<11> expected;
    expected << 5, 6, 7, 1, 2, 3, 9, 10, 11, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(FixedMeasurementVectorTest, Reassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    test_measurement.set_field<DynamicPressure>(4);
    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    test_measurement.set_field<StaticPressure>(8);
    test_measurement.set_field<Magnetometer>(UKF::FieldVector(9, 10, 11));

    EXPECT_EQ(11, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.get_field<Gyroscope>());

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(4, 5, 6));

    EXPECT_EQ(11, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.get_field<Gyroscope>());

    UKF::Vector<11> expected;
    expected << 5, 6, 7, 4, 5, 6, 9, 10, 11, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(FixedMeasurementVectorTest, MultipleReassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    test_measurement.set_field<DynamicPressure>(4);
    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    test_measurement.set_field<StaticPressure>(8);
    test_measurement.set_field<Magnetometer>(UKF::FieldVector(9, 10, 11));

    EXPECT_EQ(11, test_measurement.size());
    UKF::Vector<11> expected;
    expected << 5, 6, 7, 1, 2, 3, 9, 10, 11, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(4, 5, 6));

    EXPECT_EQ(11, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.get_field<Gyroscope>());
    expected << 5, 6, 7, 4, 5, 6, 9, 10, 11, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(7, 8, 9));

    EXPECT_EQ(11, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(7, 8, 9), test_measurement.get_field<Accelerometer>());
    expected << 7, 8, 9, 4, 5, 6, 9, 10, 11, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<DynamicPressure>(1);

    EXPECT_EQ(11, test_measurement.size());
    EXPECT_EQ(1, test_measurement.get_field<DynamicPressure>());
    expected << 7, 8, 9, 4, 5, 6, 9, 10, 11, 8, 1;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<StaticPressure>(3);

    EXPECT_EQ(11, test_measurement.size());
    EXPECT_EQ(3, test_measurement.get_field<StaticPressure>());
    expected << 7, 8, 9, 4, 5, 6, 9, 10, 11, 3, 1;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.set_field<Magnetometer>(UKF::FieldVector(1, 2, 3));

    EXPECT_EQ(11, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.get_field<Magnetometer>());
    expected << 7, 8, 9, 4, 5, 6, 1, 2, 3, 3, 1;
    EXPECT_VECTOR_EQ(expected, test_measurement);
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

namespace UKF {
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

template <> template <>
UKF::FieldVector MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * UKF::FieldVector(0.45, 0, 0);
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

template <> template <>
UKF::FieldVector MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer, UKF::Vector<3>>(const MyStateVector& state, const UKF::Vector<3>& input) {
    return state.get_field<Attitude>() * UKF::FieldVector(0.45, 0, 0) + input;
}
}

TEST(FixedMeasurementVectorTest, SigmaPointGeneration) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution((covariance *
        (MyStateVector::covariance_size() + UKF::Parameters::Lambda<MyStateVector>)).llt().matrixL());

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points, target_sigma_points;

    target_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,
                            0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,
                         -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,
                            1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                            0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                            0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                         0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45, -0.238, -0.238,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45, -0.238, -0.238,   0.45,
                            0,      0,      0,      0,      0,      0,      0,      0,      0,  0.382,      0,      0,      0,      0,      0,      0,      0,      0,      0, -0.382,      0,
                            0,      0,      0,      0,      0,      0,      0,      0, -0.382,      0,      0,      0,      0,      0,      0,      0,      0,      0,  0.382,      0,      0,
                         89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.257,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.343,
                        8.575, 20.954, 25.371, 29.788,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575, 12.121,  7.704,  3.287,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575;
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

TEST(FixedMeasurementVectorTest, SigmaPointGenerationWithInput) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution((covariance *
        (MyStateVector::covariance_size() + UKF::Parameters::Lambda<MyStateVector>)).llt().matrixL());

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points, target_sigma_points;

    target_sigma_points <<  1,      1,      1,      1,      1,      1,      1,      1, -7.313,      1,      1,      1,      1,      1,      1,      1,      1,      1,  9.313,      1,      1,
                            2,      2,      2,      2,      2,      2,      2, 10.313,      2,      2,      2,      2,      2,      2,      2,      2,      2, -6.313,      2,      2,      2,
                         -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,  8.188,  8.188,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,   -6.8,  8.188,  8.188,   -6.8,   -6.8,
                            1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                            0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                            0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                         1.45,   1.45,   1.45,   1.45,   1.45,   1.45,   1.45,   1.45,  0.762,  0.762,   1.45,   1.45,   1.45,   1.45,   1.45,   1.45,   1.45,   1.45,  0.762,  0.762,   1.45,
                            2,      2,      2,      2,      2,      2,      2,      2,      2,  2.382,      2,      2,      2,      2,      2,      2,      2,      2,      2,  1.618,      2,
                            3,      3,      3,      3,      3,      3,      3,      3,  2.618,      3,      3,      3,      3,      3,      3,      3,      3,      3,  3.382,      3,      3,
                         89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.257,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.343,
                        8.575, 20.954, 25.371, 29.788,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575, 12.121,  7.704,  3.287,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575;
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

TEST(FixedMeasurementVectorTest, SigmaPointMean) {
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points;
    MyMeasurementVector test_measurement;

    measurement_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,
                                 0,      0,      0,      0,      0,      0,      0,  8.314,      0,      0,      0,      0,      0,      0,      0,      0,      0, -8.314,      0,      0,      0,
                              -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  5.188,  5.188,   -9.8,   -9.8,
                                 1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                                 0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                                 0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                              0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45, -0.238, -0.238,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45,   0.45, -0.238, -0.238,   0.45,
                                 0,      0,      0,      0,      0,      0,      0,      0,      0,  0.382,      0,      0,      0,      0,      0,      0,      0,      0,      0, -0.382,      0,
                                 0,      0,      0,      0,      0,      0,      0,      0, -0.382,      0,      0,      0,      0,      0,      0,      0,      0,      0,  0.382,      0,      0,
                              89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.257,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3,   89.3, 89.343,
                             8.575, 20.954, 25.371, 29.788,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575, 12.121,  7.704,  3.287,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575;

    MyMeasurementVector expected_mean;

    expected_mean << 0.0, 0.0, -7.494, 1, 0, 0, 0.45, 0, 0, 89.3, 10.4125;

    EXPECT_VECTOR_EQ(expected_mean, test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points));
}

TEST(FixedMeasurementVectorTest, SigmaPointDeltas) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Magnetometer>(UKF::FieldVector(1, 0, 0));

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution((covariance *
        (MyStateVector::covariance_size() + UKF::Parameters::Lambda<MyStateVector>)).llt().matrixL());

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> sigma_point_deltas, target_sigma_point_deltas;

    target_sigma_point_deltas <<       0,       0,       0,       0,       0,       0,       0,       0,  -8.314,       0,       0,       0,       0,       0,       0,       0,       0,       0,   8.314,       0,       0,
                                       0,       0,       0,       0,       0,       0,       0,   8.314,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -8.314,       0,       0,       0,
                                  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  12.682,  12.682,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  -2.306,  12.682,  12.682,  -2.306,  -2.306,
                                       0,       0,       0,       0,   3.606,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -3.606,       0,       0,       0,       0,       0,       0,
                                       0,       0,       0,       0,       0,   3.606,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -3.606,       0,       0,       0,       0,       0,
                                       0,       0,       0,       0,       0,       0,   3.606,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -3.606,       0,       0,       0,       0,
                                       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,
                                       0,       0,       0,       0,       0,       0,       0,       0,   3.606,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -3.606,       0,       0,
                                       0,       0,       0,       0,       0,       0,       0,       0,       0,   3.606,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -3.606,       0,
                                       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,  -0.043,       0,       0,       0,       0,       0,       0,       0,       0,       0,   0.043,
                                 -1.8375,  10.542,  14.959,  19.376, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375,  1.7082, -2.7086,  -7.126, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375, -1.8375;
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

TEST(FixedMeasurementVectorTest, Innovation) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement, target_innovation;

    test_state.set_field<Velocity>(UKF::Vector<3>(10, 0, 0));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 1));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(1, 0, 0));
    test_measurement.set_field<DynamicPressure>(0.5 * 122.5);
    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 9.8));
    test_measurement.set_field<StaticPressure>(101.3);
    test_measurement.set_field<Magnetometer>(UKF::FieldVector(0, 0.45, 0));

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1e-9, 1e-9, 1e-9, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution((covariance *
        (MyStateVector::covariance_size() + UKF::Parameters::Lambda<MyStateVector>)).llt().matrixL());

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);

    target_innovation.set_field<Gyroscope>(UKF::Vector<3>(1, 0, -1));
    target_innovation.set_field<DynamicPressure>(0);
    target_innovation.set_field<Accelerometer>(UKF::Vector<3>(0, 0, 17.294));
    target_innovation.set_field<StaticPressure>(12.0);
    target_innovation.set_field<Magnetometer>(UKF::FieldVector(0, 0, 2));

    EXPECT_VECTOR_EQ(target_innovation, mean_measurement.calculate_innovation(test_measurement));
}

TEST(FixedMeasurementVectorTest, SigmaPointCovariance) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.set_field<Magnetometer>(UKF::FieldVector(1, 0, 0));

    test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
    test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 0, 0));
    test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_state.set_field<Altitude>(1000);

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution((covariance *
        (MyStateVector::covariance_size() + UKF::Parameters::Lambda<MyStateVector>)).llt().matrixL());

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::SigmaPointDeltas<MyStateVector> sigma_point_deltas =
        mean_measurement.calculate_sigma_point_deltas<MyStateVector>(measurement_sigma_points);
    MyMeasurementVector::CovarianceMatrix calculated_covariance =
        mean_measurement.calculate_sigma_point_covariance<MyStateVector>(sigma_point_deltas);

    MyMeasurementVector::CovarianceMatrix expected_covariance;

    expected_covariance << 5.31709,       0,       0,       0,       0,       0,       0, -2.3059,       0,       0,       0,
                                 0, 5.31709,       0,       0,       0,       0,       0,       0,       0,       0,       0,
                                 0,       0, 29.2440,       0,       0,       0,       0,       0,       0,       0,  -4.237,
                                 0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,
                                 0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,
                                 0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,
                                 0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,
                           -2.3059,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,
                                 0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,
                                 0,       0,       0,       0,       0,       0,       0,       0,       0, 1.44e-4,       0,
                                 0,       0,  -4.237,       0,       0,       0,       0,       0,       0,       0,  32.263;

    EXPECT_VECTOR_EQ(expected_covariance.col(0),  calculated_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_covariance.col(1),  calculated_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_covariance.col(2),  calculated_covariance.col(2));
    EXPECT_VECTOR_EQ(expected_covariance.col(3),  calculated_covariance.col(3));
    EXPECT_VECTOR_EQ(expected_covariance.col(4),  calculated_covariance.col(4));
    EXPECT_VECTOR_EQ(expected_covariance.col(5),  calculated_covariance.col(5));
    EXPECT_VECTOR_EQ(expected_covariance.col(6),  calculated_covariance.col(6));
    EXPECT_VECTOR_EQ(expected_covariance.col(7),  calculated_covariance.col(7));
    EXPECT_VECTOR_EQ(expected_covariance.col(8),  calculated_covariance.col(8));
    EXPECT_VECTOR_EQ(expected_covariance.col(9),  calculated_covariance.col(9));
    EXPECT_VECTOR_EQ(expected_covariance.col(10),  calculated_covariance.col(10));
}

TEST(FixedMeasurementVectorTest, MeasurementCovariance) {
    MyMeasurementVector test_measurement, expected_measurement;
    MyMeasurementVector::CovarianceVector measurement_covariance;

    measurement_covariance.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    measurement_covariance.set_field<DynamicPressure>(4);
    measurement_covariance.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    measurement_covariance.set_field<StaticPressure>(8);
    measurement_covariance.set_field<Magnetometer>(UKF::FieldVector(1, 2, 3));

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<DynamicPressure>(0);
    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    test_measurement.set_field<StaticPressure>(101.3);
    test_measurement.set_field<Magnetometer>(UKF::FieldVector(-1, 0, 1));

    expected_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    expected_measurement.set_field<DynamicPressure>(0);
    expected_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    expected_measurement.set_field<StaticPressure>(101.3);
    expected_measurement.set_field<Magnetometer>(UKF::FieldVector(1, 1, 0));

    MyMeasurementVector::CovarianceMatrix expected_measurement_covariance = MyMeasurementVector::CovarianceMatrix::Zero();
    expected_measurement_covariance.diagonal() << 5, 6, 7, 1, 2, 3, 0, 0, 0, 8, 4;
    expected_measurement_covariance.block<3, 3>(6, 6) <<  8, -8,  0,
                                                         -8,  8,  0,
                                                          0,  0, 16;

    MyMeasurementVector::CovarianceMatrix out_measurement_covariance =
        test_measurement.calculate_measurement_covariance(measurement_covariance, expected_measurement);

    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(0),  out_measurement_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(1),  out_measurement_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(2),  out_measurement_covariance.col(2));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(3),  out_measurement_covariance.col(3));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(4),  out_measurement_covariance.col(4));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(5),  out_measurement_covariance.col(5));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(6),  out_measurement_covariance.col(6));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(7),  out_measurement_covariance.col(7));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(8),  out_measurement_covariance.col(8));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(9),  out_measurement_covariance.col(9));
    EXPECT_VECTOR_EQ(expected_measurement_covariance.col(10),  out_measurement_covariance.col(10));
}

TEST(FixedMeasurementVectorTest, MeasurementRootCovariance) {
    MyMeasurementVector test_measurement, expected_measurement;
    MyMeasurementVector::CovarianceVector measurement_root_covariance;

    measurement_root_covariance.set_field<Gyroscope>(UKF::Vector<3>(1, 2, 3));
    measurement_root_covariance.set_field<DynamicPressure>(4);
    measurement_root_covariance.set_field<Accelerometer>(UKF::Vector<3>(5, 6, 7));
    measurement_root_covariance.set_field<StaticPressure>(8);
    measurement_root_covariance.set_field<Magnetometer>(UKF::FieldVector(1, 2, 3));

    test_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    test_measurement.set_field<DynamicPressure>(0);
    test_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    test_measurement.set_field<StaticPressure>(101.3);
    test_measurement.set_field<Magnetometer>(UKF::FieldVector(-1, 0, 1));

    expected_measurement.set_field<Gyroscope>(UKF::Vector<3>(0, 0, 0));
    expected_measurement.set_field<DynamicPressure>(0);
    expected_measurement.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    expected_measurement.set_field<StaticPressure>(101.3);
    expected_measurement.set_field<Magnetometer>(UKF::FieldVector(1, 1, 0));

    MyMeasurementVector::CovarianceMatrix expected_measurement_root_covariance = MyMeasurementVector::CovarianceMatrix::Zero();
    expected_measurement_root_covariance.diagonal() << 5, 6, 7, 1, 2, 3, 0, 0, 0, 8, 4;
    expected_measurement_root_covariance.block<3, 3>(6, 6) <<  0, -4,  0,
                                                               0,  4,  0,
                                                              -2,  0, -6;

    MyMeasurementVector::CovarianceMatrix out_measurement_root_covariance =
        test_measurement.calculate_measurement_root_covariance(measurement_root_covariance, expected_measurement);

    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(0),  out_measurement_root_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(1),  out_measurement_root_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(2),  out_measurement_root_covariance.col(2));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(3),  out_measurement_root_covariance.col(3));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(4),  out_measurement_root_covariance.col(4));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(5),  out_measurement_root_covariance.col(5));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(6),  out_measurement_root_covariance.col(6));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(7),  out_measurement_root_covariance.col(7));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(8),  out_measurement_root_covariance.col(8));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(9),  out_measurement_root_covariance.col(9));
    EXPECT_VECTOR_EQ(expected_measurement_root_covariance.col(10),  out_measurement_root_covariance.col(10));
}

TEST(FixedMeasurementVectorTest, CalculateRotationVector) {
    UKF::Vector<3> test;

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(0, 0, 0), UKF::Vector<3>(0, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(0, 0, 0), UKF::Vector<3>(1, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(1, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(-1, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, -2/std::numeric_limits<real_t>::epsilon()), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 1, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, -2), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, -1, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 2), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 0, 1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 2, 0), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 0, -1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, -2, 0), test);

    /* Some random vectors, for the a = 0, f = 2 case. */
    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(
        UKF::Vector<3>(5.0746, -2.3911, 1.3564), UKF::Vector<3>(8.3439, -4.2832, 5.1440));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0.1070, 0.2438, 0.0294), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(
        UKF::Vector<3>(5.5833, 8.6802, -7.4019), UKF::Vector<3>(-8.4829, -8.9210, 0.6160));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4.4643, -4.3661, -1.7527), test);

    test = UKF::Detail::calculate_rotation_vector<MyMeasurementVector>(
        UKF::Vector<3>(2.0396, -4.7406, 3.0816), UKF::Vector<3>(-3.7757, 0.5707, -6.6870));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(-3.9209, -0.2624, 2.1914), test);
}

TEST(FixedMeasurementVectorTest, CalculateRotationVectorJacobian) {
    UKF::Matrix<3, 3> test;

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(0, 0, 0), UKF::Vector<3>(0, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(0, 0, 0), UKF::Vector<3>(1, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 4/(std::numeric_limits<real_t>::epsilon()*std::numeric_limits<real_t>::epsilon()), 0).transpose(),
            (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 4/(std::numeric_limits<real_t>::epsilon()*std::numeric_limits<real_t>::epsilon())).transpose(),
            (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(1, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 1, 0).transpose(), (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 1).transpose(), (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(-1, 0, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 4/(std::numeric_limits<real_t>::epsilon()*std::numeric_limits<real_t>::epsilon()), 0).transpose(),
            (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 4/(std::numeric_limits<real_t>::epsilon()*std::numeric_limits<real_t>::epsilon())).transpose(),
            (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 1, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 4).transpose(), (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, -1, 0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 4).transpose(), (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 0, 1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 4, 0).transpose(), (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(UKF::Vector<3>(1, 0, 0), UKF::Vector<3>(0, 0, -1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 0, 0).transpose(), (test*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 4, 0).transpose(), (test*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(0, 0, 0).transpose(), (test*test.transpose()).row(2));

    /* Some random vectors, for the a = 0, f = 2 case. */
    Eigen::DiagonalMatrix<real_t, 3> c(UKF::Vector<3>(1, 2, 3));
    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(
        UKF::Vector<3>(5.0746, -2.3911, 1.3564), UKF::Vector<3>(8.3439, -4.2832, 5.1440));
    EXPECT_VECTOR_EQ(UKF::Vector<3>( 0.030105,  0.032038, -0.022155).transpose(), (test*c*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>( 0.032038,  0.073227,  0.009005).transpose(), (test*c*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(-0.022155,  0.009005,  0.043436).transpose(), (test*c*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(
        UKF::Vector<3>(5.5833, 8.6802, -7.4019), UKF::Vector<3>(-8.4829, -8.9210, 0.6160));
    EXPECT_VECTOR_EQ(UKF::Vector<3>( 0.790494, -0.776049, -0.353005).transpose(), (test*c*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(-0.776049,  0.768788,  0.446778).transpose(), (test*c*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(-0.353005,  0.446778,  1.609094).transpose(), (test*c*test.transpose()).row(2));

    test = UKF::Detail::calculate_rotation_vector_jacobian<MyMeasurementVector>(
        UKF::Vector<3>(2.0396, -4.7406, 3.0816), UKF::Vector<3>(-3.7757, 0.5707, -6.6870));
    EXPECT_VECTOR_EQ(UKF::Vector<3>( 1.850606, -0.474603, -1.085418).transpose(), (test*c*test.transpose()).row(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(-0.474603,  1.420474,  0.389206).transpose(), (test*c*test.transpose()).row(1));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(-1.085418,  0.389206,  0.646079).transpose(), (test*c*test.transpose()).row(2));
}
