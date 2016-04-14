#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Types.h"
#include "Integrator.h"
#include "StateVector.h"
#include "MeasurementVector.h"
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

    test_measurement.field<Gyroscope>() << 1, 2, 3;

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.field<Gyroscope>());

    test_measurement.field<DynamicPressure>() << 4;

    EXPECT_EQ(4, test_measurement.size());
    EXPECT_EQ(4, test_measurement.field<DynamicPressure>()(0));

    test_measurement.field<Accelerometer>() << 5, 6, 7;

    EXPECT_EQ(7, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(5, 6, 7), test_measurement.field<Accelerometer>());

    test_measurement.field<StaticPressure>() << 8;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(8, test_measurement.field<StaticPressure>()(0));

    UKF::Vector<8> expected;
    expected << 1, 2, 3, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(DynamicMeasurementVectorTest, Reassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.field<Gyroscope>() << 1, 2, 3;

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement);
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.field<Gyroscope>());

    test_measurement.field<Gyroscope>() << 4, 5, 6;

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement);
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.field<Gyroscope>());
}

TEST(DynamicMeasurementVectorTest, MultipleReassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.field<Gyroscope>() << 1, 2, 3;
    test_measurement.field<DynamicPressure>() << 4;
    test_measurement.field<Accelerometer>() << 5, 6, 7;
    test_measurement.field<StaticPressure>() << 8;

    EXPECT_EQ(8, test_measurement.size());
    UKF::Vector<8> expected;
    expected << 1, 2, 3, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<Gyroscope>() << 4, 5, 6;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.field<Gyroscope>());
    expected << 4, 5, 6, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<Accelerometer>() << 7, 8, 9;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(7, 8, 9), test_measurement.field<Accelerometer>());
    expected << 4, 5, 6, 4, 7, 8, 9, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<DynamicPressure>() << 1;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(1, test_measurement.field<DynamicPressure>()(0));
    expected << 4, 5, 6, 1, 7, 8, 9, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<StaticPressure>() << 3;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(3, test_measurement.field<StaticPressure>()(0));
    expected << 4, 5, 6, 1, 7, 8, 9, 3;
    EXPECT_VECTOR_EQ(expected, test_measurement);
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
    IntegratorRK4,
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
<MyStateVector, UKF::Field<Accelerometer, UKF::Vector<3>>>(const MyStateVector &state) {
    return UKF::Quaternion(state.field<Attitude>()) * UKF::Vector<3>(0, 0, -9.8);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, UKF::Field<Gyroscope, UKF::Vector<3>>>(const MyStateVector &state) {
    return state.field<AngularVelocity>();
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, UKF::Field<StaticPressure, real_t>>(const MyStateVector &state) {
    return 101300.0 - 1200.0*(state.field<Altitude>()(0) / 100.0);
}

template <> template <>
real_t MyMeasurementVector::expected_measurement
<MyStateVector, UKF::Field<DynamicPressure, real_t>>(const MyStateVector &state) {
    return 0.5 * 1.225 * state.field<Velocity>().squaredNorm();
}

TEST(DynamicMeasurementVectorTest, SigmaPointGeneration) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.field<Accelerometer>() << 0, 0, 0;
    test_measurement.field<Gyroscope>() << 0, 0, 0;
    test_measurement.field<StaticPressure>() << 0;
    test_measurement.field<DynamicPressure>() << 0;

    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<AngularVelocity>() << 1, 0, 0;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 1000;

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points, target_sigma_points;

    target_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,
                            0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,      0,
                         -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  9.590,  9.590,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  9.590,  9.590,   -9.8,   -9.8,
                            1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                            0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                            0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                        89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89257,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89343,
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
}

TEST(DynamicMeasurementVectorTest, PartialSigmaPointGeneration) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.field<Accelerometer>() << 0, 0, 0;

    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<AngularVelocity>() << 1, 0, 0;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 1000;

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points, target_sigma_points;

    target_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,
                            0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,      0,
                         -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  9.590,  9.590,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  9.590,  9.590,   -9.8,   -9.8;;
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
}

TEST(DynamicMeasurementVectorTest, SigmaPointMean) {
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points;
    MyMeasurementVector test_measurement;

    test_measurement.field<Accelerometer>() << 0, 0, 0;
    test_measurement.field<Gyroscope>() << 0, 0, 0;
    test_measurement.field<StaticPressure>() << 0;
    test_measurement.field<DynamicPressure>() << 0;

    measurement_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,
                                 0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,      0,
                              -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  9.590,  9.590,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8, -9.590, -9.590,   -9.8,   -9.8,
                                 1,      1,      1,      1,  4.606,      1,      1,      1,      1,      1,      1,      1,      1,      1, -2.606,      1,      1,      1,      1,      1,      1,
                                 0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,      0,
                                 0,      0,      0,      0,      0,      0,  3.606,      0,      0,      0,      0,      0,      0,      0,      0,      0, -3.606,      0,      0,      0,      0,
                             89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89257,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89300,  89343,
                             8.575, 20.954, 25.371, 29.788,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575, 12.121,  7.704,  3.287,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575,  8.575;

    MyMeasurementVector expected_mean;

    expected_mean << 0.0, 0.0, -9.8, 1.0, 0.0, 0.0, 89300, 10.4125;

    EXPECT_VECTOR_EQ(expected_mean, test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points));
}

TEST(DynamicMeasurementVectorTest, PartialSigmaPointMean) {
    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points;
    MyMeasurementVector test_measurement;

    test_measurement.field<Accelerometer>() << 0, 0, 0;

    measurement_sigma_points <<  0,      0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,
                                 0,      0,      0,      0,      0,      0,      0, -2.017,      0,      0,      0,      0,      0,      0,      0,      0,      0,  2.017,      0,      0,      0,
                              -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,  9.590,  9.590,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8,   -9.8, -9.590, -9.590,   -9.8,   -9.8;

    MyMeasurementVector expected_mean;

    expected_mean << 0.0, 0.0, -9.8;

    EXPECT_VECTOR_EQ(expected_mean, test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points));
}

TEST(DynamicMeasurementVectorTest, SigmaPointCovariance) {
    MyStateVector test_state;
    MyMeasurementVector test_measurement;

    test_measurement.field<Accelerometer>() << 0, 0, 0;
    test_measurement.field<Gyroscope>() << 0, 0, 0;
    test_measurement.field<StaticPressure>() << 0;
    test_measurement.field<DynamicPressure>() << 0;

    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<AngularVelocity>() << 1, 0, 0;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 1000;

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);

    MyMeasurementVector::CovarianceMatrix calculated_covariance =
        mean_measurement.calculate_sigma_point_covariance<MyStateVector>(measurement_sigma_points);

    MyMeasurementVector::CovarianceMatrix expected_covariance;

    expected_covariance << 0.31285,       0,       0,       0,       0,       0,       0,       0,
                                 0, 0.31285,       0,       0,       0,       0,       0,       0,
                                 0,       0, 48.9444,       0,       0,       0,       0,  -5.481,
                                 0,       0,       0,       1,       0,       0,       0,       0,
                                 0,       0,       0,       0,       1,       0,       0,       0,
                                 0,       0,       0,       0,       0,       1,       0,       0,
                                 0,       0,       0,       0,       0,       0,     144,       0,
                                 0,       0,  -5.481,       0,       0,       0,       0,  32.263;

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

    test_measurement.field<Accelerometer>() << 0, 0, 0;

    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<AngularVelocity>() << 1, 0, 0;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 1000;

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    MyMeasurementVector::SigmaPointDistribution<MyStateVector> measurement_sigma_points =
        test_measurement.calculate_sigma_point_distribution<MyStateVector>(sigma_points);

    MyMeasurementVector mean_measurement = test_measurement.calculate_sigma_point_mean<MyStateVector>(measurement_sigma_points);

    MyMeasurementVector::CovarianceMatrix calculated_covariance =
        mean_measurement.calculate_sigma_point_covariance<MyStateVector>(measurement_sigma_points);

    MyMeasurementVector::CovarianceMatrix expected_covariance;

    expected_covariance << 0.31285,       0,       0,
                                 0, 0.31285,       0,
                                 0,       0, 48.9444;

    EXPECT_VECTOR_EQ(expected_covariance.col(0),  calculated_covariance.col(0));
    EXPECT_VECTOR_EQ(expected_covariance.col(1),  calculated_covariance.col(1));
    EXPECT_VECTOR_EQ(expected_covariance.col(2),  calculated_covariance.col(2));
}
