#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "MeasurementVector.h"
#include "comparisons.h"

enum MyFields {
    StaticPressure,
    DynamicPressure,
    Accelerometer,
    Gyroscope
};

using MyMeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<Accelerometer, Eigen::Vector3d>,
    UKF::Field<Gyroscope, Eigen::Vector3d>,
    UKF::Field<StaticPressure, real_t>,
    UKF::Field<DynamicPressure, real_t>
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
    EXPECT_VECTOR_EQ(Eigen::Vector3d(1, 2, 3), test_measurement.field<Gyroscope>());

    test_measurement.field<DynamicPressure>() << 4;

    EXPECT_EQ(4, test_measurement.size());
    EXPECT_EQ(4, test_measurement.field<DynamicPressure>()(0));

    test_measurement.field<Accelerometer>() << 5, 6, 7;

    EXPECT_EQ(7, test_measurement.size());
    EXPECT_VECTOR_EQ(Eigen::Vector3d(5, 6, 7), test_measurement.field<Accelerometer>());

    test_measurement.field<StaticPressure>() << 8;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(8, test_measurement.field<StaticPressure>()(0));

    Eigen::Matrix<double, 8, 1> expected;
    expected << 1, 2, 3, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(DynamicMeasurementVectorTest, Reassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.field<Gyroscope>() << 1, 2, 3;

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(Eigen::Vector3d(1, 2, 3), test_measurement);
    EXPECT_VECTOR_EQ(Eigen::Vector3d(1, 2, 3), test_measurement.field<Gyroscope>());

    test_measurement.field<Gyroscope>() << 4, 5, 6;

    EXPECT_EQ(3, test_measurement.size());
    EXPECT_VECTOR_EQ(Eigen::Vector3d(4, 5, 6), test_measurement);
    EXPECT_VECTOR_EQ(Eigen::Vector3d(4, 5, 6), test_measurement.field<Gyroscope>());
}

TEST(DynamicMeasurementVectorTest, MultipleReassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.field<Gyroscope>() << 1, 2, 3;
    test_measurement.field<DynamicPressure>() << 4;
    test_measurement.field<Accelerometer>() << 5, 6, 7;
    test_measurement.field<StaticPressure>() << 8;

    EXPECT_EQ(8, test_measurement.size());
    Eigen::Matrix<double, 8, 1> expected;
    expected << 1, 2, 3, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<Gyroscope>() << 4, 5, 6;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(Eigen::Vector3d(4, 5, 6), test_measurement.field<Gyroscope>());
    expected << 4, 5, 6, 4, 5, 6, 7, 8;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<Accelerometer>() << 7, 8, 9;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(Eigen::Vector3d(7, 8, 9), test_measurement.field<Accelerometer>());
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
