#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Types.h"
#include "MeasurementVector.h"
#include "comparisons.h"

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

TEST(FixedMeasurementVectorTest, Instantiation) {
    MyMeasurementVector test_measurement;

    EXPECT_EQ(8, MyMeasurementVector::MaxRowsAtCompileTime);
    EXPECT_EQ(8, test_measurement.size());
}

TEST(FixedMeasurementVectorTest, Assignment) {
    MyMeasurementVector test_measurement;

    test_measurement.field<Gyroscope>() << 1, 2, 3;
    test_measurement.field<DynamicPressure>() << 4;
    test_measurement.field<Accelerometer>() << 5, 6, 7;
    test_measurement.field<StaticPressure>() << 8;

    EXPECT_EQ(8, test_measurement.size());

    EXPECT_EQ(8, test_measurement.field<StaticPressure>()(0));
    EXPECT_EQ(4, test_measurement.field<DynamicPressure>()(0));
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.field<Gyroscope>());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(5, 6, 7), test_measurement.field<Accelerometer>());

    UKF::Vector<8> expected;
    expected << 5, 6, 7, 1, 2, 3, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(FixedMeasurementVectorTest, Reassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.field<Gyroscope>() << 1, 2, 3;
    test_measurement.field<DynamicPressure>() << 4;
    test_measurement.field<Accelerometer>() << 5, 6, 7;
    test_measurement.field<StaticPressure>() << 8;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_measurement.field<Gyroscope>());

    test_measurement.field<Gyroscope>() << 4, 5, 6;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.field<Gyroscope>());

    UKF::Vector<8> expected;
    expected << 5, 6, 7, 4, 5, 6, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(FixedMeasurementVectorTest, MultipleReassignment) {
    MyMeasurementVector test_measurement;

    test_measurement.field<Gyroscope>() << 1, 2, 3;
    test_measurement.field<DynamicPressure>() << 4;
    test_measurement.field<Accelerometer>() << 5, 6, 7;
    test_measurement.field<StaticPressure>() << 8;

    EXPECT_EQ(8, test_measurement.size());
    UKF::Vector<8> expected;
    expected << 5, 6, 7, 1, 2, 3, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<Gyroscope>() << 4, 5, 6;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(4, 5, 6), test_measurement.field<Gyroscope>());
    expected << 5, 6, 7, 4, 5, 6, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<Accelerometer>() << 7, 8, 9;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(7, 8, 9), test_measurement.field<Accelerometer>());
    expected << 7, 8, 9, 4, 5, 6, 8, 4;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<DynamicPressure>() << 1;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(1, test_measurement.field<DynamicPressure>()(0));
    expected << 7, 8, 9, 4, 5, 6, 8, 1;
    EXPECT_VECTOR_EQ(expected, test_measurement);

    test_measurement.field<StaticPressure>() << 3;

    EXPECT_EQ(8, test_measurement.size());
    EXPECT_EQ(3, test_measurement.field<StaticPressure>()(0));
    expected << 7, 8, 9, 4, 5, 6, 3, 1;
    EXPECT_VECTOR_EQ(expected, test_measurement);
}

TEST(FixedMeasurementVectorTest, Arithmetic) {
    
}
