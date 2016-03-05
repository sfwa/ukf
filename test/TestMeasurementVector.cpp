#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "MeasurementVector.h"
#include "comparisons.h"

TEST(MeasurementTest, Instantiation) {
    enum MyFields {
        StaticPressure,
        DynamicPressure,
        Accelerometer,
        Gyroscope
    };

    using MyMeasurementVector = UKF::MeasurementVector<
        UKF::Field<Accelerometer, Eigen::Vector3d>,
        UKF::Field<Gyroscope, Eigen::Vector3d>,
        UKF::Field<StaticPressure, real_t>,
        UKF::Field<DynamicPressure, real_t>
    >;

    MyMeasurementVector test_measurement;

    EXPECT_EQ(8, MyMeasurementVector::MaxRowsAtCompileTime);
    EXPECT_EQ(8, test_measurement.max_size());
}
