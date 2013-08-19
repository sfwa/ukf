#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "types.h"
#include "state.h"
#include "sensors.h"

TEST(IOBoardModelTest, Instantiation) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));

    EXPECT_EQ(test.collate(), MeasurementVector());
}

TEST(IOBoardModelTest, SetAccelerometer) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(3);

    test.set_accelerometer(Vector3r(1, 2, 3));
    target << 1, 2, 3;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetGyroscope) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(3);

    test.set_gyroscope(Vector3r(1, 2, 3));
    target << 1, 2, 3;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetMagnetometer) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(3);

    test.set_magnetometer(Vector3r(1, 2, 3));
    target << 1, 2, 3;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetGPSPosition) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(3);

    test.set_gps_position(Vector3r(1, 2, 3));
    target << 1, 2, 3;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetGPSVelocity) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(3);

    test.set_gps_velocity(Vector3r(1, 2, 3));
    target << 1, 2, 3;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetPitotTAS) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(1);

    test.set_pitot_tas(1);
    target << 1;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetBarometerAMSL) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(1);

    test.set_barometer_amsl(1);
    target << 1;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetBarometerROC) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(1);

    test.set_barometer_roc(1);
    target << 1;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetMultipleInOrder) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(7);

    test.set_accelerometer(Vector3r(1, 2, 3));
    test.set_gyroscope(Vector3r(4, 5, 6));
    test.set_barometer_roc(7);
    target << 1, 2, 3, 4, 5, 6, 7;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetMultipleOutOfOrder) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    MeasurementVector target(10);

    test.set_magnetometer(Vector3r(4, 5, 6));
    test.set_accelerometer(Vector3r(1, 2, 3));
    test.set_barometer_roc(10);
    test.set_gps_velocity(Vector3r(7, 8, 9));
    target << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    EXPECT_EQ(test.collate(), target);
}

TEST(IOBoardModelTest, SetThenClear) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));

    test.set_accelerometer(Vector3r(1, 2, 3));
    test.set_gyroscope(Vector3r(4, 5, 6));
    test.set_barometer_roc(7);
    test.clear();

    EXPECT_EQ(test.collate(), MeasurementVector());
}

TEST(IOBoardModelTest, PredictAccelerometer) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    State test_state;
    test_state << 0, 0, 0,
            0, 0, 0,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0;
    Eigen::Matrix<real_t, 6, 1> target;

    test.set_accelerometer(Vector3r(0, 0, 0));

    target << 1, 2, 3, 0, G_ACCEL, 0;
    EXPECT_TRUE(test.predict(test_state).isApprox(target, 0.001));

    test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(1, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    test_state << 0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0, 1,
            0, 0, 1,
            0, 0, 1,
            0, 0, 0,
            0, 0, 0;

    test.set_accelerometer(Vector3r(0, 0, 0));

    target << -1, 1, 0, 0, 0, -G_ACCEL;
    EXPECT_TRUE(test.predict(test_state).isApprox(target, 0.001));
}

TEST(IOBoardModelTest, PredictGyroscope) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    State test_state;
    test_state << 0, 0, 0,
            0, 0, 0,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            4, 5, 6,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0;

    test.set_gyroscope(Vector3r(0, 0, 0));

    EXPECT_EQ(test.predict(test_state), Vector3r(4, 5, 6));
}

TEST(IOBoardModelTest, PredictMagnetometer) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(7, 8, 9));
    State test_state;
    test_state << 0, 0, 0,
            0, 0, 0,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            4, 5, 6,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0;

    test.set_magnetometer(Vector3r(0, 0, 0));

    EXPECT_TRUE(test.predict(test_state).isApprox(
        Vector3r(7, -9, 8), 0.001));
}

TEST(IOBoardModelTest, PredictGPSPosition) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(7, 8, 9));
    State test_state;
    test_state << -37, 145, 0,
            0, 0, 0,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            4, 5, 6,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0;

    test.set_gps_position(Vector3r(0, 0, 0));

    EXPECT_TRUE(test.predict(test_state).isApprox(
        Vector3r(-37, 145, 0), 0.001));
}

TEST(IOBoardModelTest, PredictGPSVelocity) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(7, 8, 9));
    State test_state;
    test_state << 0, 0, 0,
            7, 8, 9,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            4, 5, 6,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0;

    test.set_gps_velocity(Vector3r(0, 0, 0));

    EXPECT_TRUE(test.predict(test_state).isApprox(
        Vector3r(7, 8, 9), 0.001));
}

TEST(IOBoardModelTest, PredictPitotTAS) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(7, 8, 9));
    State test_state;
    test_state << 0, 0, 0,
            7, 8, 9,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            4, 5, 6,
            0, 0, 0,
            1, 2, 3,
            0, 0, 0;

    test.set_pitot_tas(0);

    EXPECT_EQ(test.predict(test_state)(0), 6);
}

TEST(IOBoardModelTest, PredictBarometerAMSL) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(7, 8, 9));
    State test_state;
    test_state << -37, 145, 150,
            7, 8, 9,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            4, 5, 6,
            0, 0, 0,
            1, 2, 3,
            0, 0, 0;

    test.set_barometer_amsl(0);

    EXPECT_TRUE(abs(test.predict(test_state)(0) - 150) < 0.1);
}

TEST(IOBoardModelTest, PredictBarometerROC) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(7, 8, 9));
    State test_state;
    test_state << 0, 0, 0,
            0, 0, 10,
            1, 2, 3,
            0.7071, 0, 0, 0.7071,
            4, 5, 6,
            0, 0, 0,
            1, 2, 3,
            0, 0, 0;

    test.set_barometer_roc(0);

    EXPECT_EQ(test.predict(test_state)[0], -10);
}
