#include <gtest/gtest.h>
#include "types.h"
#include "state.h"
#include "sensors.h"
#include "comparisons.h"

TEST(IOBoardModelTest, Instantiation) {
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));

    EXPECT_EQ(MeasurementVector(), test.collate());
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

    EXPECT_EQ(target, test.collate());
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

    EXPECT_EQ(target, test.collate());
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

    EXPECT_EQ(target, test.collate());
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

    EXPECT_EQ(target, test.collate());
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

    EXPECT_EQ(target, test.collate());
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

    EXPECT_EQ(target, test.collate());
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

    EXPECT_EQ(target, test.collate());
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
    test.set_barometer_amsl(7);
    target << 1, 2, 3, 4, 5, 6, 7;

    EXPECT_EQ(target, test.collate());
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
    test.set_barometer_amsl(10);
    test.set_gps_velocity(Vector3r(7, 8, 9));
    target << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    EXPECT_EQ(target, test.collate());
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
    test.set_barometer_amsl(7);
    test.clear();

    EXPECT_EQ(MeasurementVector(), test.collate());
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
    MeasurementVector target(6);

    test.set_accelerometer(Vector3r(0, 0, 0));

    target << 1, 2, 3, 0, G_ACCEL, 0;
    EXPECT_MEASUREMENT_EQ(target, test.predict(test_state));

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
    EXPECT_MEASUREMENT_EQ(target, test.predict(test_state));
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

    EXPECT_MEASUREMENT_EQ(Vector3r(4, 5, 6), test.predict(test_state));
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

    EXPECT_MEASUREMENT_EQ(Vector3r(7, -9, 8), test.predict(test_state));
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

    EXPECT_MEASUREMENT_EQ(Vector3r(-37, 145, 0), test.predict(test_state));
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

    EXPECT_MEASUREMENT_EQ(Vector3r(7, 8, 9), test.predict(test_state));
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

    EXPECT_EQ(6, test.predict(test_state)(0));
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

    EXPECT_GE(0.1, std::abs(test.predict(test_state)(0) - 150));
}
