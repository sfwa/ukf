#include <gtest/gtest.h>
#include "types.h"
#include "state.h"
#include "comparisons.h"

TEST(StateTest, Instantiation) {
    State test = State();
}

TEST(StateTest, Initialisation) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;
    EXPECT_EQ(0, test[0]);
    EXPECT_EQ(1, test[1]);
    EXPECT_EQ(2, test[2]);
    EXPECT_EQ(3, test[3]);
    EXPECT_EQ(4, test[4]);
    EXPECT_EQ(5, test[5]);
    EXPECT_EQ(6, test[6]);
    EXPECT_EQ(7, test[7]);
    EXPECT_EQ(8, test[8]);
    EXPECT_EQ(9, test[9]);
    EXPECT_EQ(10, test[10]);
    EXPECT_EQ(11, test[11]);
    EXPECT_EQ(12, test[12]);
    EXPECT_EQ(13, test[13]);
    EXPECT_EQ(14, test[14]);
    EXPECT_EQ(15, test[15]);
    EXPECT_EQ(16, test[16]);
    EXPECT_EQ(17, test[17]);
    EXPECT_EQ(18, test[18]);
    EXPECT_EQ(19, test[19]);
    EXPECT_EQ(20, test[20]);
    EXPECT_EQ(21, test[21]);
    EXPECT_EQ(22, test[22]);
    EXPECT_EQ(23, test[23]);
    EXPECT_EQ(24, test[24]);
}

TEST(StateTest, PositionAccessor) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_VECTOR_EQ(Vector3r(0, 1, 2), test.position());

    test.position() *= 2;
    EXPECT_VECTOR_EQ(Vector3r(0, 2, 4), Vector3r(test[0], test[1], test[2]));
}

TEST(StateTest, VelocityAccessor) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_VECTOR_EQ(Vector3r(3, 4, 5), test.velocity());

    test.velocity() *= 2;
    EXPECT_VECTOR_EQ(Vector3r(6, 8, 10), Vector3r(test[3], test[4], test[5]));
}

TEST(StateTest, AccelerationAccessor) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_VECTOR_EQ(Vector3r(6, 7, 8), test.acceleration());

    test.acceleration() *= 2;
    EXPECT_VECTOR_EQ(Vector3r(12, 14, 16),
        Vector3r(test[6], test[7], test[8]));
}

TEST(StateTest, AttitudeAccessor) {
    State test;
    Quaternionr a = Quaternionr(1, 0, 0, 0);
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_EQ(9, test.attitude()[0]);
    EXPECT_EQ(10, test.attitude()[1]);
    EXPECT_EQ(11, test.attitude()[2]);
    EXPECT_EQ(12, test.attitude()[3]);

    test.attitude() << a.vec(), a.w();
    EXPECT_EQ(0, test[9]);
    EXPECT_EQ(0, test[10]);
    EXPECT_EQ(0, test[11]);
    EXPECT_EQ(1, test[12]);
}

TEST(StateTest, AngularVelocityAccessor) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_VECTOR_EQ(Vector3r(13, 14, 15), test.angular_velocity());

    test.angular_velocity() *= 2;
    EXPECT_VECTOR_EQ(Vector3r(26, 28, 30),
        Vector3r(test[13], test[14], test[15]));
}

TEST(StateTest, AngularAccelerationAccessor) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_VECTOR_EQ(Vector3r(16, 17, 18), test.angular_acceleration());

    test.angular_acceleration() *= 2;
    EXPECT_VECTOR_EQ(Vector3r(32, 34, 36),
        Vector3r(test[16], test[17], test[18]));
}

TEST(StateTest, WindVelocityAccessor) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_VECTOR_EQ(Vector3r(19, 20, 21), test.wind_velocity());

    test.wind_velocity() *= 2;
    EXPECT_VECTOR_EQ(Vector3r(38, 40, 42),
        Vector3r(test[19], test[20], test[21]));
}

TEST(StateTest, GyroBiasAccessor) {
    State test;
    test << 0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15,
            16, 17, 18,
            19, 20, 21,
            22, 23, 24;

    EXPECT_VECTOR_EQ(Vector3r(22, 23, 24), test.gyro_bias());

    test.gyro_bias() *= 2;
    EXPECT_VECTOR_EQ(Vector3r(44, 46, 48),
        Vector3r(test[22], test[23], test[24]));
}


TEST(StateTest, KinematicsModelBasic) {
    State test;
    StateVectorDerivative derivatives;
    test << 0, 0, 0,
            0, 1, 2,
            3, 4, 5,
            0, 0, 0, 1,
            1, 0, 0,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    derivatives = test.model();

    EXPECT_VECTOR_EQ(Vector3r(0, 1.5678559e-07, -2),
        Vector3r(derivatives[0], derivatives[1], derivatives[2]));
    EXPECT_VECTOR_EQ(Vector3r(3, 4, 5),
        Vector3r(derivatives[3], derivatives[4], derivatives[5]));
    EXPECT_VECTOR_EQ(Vector3r(6, 7, 8),
        Vector3r(derivatives[13], derivatives[14], derivatives[15]));

    EXPECT_FLOAT_EQ(-0.5, derivatives[9]);
    EXPECT_FLOAT_EQ(0, derivatives[10]);
    EXPECT_FLOAT_EQ(0, derivatives[11]);
    EXPECT_FLOAT_EQ(0, derivatives[12]);
}

TEST(StateTest, KinematicsModelRotated) {
    State test;
    StateVectorDerivative derivatives;
    test << 0, 0, 0,
            0, 1, 2,
            3, 4, 5,
            0.707, 0, 0, 0.707,
            1, 0, 0,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    derivatives = test.model();

    EXPECT_VECTOR_EQ(Vector3r(0, 1.5678559e-07, -2),
        Vector3r(derivatives[0], derivatives[1], derivatives[2]));
    EXPECT_VECTOR_EQ(Vector3r(3, 5, -4),
        Vector3r(derivatives[3], derivatives[4], derivatives[5]));

    EXPECT_QUATERNION_EQ(
        Quaternionr(derivatives[12], derivatives[9], derivatives[10],
            derivatives[11]),
        Quaternionr(0.3535, -0.3535, 0, 0));

    EXPECT_VECTOR_EQ(Vector3r(6, 7, 8),
        Vector3r(derivatives[13], derivatives[14], derivatives[15]));
}
