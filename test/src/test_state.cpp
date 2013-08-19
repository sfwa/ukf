#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "types.h"
#include "state.h"

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
    EXPECT_EQ(test(0), 0);
    EXPECT_EQ(test(1), 1);
    EXPECT_EQ(test(2), 2);
    EXPECT_EQ(test(3), 3);
    EXPECT_EQ(test(4), 4);
    EXPECT_EQ(test(5), 5);
    EXPECT_EQ(test(6), 6);
    EXPECT_EQ(test(7), 7);
    EXPECT_EQ(test(8), 8);
    EXPECT_EQ(test(9), 9);
    EXPECT_EQ(test(10), 10);
    EXPECT_EQ(test(11), 11);
    EXPECT_EQ(test(12), 12);
    EXPECT_EQ(test(13), 13);
    EXPECT_EQ(test(14), 14);
    EXPECT_EQ(test(15), 15);
    EXPECT_EQ(test(16), 16);
    EXPECT_EQ(test(17), 17);
    EXPECT_EQ(test(18), 18);
    EXPECT_EQ(test(19), 19);
    EXPECT_EQ(test(20), 20);
    EXPECT_EQ(test(21), 21);
    EXPECT_EQ(test(22), 22);
    EXPECT_EQ(test(23), 23);
    EXPECT_EQ(test(24), 24);
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

    EXPECT_EQ(test.position()(0), 0);
    EXPECT_EQ(test.position()(1), 1);
    EXPECT_EQ(test.position()(2), 2);

    test.position() *= 2;
    EXPECT_EQ(test(0), 0);
    EXPECT_EQ(test(1), 2);
    EXPECT_EQ(test(2), 4);
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

    EXPECT_EQ(test.velocity()(0), 3);
    EXPECT_EQ(test.velocity()(1), 4);
    EXPECT_EQ(test.velocity()(2), 5);

    test.velocity() *= 2;
    EXPECT_EQ(test(3), 6);
    EXPECT_EQ(test(4), 8);
    EXPECT_EQ(test(5), 10);
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

    EXPECT_EQ(test.acceleration()(0), 6);
    EXPECT_EQ(test.acceleration()(1), 7);
    EXPECT_EQ(test.acceleration()(2), 8);

    test.acceleration() *= 2;
    EXPECT_EQ(test(6), 12);
    EXPECT_EQ(test(7), 14);
    EXPECT_EQ(test(8), 16);
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

    EXPECT_EQ(test.attitude()(0), 9);
    EXPECT_EQ(test.attitude()(1), 10);
    EXPECT_EQ(test.attitude()(2), 11);
    EXPECT_EQ(test.attitude()(3), 12);

    test.attitude() << a.vec(), a.w();
    EXPECT_EQ(test(9), 0);
    EXPECT_EQ(test(10), 0);
    EXPECT_EQ(test(11), 0);
    EXPECT_EQ(test(12), 1);
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

    EXPECT_EQ(test.angular_velocity()(0), 13);
    EXPECT_EQ(test.angular_velocity()(1), 14);
    EXPECT_EQ(test.angular_velocity()(2), 15);

    test.angular_velocity() *= 2;
    EXPECT_EQ(test(13), 26);
    EXPECT_EQ(test(14), 28);
    EXPECT_EQ(test(15), 30);
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

    EXPECT_EQ(test.angular_acceleration()(0), 16);
    EXPECT_EQ(test.angular_acceleration()(1), 17);
    EXPECT_EQ(test.angular_acceleration()(2), 18);

    test.angular_acceleration() *= 2;
    EXPECT_EQ(test(16), 32);
    EXPECT_EQ(test(17), 34);
    EXPECT_EQ(test(18), 36);
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

    EXPECT_EQ(test.wind_velocity()(0), 19);
    EXPECT_EQ(test.wind_velocity()(1), 20);
    EXPECT_EQ(test.wind_velocity()(2), 21);

    test.wind_velocity() *= 2;
    EXPECT_EQ(test(19), 38);
    EXPECT_EQ(test(20), 40);
    EXPECT_EQ(test(21), 42);
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

    EXPECT_EQ(test.gyro_bias()(0), 22);
    EXPECT_EQ(test.gyro_bias()(1), 23);
    EXPECT_EQ(test.gyro_bias()(2), 24);

    test.gyro_bias() *= 2;
    EXPECT_EQ(test(22), 44);
    EXPECT_EQ(test(23), 46);
    EXPECT_EQ(test(24), 48);
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

    EXPECT_FLOAT_EQ(derivatives(0), 0);
    EXPECT_FLOAT_EQ(derivatives(1), 1.5678559e-07);
    EXPECT_EQ(derivatives(2), -2);

    EXPECT_EQ(derivatives(3), 3);
    EXPECT_EQ(derivatives(4), 4);
    EXPECT_EQ(derivatives(5), 5);

    EXPECT_EQ(derivatives(9), -0.5);
    EXPECT_EQ(derivatives(10), 0);
    EXPECT_EQ(derivatives(11), 0);
    EXPECT_EQ(derivatives(12), 0);

    EXPECT_EQ(derivatives(13), 6);
    EXPECT_EQ(derivatives(14), 7);
    EXPECT_EQ(derivatives(15), 8);
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

    EXPECT_FLOAT_EQ(derivatives(0), 0);
    EXPECT_FLOAT_EQ(derivatives(1), 1.5678559e-07);
    EXPECT_EQ(derivatives(2), -2);

    EXPECT_TRUE(derivatives.segment<3>(3).isApprox(
        Eigen::Matrix<real_t, 3, 1>(3, 5, -4), 0.01));

    EXPECT_TRUE(derivatives.segment<4>(9).isApprox(
        Eigen::Matrix<real_t, 4, 1>(-0.3535, 0, 0, 0.3535)));

    EXPECT_EQ(derivatives(13), 6);
    EXPECT_EQ(derivatives(14), 7);
    EXPECT_EQ(derivatives(15), 8);
}
