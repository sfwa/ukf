#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "types.h"
#include "state.h"
#include "integrator.h"
#include <iostream>

TEST(IntegratorRK4Test, Instantiation) {
    IntegratorRK4 test = IntegratorRK4();
}

TEST(IntegratorRK4Test, ConstantAcceleration) {
    IntegratorRK4 test = IntegratorRK4();
    State test_state;
    test_state << 0, 0, 0,
                  3, 4, 5,
                  6, 7, 8,
                  0, 0, 0, 1,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.01) {
        test_state = test.integrate(test_state, 0.01);
    }

    EXPECT_TRUE(test_state.position().isApprox(
        Eigen::Matrix<real_t, 3, 1>(5.20898e-05, 6.11485e-05, -450), 1e-8));
    EXPECT_TRUE(test_state.velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(63, 74, 85), 0.001));
}

TEST(IntegratorRK4Test, ConstantAngularAcceleration) {
    IntegratorRK4 test = IntegratorRK4();
    State test_state;
    test_state << 0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 0, 0,
                  0, 0, 1,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.01) {
        test_state = test.integrate(test_state, 0.01);
    }

    EXPECT_TRUE(test_state.attitude().isApprox(
        Eigen::Matrix<real_t, 4, 1>(0, 0, 0.1392, 0.9903), 0.01));
    

    EXPECT_TRUE(test_state.angular_velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(0, 0, 10), 0.001));
}

TEST(IntegratorRK4Test, CircularMotion) {
    IntegratorRK4 test = IntegratorRK4();
    State test_state;
    test_state << 0, 0, 0,
                  10, 0, 0,
                  0, 6.283, 0,
                  0, 0, 0, 1,
                  0, 0, 0.6283,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.01) {
        test_state = test.integrate(test_state, 0.01);
    }

    EXPECT_TRUE((test_state.position() - Vector3r(0, 0, 0)).isZero(1e-8));
    EXPECT_TRUE(test_state.velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(10, 0, 0), 0.001));
    EXPECT_TRUE(test_state.acceleration().isApprox(
        Eigen::Matrix<real_t, 3, 1>(0, 6.283, 0), 0.001));
    Quaternionr attitude = Quaternionr(test_state.attitude());
    EXPECT_TRUE(attitude.matrix().isApprox(
        Quaternionr(1, 0, 0, 0).matrix(), 0.001));
    EXPECT_TRUE(test_state.angular_velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(0, 0, 0.6283), 0.001));
    EXPECT_TRUE(test_state.angular_acceleration().isZero(0.001));
}

TEST(IntegratorHeunTest, Instantiation) {
    IntegratorHeun test = IntegratorHeun();
}

TEST(IntegratorHeunTest, ConstantAcceleration) {
    IntegratorHeun test = IntegratorHeun();
    State test_state;
    test_state << 0, 0, 0,
                  3, 4, 5,
                  6, 7, 8,
                  0, 0, 0, 1,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.01) {
        test_state = test.integrate(test_state, 0.01);
    }

    EXPECT_TRUE(test_state.position().isApprox(
        Eigen::Matrix<real_t, 3, 1>(5.20898e-05, 6.11485e-05, -450), 1e-8));
    EXPECT_TRUE(test_state.velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(63, 74, 85), 0.001));
}

TEST(IntegratorHeunTest, ConstantAngularAcceleration) {
    IntegratorHeun test = IntegratorHeun();
    State test_state;
    test_state << 0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 0, 0,
                  0, 0, 1,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.01) {
        test_state = test.integrate(test_state, 0.01);
    }

    EXPECT_TRUE(test_state.attitude().isApprox(
        Eigen::Matrix<real_t, 4, 1>(0, 0, 0.1392, 0.9903), 0.1));

    EXPECT_TRUE(test_state.angular_velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(0, 0, 10), 0.001));
}

TEST(IntegratorEulerTest, Instantiation) {
    IntegratorEuler test = IntegratorEuler();
}

TEST(IntegratorEulerTest, ConstantAcceleration) {
    IntegratorEuler test = IntegratorEuler();
    State test_state;
    test_state << 0, 0, 0,
                  3, 4, 5,
                  6, 7, 8,
                  0, 0, 0, 1,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.01) {
        test_state = test.integrate(test_state, 0.01);
    }

    EXPECT_TRUE(test_state.position().isApprox(
        Eigen::Matrix<real_t, 3, 1>(5.20898e-05, 6.11485e-05, -450), 0.001));
    EXPECT_TRUE(test_state.velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(63, 74, 85), 0.001));
}

TEST(IntegratorEulerTest, ConstantAngularAcceleration) {
    IntegratorEuler test = IntegratorEuler();
    State test_state;
    test_state << 0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 0, 0,
                  0, 0, 1,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.001) {
        test_state = test.integrate(test_state, 0.001);
    }

    EXPECT_TRUE(test_state.attitude().isApprox(
        Eigen::Matrix<real_t, 4, 1>(0, 0, 0.1392, 0.9903), 0.1));

    EXPECT_TRUE(test_state.angular_velocity().isApprox(
        Eigen::Matrix<real_t, 3, 1>(0, 0, 10), 0.001));
}
