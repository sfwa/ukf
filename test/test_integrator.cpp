#include <gtest/gtest.h>
#include "types.h"
#include "state.h"
#include "integrator.h"
#include "comparisons.h"

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

    EXPECT_VECTOR_EQ(Vector3r(5.20898e-05, 6.11485e-05, -450),
        test_state.position());
    EXPECT_VECTOR_EQ(Vector3r(63, 74, 85), test_state.velocity());
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

    EXPECT_QUATERNION_EQ(Quaternionr(0.9903, 0, 0, 0.1392),
        Quaternionr(test_state.attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 10), test_state.angular_velocity());
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

    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), test_state.position());
    EXPECT_VECTOR_EQ(Vector3r(10, 0, 0), test_state.velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 6.283, 0), test_state.acceleration());

    EXPECT_QUATERNION_EQ(Quaternionr(1, 0, 0, 0),
        Quaternionr(test_state.attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0.6283), test_state.angular_velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), test_state.angular_acceleration());
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

    EXPECT_VECTOR_EQ(Vector3r(5.20898e-05, 6.11485e-05, -450),
        test_state.position());
    EXPECT_VECTOR_EQ(Vector3r(63, 74, 85), test_state.velocity());
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

    for(float i = 0; i < 10; i += 0.001) {
        test_state = test.integrate(test_state, 0.001);
    }

    EXPECT_QUATERNION_EQ(Quaternionr(0.9903, 0, 0, 0.1392),
        Quaternionr(test_state.attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 10), test_state.angular_velocity());
}

TEST(IntegratorHeunTest, CircularMotion) {
    IntegratorHeun test = IntegratorHeun();
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

    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), test_state.position());
    EXPECT_VECTOR_EQ(Vector3r(10, 0, 0), test_state.velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 6.283, 0), test_state.acceleration());

    EXPECT_QUATERNION_EQ(Quaternionr(1, 0, 0, 0),
        Quaternionr(test_state.attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0.6283), test_state.angular_velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), test_state.angular_acceleration());
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

    for(float i = 0; i < 10; i += 0.001) {
        test_state = test.integrate(test_state, 0.001);
    }

    EXPECT_VECTOR_EQ(Vector3r(5.20898e-05, 6.11485e-05, -450),
        test_state.position());
    EXPECT_VECTOR_EQ(Vector3r(63, 74, 85), test_state.velocity());
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

    EXPECT_QUATERNION_EQ(Quaternionr(0.9903, 0, 0, 0.1392),
        Quaternionr(test_state.attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 10), test_state.angular_velocity());
}

TEST(IntegratorEulerTest, CircularMotion) {
    IntegratorEuler test = IntegratorEuler();
    State test_state;
    test_state << 0, 0, 0,
                  10, 0, 0,
                  0, 6.283, 0,
                  0, 0, 0, 1,
                  0, 0, 0.6283,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;

    for(float i = 0; i < 10; i += 0.001) {
        test_state = test.integrate(test_state, 0.001);
    }

    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), test_state.position());
    EXPECT_VECTOR_EQ(Vector3r(10, 0, 0), test_state.velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 6.283, 0), test_state.acceleration());

    EXPECT_QUATERNION_EQ(Quaternionr(1, 0, 0, 0),
        Quaternionr(test_state.attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0.6283), test_state.angular_velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), test_state.angular_acceleration());
}
