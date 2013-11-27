#include <gtest/gtest.h>
#include "types.h"
#include "state.h"
#include "dynamics.h"
#include "comparisons.h"

TEST(CentripetalModelTest, Instantiation) {
    CentripetalModel test = CentripetalModel();
}

TEST(CentripetalModelTest, XAxisPositive) {
    State test;
    CentripetalModel model = CentripetalModel();
    AccelerationVector accelerations;
    test << 0, 0, 0,
            10, 0, 0,
            3, 4, 5,
            0, 0, 0, 1,
            0, 1, 0,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    accelerations = model.evaluate(test, ControlVector());

    EXPECT_EQ(0, accelerations[0]);
    EXPECT_EQ(0, accelerations[1]);
    EXPECT_EQ(-10, accelerations[2]);
    EXPECT_EQ(0, accelerations[3]);
    EXPECT_EQ(0, accelerations[4]);
    EXPECT_EQ(0, accelerations[5]);
}

TEST(CentripetalModelTest, XAxisNegative) {
    State test;
    CentripetalModel model = CentripetalModel();
    AccelerationVector accelerations;
    test << 0, 0, 0,
            10, 0, 0,
            3, 4, 5,
            0, 0, 0, 1,
            0, -1, 0,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    accelerations = model.evaluate(test, ControlVector());

    EXPECT_EQ(0, accelerations[0]);
    EXPECT_EQ(0, accelerations[1]);
    EXPECT_EQ(10, accelerations[2]);
    EXPECT_EQ(0, accelerations[3]);
    EXPECT_EQ(0, accelerations[4]);
    EXPECT_EQ(0, accelerations[5]);
}

TEST(CentripetalModelTest, YAxisPositive) {
    State test;
    CentripetalModel model = CentripetalModel();
    AccelerationVector accelerations;
    test << 0, 0, 0,
            0, 10, 0,
            3, 4, 5,
            0, 0, 0, 1,
            0, 0, -1,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    accelerations = model.evaluate(test, ControlVector());

    EXPECT_EQ(10, accelerations[0]);
    EXPECT_EQ(0, accelerations[1]);
    EXPECT_EQ(0, accelerations[2]);
    EXPECT_EQ(0, accelerations[3]);
    EXPECT_EQ(0, accelerations[4]);
    EXPECT_EQ(0, accelerations[5]);
}

TEST(CentripetalModelTest, YAxisNegative) {
    State test;
    CentripetalModel model = CentripetalModel();
    AccelerationVector accelerations;
    test << 0, 0, 0,
            0, 10, 0,
            3, 4, 5,
            0, 0, 0, 1,
            0, 0, 1,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    accelerations = model.evaluate(test, ControlVector());

    EXPECT_EQ(-10, accelerations[0]);
    EXPECT_EQ(0, accelerations[1]);
    EXPECT_EQ(0, accelerations[2]);
    EXPECT_EQ(0, accelerations[3]);
    EXPECT_EQ(0, accelerations[4]);
    EXPECT_EQ(0, accelerations[5]);
}

TEST(CentripetalModelTest, ZAxisPositive) {
    State test;
    CentripetalModel model = CentripetalModel();
    AccelerationVector accelerations;
    test << 0, 0, 0,
            0, 0, 10,
            3, 4, 5,
            0, 0, 0, 1,
            1, 0, 0,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    accelerations = model.evaluate(test, ControlVector());

    EXPECT_EQ(0, accelerations[0]);
    EXPECT_EQ(-10, accelerations[1]);
    EXPECT_EQ(0, accelerations[2]);
    EXPECT_EQ(0, accelerations[3]);
    EXPECT_EQ(0, accelerations[4]);
    EXPECT_EQ(0, accelerations[5]);
}

TEST(CentripetalModelTest, ZAxisNegative) {
    State test;
    CentripetalModel model = CentripetalModel();
    AccelerationVector accelerations;
    test << 0, 0, 0,
            0, 0, 10,
            3, 4, 5,
            0, 0, 0, 1,
            -1, 0, 0,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    accelerations = model.evaluate(test, ControlVector());

    EXPECT_EQ(0, accelerations[0]);
    EXPECT_EQ(10, accelerations[1]);
    EXPECT_EQ(0, accelerations[2]);
    EXPECT_EQ(0, accelerations[3]);
    EXPECT_EQ(0, accelerations[4]);
    EXPECT_EQ(0, accelerations[5]);
}

void _constant_dynamics(const real_t *state,
const real_t *control, real_t *output) {
    output[0] = 1.0;
    output[1] = 2.0;
    output[2] = 3.0;
    output[3] = 4.0;
    output[4] = 5.0;
    output[5] = 6.0;
}

TEST(CustomModelTest, Instantiation) {
    CustomDynamicsModel model = CustomDynamicsModel();
    model.set_function(_constant_dynamics);
}

TEST(CustomModelTest, Output) {
    CustomDynamicsModel model = CustomDynamicsModel();
    model.set_function(_constant_dynamics);

    State test;
    AccelerationVector accelerations;
    test << 0, 0, 0,
            0, 0, 10,
            3, 4, 5,
            0, 0, 0, 1,
            -1, 0, 0,
            6, 7, 8,
            9, 10, 11,
            12, 13, 14;

    accelerations = model.evaluate(test, ControlVector());

    EXPECT_EQ(1.0, accelerations[0]);
    EXPECT_EQ(2.0, accelerations[1]);
    EXPECT_EQ(3.0, accelerations[2]);
    EXPECT_EQ(4.0, accelerations[3]);
    EXPECT_EQ(5.0, accelerations[4]);
    EXPECT_EQ(6.0, accelerations[5]);
}
