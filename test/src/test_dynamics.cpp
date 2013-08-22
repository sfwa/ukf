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

TEST(FixedWingFlightDynamicsModelTest, Instantiation) {
    FixedWingFlightDynamicsModel test = FixedWingFlightDynamicsModel();
}

TEST(FixedWingFlightDynamicsModelTest, AirframeProperties) {
    FixedWingFlightDynamicsModel model = FixedWingFlightDynamicsModel();
    model.set_mass(5.0);

    Matrix3x3r tensor;
    tensor << 2, 0, 0,
              0, 1, 0,
              0, 0, 3;
    model.set_inertia_tensor(tensor);
}
