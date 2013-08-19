#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "types.h"
#include "state.h"
#include "dynamics.h"

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

    EXPECT_EQ(accelerations(0), 0);
    EXPECT_EQ(accelerations(1), 0);
    EXPECT_EQ(accelerations(2), -10);
    EXPECT_EQ(accelerations(3), 0);
    EXPECT_EQ(accelerations(4), 0);
    EXPECT_EQ(accelerations(5), 0);
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

    EXPECT_EQ(accelerations(0), 0);
    EXPECT_EQ(accelerations(1), 0);
    EXPECT_EQ(accelerations(2), 10);
    EXPECT_EQ(accelerations(3), 0);
    EXPECT_EQ(accelerations(4), 0);
    EXPECT_EQ(accelerations(5), 0);
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

    EXPECT_EQ(accelerations(0), 10);
    EXPECT_EQ(accelerations(1), 0);
    EXPECT_EQ(accelerations(2), 0);
    EXPECT_EQ(accelerations(3), 0);
    EXPECT_EQ(accelerations(4), 0);
    EXPECT_EQ(accelerations(5), 0);
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

    EXPECT_EQ(accelerations(0), -10);
    EXPECT_EQ(accelerations(1), 0);
    EXPECT_EQ(accelerations(2), 0);
    EXPECT_EQ(accelerations(3), 0);
    EXPECT_EQ(accelerations(4), 0);
    EXPECT_EQ(accelerations(5), 0);
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

    EXPECT_EQ(accelerations(0), 0);
    EXPECT_EQ(accelerations(1), -10);
    EXPECT_EQ(accelerations(2), 0);
    EXPECT_EQ(accelerations(3), 0);
    EXPECT_EQ(accelerations(4), 0);
    EXPECT_EQ(accelerations(5), 0);
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

    EXPECT_EQ(accelerations(0), 0);
    EXPECT_EQ(accelerations(1), 10);
    EXPECT_EQ(accelerations(2), 0);
    EXPECT_EQ(accelerations(3), 0);
    EXPECT_EQ(accelerations(4), 0);
    EXPECT_EQ(accelerations(5), 0);
}
