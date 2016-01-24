#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Integrator.h"
#include "StateVector.h"

TEST(StateTest, Instantiation) {
    UKF::StateVector<
        IntegratorRK4,
        Eigen::Vector2f,
        Eigen::Vector3f,
        Eigen::Quaternionf,
        float
    > test_state;
    EXPECT_EQ(10, test_state.GetDimension());
}
