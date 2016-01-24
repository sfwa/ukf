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
        Eigen::Quaternionf
    > test_state;
    EXPECT_EQ(9, test_state.GetDimension());
}
