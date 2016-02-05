#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Integrator.h"
#include "StateVector.h"

TEST(StateTest, Instantiation) {
    using MyStateVector = UKF::StateVector<
        IntegratorRK4,
        Eigen::Vector2f,
        Eigen::Vector3f,
        Eigen::Quaternionf,
        float
    >;

    MyStateVector test_state;

    EXPECT_EQ(10, test_state.GetDimension());
}
