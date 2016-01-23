#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Integrator.h"
#include "StateModel.h"
#include "comparisons.h"

TEST(StateTest, Instantiation) {
    StateModel<
        IntegratorRK4,
        Eigne::Vector2f,
        Eigen::Vector3f,
        Eigen::Quaternionf
    > test_state;
    EXPECT_EQ(9, test_state.GetDimension());
}
