#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Integrator.h"
#include "StateVector.h"

TEST(StateTest, Instantiation) {
    using MyStateVector = UKF::StateVector<
        IntegratorRK4,
        UKF::Field<1, Eigen::Vector2f>,
        UKF::Field<2, Eigen::Vector3f>,
        UKF::Field<3, Eigen::Quaternionf>,
        UKF::Field<4, float>
    >;

    MyStateVector test_state;

    EXPECT_EQ(10, MyStateVector::MaxRowsAtCompileTime);
    EXPECT_EQ(10, test_state.size());
}
