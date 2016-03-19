#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Integrator.h"
#include "StateVector.h"
#include "comparisons.h"

enum MyFields {
    LatLon,
    Altitude,
    Velocity,
    Attitude
};

using MyStateVector = UKF::StateVector<
    IntegratorRK4,
    UKF::Field<LatLon, Eigen::Vector2d>,
    UKF::Field<Velocity, Eigen::Vector3d>,
    UKF::Field<Attitude, Eigen::Quaterniond>,
    UKF::Field<Altitude, real_t>
>;

TEST(StateVectorTest, Instantiation) {
    MyStateVector test_state;

    EXPECT_EQ(10, MyStateVector::MaxRowsAtCompileTime);
    EXPECT_EQ(10, test_state.size());
}

TEST(StateVectorTest, Assignment) {
    MyStateVector test_state;

    test_state.field<LatLon>() << -37.8136, 144.9631;
    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 10;

    EXPECT_VECTOR_EQ(Eigen::Vector2d(-37.8136, 144.9631), test_state.field<LatLon>());
    EXPECT_VECTOR_EQ(Eigen::Vector3d(1, 2, 3), test_state.field<Velocity>());
    EXPECT_QUATERNION_EQ(Eigen::Quaterniond(1, 0, 0, 0), Eigen::Quaterniond(test_state.field<Attitude>()));
    EXPECT_EQ(10, test_state.field<Altitude>()(0));
}