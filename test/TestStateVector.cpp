#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Integrator.h"
#include "StateVector.h"

TEST(StateTest, Instantiation) {
    enum MyFields {
        LatLon,
        Altitude,
        Velocity,
        Attitude
    };

    using MyStateVector = UKF::StateVector<
        IntegratorRK4,
        UKF::Field<LatLon, Eigen::Vector2f>,
        UKF::Field<Velocity, Eigen::Vector3f>,
        UKF::Field<Attitude, Eigen::Quaternionf>,
        UKF::Field<Altitude, real_t>
    >;

    MyStateVector test_state;

    EXPECT_EQ(10, MyStateVector::MaxRowsAtCompileTime);
    EXPECT_EQ(10, test_state.size());
}

TEST(StateTest, Assignment) {
    enum MyFields {
        LatLon,
        Altitude,
        Velocity,
        Attitude
    };

    using MyStateVector = UKF::StateVector<
        IntegratorRK4,
        UKF::Field<LatLon, Eigen::Vector2f>,
        UKF::Field<Velocity, Eigen::Vector3f>,
        UKF::Field<Attitude, Eigen::Quaternionf>,
        UKF::Field<Altitude, real_t>
    >;

    MyStateVector test_state;

    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<Altitude>() << 10;

    EXPECT_EQ(1, test_state.field<Velocity>()[0]);
    EXPECT_EQ(2, test_state.field<Velocity>()[1]);
    EXPECT_EQ(3, test_state.field<Velocity>()[2]);
    EXPECT_EQ(10, test_state.field<Altitude>()[0]);
}