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

    EXPECT_EQ(0, test_state.offset<LatLon>());
    EXPECT_EQ(2, test_state.offset<Velocity>());
    EXPECT_EQ(5, test_state.offset<Attitude>());
    EXPECT_EQ(9, test_state.offset<Altitude>());

    EXPECT_EQ(10, MyStateVector::MaxRowsAtCompileTime);
    EXPECT_EQ(10, test_state.size());
}
