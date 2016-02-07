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
        UKF::Field<Altitude, Eigen::Vector3f>,
        UKF::Field<Velocity, Eigen::Quaternionf>,
        UKF::Field<Attitude, float>
    >;

    MyStateVector test_state;

    EXPECT_EQ(10, MyStateVector::MaxRowsAtCompileTime);
    EXPECT_EQ(10, test_state.size());
}
