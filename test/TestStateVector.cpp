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
    EXPECT_EQ(9, MyStateVector::CovarianceMatrix::MaxRowsAtCompileTime);
    EXPECT_EQ(9, MyStateVector::CovarianceMatrix::MaxColsAtCompileTime);
    EXPECT_EQ(10, MyStateVector::SigmaPointDistribution::MaxRowsAtCompileTime);
    EXPECT_EQ(19, MyStateVector::SigmaPointDistribution::MaxColsAtCompileTime);
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

TEST(StateVectorTest, Arithmetic) {

}

TEST(StateVectorTest, DefaultParameters) {
    MyStateVector test_state;

    EXPECT_EQ(1.0, UKF::Parameters::AlphaSquared<MyStateVector>);
    EXPECT_EQ(0.0, UKF::Parameters::Beta<MyStateVector>);
    EXPECT_EQ(3.0, UKF::Parameters::Kappa<MyStateVector>);
    EXPECT_EQ(3.0, UKF::Parameters::Lambda<MyStateVector>);
    EXPECT_EQ(1.0, UKF::Parameters::MRP_A<MyStateVector>);
    EXPECT_EQ(4.0, UKF::Parameters::MRP_F<MyStateVector>);
    EXPECT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WM0<MyStateVector>);
    EXPECT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WC0<MyStateVector>);
    EXPECT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WMI<MyStateVector>);
    EXPECT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WCI<MyStateVector>);
}

using AlternateStateVector = UKF::StateVector<
    IntegratorRK4,
    UKF::Field<LatLon, Eigen::Vector2d>,
    UKF::Field<Velocity, Eigen::Vector3d>,
    UKF::Field<Attitude, Eigen::Quaterniond>
>;

template <> constexpr real_t UKF::Parameters::AlphaSquared<AlternateStateVector> = 2.0;
template <> constexpr real_t UKF::Parameters::Beta<AlternateStateVector> = 1.0;
template <> constexpr real_t UKF::Parameters::Kappa<AlternateStateVector> = 4.0;
template <> constexpr real_t UKF::Parameters::MRP_F<AlternateStateVector> =
    2.0 * (UKF::Parameters::AlphaSquared<AlternateStateVector> + 3.0);


TEST(StateVectorTest, CustomParameters) {
    AlternateStateVector test_state;

    EXPECT_EQ(2.0, UKF::Parameters::AlphaSquared<AlternateStateVector>);
    EXPECT_EQ(1.0, UKF::Parameters::Beta<AlternateStateVector>);
    EXPECT_EQ(4.0, UKF::Parameters::Kappa<AlternateStateVector>);
    EXPECT_EQ(16.0, UKF::Parameters::Lambda<AlternateStateVector>);
    EXPECT_EQ(1.0, UKF::Parameters::MRP_A<AlternateStateVector>);
    EXPECT_EQ(10.0, UKF::Parameters::MRP_F<AlternateStateVector>);
    EXPECT_EQ(16.0 / 24.0, UKF::Parameters::Sigma_WM0<AlternateStateVector>);
    EXPECT_EQ(16.0 / 24.0, UKF::Parameters::Sigma_WC0<AlternateStateVector>);
    EXPECT_EQ(1.0 / 48.0, UKF::Parameters::Sigma_WMI<AlternateStateVector>);
    EXPECT_EQ(1.0 / 48.0, UKF::Parameters::Sigma_WCI<AlternateStateVector>);

    EXPECT_EQ(1.0, UKF::Parameters::AlphaSquared<MyStateVector>);
    EXPECT_EQ(0.0, UKF::Parameters::Beta<MyStateVector>);
    EXPECT_EQ(3.0, UKF::Parameters::Kappa<MyStateVector>);
    EXPECT_EQ(3.0, UKF::Parameters::Lambda<MyStateVector>);
    EXPECT_EQ(1.0, UKF::Parameters::MRP_A<MyStateVector>);
    EXPECT_EQ(4.0, UKF::Parameters::MRP_F<MyStateVector>);
    EXPECT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WM0<MyStateVector>);
    EXPECT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WC0<MyStateVector>);
    EXPECT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WMI<MyStateVector>);
    EXPECT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WCI<MyStateVector>);
}

TEST(StateVectorTest, SigmaPointGeneration) {
    MyStateVector test_state;

    test_state.field<LatLon>() << -37.8136, 144.9631;
    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 10;

    MyStateVector::CovarianceMatrix covariance;
    MyStateVector::SigmaPointDistribution sigma_points;


}
