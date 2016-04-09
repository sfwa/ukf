#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Types.h"
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
    UKF::Field<LatLon, UKF::Vector<2>>,
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>,
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

    EXPECT_VECTOR_EQ(UKF::Vector<2>(-37.8136, 144.9631), test_state.field<LatLon>());
    EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_state.field<Velocity>());
    EXPECT_QUATERNION_EQ(UKF::Quaternion(1, 0, 0, 0), UKF::Quaternion(test_state.field<Attitude>()));
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
    UKF::Field<LatLon, UKF::Vector<2>>,
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>
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

    test_state.field<LatLon>() << 45, 135;
    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 10;

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points, target_sigma_points;

    target_sigma_points << 45, 48.464,     45,     45,     45,     45,     45,     45,     45,     45, 41.536,     45,     45,     45,     45,     45,     45,     45,     45,
                          135,    135, 138.46,    135,    135,    135,    135,    135,    135,    135,    135, 131.54,    135,    135,    135,    135,    135,    135,    135,
                            1,      1,      1,  4.464,      1,      1,      1,      1,      1,      1,      1,      1, -2.464,      1,      1,      1,      1,      1,      1,
                            2,      2,      2,      2,  5.464,      2,      2,      2,      2,      2,      2,      2,      2, -1.464,      2,      2,      2,      2,      2,
                            3,      3,      3,      3,      3,  6.464,      3,      3,      3,      3,      3,      3,      3,      3, -0.464,      3,      3,      3,      3,
                            0,      0,      0,      0,      0,      0, 0.9897,      0,      0,      0,      0,      0,      0,      0,      0, -0.990,      0,      0,      0,
                            0,      0,      0,      0,      0,      0,      0, 0.9897,      0,      0,      0,      0,      0,      0,      0,      0, -0.990,      0,      0,
                            0,      0,      0,      0,      0,      0,      0,      0, 0.9897,      0,      0,      0,      0,      0,      0,      0,      0, -0.990,      0,
                            1,      1,      1,      1,      1,      1, 0.1428, 0.1428, 0.1428,      1,      1,      1,      1,      1,      1, 0.1428, 0.1428, 0.1428,      1,
                           10,     10,     10,     10,     10,     10,     10,     10,     10, 13.464,     10,     10,     10,     10,     10,     10,     10,     10,  6.536;
    sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    EXPECT_VECTOR_EQ(target_sigma_points.col(0),  sigma_points.col(0));
    EXPECT_VECTOR_EQ(target_sigma_points.col(1),  sigma_points.col(1));
    EXPECT_VECTOR_EQ(target_sigma_points.col(2),  sigma_points.col(2));
    EXPECT_VECTOR_EQ(target_sigma_points.col(3),  sigma_points.col(3));
    EXPECT_VECTOR_EQ(target_sigma_points.col(4),  sigma_points.col(4));
    EXPECT_VECTOR_EQ(target_sigma_points.col(5),  sigma_points.col(5));
    EXPECT_VECTOR_EQ(target_sigma_points.col(6),  sigma_points.col(6));
    EXPECT_VECTOR_EQ(target_sigma_points.col(7),  sigma_points.col(7));
    EXPECT_VECTOR_EQ(target_sigma_points.col(8),  sigma_points.col(8));
    EXPECT_VECTOR_EQ(target_sigma_points.col(9),  sigma_points.col(9));
    EXPECT_VECTOR_EQ(target_sigma_points.col(10), sigma_points.col(10));
    EXPECT_VECTOR_EQ(target_sigma_points.col(11), sigma_points.col(11));
    EXPECT_VECTOR_EQ(target_sigma_points.col(12), sigma_points.col(12));
    EXPECT_VECTOR_EQ(target_sigma_points.col(13), sigma_points.col(13));
    EXPECT_VECTOR_EQ(target_sigma_points.col(14), sigma_points.col(14));
    EXPECT_VECTOR_EQ(target_sigma_points.col(15), sigma_points.col(15));
    EXPECT_VECTOR_EQ(target_sigma_points.col(16), sigma_points.col(16));
    EXPECT_VECTOR_EQ(target_sigma_points.col(17), sigma_points.col(17));
    EXPECT_VECTOR_EQ(target_sigma_points.col(18), sigma_points.col(18));
}

TEST(StateVectorTest, SigmaPointMean) {
    MyStateVector test_state;

    test_state.field<LatLon>() << 45, 135;
    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 10;

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);

    EXPECT_VECTOR_EQ(test_state, MyStateVector::calculate_sigma_point_mean(sigma_points));
}

TEST(StateVectorTest, SigmaPointCovariance) {
    MyStateVector test_state;

    test_state.field<LatLon>() << 45, 135;
    test_state.field<Velocity>() << 1, 2, 3;
    test_state.field<Attitude>() << 0, 0, 0, 1;
    test_state.field<Altitude>() << 10;

    MyStateVector::CovarianceMatrix covariance = MyStateVector::CovarianceMatrix::Zero();
    covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

    MyStateVector::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);
    MyStateVector test_mean = MyStateVector::calculate_sigma_point_mean(sigma_points);
    MyStateVector::CovarianceMatrix calculated_covariance = test_mean.calculate_sigma_point_covariance(sigma_points);

    EXPECT_VECTOR_EQ(covariance.col(0),  calculated_covariance.col(0));
    EXPECT_VECTOR_EQ(covariance.col(1),  calculated_covariance.col(1));
    EXPECT_VECTOR_EQ(covariance.col(2),  calculated_covariance.col(2));
    EXPECT_VECTOR_EQ(covariance.col(3),  calculated_covariance.col(3));
    EXPECT_VECTOR_EQ(covariance.col(4),  calculated_covariance.col(4));
    EXPECT_VECTOR_EQ(covariance.col(5),  calculated_covariance.col(5));
    EXPECT_VECTOR_EQ(covariance.col(6),  calculated_covariance.col(6));
    EXPECT_VECTOR_EQ(covariance.col(7),  calculated_covariance.col(7));
    EXPECT_VECTOR_EQ(covariance.col(8),  calculated_covariance.col(8));
}
