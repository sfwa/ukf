#include <gtest/gtest.h>

#include "UKF/Core.h"
#include "UKF/MeasurementVector.h"

using ScaledState = UKF::StateVector<UKF::Field<0, UKF::Vector<3>>>;
using ScaledMeasurement =
    UKF::FixedMeasurementVector<UKF::Field<1, UKF::Vector<3>>>;

namespace UKF {
template <>
template <>
ScaledState ScaledState::process_model<>(real_t, const real_t&) const {
  return *this;
}

template <>
template <>
UKF::Vector<3> ScaledMeasurement::expected_measurement<ScaledState, 1, real_t>(
    const ScaledState& state, const real_t& scale) {
  return UKF::Vector<3>(state(0), scale * state(1), state(0));
}
}  // namespace UKF

TEST(InnovationScalingTest, IndependentChannelsRetainTheirGainInAnyUnits) {
  for (bool singular : {false, true}) {
    UKF::Core<ScaledState, ScaledMeasurement> filter;
    const real_t scale = real_t(1e-10);
    filter.state.setZero();
    filter.process_noise_covariance.setZero();
    filter.covariance.setIdentity();
    // The first and last measurements observe the same state. They make
    // the innovation singular when both observations are exact.
    const real_t variance = singular ? real_t(0) : real_t(1);
    filter.measurement_covariance << variance, scale * scale, variance;
    ScaledMeasurement z;
    z << 1, scale, 1;
    ASSERT_TRUE(filter.step(real_t(0), z, scale));
    EXPECT_NEAR(singular ? real_t(1) : real_t(2) / real_t(3), filter.state(0),
                real_t(1e-5));
    EXPECT_NEAR(real_t(0.5), filter.state(1), real_t(1e-5));
    EXPECT_NEAR(real_t(0.5), filter.covariance(1, 1), real_t(1e-5));
  }
}

TEST(InnovationScalingTest, ZeroVarianceChannelKeepsAnExactlyZeroRootRow) {
  using Matrix = UKF::Matrix<5, 5>;
  Matrix factor;
  factor << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 0, 0, -1, 0, 2, 4, 0, 3, 0,
      -1, 2, 1;
  const Matrix covariance = factor * factor.transpose();
  Matrix root;
  ASSERT_TRUE(UKF::Detail::innovation_root(covariance, root));
  EXPECT_EQ(real_t(0), root.row(1).cwiseAbs().maxCoeff());
  UKF::Vector<5> gain;
  gain << 1, 0, 2, -1, 3;
  auto cross = (gain.transpose() * covariance).eval();
  UKF::Vector<5> residual = factor * UKF::Vector<5>::Ones();
  const real_t expected = gain.dot(residual);
  UKF::Detail::whiten(cross, residual, root);
  EXPECT_NEAR(expected, (cross * residual)(0), real_t(1e-4));
}

TEST(InnovationScalingTest, SplitPosteriorUsesThePublicInnovationCovariance) {
  UKF::Core<ScaledState, ScaledMeasurement> filter;
  filter.state.setZero();
  filter.covariance.setIdentity();
  filter.process_noise_covariance.setZero();
  filter.measurement_covariance.setOnes();
  ScaledMeasurement z;
  z << 0, 1, 0;
  const real_t scale = 1;
  ASSERT_TRUE(filter.a_priori_step(real_t(0), scale));
  ASSERT_TRUE(filter.innovation_step(z, scale));
  filter.innovation_covariance(1, 1) = 4;
  ASSERT_TRUE(filter.a_posteriori_step());
  EXPECT_NEAR(real_t(0.25), filter.state(1), real_t(1e-5));
  EXPECT_NEAR(real_t(0.75), filter.covariance(1, 1), real_t(1e-5));
}
