#include <gtest/gtest.h>

#include <limits>

#include "UKF/Core.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Smoother.h"

namespace {
using State = UKF::StateVector<UKF::Field<910, real_t>>;
using Sensor = UKF::FixedMeasurementVector<UKF::Field<911, real_t>>;
using Estimate = UKF::SmoothingEstimate<State>;
constexpr real_t tolerance = sizeof(real_t) == 4 ? real_t(1e-5) : real_t(1e-12);

Estimate scalar(real_t x, real_t p) {
  Estimate result;
  result.state(0) = x;
  result.covariance(0) = p;
  return result;
}
}  // namespace

namespace UKF {
template <>
template <>
State State::process_model<>(real_t dt, const real_t& acceleration) const {
  State result;
  result(0) = (*this)(0) + dt * acceleration * (*this)(0) * (*this)(0);
  return result;
}
}  // namespace UKF

TEST(UnscentedSmoother, ForwardCrossPairsNonlinearTrajectories) {
  UKF::Core<State, Sensor> standard;
  UKF::SquareRootCore<State, Sensor> root;
  standard.state(0) = 2;
  standard.covariance(0) = real_t(0.25);
  standard.process_noise_covariance(0) = 4;
  root.state = standard.state;
  root.root_covariance(0) = real_t(0.5);
  root.process_noise_root_covariance(0) = 2;
  State::CovarianceMatrix cross;
  ASSERT_TRUE(
      standard.a_priori_step_with_cross_covariance(cross, 1, real_t(1)));
  // For a symmetric input with mean 2, Cov(x, x+x²) = (1+2*2) P.
  EXPECT_NEAR(real_t(1.25), cross(0), tolerance);
  const auto expected = standard.covariance(0);
  EXPECT_GT(expected, real_t(4));
  ASSERT_TRUE(root.a_priori_step_with_cross_covariance(cross, 1, real_t(1)));
  EXPECT_NEAR(real_t(1.25), cross(0), tolerance);
  EXPECT_NEAR(standard.state(0), root.state(0), tolerance);
  EXPECT_NEAR(expected, root.root_covariance(0) * root.root_covariance(0),
              tolerance);
  standard.process_noise_covariance(0) =
      std::numeric_limits<real_t>::quiet_NaN();
  cross(0) = -42;
  EXPECT_FALSE(
      standard.a_priori_step_with_cross_covariance(cross, 1, real_t(1)));
  EXPECT_EQ(real_t(-42), cross(0));
}

TEST(UnscentedSmoother, FixedLagMatchesLinearGaussianConditioningAndWraps) {
  UKF::FixedLagSmoother<State, 2> window;
  auto output = scalar(-10, 3);
  EXPECT_FALSE(window.estimate(0, output));
  ASSERT_TRUE(window.reset(scalar(0, 1)));
  // Random walk Q=1, scalar sensor R=1. The two readings are 1 and 2.
  const auto first = scalar(real_t(2) / 3, real_t(2) / 3);
  auto prior = scalar(0, 2);
  auto cross = scalar(0, 1).covariance;
  ASSERT_TRUE(window.push(prior, cross, first));
  ASSERT_TRUE(window.estimate(1, output));
  EXPECT_NEAR(real_t(1) / 3, output.state(0), tolerance);
  EXPECT_NEAR(real_t(2) / 3, output.covariance(0), tolerance);
  prior = scalar(real_t(2) / 3, real_t(5) / 3);
  cross(0) = real_t(2) / 3;
  const auto second = scalar(real_t(1.5), real_t(5) / 8);
  ASSERT_TRUE(window.push(prior, cross, second));
  ASSERT_TRUE(window.estimate(2, output));
  EXPECT_NEAR(real_t(0.5), output.state(0), tolerance);
  EXPECT_NEAR(real_t(5) / 8, output.covariance(0), tolerance);
  ASSERT_TRUE(window.estimate(1, output));
  EXPECT_NEAR(1, output.state(0), tolerance);
  EXPECT_NEAR(real_t(0.5), output.covariance(0), tolerance);
  ASSERT_TRUE(window.estimate(0, output));
  EXPECT_EQ(second.state, output.state);
  ASSERT_TRUE(window.push(second, second.covariance, second));
  EXPECT_EQ(3u, window.size());
  ASSERT_TRUE(window.estimate(2, output));
  EXPECT_NEAR(1, output.state(0), tolerance);
  EXPECT_FALSE(window.estimate(3, output));
  EXPECT_NEAR(1, output.state(0), tolerance);
}

TEST(UnscentedSmoother, RejectsBadInputAndAcceptsExactConstraints) {
  UKF::FixedLagSmoother<State, 0> zero;
  Estimate out = scalar(9, 4);
  const auto exact = scalar(2, 0);
  EXPECT_FALSE(zero.push(exact, exact.covariance, exact));
  ASSERT_TRUE(zero.reset(exact));
  ASSERT_TRUE(zero.push(exact, exact.covariance, exact));
  EXPECT_EQ(1u, zero.size());
  EXPECT_FALSE(zero.estimate(1, out));
  EXPECT_EQ(9, out.state(0));
  auto bad = exact;
  bad.covariance(0) = -1;
  EXPECT_FALSE(zero.reset(bad));
  EXPECT_FALSE(zero.push(bad, exact.covariance, exact));
  auto invalid_cross = exact.covariance;
  invalid_cross(0) = std::numeric_limits<real_t>::quiet_NaN();
  EXPECT_FALSE(zero.push(exact, invalid_cross, exact));
  EXPECT_FALSE(
      UKF::unscented_rts_step(exact, exact, invalid_cross, exact, out));
  EXPECT_EQ(9, out.state(0));
  ASSERT_TRUE(
      UKF::unscented_rts_step(exact, exact, exact.covariance, exact, out));
  EXPECT_EQ(exact.state, out.state);
  EXPECT_EQ(exact.covariance, out.covariance);
  ASSERT_TRUE(zero.estimate(0, out));
  EXPECT_EQ(exact.state, out.state);
}

TEST(UnscentedSmoother, QuaternionIdentityTransitionTransportsFullCovariance) {
  using Q = UKF::StateVector<UKF::Field<912, UKF::Quaternion>,
                             UKF::Field<913, real_t>>;
  using Pose = UKF::SmoothingEstimate<Q>;
  Pose filtered;
  filtered.state.set_field<912>(UKF::Quaternion::Identity());
  filtered.state.set_field<913>(real_t(1));
  filtered.covariance.setIdentity();
  Pose next = filtered;
  Q::StateVectorDelta delta;
  delta << real_t(.4), real_t(-.2), real_t(.3), real_t(.6);
  next.state.apply_delta(delta);
  next.covariance.diagonal() << real_t(.4), real_t(.3), real_t(.2), real_t(.1);
  next.covariance(1, 3) = next.covariance(3, 1) = real_t(.05);
  // An opposite-sign representative describes the same future orientation.
  next.state.template head<4>() *= real_t(-1);
  Pose out;
  ASSERT_TRUE(UKF::unscented_rts_step(filtered, filtered, filtered.covariance,
                                      next, out));
  EXPECT_LT(out.state.delta_from(next.state).norm(), tolerance);
  EXPECT_TRUE(out.covariance.isApprox(next.covariance, tolerance));
}

TEST(UnscentedSmoother,
     FractionalGainMatchesDifferentiatedQuaternionChartChanges) {
  using Q = UKF::StateVector<UKF::Field<912, UKF::Quaternion>,
                             UKF::Field<913, real_t>>;
  using Pose = UKF::SmoothingEstimate<Q>;
  Pose filtered;
  filtered.state.set_field<912>(UKF::Quaternion::Identity());
  filtered.state.set_field<913>(real_t(1));
  filtered.covariance.setIdentity();
  Pose predicted = filtered;
  predicted.covariance *= real_t(2);
  Pose next = predicted;
  Q::StateVectorDelta correction;
  correction << real_t(.9), real_t(-.5), real_t(.4), real_t(.3);
  next.state.apply_delta(correction);
  next.covariance.diagonal() << real_t(.2), real_t(.6), real_t(.3), real_t(.4);
  next.covariance(0, 3) = next.covariance(3, 0) = real_t(.05);
  const Q::StateVectorDelta half = correction * real_t(.5);
  Q expected_state = filtered.state;
  expected_state.apply_delta(half);

  // Differentiate the coordinate maps independently of the analytic Jacobian.
  Q::CovarianceMatrix back, forward;
  const real_t h = sizeof(real_t) == 4 ? real_t(.001) : real_t(.00001);
  for (int i = 0; i < 4; ++i) {
    Q::StateVectorDelta step = Q::StateVectorDelta::Zero();
    step(i) = h;
    Q plus = next.state, minus = next.state;
    plus.apply_delta(step);
    minus.apply_delta(-step);
    back.col(i) =
        (plus.delta_from(predicted.state) - minus.delta_from(predicted.state)) /
        (real_t(2) * h);
    plus = minus = filtered.state;
    plus.apply_delta(half + step);
    minus.apply_delta(half - step);
    forward.col(i) =
        (plus.delta_from(expected_state) - minus.delta_from(expected_state)) /
        (real_t(2) * h);
  }
  const Q::CovarianceMatrix local =
      filtered.covariance +
      real_t(.25) *
          (back * next.covariance * back.transpose() - predicted.covariance);
  const Q::CovarianceMatrix expected = forward * local * forward.transpose();
  Pose result;
  ASSERT_TRUE(UKF::unscented_rts_step(filtered, predicted, filtered.covariance,
                                      next, result));
  EXPECT_LT(result.state.delta_from(expected_state).norm(), tolerance);
  EXPECT_TRUE(result.covariance.isApprox(
      expected, sizeof(real_t) == 4 ? real_t(.001) : real_t(1e-8)));
}
