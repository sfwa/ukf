// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>

#include <limits>
#include <type_traits>

#include "UKF/Core.h"
#include "UKF/MeasurementVector.h"

using TimedState = UKF::StateVector<UKF::Field<410, real_t>>;
using TimedObservation = UKF::DynamicMeasurementVector<UKF::Field<411, real_t>>;

namespace UKF {
// Deliberately define only the measurement model. A measurement-only step can
// be used without a transition model, even for the two state estimators.
template <>
template <>
real_t TimedObservation::expected_measurement<TimedState, 411, real_t>(
    const TimedState& state, const real_t& bias) {
  return state(0) + bias;
}
}  // namespace UKF

namespace {
using Standard = UKF::Core<TimedState, TimedObservation>;
using Root = UKF::SquareRootCore<TimedState, TimedObservation>;
using Parameter =
    UKF::SquareRootParameterEstimationCore<TimedState, TimedObservation>;

template <typename Filter>
void check_measurements() {
  Filter f;
  f.state << 0;
  const real_t bias = 3;
  if constexpr (std::is_same_v<Filter, Standard>) {
    f.covariance << 1;
    f.measurement_covariance << 1;
    f.process_noise_covariance << 100;
  } else {
    f.root_covariance << 1;
    f.measurement_root_covariance << 1;
    f.process_noise_root_covariance << 10;
  }
  TimedObservation sample;
  sample.set_field<411>(1 + bias);
  ASSERT_TRUE(f.measurement_step(sample, bias));
  EXPECT_NEAR(real_t(0.5), f.state(0), real_t(1e-5));
  sample.set_field<411>(2 + bias);
  ASSERT_TRUE(f.measurement_step(sample, bias));
  EXPECT_NEAR(real_t(1), f.state(0), real_t(1e-5));
  const auto variance = [&]() {
    if constexpr (std::is_same_v<Filter, Standard>) {
      return f.covariance(0, 0);
    } else {
      return f.root_covariance.squaredNorm();
    }
  };
  EXPECT_NEAR(real_t(1) / 3, variance(), real_t(1e-5));

  const real_t x = f.state(0);
  const real_t p = variance();
  const real_t previous_innovation = f.innovation(0);
  sample.set_field<411>(std::numeric_limits<real_t>::quiet_NaN());
  EXPECT_FALSE(f.measurement_step(sample, bias));
  EXPECT_EQ(previous_innovation, f.innovation(0));
  EXPECT_EQ(x, f.state(0));
  EXPECT_EQ(p, variance());
  EXPECT_FALSE(f.a_posteriori_step());

  ASSERT_TRUE(f.measurement_step(TimedObservation(), bias));
  EXPECT_EQ(0, f.innovation.size());
  EXPECT_EQ(x, f.state(0));
  EXPECT_EQ(p, variance());
  EXPECT_FALSE(f.innovation_step(sample, bias));

  sample.set_field<411>(3 + bias);
  ASSERT_TRUE(f.prepare_measurement_step());
  ASSERT_TRUE(f.innovation_step(sample, bias));
  EXPECT_FALSE(f.a_posteriori_step([](const Filter&) { return false; }));
  EXPECT_EQ(x, f.state(0));
  EXPECT_EQ(p, variance());
}
}  // namespace

TEST(MeasurementStepTest, IndependentAtSameTime) {
  check_measurements<Standard>();
  check_measurements<Root>();
  check_measurements<Parameter>();
}
