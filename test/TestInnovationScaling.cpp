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
