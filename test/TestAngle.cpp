#include <gtest/gtest.h>

#include <cmath>
#include <type_traits>

#include "UKF/Core.h"
#include "UKF/MeasurementVector.h"

namespace {
enum Keys { Heading, Bearing, Range };
using AngleState = UKF::StateVector<UKF::Field<Heading, UKF::Angle>>;
using FixedBearing =
    UKF::FixedMeasurementVector<UKF::Field<Bearing, UKF::Angle>>;
using DynamicBearing =
    UKF::DynamicMeasurementVector<UKF::Field<Bearing, UKF::Angle>,
                                  UKF::Field<Range, real_t>>;
constexpr real_t radians(real_t degrees) {
  return degrees * (real_t(3.14159265358979323846) / real_t(180));
}
constexpr real_t tolerance = real_t(2e-6);
}  // namespace

namespace UKF {
template <>
template <>
AngleState AngleState::process_model<>(real_t) const {
  return *this;
}

template <>
template <>
Angle FixedBearing::expected_measurement<AngleState, Bearing>(
    const AngleState& state) {
  return state.get_field<Heading>();
}

template <>
template <>
Angle DynamicBearing::expected_measurement<AngleState, Bearing>(
    const AngleState& state) {
  return state.get_field<Heading>();
}

template <>
template <>
real_t DynamicBearing::expected_measurement<AngleState, Range>(
    const AngleState&) {
  return real_t(0);
}
}  // namespace UKF

TEST(AngleTest, StateSigmaPointsMeanDeltasAndRetractionCrossTheCut) {
  AngleState state;
  state.set_field<Heading>(UKF::Angle(radians(179)));
  const auto sigma = state.calculate_sigma_point_distribution(
      AngleState::CovarianceMatrix::Constant(radians(2)));
  EXPECT_NEAR(radians(-179), sigma(0, 1), tolerance);
  EXPECT_NEAR(radians(177), sigma(0, 2), tolerance);
  const auto mean = AngleState::calculate_sigma_point_mean(sigma);
  EXPECT_NEAR(radians(179), mean.get_field<Heading>().radians(), tolerance);
  const auto deltas = mean.calculate_sigma_point_deltas(sigma);
  EXPECT_NEAR(radians(2), deltas(0, 1), tolerance);
  EXPECT_NEAR(radians(-2), deltas(0, 2), tolerance);
  state.apply_delta(AngleState::StateVectorDelta::Constant(radians(3)));
  EXPECT_NEAR(radians(-178), state.get_field<Heading>().radians(), tolerance);
}

TEST(AngleTest, PeriodicInputsAndHalfTurnTieHaveConsistentInnovations) {
  FixedBearing predicted, actual;
  predicted.set_field<Bearing>(UKF::Angle(radians(179)));
  actual.set_field<Bearing>(UKF::Angle(radians(541)));
  EXPECT_NEAR(radians(2),
              predicted.calculate_innovation(actual).get_field<Bearing>(),
              tolerance);
  predicted.set_field<Bearing>(UKF::Angle(0));
  for (real_t sign : {real_t(-1), real_t(1)}) {
    actual.set_field<Bearing>(UKF::Angle(sign * radians(180)));
    EXPECT_NEAR(radians(-180),
                predicted.calculate_innovation(actual).get_field<Bearing>(),
                tolerance);
  }
}

TEST(AngleTest, DynamicInnovationRetainsFieldOrderAndLinearNoiseUnits) {
  DynamicBearing predicted, actual;
  predicted.set_field<Range>(real_t(100));
  predicted.set_field<Bearing>(UKF::Angle(radians(179)));
  actual.set_field<Bearing>(UKF::Angle(radians(-179)));
  actual.set_field<Range>(real_t(99));
  const auto innovation = predicted.calculate_innovation(actual);
  ASSERT_EQ(2, innovation.size());
  EXPECT_EQ(-1, innovation(0));
  EXPECT_NEAR(radians(2), innovation(1), tolerance);
  static_assert(
      std::is_same_v<decltype(innovation.get_field<Bearing>()), real_t>);
  DynamicBearing::CovarianceVector noise;
  noise.set_field<Bearing>(real_t(0.02));
  noise.set_field<Range>(real_t(4));
  const auto covariance =
      actual.calculate_measurement_covariance(noise, predicted);
  EXPECT_EQ(real_t(0.02), covariance(0, 0));
  EXPECT_EQ(real_t(4), covariance(1, 1));
}

template <typename Filter, typename Measurement>
void check_heading_filter() {
  Filter filter;
  filter.state.template set_field<Heading>(UKF::Angle(radians(179)));
  const real_t sd = radians(2);
  constexpr bool standard =
      std::is_same_v<Filter, UKF::Core<AngleState, Measurement>>;
  if constexpr (standard) {
    filter.covariance.setConstant(sd * sd);
    filter.process_noise_covariance.setZero();
    filter.measurement_covariance.setConstant(sd * sd);
  } else {
    filter.root_covariance.setConstant(sd);
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance.setConstant(sd);
  }
  Measurement measured;
  measured.template set_field<Bearing>(UKF::Angle(radians(-177)));
  if constexpr (std::is_same_v<Filter, UKF::SquareRootParameterEstimationCore<
                                           AngleState, Measurement>>) {
    ASSERT_TRUE(filter.step(measured));
  } else {
    ASSERT_TRUE(filter.step(0, measured));
  }
  EXPECT_NEAR(radians(4), filter.innovation.template get_field<Bearing>(),
              tolerance);
  EXPECT_NEAR(radians(-179),
              filter.state.template get_field<Heading>().radians(), tolerance);
  if constexpr (standard) {
    EXPECT_NEAR(sd * sd / 2, filter.covariance(0, 0), tolerance);
  } else {
    EXPECT_NEAR(sd * sd / 2,
                filter.root_covariance(0, 0) * filter.root_covariance(0, 0),
                tolerance);
  }
}

TEST(AngleTest, AllFiltersUpdateAcrossTheBranchCut) {
  check_heading_filter<UKF::Core<AngleState, FixedBearing>, FixedBearing>();
  check_heading_filter<UKF::SquareRootCore<AngleState, DynamicBearing>,
                       DynamicBearing>();
  check_heading_filter<
      UKF::SquareRootParameterEstimationCore<AngleState, FixedBearing>,
      FixedBearing>();
}
