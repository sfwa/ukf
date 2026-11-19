#include <gtest/gtest.h>

#include <cmath>
#include <limits>
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

TEST(AngleTest, DefaultsToZero) {
  static_assert(std::is_default_constructible<UKF::Angle>::value,
                "An angle can be declared before it is set");
  const UKF::Angle angles[2];
  EXPECT_EQ(real_t(0), angles[0].radians());
  EXPECT_EQ(real_t(0), angles[1].radians());
  EXPECT_EQ(real_t(1.5), UKF::Angle(real_t(1.5)).radians());
}

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

/*
Headings many turns out, each exact in both precisions, and their principal
values, from x - 2 pi round(x / 2 pi) in 120-digit decimal arithmetic. They
do not come from the library's own reduction.
*/
struct DistantHeading {
  real_t radians;
  double principal;
};
constexpr DistantHeading distant_headings[] = {
    {real_t(1e6), -0.357564167085735044015},
    {real_t(-1e6), 0.357564167085735044015},
    {real_t(1e10), -0.509231072165734782829},
    {real_t(0x1p60), -2.16131999313972709154},
    {real_t(-0x1p100), 1.05964853168015739226},
    {real_t(0x1.8p126), 0.504798781762045098518},
};

/* A few roundings of a half-turn angle in the stored precision. */
constexpr real_t phase_tolerance =
    real_t(8) * std::numeric_limits<real_t>::epsilon() * real_t(3.14159265);

template <typename Filter>
Filter heading_filter(real_t start) {
  Filter filter;
  filter.state.template set_field<Heading>(UKF::Angle(start));
  if constexpr (std::is_same_v<Filter, UKF::Core<AngleState, FixedBearing>>) {
    filter.covariance.setConstant(real_t(0.01));
    filter.process_noise_covariance.setZero();
    filter.measurement_covariance.setConstant(real_t(0.01));
  } else {
    filter.root_covariance.setConstant(real_t(0.1));
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance.setConstant(real_t(0.1));
  }
  return filter;
}

template <typename Filter>
void step_heading(Filter& filter, real_t bearing) {
  FixedBearing measured;
  measured.set_field<Bearing>(UKF::Angle(bearing));
  if constexpr (std::is_same_v<Filter, UKF::SquareRootParameterEstimationCore<
                                           AngleState, FixedBearing>>) {
    ASSERT_TRUE(filter.step(measured));
  } else {
    ASSERT_TRUE(filter.step(0, measured));
  }
}

/*
A heading many turns out must read as its principal value: observing it from
that value gives no innovation, in either precision.
*/
template <typename Filter>
void check_distant_reading(const DistantHeading& heading) {
  const real_t principal = static_cast<real_t>(heading.principal);
  Filter filter = heading_filter<Filter>(principal);
  step_heading(filter, heading.radians);
  EXPECT_NEAR(real_t(0), filter.innovation.template get_field<Bearing>(),
              phase_tolerance);
  EXPECT_NEAR(principal, filter.state.template get_field<Heading>().radians(),
              phase_tolerance);
}

TEST(AngleTest, DistantHeadingReadsAsItsPrincipalValue) {
  for (const DistantHeading& heading : distant_headings) {
    SCOPED_TRACE(heading.radians);
    check_distant_reading<UKF::Core<AngleState, FixedBearing>>(heading);
    check_distant_reading<UKF::SquareRootCore<AngleState, FixedBearing>>(
        heading);
    check_distant_reading<
        UKF::SquareRootParameterEstimationCore<AngleState, FixedBearing>>(
        heading);
  }
}

/*
A stored heading need not be wrapped, and an unwrapped one must take the
same correction as its principal value. In single precision a correction of
a hundredth of a radian added to a heading of a million radians would round
away.
*/
template <typename Filter>
void check_unwrapped_heading_correction(real_t start, real_t principal) {
  Filter filter = heading_filter<Filter>(start);
  step_heading(filter, principal + real_t(0.02));
  EXPECT_NEAR(real_t(0.02), filter.innovation.template get_field<Bearing>(),
              tolerance);
  EXPECT_NEAR(principal + real_t(0.01),
              filter.state.template get_field<Heading>().radians(), tolerance);
}

TEST(AngleTest, UnwrappedHeadingTakesTheWholeCorrection) {
  for (const DistantHeading& heading : distant_headings) {
    const real_t principal = static_cast<real_t>(heading.principal);
    for (real_t start : {heading.radians, principal}) {
      SCOPED_TRACE(start);
      check_unwrapped_heading_correction<UKF::Core<AngleState, FixedBearing>>(
          start, principal);
      check_unwrapped_heading_correction<
          UKF::SquareRootCore<AngleState, FixedBearing>>(start, principal);
      check_unwrapped_heading_correction<
          UKF::SquareRootParameterEstimationCore<AngleState, FixedBearing>>(
          start, principal);
    }
  }
}
