#include <cmath>
#include <cstdio>

#include "UKF/Core.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Smoother.h"

using Position = UKF::StateVector<UKF::Field<0, real_t>>;
using Sensor = UKF::FixedMeasurementVector<UKF::Field<1, real_t>>;

namespace UKF {
template <>
template <>
Position Position::process_model<>(real_t) const {
  return *this;  // A discrete random walk, with Q supplied by the filter.
}
template <>
template <>
real_t Sensor::expected_measurement<Position, 1>(const Position& state) {
  return state.get_field<0>();
}
}  // namespace UKF

int main() {
  UKF::Core<Position, Sensor> filter;
  filter.state(0) = 0;
  filter.covariance(0) = 1;
  filter.process_noise_covariance(0) = 1;
  filter.measurement_covariance(0) = 1;
  UKF::FixedLagSmoother<Position, 2> smoother;
  if (!smoother.reset({filter.state, filter.covariance})) return 1;
  for (real_t reading : {real_t(1), real_t(2)}) {
    Position::CovarianceMatrix cross;
    if (!filter.a_priori_step_with_cross_covariance(cross, real_t(1))) return 1;
    const UKF::SmoothingEstimate<Position> predicted{filter.state,
                                                     filter.covariance};
    Sensor measurement;
    measurement.set_field<1>(reading);
    if (!filter.innovation_step(measurement) || !filter.a_posteriori_step() ||
        !smoother.push(predicted, cross, {filter.state, filter.covariance}))
      return 1;
  }
  UKF::SmoothingEstimate<Position> initial;
  if (!smoother.estimate(2, initial)) return 1;
  std::printf(
      "Initial position given two later observations: %.3f, variance %.3f\n",
      static_cast<double>(initial.state(0)),
      static_cast<double>(initial.covariance(0)));
  return std::abs(initial.state(0) - real_t(.5)) < real_t(1e-5) &&
                 std::abs(initial.covariance(0) - real_t(.625)) < real_t(1e-5)
             ? 0
             : 1;
}
