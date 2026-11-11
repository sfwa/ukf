#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"

enum StateFields { Position, Velocity };
enum MeasurementFields { ObservedPosition = 100, ObservedVelocity };

using TrackingState = UKF::StateVector<
    UKF::Field<Position, real_t>, UKF::Field<Velocity, real_t>>;
using TrackingMeasurement = UKF::DynamicMeasurementVector<
    UKF::Field<ObservedPosition, real_t>, UKF::Field<ObservedVelocity, real_t>>;
using TrackingFilter = UKF::Core<TrackingState, TrackingMeasurement,
    UKF::NotUsedIntegrator, TrackingMeasurement::DenseCovarianceMatrix>;

namespace UKF {

// A direct discrete-time model needs no derivative or numerical integrator.
template <> template <>
TrackingState TrackingState::process_model<>(real_t dt) const {
    TrackingState predicted = *this;
    predicted.set_field<Position>(get_field<Position>() + dt * get_field<Velocity>());
    return predicted;
}

// State and measurement keys deliberately differ. These functions connect them.
template <> template <>
real_t TrackingMeasurement::expected_measurement<TrackingState, ObservedPosition>(
        const TrackingState& state) {
    return state.get_field<Position>();
}

template <> template <>
real_t TrackingMeasurement::expected_measurement<TrackingState, ObservedVelocity>(
        const TrackingState& state) {
    return state.get_field<Velocity>();
}

}

namespace {

TrackingState::CovarianceMatrix process_noise(real_t dt) {
    // Discrete covariance for a constant-velocity model driven by continuous
    // white acceleration with spectral density q. Recompute it when dt changes.
    const real_t q = real_t(0.005);
    TrackingState::CovarianceMatrix covariance;
    covariance << dt * dt * dt / real_t(3), dt * dt / real_t(2),
                  dt * dt / real_t(2), dt;
    return q * covariance;
}

bool valid_score(real_t score) {
    return std::isfinite(score) && score > real_t(0) && score <= real_t(1);
}

int run_example() {
    TrackingFilter filter;
    filter.state << real_t(0), real_t(0);
    filter.covariance.setZero();
    filter.covariance.diagonal() << real_t(10), real_t(4);
    filter.measurement_covariance.setIdentity();
    const real_t dt = real_t(0.1);
    filter.process_noise_covariance = process_noise(dt);

    int gate_calls = 0;
    auto gate = [&gate_calls](const TrackingFilter& current, const TrackingMeasurement&) {
        ++gate_calls;
        // R is positive definite here, so the innovation covariance is too.
        const auto decomposition = current.innovation_covariance.llt();
        if (decomposition.info() != Eigen::Success) {
            return false;
        }
        const auto weighted = decomposition.solve(current.innovation).eval();
        const real_t nis = current.innovation.dot(weighted);
        // Approximate 99% chi-square thresholds under a Gaussian noise model,
        // for one or two observations. Our deterministic signals only illustrate it.
        const real_t threshold = current.innovation.size() == 1 ? real_t(6.63) : real_t(9.21);
        return std::isfinite(nis) && nis >= real_t(0) && nis <= threshold;
    };

    const std::array<real_t, 4> scores{{real_t(1), real_t(0.75), real_t(0.5), real_t(0.25)}};
    const real_t true_velocity = real_t(0.7);
    const real_t correlation = real_t(0.4);
    real_t raw_squared_error = 0, filtered_squared_error = 0;
    int measured = 0, accepted = 0, rejected = 0, missing = 0, invalid = 0;
    int injected_outliers = 0, rejected_outliers = 0;

    for (int i = 1; i <= 600; ++i) {
        const real_t true_position = real_t(1) + true_velocity * dt * real_t(i);
        real_t score = scores[(i / 30) % scores.size()];
        if (i % 31 == 0) score = real_t(0);
        if (i % 97 == 0) score = std::numeric_limits<real_t>::quiet_NaN();

        TrackingMeasurement measurement;
        const bool outlier = i % 53 == 0;
        real_t observed_position = 0;
        if (i % 17 == 0) {
            ++missing;
        } else if (!valid_score(score)) {
            // Validate confidence before constructing noise or calling the gate.
            ++invalid;
        } else {
            // This calibration belongs to our synthetic sensor. Real detector
            // confidence is not automatically an inverse measurement variance.
            const real_t position_sigma = real_t(0.25) + real_t(0.75) * (real_t(1) - score);
            const real_t velocity_sigma = real_t(0.10) + real_t(0.25) * (real_t(1) - score);
            // Deterministic unit-variance sinusoidal signals keep this demo repeatable.
            const real_t u = std::sqrt(real_t(2)) * std::sin(real_t(1.731) * real_t(i));
            const real_t v = std::sqrt(real_t(2)) * std::cos(real_t(2.117) * real_t(i));
            observed_position = true_position + position_sigma * u + (outlier ? real_t(8) : real_t(0));
            const real_t observed_velocity = true_velocity + velocity_sigma *
                (correlation * u + std::sqrt(real_t(1) - correlation * correlation) * v);

            // Dense noise always follows declaration order: position, velocity.
            // These entries are variances/covariances, not standard deviations.
            filter.measurement_covariance << position_sigma * position_sigma,
                correlation * position_sigma * velocity_sigma,
                correlation * position_sigma * velocity_sigma, velocity_sigma * velocity_sigma;
            if (i % 3 == 0) {
                // Insert velocity first to demonstrate that dynamic field order
                // does not change the declared ordering of dense noise.
                measurement.set_field<ObservedVelocity>(observed_velocity);
            }
            measurement.set_field<ObservedPosition>(observed_position);
            if (outlier) ++injected_outliers;
        }

        TrackingState predicted_state = filter.state;
        TrackingState::CovarianceMatrix predicted_covariance = filter.covariance;
        if (outlier) {
            TrackingFilter predicted = filter;
            predicted.a_priori_step(dt);
            predicted_state = predicted.state;
            predicted_covariance = predicted.covariance;
        }
        const bool applied = filter.step_with_gate(dt, measurement, gate);
        if (measurement.size() != 0) {
            applied ? ++accepted : ++rejected;
            if (outlier && !applied) {
                ++rejected_outliers;
                if ((filter.state - predicted_state).norm() > real_t(1e-10) ||
                        (filter.covariance - predicted_covariance).norm() > real_t(1e-10)) {
                    std::cerr << "Rejected measurement changed the prediction\n";
                    return 1;
                }
            }
            if (i > 20) { // Exclude the deliberately uncertain initialization.
                const real_t raw_error = observed_position - true_position;
                const real_t estimated_error = filter.state.get_field<Position>() - true_position;
                raw_squared_error += raw_error * raw_error;
                filtered_squared_error += estimated_error * estimated_error;
                ++measured;
            }
        }
        if (!filter.state.allFinite() || !filter.covariance.allFinite()) {
            std::cerr << "Nonfinite filter estimate\n";
            return 1;
        }
    }

    // Forecast without advancing the live filter. Q must match the new interval.
    const TrackingState saved_state = filter.state;
    const TrackingState::CovarianceMatrix saved_covariance = filter.covariance;
    TrackingFilter forecast = filter;
    const real_t horizon = real_t(2);
    forecast.process_noise_covariance = process_noise(horizon);
    forecast.a_priori_step(horizon);
    const real_t expected_forecast = saved_state.get_field<Position>() +
        horizon * saved_state.get_field<Velocity>();
    const real_t raw_rmse = std::sqrt(raw_squared_error / real_t(measured));
    const real_t filtered_rmse = std::sqrt(filtered_squared_error / real_t(measured));

    std::cout << std::fixed << std::setprecision(4)
        << "Position RMSE: raw=" << raw_rmse << ", estimated=" << filtered_rmse << '\n'
        << "Updates: accepted=" << accepted << ", gated=" << rejected
        << ", missing=" << missing << ", invalid-confidence=" << invalid << '\n'
        << "Injected outliers rejected: " << rejected_outliers << '/' << injected_outliers << '\n'
        << "Estimated state: position=" << filter.state.get_field<Position>()
        << ", velocity=" << filter.state.get_field<Velocity>() << '\n'
        << "2-second forecast position=" << forecast.state.get_field<Position>() << '\n';

    // Check useful behavior, not exact floating-point output across toolchains.
    if (measured == 0 || !(filtered_rmse < real_t(0.4) * raw_rmse) ||
            injected_outliers == 0 || rejected_outliers != injected_outliers ||
            gate_calls != accepted + rejected ||
            std::abs(forecast.state.get_field<Position>() - expected_forecast) > real_t(1e-9) ||
            (filter.state - saved_state).norm() != real_t(0) ||
            (filter.covariance - saved_covariance).norm() != real_t(0)) {
        std::cerr << "Tracking consistency check failed\n";
        return 1;
    }
    return 0;
}

}

int main() {
    try {
        return run_example();
    } catch (const std::exception& error) {
        std::cerr << "Tracking failed: " << error.what() << '\n';
        return 1;
    }
}
