#include <cmath>
#include <iomanip>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"

// State and observation keys are deliberately different. The measurement
// specializations below, rather than matching names or numbers, connect them.
enum StateKeys { Orientation = 20, BodyRate = 21 };
enum ObservationKeys { Gyroscope = 100, GravityDirection = 101, MagneticDirection = 102 };

using AttitudeState = UKF::StateVector<
    UKF::Field<Orientation, UKF::Quaternion>,
    UKF::Field<BodyRate, UKF::Vector<3>>>;

using AttitudeMeasurement = UKF::FixedMeasurementVector<
    UKF::Field<Gyroscope, UKF::Vector<3>>,
    UKF::Field<GravityDirection, UKF::Vector<3>>,
    UKF::Field<MagneticDirection, UKF::Vector<3>>>;

using AttitudeFilter = UKF::SquareRootCore<
    AttitudeState, AttitudeMeasurement, UKF::IntegratorRK4>;

namespace UKF {

template <> template <>
AttitudeState AttitudeState::derivative<>() const {
    AttitudeState derivative;
    // q rotates navigation-frame vectors into the body frame. Angular
    // velocity is expressed in body coordinates, in radians per second:
    // q_dot = -0.5 * [0, omega_body] * q.
    Quaternion half_rate;
    half_rate.w() = real_t(0);
    half_rate.vec() = real_t(0.5) * get_field<BodyRate>();
    derivative.set_field<Orientation>(half_rate.conjugate() * get_field<Orientation>());
    // Constant angular velocity is this example's physical process model.
    derivative.set_field<BodyRate>(Vector<3>::Zero());
    return derivative;
}

template <> template <>
Vector<3> AttitudeMeasurement::expected_measurement<AttitudeState, Gyroscope>(
        const AttitudeState& state) {
    return state.get_field<BodyRate>();
}

template <> template <>
Vector<3> AttitudeMeasurement::expected_measurement<AttitudeState, GravityDirection>(
        const AttitudeState& state) {
    return state.get_field<Orientation>() * Vector<3>(0, 0, -1);
}

template <> template <>
Vector<3> AttitudeMeasurement::expected_measurement<AttitudeState, MagneticDirection>(
        const AttitudeState& state) {
    // A second, nonparallel reference direction makes heading observable.
    return state.get_field<Orientation>() * Vector<3>(1, 0, 0);
}

} // namespace UKF

namespace {

// Bounded deterministic disturbances keep the example reproducible across
// standard libraries. These are illustrative sensor errors, not a calibration.
UKF::Vector<3> disturbance(real_t time, real_t phase) {
    return UKF::Vector<3>(std::sin(real_t(13) * time + phase),
        std::cos(real_t(17) * time + phase), std::sin(real_t(19) * time - phase));
}

real_t angular_error(const UKF::Quaternion& estimate, const UKF::Quaternion& truth) {
    const UKF::Quaternion difference = estimate * truth.conjugate();
    // q and -q represent the same orientation.
    return real_t(2) * std::atan2(difference.vec().norm(), std::abs(difference.w()));
}

} // namespace

int main() {
    constexpr real_t dt = real_t(0.01);
    constexpr int steps = 400;
    const UKF::Vector<3> true_rate(real_t(0.1), real_t(-0.05), real_t(0.2));

    AttitudeFilter filter;
    filter.state.set_field<Orientation>(UKF::Quaternion(
        Eigen::AngleAxis<real_t>(real_t(0.15), UKF::Vector<3>::UnitX())));
    filter.state.set_field<BodyRate>(UKF::Vector<3>::Zero());

    // The quaternion occupies four state coefficients, but has a three-
    // dimensional rotation-error covariance. The full state covariance is 6x6.
    filter.root_covariance.setZero();
    filter.root_covariance.diagonal() << real_t(0.2), real_t(0.2), real_t(0.2),
        real_t(0.1), real_t(0.1), real_t(0.1);

    // Explicit diagonal random-walk model: Q(dt) = diag(q_attitude, q_rate)*dt.
    // These values are noise intensities; the filter receives sqrt(Q(dt)).
    // Applications with coupled dynamics should discretize their full Q(dt).
    const real_t attitude_noise = std::sqrt(real_t(1e-5) * dt);
    const real_t rate_noise = std::sqrt(real_t(1e-4) * dt);
    filter.process_noise_root_covariance.setZero();
    filter.process_noise_root_covariance.diagonal() <<
        attitude_noise, attitude_noise, attitude_noise, rate_noise, rate_noise, rate_noise;

    // Square-root noise entries are standard deviations, not variances.
    // They follow the declared gyro, gravity, magnetic measurement order.
    filter.measurement_root_covariance << real_t(0.004), real_t(0.004), real_t(0.004),
        real_t(0.02), real_t(0.02), real_t(0.02), real_t(0.02), real_t(0.02), real_t(0.02);

    UKF::Quaternion truth = UKF::Quaternion::Identity();
    real_t worst_settled_error = real_t(0);
    for (int step = 1; step <= steps; ++step) {
        const real_t time = dt * real_t(step);
        // Exact solution for constant body angular velocity, starting at identity.
        truth = UKF::Quaternion(Eigen::AngleAxis<real_t>(-true_rate.norm() * time,
            true_rate.normalized()));
        AttitudeMeasurement measurement;
        measurement.set_field<Gyroscope>(true_rate + real_t(0.003) * disturbance(time, real_t(0)));
        measurement.set_field<GravityDirection>(truth * UKF::Vector<3>(0, 0, -1) +
            real_t(0.01) * disturbance(time, real_t(1)));
        measurement.set_field<MagneticDirection>(truth * UKF::Vector<3>(1, 0, 0) +
            real_t(0.01) * disturbance(time, real_t(2)));

        filter.step(dt, measurement);
        const UKF::Quaternion estimate = filter.state.get_field<Orientation>();
        if (!filter.state.allFinite() || !filter.root_covariance.allFinite() ||
                std::abs(estimate.norm() - real_t(1)) > real_t(1e-8)) {
            std::cerr << "Nonfinite filter output or non-unit quaternion\n";
            return 1;
        }
        if (step > 100) {
            const real_t error = angular_error(estimate, truth);
            if (error > worst_settled_error) {
                worst_settled_error = error;
            }
        }
    }

    // Read estimates from the filter state, not from the measurement vector.
    const UKF::Quaternion estimate = filter.state.get_field<Orientation>();
    const real_t final_error = angular_error(estimate, truth);
    const real_t rate_error = (filter.state.get_field<BodyRate>() - true_rate).norm();
    std::cout << std::fixed << std::setprecision(6)
        << "Estimated quaternion (w, x, y, z): " << estimate.w() << ", "
        << estimate.x() << ", " << estimate.y() << ", " << estimate.z() << '\n'
        << "Estimated body angular velocity (rad/s): "
        << filter.state.get_field<BodyRate>().transpose() << '\n'
        << "Final orientation error (rad): " << final_error << '\n'
        << "Worst orientation error after settling (rad): " << worst_settled_error << '\n';

    // Conservative bounds validate the example's model, frame conventions and
    // correction path without relying on exact platform-specific rounding.
    if (final_error > real_t(0.03) || worst_settled_error > real_t(0.04) ||
            rate_error > real_t(0.02)) {
        std::cerr << "Attitude tracking accuracy check failed\n";
        return 1;
    }
    return 0;
}
