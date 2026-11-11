// MIT license; see the repository LICENSE.
// Robust sensor interpretation around sfwa/ukf's square-root filters.
// The gyro-frame acceleration filtering and separation of magnetic heading
// follow the design described by Laidig and Seel, VQF, Information Fusion 2023,
// https://doi.org/10.1016/j.inffus.2022.10.014 . This is an independent
// implementation; it does not call or vendor the VQF implementation.
#include "robust_estimator.h"
#include "robust_geometry.h"
#include "robust_correction.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "UKF/Types.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Integrator.h"
#include "UKF/Core.h"

namespace Robust {
using V3 = UKF::Vector<3>;
using V2 = UKF::Vector<2>;
using Q = UKF::Quaternion;
constexpr double Pi = 3.14159265358979323846;
constexpr double Rad = Pi / 180;
enum Keys { Attitude = 100, Bias, Gravity, RestGyro };
using AttitudeState = UKF::StateVector<UKF::Field<Attitude, Q>>;
using BiasState = UKF::StateVector<UKF::Field<Bias, V3>>;
using GravityMeasurement = UKF::DynamicMeasurementVector<UKF::Field<Gravity, V2>>;
using BiasMeasurement = UKF::DynamicMeasurementVector<UKF::Field<RestGyro, V3>>;
struct AttitudeInput { V3 omega; Eigen::Matrix<real_t, 3, 2> tangent; };

Q increment(const V3& rotation) {
    const double angle = rotation.norm();
    return angle < 1e-12 ? Q::Identity() : Q(Eigen::AngleAxis<real_t>(angle, rotation / angle));
}
double wrap(double angle) { return std::remainder(angle, 2 * Pi); }
V3 vector(const std::array<double, 3>& x) { return V3(x[0], x[1], x[2]); }
}

namespace UKF {
namespace Parameters {
template <> constexpr real_t AlphaSquared<Robust::AttitudeState> = 1;
template <> constexpr real_t Beta<Robust::AttitudeState> = 2;
template <> constexpr real_t Kappa<Robust::AttitudeState> = 0;
template <> constexpr real_t AlphaSquared<Robust::BiasState> = 1;
template <> constexpr real_t Beta<Robust::BiasState> = 2;
template <> constexpr real_t Kappa<Robust::BiasState> = 0;
}
template <> template <>
Robust::AttitudeState Robust::AttitudeState::process_model<>(
        real_t dt, const Robust::AttitudeInput& input) const {
    Robust::AttitudeState result = *this;
    result.set_field<Robust::Attitude>((get_field<Robust::Attitude>() * Robust::increment(input.omega * dt)).normalized());
    return result;
}
template <> template <>
Robust::V2 Robust::GravityMeasurement::expected_measurement<Robust::AttitudeState, Robust::Gravity, Robust::AttitudeInput>(
        const Robust::AttitudeState& state, const Robust::AttitudeInput& input) {
    return input.tangent.transpose() * (state.get_field<Robust::Attitude>().conjugate() * Robust::V3(0, 0, 1));
}
template <> template <>
Robust::V3 Robust::BiasMeasurement::expected_measurement<Robust::BiasState, Robust::RestGyro>(
        const Robust::BiasState& state) { return state.get_field<Robust::Bias>(); }
}

namespace Robust {
struct AttitudeIntegrator {
    static AttitudeState integrate(real_t dt, const AttitudeState& state, const AttitudeInput& input) {
        return state.process_model(dt, input);
    }
};
struct BiasIntegrator {
    static BiasState integrate(real_t, const BiasState& state) { return state; }
};
class AttitudeFilter : public UKF::SquareRootCore<AttitudeState, GravityMeasurement, AttitudeIntegrator> {
public:
    void a_posteriori_step() {
        if (innovation.size() == 0) return;
        // Keep the UKF's nonlinear prediction, measurement moments and cross
        // covariance. Gravity can correct only the two horizontal components
        // of a LEFT/world attitude error, so constrain the gain's world-Z row.
        UKF::Matrix<3, 2> cross_covariance = UKF::Matrix<3, 2>::Zero();
        for (std::size_t i = 1; i < AttitudeState::num_sigma(); ++i)
            cross_covariance.noalias() += UKF::Parameters::Sigma_WCI<AttitudeState>
                * (w_prime.col(i) * z_prime.col(i).transpose());
        cross_covariance.noalias() += UKF::Parameters::Sigma_WC0<AttitudeState>
            * (w_prime.col(0) * z_prime.col(0).transpose());
        const UKF::Matrix<3, 3> prior = root_covariance * root_covariance.transpose();
        const UKF::Matrix<2, 2> observation = innovation_root_covariance * innovation_root_covariance.transpose();
        const auto update = robust_detail::constrained_gravity_update(prior, cross_covariance, observation, innovation);
        const auto posterior_root = UKF::Detail::covariance_root(update.posteriorP, prior.cwiseAbs().maxCoeff());
        auto corrected_state = state;
        corrected_state.apply_delta(update.delta);
        state = corrected_state;
        root_covariance = posterior_root;
    }
};
using BiasFilter = UKF::SquareRootCore<BiasState, BiasMeasurement, BiasIntegrator>;

// Exact zero-order-hold solution of y'' + 2/tau*y' + 2/tau^2*y = 2/tau^2*x.
// This is a second-order Butterworth low-pass, with cutoff sqrt(2)/(2*pi*tau).
// Continuous-state integration preserves its time constant for variable dt.
class LowPass {
    V3 value = V3::Zero(), rate = V3::Zero();
    bool initialized = false;
public:
    V3 update(const V3& sample, double dt, double tau) {
        if (!initialized) { value = sample; initialized = true; return value; }
        const double h = dt / tau, decay = std::exp(-h), c = std::cos(h), s = std::sin(h);
        const V3 error = value - sample;
        value = sample + decay * ((c + s) * error + tau * s * rate);
        rate = decay * ((-2 * s / tau) * error + (c - s) * rate);
        return value;
    }
};

struct RestDetector {
    LowPass gyro_filter, acc_filter;
    V3 gyro_average = V3::Zero();
    double duration = 0;
    bool update(const V3& gyro, const V3& acc, double dt) {
        gyro_average = gyro_filter.update(gyro, dt, .5);
        const V3 acc_average = acc_filter.update(acc, dt, .5);
        const bool quiet = (gyro - gyro_average).norm() <= 2 * Rad &&
                           (acc - acc_average).norm() <= .5 &&
                           gyro_average.cwiseAbs().maxCoeff() <= 2 * Rad;
        duration = quiet ? duration + dt : 0;
        return duration >= 1.5;
    }
};

// Magnetic data only updates this scalar world-yaw correction; it has no path
// into the six-axis attitude, acceleration filtering, or bias state.
struct HeadingFilter {
    double delta = 0, norm = 0, dip = 0, ref_norm = 0, ref_dip = 0;
    double good_time = 0, candidate_norm = 0, candidate_dip = 0, candidate_time = 0;
    bool initialized = false, rejected = false;

    void update(const V3& magnetic_world, double dt, bool protect, double angular_speed) {
        const double observed_norm = magnetic_world.norm();
        if (observed_norm < 1e-12) throw std::runtime_error("Zero magnetic field");
        const double observed_dip = std::asin(std::clamp(magnetic_world.z() / observed_norm, -1.0, 1.0));
        const double horizontal = magnetic_world.head<2>().norm();
        if (!initialized) {
            if (horizontal < 1e-9 * observed_norm) throw std::runtime_error("Initial magnetic field is vertical");
            norm = ref_norm = candidate_norm = observed_norm;
            dip = ref_dip = candidate_dip = observed_dip;
            delta = -std::atan2(magnetic_world.y(), magnetic_world.x());
            initialized = true;
            return;
        }
        const double smooth = -std::expm1(-dt / .05);
        norm += smooth * (observed_norm - norm);
        dip += smooth * (observed_dip - dip);
        bool acceptable = horizontal > 1e-6 * observed_norm;
        if (protect) {
            acceptable = acceptable && std::abs(norm - ref_norm) <= .1 * ref_norm && std::abs(dip - ref_dip) <= 10 * Rad;
            good_time = acceptable ? good_time + dt : 0;
            if (!acceptable) rejected = true;
            else if (good_time >= .5) {
                rejected = false;
                const double gain = -std::expm1(-dt / 20);
                ref_norm += gain * (norm - ref_norm);
                ref_dip += gain * (dip - ref_dip);
            }
            // A new spatially homogeneous field can be accepted after sustained
            // rotation. Stationary readings alone cannot certify a replacement.
            if (std::abs(norm - candidate_norm) < .1 * candidate_norm && std::abs(dip - candidate_dip) < 10 * Rad) {
                if (angular_speed >= 20 * Rad) candidate_time += dt;
                const double gain = -std::expm1(-dt / 20);
                candidate_norm += gain * (norm - candidate_norm);
                candidate_dip += gain * (dip - candidate_dip);
            } else {
                candidate_norm = norm; candidate_dip = dip; candidate_time = 0;
            }
            if (candidate_time >= 20) {
                ref_norm = candidate_norm; ref_dip = candidate_dip;
                rejected = false; good_time = .5; candidate_time = 0;
            }
        } else rejected = !acceptable;
        if (!rejected && horizontal > 1e-6 * observed_norm) {
            const double residual = wrap(-std::atan2(magnetic_world.y(), magnetic_world.x()) - delta);
            delta = wrap(delta + (-std::expm1(-dt / 9)) * residual);
        }
    }
};

Eigen::Matrix<real_t, 3, 2> tangent_basis(const V3& up) {
    Eigen::Matrix<real_t, 3, 2> basis;
    const V3 axis = std::abs(up.z()) < .9 ? V3(0, 0, 1) : V3(1, 0, 0);
    basis.col(0) = up.cross(axis).normalized();
    basis.col(1) = up.cross(basis.col(0)).normalized();
    return basis;
}

void recenter(AttitudeFilter& filter, const Q& predicted) {
    // Transport the covariance into the corrected local Gibbs chart. The core
    // update's covariance is still in the pre-correction tangent coordinates.
    filter.root_covariance = robust_detail::recenter_attitude_root(
        predicted, filter.state.get_field<Attitude>(), filter.root_covariance);
}
}

RobustResult estimate_robust(const std::vector<RobustSample>& samples, int axes,
                             const std::string& stage, const std::string& diagnostics) {
    using namespace Robust;
    const bool filtered_acceleration = stage != "robust-rest";
    const bool protected_magnetometer = stage == "robust";
    const bool protected_acceleration = stage == "robust";
    if (stage != "robust-rest" && stage != "robust-acc" && stage != "robust") throw std::runtime_error("Unknown robust stage");
    if (samples.size() < 2) throw std::runtime_error("At least two samples required");
    AttitudeFilter attitude;
    BiasFilter bias;
    constexpr double Gravity = 9.80665;
    auto plausible_acceleration = [](double norm) { return norm >= .7 * Gravity && norm <= 1.3 * Gravity; };
    bool gravity_initialized = !protected_acceleration || plausible_acceleration(vector(samples.front().acc).norm());
    const Q initial = gravity_initialized ? Q(Q::FromTwoVectors(vector(samples.front().acc).normalized(), V3(0, 0, 1))) : Q::Identity();
    attitude.state.set_field<Attitude>(initial);
    attitude.root_covariance.setIdentity(); attitude.root_covariance *= .05;
    bias.state.set_field<Bias>(V3::Zero());
    bias.root_covariance.setIdentity(); bias.root_covariance *= .5 * Rad;
    // Effective rest-observation noise is set from the desired converged
    // posterior uncertainty at each update below. It must not be interpreted
    // as independent raw-sample noise for the correlated gyro averages.
    bias.measurement_root_covariance.setConstant(1);
    Q gyro_frame = Q::Identity();
    LowPass acceleration_filter;
    V3 last_accepted_acceleration = V3::Zero();
    bool has_accepted_acceleration = false, direction_rejected = false;
    RestDetector rest_detector;
    HeadingFilter heading;
    double correction_clock = 0, acc_innovation = 0;
    constexpr double UpdatePeriod = .05;
    std::ofstream diagnostic;
    if (!diagnostics.empty()) {
        diagnostic.open(diagnostics);
        if (!diagnostic) throw std::runtime_error("Cannot open diagnostics: " + diagnostics);
        diagnostic << "time,bx,by,bz,rest,mag_rejected,acc_innovation_rad,attitude_covariance_trace,bias_covariance_trace,mag_norm,mag_dip_rad,heading_correction_rad,acc_rejected,acc_norm_m_s2\n" << std::setprecision(17);
    }
    RobustResult result;
    result.quaternions.reserve(samples.size());
    const auto begin = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        const double dt = i ? sample.time - samples[i - 1].time : 0;
        const V3 gyro = vector(sample.gyro), acc = vector(sample.acc);
        const double acceleration_norm = acc.norm();
        const bool norm_accepted = !protected_acceleration || plausible_acceleration(acceleration_norm);
        bool rest = rest_detector.update(gyro, acc, dt);
        if (!norm_accepted) { rest_detector.duration = 0; rest = false; }
        const V3 previous_bias = bias.state.get_field<Bias>();
        AttitudeInput input;
        input.omega = V3::Zero();
        if (i) input.omega = .5 * (gyro + vector(samples[i - 1].gyro)) - previous_bias;
        input.tangent.setZero();
        if (i) {
            gyro_frame = (gyro_frame * increment(input.omega * dt)).normalized();
            attitude.process_noise_root_covariance.setIdentity();
            attitude.process_noise_root_covariance *= .02 * std::sqrt(dt);
            attitude.a_priori_step(dt, input);
            bias.process_noise_root_covariance.setIdentity();
            bias.process_noise_root_covariance *= (.01 * Rad) * std::sqrt(dt);
            bias.a_priori_step(dt);
        }
        bool just_initialized = false;
        if (!gravity_initialized && norm_accepted) {
            // An invalid startup reading cannot define gravity. Initialize on
            // the first usable current reading, without consulting future data.
            const Q predicted = attitude.state.get_field<Attitude>();
            const Q correction = Q::FromTwoVectors(predicted * acc.normalized(), V3(0, 0, 1));
            attitude.state.set_field<Attitude>((correction * predicted).normalized());
            attitude.root_covariance.setIdentity(); attitude.root_covariance *= .05;
            acceleration_filter = LowPass(); has_accepted_acceleration = false;
            gravity_initialized = true; just_initialized = true; correction_clock = 0;
        }
        V3 measured_up = acc;
        if (filtered_acceleration) {
            if (acceleration_norm > 1e-9) {
                last_accepted_acceleration = gyro_frame * acc;
                has_accepted_acceleration = true;
            }
            // Preserve finite nonzero forces in the linear low-pass: selectively
            // dropping large acceleration samples biases cancellation over a
            // motion cycle. Norm gating controls the measurement UPDATE below.
            // A zero vector carries no directional information, so hold the
            // previous input while advancing the filter clock during freefall.
            if (has_accepted_acceleration)
                measured_up = gyro_frame.conjugate() * acceleration_filter.update(last_accepted_acceleration, dt, 3);
        }
        correction_clock += dt;
        if (i && correction_clock + 1e-12 >= UpdatePeriod) {
            const double observation_dt = correction_clock;
            correction_clock = 0;
            if (rest) {
                const double target_variance = std::pow(.03 * Rad, 2);
                const double process_variance = std::pow(.01 * Rad, 2) * observation_dt;
                // Scalar steady-state Riccati solution: R=P_target^2/Q+P_target.
                // This sets a 0.03 deg/s converged bias uncertainty independently
                // of the actual sampling/correction interval.
                bias.measurement_root_covariance.setConstant(std::sqrt(
                    target_variance * target_variance / process_variance + target_variance));
                BiasMeasurement measurement;
                measurement.set_field<RestGyro>(rest_detector.gyro_average);
                bias.innovation_step(measurement);
                bias.a_posteriori_step();
                // Rest measurements and accepted low-rate observations bound
                // this linear estimate naturally; reject impossible states.
                if (bias.state.get_field<Bias>().cwiseAbs().maxCoeff() > 2.1 * Rad)
                    throw std::runtime_error("Rest bias estimate exceeded physical bound");
            }
            if (norm_accepted && !just_initialized && measured_up.norm() > 1e-6) {
                measured_up.normalize();
                const Q predicted = attitude.state.get_field<Attitude>();
                const V3 expected_up = predicted.conjugate() * V3(0, 0, 1);
                input.tangent = tangent_basis(expected_up);
                acc_innovation = std::acos(std::clamp(expected_up.dot(measured_up), -1.0, 1.0));
                direction_rejected = protected_acceleration && acc_innovation > 60 * Rad;
                // Tangent-plane observations have two independent components.
                // Inflate R for overlapping filtered outputs; this is a
                // conservative effective-noise model, not a whitened sensor.
                const double sigma = filtered_acceleration ? .002 * std::sqrt(std::max(1.0, 6 / observation_dt)) : .1;
                attitude.measurement_root_covariance.setConstant(sigma);
                if (!direction_rejected) {
                    GravityMeasurement measurement;
                    measurement.set_field<Robust::Gravity>(input.tangent.transpose() * measured_up);
                    attitude.innovation_step(measurement, input);
                    attitude.a_posteriori_step();
                    recenter(attitude, predicted);
                }
            }
        }
        const Q six_axis = attitude.state.get_field<Attitude>();
        Q output = six_axis;
        if (axes == 9 && gravity_initialized) {
            heading.update(six_axis * vector(sample.mag), dt, protected_magnetometer, input.omega.norm());
            output = Q(Eigen::AngleAxis<real_t>(heading.delta, V3(0, 0, 1))) * six_axis;
        }
        const V3 b = bias.state.get_field<Bias>();
        if (!output.coeffs().allFinite() || !b.allFinite() || !attitude.root_covariance.allFinite() || !bias.root_covariance.allFinite())
            throw std::runtime_error("Non-finite robust filter state");
        output.normalize();
        result.quaternions.push_back({output.w(), output.x(), output.y(), output.z()});
        if (diagnostic.is_open()) diagnostic << sample.time << ',' << b.x() << ',' << b.y() << ',' << b.z() << ',' << rest << ','
            << (axes == 9 && heading.rejected) << ',' << acc_innovation << ',' << attitude.root_covariance.squaredNorm() << ','
            << bias.root_covariance.squaredNorm() << ',' << heading.norm << ',' << heading.dip << ',' << heading.delta << ','
            << (!norm_accepted || direction_rejected) << ',' << acceleration_norm << '\n';
    }
    result.elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    if (diagnostic.is_open()) {
        diagnostic.close();
        if (!diagnostic) throw std::runtime_error("Failed writing diagnostics");
    }
    return result;
}

std::string robust_parameters(const std::string& stage) {
    if (stage != "robust-rest" && stage != "robust-acc" && stage != "robust")
        throw std::runtime_error("Unknown robust stage");
    return "{\"stage\":\"" + stage + "\",\"attitude\":\"quaternion-only SR-UKF\",\"bias\":\"separate rest-only SR-UKF\","
        "\"gravity_gain_constraint\":\"horizontal left/world Gibbs error\",\"covariance_reset\":\"analytic left-Gibbs transport\","
        "\"unscented_alpha_squared\":1,\"unscented_beta\":2,\"unscented_kappa\":0,\"mrp_a\":0,\"mrp_f\":2,"
        "\"gyro_integration\":\"exact quaternion increment, midpoint gyro minus previous bias estimate\","
        "\"initial_attitude_sigma_rad\":0.05,\"initial_bias_sigma_deg_s\":0.5,"
        "\"gravity_observation\":\"two-dimensional predicted-gravity tangent plane\","
        "\"raw_acc_direction_sigma\":0.1,\"filtered_acc_direction_sigma_base\":0.002,"
        "\"filtered_acc_sigma_formula\":\"0.002*sqrt(max(1,6/observation_dt))\","
        "\"attitude_covariance_scope\":\"conditional on separate bias estimate; filtered-observation R is an effective model\","
        "\"acceleration_filter_enabled\":" + (stage != "robust-rest" ? "true" : "false") + ","
        "\"magnetic_rejection_enabled\":" + (stage == "robust" ? "true" : "false") + ","
        "\"acceleration_rejection_enabled\":" + (stage == "robust" ? "true" : "false") + ","
        "\"acc_tau_seconds\":3,\"correction_period_seconds\":0.05,\"rest_min_seconds\":1.5,\"rest_filter_tau_seconds\":0.5,"
        "\"rest_gyro_threshold_deg_s\":2,\"rest_absolute_filtered_gyro_limit_deg_s\":2,"
        "\"rest_acc_threshold_m_s2\":0.5,\"gyro_noise_rad_sqrt_s\":0.02,"
        "\"bias_random_walk_deg_s_sqrt_s\":0.01,\"mag_tau_seconds\":9,\"mag_norm_threshold_fraction\":0.1,"
        "\"rest_bias_target_sigma_deg_s\":0.03,\"rest_bias_effective_R\":\"P_target^2/Q_interval+P_target\","
        "\"bias_health_bound_deg_s\":2.1,\"mag_norm_dip_lowpass_tau_seconds\":0.05,\"mag_reference_tau_seconds\":20,"
        "\"mag_dip_threshold_degrees\":10,\"mag_good_dwell_seconds\":0.5,\"mag_new_field_rotation_seconds\":20,"
        "\"mag_new_field_min_rotation_deg_s\":20,\"heading_startup\":\"first usable magnetic sample\","
        "\"gravity_m_s2\":9.80665,\"acc_norm_reject_fraction\":0.3,\"acc_max_innovation_degrees\":60,"
        "\"acc_norm_gating_location\":\"measurement and rest updates; nonzero LPF input retained\"}";
}
