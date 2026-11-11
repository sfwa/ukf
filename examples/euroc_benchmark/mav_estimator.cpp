// MIT license; see the repository LICENSE.
// Nominal-thrust multirotor rotor-drag model, independently implemented from
// Leishman et al., Quadrotors and Accelerometers, IEEE CSM 34(1), 2014.
#include "mav_estimator.h"
#include "UKF/Core.h"
#include "UKF/Covariance.h"
#include "UKF/Integrator.h"
#include "UKF/MeasurementVector.h"
#include "UKF/StateVector.h"
#include "UKF/Types.h"
#include "mav_correction.h"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
namespace Mav {
using V3 = UKF::Vector<3>;
using V2 = UKF::Vector<2>;
using Q = UKF::Quaternion;
constexpr double Pi = 3.14159265358979323846, Rad = Pi / 180, G = 9.80665;
enum Keys { Attitude = 300, Velocity, Bias, Drag, Gravity, ZeroVelocity, RestGyro };
using State = UKF::StateVector<UKF::Field<Attitude, Q>, UKF::Field<Velocity, V3>, UKF::Field<Bias, V3>>;
using Measurement =
    UKF::DynamicMeasurementVector<UKF::Field<Drag, V2>, UKF::Field<Gravity, V2>, UKF::Field<RestGyro, V3>>;
struct Input {
    V3 omega, thrust_axis;
    Eigen::Matrix<real_t, 3, 2> lateral, tangent;
    double collective, drag;
    bool model_active, planar;
};
Q increment(const V3 &r) {
    double a = r.norm();
    return a < 1e-12 ? Q::Identity() : Q(Eigen::AngleAxis<real_t>(a, r / a));
}
V3 vec(const std::array<double, 3> &a) {
    return V3(a[0], a[1], a[2]);
}
V3 acceleration(const Q &q, const V3 &v, const Input &in) {
    if (!in.model_active)
        return V3::Zero();
    V3 derivative = q.conjugate() * V3(0, 0, -G) - in.omega.cross(v) + in.collective * in.thrust_axis -
                    in.drag * in.lateral * (in.lateral.transpose() * v);
    if (in.planar)
        derivative -= in.thrust_axis * in.thrust_axis.dot(derivative);
    return derivative;
}
} // namespace Mav
namespace UKF {
namespace Parameters {
template <> constexpr real_t AlphaSquared<Mav::State> = 1;
template <> constexpr real_t Beta<Mav::State> = 2;
template <> constexpr real_t Kappa<Mav::State> = 0;
} // namespace Parameters
template <> template <> Mav::State Mav::State::process_model<>(real_t dt, const Mav::Input &in) const {
    auto result = *this;
    const auto q = get_field<Mav::Attitude>();
    const auto v = get_field<Mav::Velocity>();
    Mav::Input actual = in;
    actual.omega -= get_field<Mav::Bias>();
    const auto midpoint = q * Mav::increment(actual.omega * (dt * .5));
    const auto v_midpoint = v + Mav::acceleration(q, v, actual) * (dt * .5);
    result.set_field<Mav::Attitude>((q * Mav::increment(actual.omega * dt)).normalized());
    result.set_field<Mav::Velocity>(v + Mav::acceleration(midpoint, v_midpoint, actual) * dt);
    return result;
}
template <>
template <>
Mav::V2 Mav::Measurement::expected_measurement<Mav::State, Mav::Drag, Mav::Input>(const Mav::State &s,
                                                                                  const Mav::Input &in) {
    return -in.drag * in.lateral.transpose() * s.get_field<Mav::Velocity>();
}
template <>
template <>
Mav::V2 Mav::Measurement::expected_measurement<Mav::State, Mav::Gravity, Mav::Input>(const Mav::State &s,
                                                                                     const Mav::Input &in) {
    return in.tangent.transpose() * (s.get_field<Mav::Attitude>().conjugate() * Mav::V3(0, 0, 1));
}
template <>
template <>
Mav::V3 Mav::Measurement::expected_measurement<Mav::State, Mav::RestGyro, Mav::Input>(const Mav::State &s,
                                                                                      const Mav::Input &) {
    return s.get_field<Mav::Bias>();
}
} // namespace UKF
namespace Mav {
struct Integrator {
    static State integrate(real_t dt, const State &s, const Input &in) {
        return s.process_model(dt, in);
    }
};
class Filter : public UKF::SquareRootCore<State, Measurement, Integrator> {
  public:
    void correct(bool freeze_bias) {
        const int m = innovation.size();
        if (!m)
            return;
        UKF::Matrix<9, Eigen::Dynamic> C(9, m);
        C.setZero();
        for (std::size_t i = 1; i < State::num_sigma(); ++i)
            C.noalias() += UKF::Parameters::Sigma_WCI<State> * w_prime.col(i) * z_prime.col(i).transpose();
        C.noalias() += UKF::Parameters::Sigma_WC0<State> * w_prime.col(0) * z_prime.col(0).transpose();
        const UKF::Matrix<9, 9> P = root_covariance * root_covariance.transpose();
        const Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> S =
            innovation_root_covariance * innovation_root_covariance.transpose();
        const auto update = mav_detail::constrained_update<9>(P, C, S, innovation, freeze_bias);
        const Q old_q = state.get_field<Attitude>();
        state.apply_delta(update.delta);
        root_covariance = mav_detail::recenter_root<9>(old_q, state.get_field<Attitude>(), update.posterior,
                                                       P.cwiseAbs().maxCoeff());
    }
};
class LowPass {
    V3 value = V3::Zero(), rate = V3::Zero();
    bool initialized = false;

  public:
    V3 update(const V3 &sample, double dt, double tau) {
        if (!initialized) {
            value = sample;
            initialized = true;
            return value;
        }
        const double h = dt / tau, d = std::exp(-h), c = std::cos(h), s = std::sin(h);
        const V3 e = value - sample;
        value = sample + d * ((c + s) * e + tau * s * rate);
        rate = d * ((-2 * s / tau) * e + (c - s) * rate);
        return value;
    }
};
struct RestDetector {
    LowPass gyro, acc, gyro_prefilter, acc_prefilter;
    V3 mean_gyro = V3::Zero(), mean_acc = V3::Zero();
    double duration = 0;
    bool update(const V3 &w, const V3 &a, double dt) {
        const V3 wf = gyro_prefilter.update(w, dt, std::sqrt(2.) / (2 * Pi * 2));
        const V3 af = acc_prefilter.update(a, dt, std::sqrt(2.) / (2 * Pi * 2));
        mean_gyro = gyro.update(wf, dt, .5);
        mean_acc = acc.update(af, dt, .5);
        bool quiet = (wf - mean_gyro).norm() < 2 * Rad && (af - mean_acc).norm() < .5 &&
                     mean_gyro.cwiseAbs().maxCoeff() < 10 * Rad;
        duration = quiet ? duration + dt : 0;
        return duration >= 1.5;
    }
};
Eigen::Matrix<real_t, 3, 2> tangent(const V3 &up) {
    Eigen::Matrix<real_t, 3, 2> b;
    const V3 axis = std::abs(up.z()) < .9 ? V3::UnitZ() : V3::UnitX();
    b.col(0) = up.cross(axis).normalized();
    b.col(1) = up.cross(b.col(0)).normalized();
    return b;
}
} // namespace Mav
void validate_mav_parameters(const MavParameters &p) {
    for (double value : {p.drag_rate, p.thrust_tilt_deg, p.gyro_noise, p.velocity_noise, p.drag_noise})
        if (!std::isfinite(value))
            throw std::runtime_error("MAV parameters must be finite");
    if (!(p.drag_rate > 0 && p.gyro_noise > 0 && p.velocity_noise > 0 && p.drag_noise > 0))
        throw std::runtime_error("MAV noise and drag parameters must be positive");
}
MavResult estimate_mav(const std::vector<MavSample> &samples, const MavParameters &p,
                       const std::string &diagnostics) {
    using namespace Mav;
    if (samples.size() < 2)
        throw std::runtime_error("At least two samples required");
    validate_mav_parameters(p);
    Filter filter;
    RestDetector rest_detector;
    LowPass acceleration_filter, vibration_lowpass;
    Q gyro_frame = Q::Identity();
    const double initial_acc_norm = vec(samples.front().acc).norm();
    if (!std::isfinite(initial_acc_norm) || initial_acc_norm <= 1e-6)
        throw std::runtime_error("A nonzero finite initial acceleration is required to initialize tilt");
    filter.state.set_field<Attitude>(
        Q::FromTwoVectors(vec(samples.front().acc) / initial_acc_norm, V3::UnitZ()));
    filter.state.set_field<Velocity>(V3::Zero());
    filter.root_covariance.setIdentity();
    filter.root_covariance.diagonal() << .05, .05, .05, 1, 1, 1, 5 * Rad, 5 * Rad, 5 * Rad;
    filter.state.set_field<Bias>(V3::Zero());
    Input in;
    in.thrust_axis = V3(std::cos(p.thrust_tilt_deg * Rad), 0, -std::sin(p.thrust_tilt_deg * Rad));
    in.lateral.col(0) = V3::UnitY();
    in.lateral.col(1) = in.thrust_axis.cross(V3::UnitY());
    in.tangent.setZero();
    in.drag = p.drag_rate;
    in.planar = p.planar;
    in.omega.setZero();
    in.model_active = false;
    in.collective = 0;
    filter.measurement_root_covariance.set_field<Drag>(V2::Constant(p.drag_noise));
    filter.measurement_root_covariance.set_field<Gravity>(V2::Constant(.02));
    std::ofstream diagnostic;
    if (!diagnostics.empty()) {
        diagnostic.open(diagnostics);
        if (!diagnostic)
            throw std::runtime_error("Cannot open diagnostics");
        diagnostic << "time,bx,by,bz,rest,vx,vy,vz,model_active,innovation_norm,covariance_trace,rotor_"
                      "vibration_rms\n"
                   << std::setprecision(17);
    }
    MavResult result;
    result.quaternions.reserve(samples.size());
    double clock = 0, innovation_norm = 0, vibration_energy = 0;
    bool rotors = false;
    V3 acc_sum = V3::Zero(), world_acc_sum = V3::Zero();
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto &sample = samples[i];
        const double dt = i ? sample.time - samples[i - 1].time : 0;
        const V3 w = vec(sample.gyro), a = vec(sample.acc);
        if (!std::isfinite(dt) || (i && dt <= 0) || !std::isfinite(a.norm()) || !std::isfinite(w.norm()))
            throw std::runtime_error("Invalid MAV input magnitude or interval");
        const V3 vibration = a - vibration_lowpass.update(a, dt, std::sqrt(2.) / (2 * Pi * 20));
        const double energy_gain = -std::expm1(-dt / .5);
        vibration_energy += energy_gain * (vibration.squaredNorm() - vibration_energy);
        const double vibration_rms = std::sqrt(std::max(0., vibration_energy));
        if (vibration_rms > .8)
            rotors = true;
        else if (vibration_rms < .5)
            rotors = false;
        bool rest = rest_detector.update(w, a, dt) && a.norm() > .7 * G && a.norm() < 1.3 * G;
        if (rotors) {
            rest = false;
            rest_detector.duration = 0;
        }
        const V3 b = filter.state.get_field<Bias>();
        in.omega = V3::Zero();
        if (i)
            in.omega = .5 * (w + vec(samples[i - 1].gyro));
        in.collective = a.dot(in.thrust_axis);
        bool active = rotors && in.collective > .5 * G && in.collective < 2 * G;
        in.drag = active ? p.drag_rate : 0;
        in.model_active = active && !p.gravity_baseline;
        if (i)
            gyro_frame = (gyro_frame * increment((in.omega - b) * dt)).normalized();
        if (rest)
            acceleration_filter = LowPass();
        const V3 acceleration_input = rest ? rest_detector.mean_acc : a;
        const V3 filtered_acc =
            gyro_frame.conjugate() * acceleration_filter.update(gyro_frame * acceleration_input, dt, 3);
        if (i) {
            filter.process_noise_root_covariance.setZero();
            filter.process_noise_root_covariance.diagonal() << p.gyro_noise, p.gyro_noise, p.gyro_noise,
                p.velocity_noise, p.velocity_noise, p.velocity_noise, .01 * Rad, .01 * Rad, .01 * Rad;
            filter.process_noise_root_covariance *= std::sqrt(dt);
            filter.a_priori_step(dt, in);
        }
        clock += dt;
        acc_sum += a * dt;
        world_acc_sum += (gyro_frame * a) * dt;
        if (i && clock + 1e-12 >= .05) {
            const double observation_dt = clock;
            clock = 0;
            const V3 observation_acc = acc_sum / observation_dt;
            const V3 deskewed_acc = gyro_frame.conjugate() * (world_acc_sum / observation_dt);
            acc_sum.setZero();
            world_acc_sum.setZero();
            if (rest) {
                const double target = std::pow(.03 * Rad, 2),
                             process = std::pow(.01 * Rad, 2) * observation_dt;
                filter.measurement_root_covariance.set_field<RestGyro>(
                    V3::Constant(std::sqrt(target * target / process + target)));
            }
            Measurement z;
            if (rest)
                z.set_field<RestGyro>(rest_detector.mean_gyro);
            if ((!active || p.gravity_baseline) && a.norm() > .7 * G && a.norm() < 1.3 * G) {
                const V3 up = filter.state.get_field<Attitude>().conjugate() * V3::UnitZ();
                in.tangent = tangent(up);
                const V3 gravity_observation = rotors ? filtered_acc : deskewed_acc;
                filter.measurement_root_covariance.set_field<Gravity>(V2::Constant(rotors ? .02 : .1));
                if (gravity_observation.norm() > 1e-6)
                    z.set_field<Gravity>(in.tangent.transpose() * gravity_observation.normalized());
            } else if (active && !p.gravity_baseline)
                z.set_field<Drag>(in.lateral.transpose() * observation_acc);
            if (z.size()) {
                filter.innovation_step(z, in);
                innovation_norm = filter.innovation.norm();
                filter.correct(false);
            }
        }
        const Q q = filter.state.get_field<Attitude>();
        const V3 current_bias = filter.state.get_field<Bias>(), v = filter.state.get_field<Velocity>();
        if (!q.coeffs().allFinite() || !v.allFinite() || !current_bias.allFinite() ||
            !filter.root_covariance.allFinite())
            throw std::runtime_error("Non-finite MAV filter state");
        result.quaternions.push_back({q.w(), q.x(), q.y(), q.z()});
        if (diagnostic)
            diagnostic << sample.time << ',' << current_bias.x() << ',' << current_bias.y() << ','
                       << current_bias.z() << ',' << rest << ',' << v.x() << ',' << v.y() << ',' << v.z()
                       << ',' << in.model_active << ',' << innovation_norm << ','
                       << filter.root_covariance.squaredNorm() << ',' << vibration_rms << '\n';
    }
    result.elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (diagnostic.is_open()) {
        diagnostic.close();
        if (!diagnostic)
            throw std::runtime_error("Failed writing diagnostics");
    }
    return result;
}
std::string mav_parameters(const MavParameters &p) {
    validate_mav_parameters(p);
    std::ostringstream s;
    s << std::setprecision(17);
    s << "{\"stage\":\"" << (p.gravity_baseline ? "mav-srukf-gravity" : "mav-srukf") << "\","
      << "\"state\":\"quaternion, body-frame velocity, gyro bias; nine-dimensional error covariance\","
      << "\"model\":\""
      << (p.gravity_baseline ? "gyro-frame filtered gravity control"
                             : "multirotor rotor drag and measured axial collective")
      << "\","
      << "\"drag_rate_s_inv\":" << p.drag_rate << ",\"thrust_axis_imu\":["
      << std::cos(p.thrust_tilt_deg * Mav::Rad) << ",0," << -std::sin(p.thrust_tilt_deg * Mav::Rad) << "],"
      << "\"thrust_tilt_deg\":" << p.thrust_tilt_deg
      << ",\"thrust_mount_source\":\"sensor-only development calibration rounded to 20 degrees; not "
         "published factory calibration\","
      << "\"gyro_noise_rad_sqrt_s\":" << p.gyro_noise << ",\"velocity_noise_m_s_sqrt_s\":" << p.velocity_noise
      << ",\"drag_observation_sigma_m_s2\":" << p.drag_noise
      << ",\"planar\":" << (p.planar ? "true" : "false") << ","
      << "\"initial_attitude_sigma_rad\":0.05,\"initial_velocity_sigma_m_s\":1,\"initial_bias_sigma_deg_s\":"
         "5,"
      << "\"gyro_bias_random_walk_deg_s_sqrt_s\":0.01,\"rest_bias_limit_deg_s\":10,\"rest_gyro_deviation_deg_"
         "s\":2,"
      << "\"rest_acc_deviation_m_s2\":0.5,\"rest_prefilter_hz\":2,\"rest_reference_tau_s\":0.5,\"rest_dwell_"
         "s\":1.5,"
      << "\"rest_bias_target_sigma_deg_s\":0.03,\"rest_bias_observation_variance\":\"P_target^2/"
         "Q_interval+P_target\","
      << "\"correction_period_s\":0.05,\"vibration_lowpass_hz\":20,\"vibration_energy_tau_s\":0.5,"
      << "\"rotor_on_rms_m_s2\":0.8,\"rotor_off_rms_m_s2\":0.5,\"model_collective_min_g\":0.5,\"model_"
         "collective_max_g\":2,"
      << "\"off_rotor_gravity_sigma\":0.1,\"control_filtered_gravity_sigma\":0.02,\"gravity_filter_tau_s\":3,"
      << "\"gravity_norm_min_g\":0.7,\"gravity_norm_max_g\":1.3,\"gravity_m_s2\":9.80665,"
      << "\"gyro_integration\":\"midpoint raw gyro minus each sigma-point gyro bias; exact quaternion "
         "exponential\","
      << "\"velocity_integration\":\"explicit midpoint body-frame Newtonian dynamics\","
      << "\"gravity_observation\":\"gyro-deskewed interval mean off rotor; gyro-frame low-pass for control "
         "in powered flight\","
      << "\"motion_bias_updates\":\"enabled in both model and gravity "
         "control\",\"zero_velocity_updates\":false,"
      << "\"gain_constraint\":\"zero world-yaw attitude gain row\",\"covariance_reset\":\"full nine-state "
         "block-diagonal Gibbs transport\","
      << "\"unscented_alpha_squared\":1,\"unscented_beta\":2,\"unscented_kappa\":0,\"mrp_a\":0,\"mrp_f\":2}";
    return s.str();
}
