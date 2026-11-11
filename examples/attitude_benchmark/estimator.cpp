// MIT license; see the repository LICENSE.
// A deliberately small, causal attitude estimator for dataset comparisons.
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"

enum Fields { Attitude, GyroBias, AccDirection, MagDirection };
using Vec3 = UKF::Vector<3>;
using Quaternion = UKF::Quaternion;
using State = UKF::StateVector<UKF::Field<Attitude, Quaternion>,
    UKF::Field<GyroBias, Vec3>>;
using Measurement = UKF::DynamicMeasurementVector<UKF::Field<AccDirection, Vec3>,
    UKF::Field<MagDirection, Vec3>>;

struct Input {
    Vec3 gyro;
    Vec3 magnetic_world;
};

namespace UKF {
namespace Parameters {
template <> constexpr real_t AlphaSquared<State> = 1.0;
template <> constexpr real_t Beta<State> = 2.0;
template <> constexpr real_t Kappa<State> = 0.0;
}

// q rotates body into a right-handed world with gravity along -Z.
// The exact exponential integrates the constant midpoint angular rate.
template <> template <>
State State::process_model<>(real_t dt, const Input& input) const {
    State result = *this;
    const Vec3 rotation = (input.gyro - get_field<GyroBias>()) * dt;
    const real_t angle = rotation.norm();
    Quaternion increment = Quaternion::Identity();
    if (angle > 1e-12) {
        increment = Quaternion(Eigen::AngleAxis<real_t>(angle, rotation / angle));
    }
    result.set_field<Attitude>((get_field<Attitude>() * increment).normalized());
    return result;
}

template <> template <>
Vec3 Measurement::expected_measurement<State, AccDirection, Input>(
    const State& state, const Input&) {
    return state.get_field<Attitude>().conjugate() * Vec3(0, 0, 1);
}

template <> template <>
Vec3 Measurement::expected_measurement<State, MagDirection, Input>(
    const State& state, const Input& input) {
    return state.get_field<Attitude>().conjugate() * input.magnetic_world;
}
}

// SquareRootCore always dispatches through an integrator. This adapter gives
// both cores the same exact discrete process model without changing the library.
struct ExponentialIntegrator {
    static State integrate(real_t dt, const State& state, const Input& input) {
        return state.process_model(dt, input);
    }
};
using CovarianceFilter = UKF::Core<State, Measurement, ExponentialIntegrator>;
using SquareRootFilter = UKF::SquareRootCore<State, Measurement, ExponentialIntegrator>;

struct Sample {
    double time;
    Vec3 gyro, acc, mag;
};
using Samples = std::vector<Sample, Eigen::aligned_allocator<Sample>>;
using Quaternions = std::vector<Quaternion, Eigen::aligned_allocator<Quaternion>>;

Samples read_samples(const std::string& path, int axes) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open input: " + path);
    std::string line;
    if (!std::getline(file, line)) throw std::runtime_error("Empty input");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != "time,gx,gy,gz,ax,ay,az,mx,my,mz") {
        throw std::runtime_error("Expected header time,gx,gy,gz,ax,ay,az,mx,my,mz");
    }
    Samples samples;
    std::size_t line_number = 1;
    while (std::getline(file, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) throw std::runtime_error("Blank input row at line " + std::to_string(line_number));
        std::istringstream row(line);
        double values[10];
        for (int i = 0; i < 10; ++i) {
            std::string token;
            if (!std::getline(row, token, ',')) throw std::runtime_error("Missing CSV value at line " + std::to_string(line_number));
            std::size_t used = 0;
            values[i] = std::stod(token, &used);
            if (used != token.size() || !std::isfinite(values[i])) {
                throw std::runtime_error("Invalid CSV number at line " + std::to_string(line_number));
            }
        }
        std::string extra;
        if (std::getline(row, extra, ',') || line.back() == ',') {
            throw std::runtime_error("Extra CSV column at line " + std::to_string(line_number));
        }
        Sample sample{values[0], Vec3(values[1], values[2], values[3]),
            Vec3(values[4], values[5], values[6]), Vec3(values[7], values[8], values[9])};
        if ((!samples.empty() && sample.time <= samples.back().time) ||
            sample.acc.norm() < 1e-12 || (axes == 9 && sample.mag.norm() < 1e-12)) {
            throw std::runtime_error("Non-increasing time or zero sensor vector at line " + std::to_string(line_number));
        }
        samples.push_back(sample);
    }
    if (samples.size() < 2) throw std::runtime_error("At least two samples are required");
    return samples;
}

Quaternion initial_attitude(const Sample& first, int axes) {
    Quaternion q = Quaternion::FromTwoVectors(first.acc.normalized(), Vec3(0, 0, 1));
    if (axes == 9) {
        const Vec3 magnetic = q * first.mag.normalized();
        if (magnetic.head<2>().norm() < 1e-6) throw std::runtime_error("Initial magnetic field has no horizontal component");
        q = Quaternion(Eigen::AngleAxis<real_t>(-std::atan2(magnetic.y(), magnetic.x()), Vec3(0, 0, 1))) * q;
    }
    return q.normalized();
}

template <typename Filter>
Quaternions estimate(const Samples& samples, int axes, double& elapsed_seconds) {
    Filter filter;
    const Quaternion initial = initial_attitude(samples.front(), axes);
    filter.state.template set_field<Attitude>(initial);
    filter.state.template set_field<GyroBias>(Vec3::Zero());
    Input input{Vec3::Zero(), axes == 9 ? Vec3(initial * samples.front().mag.normalized()) : Vec3(1, 0, 0)};

    // Fixed tuning for every dataset and both filter implementations.
    // Initial attitude SD 0.1 rad, gyro bias SD 0.01 rad/s.
    State::CovarianceMatrix initial_root = State::CovarianceMatrix::Zero();
    initial_root.diagonal() << 0.1, 0.1, 0.1, 0.01, 0.01, 0.01;
    // Normalized direction measurement SDs: acc 0.1, magnetic field 0.1.
    if constexpr (std::is_same_v<Filter, CovarianceFilter>) {
        filter.covariance = initial_root * initial_root.transpose();
        filter.measurement_covariance.setConstant(0.01);
    } else {
        filter.root_covariance = initial_root;
        filter.measurement_root_covariance.setConstant(0.1);
    }

    Quaternions output;
    output.reserve(samples.size());
    output.push_back(initial);
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 1; i < samples.size(); ++i) {
        const double dt = samples[i].time - samples[i - 1].time;
        input.gyro = (samples[i].gyro + samples[i - 1].gyro) * 0.5;
        // Discrete process covariance Q = diag(sigma_g^2, sigma_b^2) * dt.
        // sigma_g = 0.02 rad/sqrt(s); sigma_b = 0.0002 rad/s/sqrt(s).
        // Isotropic attitude noise is invariant under the body/world rotation.
        // Bias coupling is propagated through the sigma points; this is a
        // first-order discretization of the additive driving noise.
        State::CovarianceMatrix root_q = State::CovarianceMatrix::Zero();
        root_q.diagonal() << 0.02, 0.02, 0.02, 0.0002, 0.0002, 0.0002;
        root_q *= std::sqrt(dt);
        if constexpr (std::is_same_v<Filter, CovarianceFilter>) {
            filter.process_noise_covariance = root_q * root_q.transpose();
        } else {
            filter.process_noise_root_covariance = root_q;
        }
        Measurement measurement;
        measurement.set_field<AccDirection>(samples[i].acc.normalized());
        if (axes == 9) measurement.set_field<MagDirection>(samples[i].mag.normalized());
        filter.step(dt, measurement, input);
        Quaternion q = filter.state.template get_field<Attitude>();
        if (!q.coeffs().allFinite() || std::abs(q.norm() - 1.0) > 1e-5 ||
            !filter.state.template get_field<GyroBias>().allFinite()) {
            throw std::runtime_error("Non-finite state or invalid quaternion at sample " + std::to_string(i));
        }
        output.push_back(q.normalized());
    }
    elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return output;
}

int main(int argc, char** argv) {
    try {
        std::string input_path, output_path, filter_name;
        int axes = 0;
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--help") {
                std::cout << "Usage: attitude_benchmark --input sensors.csv --output attitude.csv --filter ukf|srukf --axes 6|9\n"
                    "Input: time (s), gyro (rad/s), accelerometer (specific force, +Z at rest), magnetic field.\n"
                    "Output: time,qw,qx,qy,qz; body-to-world, Z up, magnetic north +X.\n";
                return 0;
            }
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + option);
            const std::string value = argv[++i];
            if (option == "--input") input_path = value;
            else if (option == "--output") output_path = value;
            else if (option == "--filter") filter_name = value;
            else if (option == "--axes" && (value == "6" || value == "9")) axes = std::stoi(value);
            else throw std::runtime_error("Unknown option or invalid value: " + option + " " + value);
        }
        if (input_path.empty() || output_path.empty() || axes == 0 ||
            (filter_name != "ukf" && filter_name != "srukf")) {
            throw std::runtime_error("Required: --input PATH --output PATH --filter ukf|srukf --axes 6|9");
        }
        if (input_path == output_path) throw std::runtime_error("Input and output paths must differ");
        const Samples samples = read_samples(input_path, axes);
        double elapsed = 0;
        const Quaternions quats = filter_name == "ukf"
            ? estimate<CovarianceFilter>(samples, axes, elapsed)
            : estimate<SquareRootFilter>(samples, axes, elapsed);
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("Cannot open output: " + output_path);
        output << "time,qw,qx,qy,qz\n" << std::setprecision(17);
        for (std::size_t i = 0; i < samples.size(); ++i) {
            output << samples[i].time << ',' << quats[i].w() << ',' << quats[i].x() << ','
                << quats[i].y() << ',' << quats[i].z() << '\n';
        }
        output.close();
        if (!output) throw std::runtime_error("Failed writing output: " + output_path);
        std::cerr << std::setprecision(17) << "{\"filter\":\"" << filter_name << "\",\"axes\":" << axes
            << ",\"samples\":" << samples.size() << ",\"updates\":" << samples.size() - 1
            << ",\"elapsed_seconds\":" << elapsed << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "attitude_benchmark: " << error.what() << '\n';
        return 1;
    }
}
