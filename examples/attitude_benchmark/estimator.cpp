/*
Command-line program for the attitude benchmark. It reads recorded IMU
samples from a CSV file, runs the attitude estimator over them, and writes
one attitude quaternion per sample. The filter itself is in
attitude_estimator.h.

Usage:
    attitude_benchmark --input sensors.csv --output attitude.csv
        --filter ukf|srukf --axes 6|9 [--diagnostics states.csv]
    attitude_benchmark --describe-filter ukf|srukf

Failures are reported on standard error with a non-zero exit code.
*/

#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <vector>
#include "attitude_estimator.h"
#include "imu_csv.h"

static int fail(const std::string& message) {
    std::cerr << "attitude_benchmark: " << message << '\n';
    return 1;
}

/* Print the filter configuration as JSON, for the benchmark reports. */
static void describe_filter(const std::string& filter) {
    Attitude_Tuning tuning;
    std::cout << std::setprecision(17)
        << "{\"stage\":\"" << filter << "\""
        << ",\"core\":\"" << (filter == "ukf" ? "UKF::Core" : "UKF::SquareRootCore") << "\""
        << ",\"state\":\"attitude quaternion, world-frame velocity, gyroscope bias\""
        << ",\"measurements\":\"zero velocity; gyroscope at rest; magnetic heading (separate filter)\""
        << ",\"integrator\":\"UKF::IntegratorRK4\""
        << ",\"initial_tilt_rad\":" << tuning.initial_tilt
        << ",\"initial_velocity_m_s\":" << tuning.initial_velocity
        << ",\"initial_gyro_bias_rad_s\":" << tuning.initial_gyro_bias
        << ",\"gyroscope_noise_rad_sqrt_s\":" << tuning.gyroscope_noise
        << ",\"gyroscope_scale_error\":" << tuning.gyroscope_scale_error
        << ",\"accelerometer_noise_m_s_sqrt_s\":" << tuning.accelerometer_noise
        << ",\"gyro_bias_walk_rad_s_sqrt_s\":" << tuning.gyro_bias_walk
        << ",\"zero_velocity_noise_m_s_sqrt_s\":" << tuning.zero_velocity_noise
        << ",\"rest_gyroscope_noise_rad_s_sqrt_s\":" << tuning.rest_gyroscope_noise
        << ",\"heading_sigma_limit_rad\":" << tuning.heading_sigma_limit
        << ",\"rest_averaging_time_s\":" << tuning.rest.averaging_time
        << ",\"rest_gyroscope_change_rad_s\":" << tuning.rest.gyroscope_change
        << ",\"rest_accelerometer_change_m_s2\":" << tuning.rest.accelerometer_change
        << ",\"rest_accelerometer_minimum_m_s2\":" << tuning.rest.accelerometer_minimum
        << ",\"rest_accelerometer_maximum_m_s2\":" << tuning.rest.accelerometer_maximum
        << ",\"rest_gyroscope_bias_rad_s\":" << tuning.rest.gyroscope_bias
        << ",\"rest_minimum_time_s\":" << tuning.rest.minimum_time
        << ",\"initial_heading_rad\":" << tuning.initial_heading
        << ",\"heading_walk_rad_sqrt_s\":" << tuning.heading_walk
        << ",\"magnetic_heading_noise_rad_sqrt_s\":" << tuning.magnetic_heading_noise
        << ",\"magnetic_strength_change\":" << tuning.magnetic_field.strength_change
        << ",\"magnetic_dip_change_rad\":" << tuning.magnetic_field.dip_change
        << ",\"magnetic_smoothing_time_s\":" << tuning.magnetic_field.smoothing_time
        << ",\"magnetic_settle_time_s\":" << tuning.magnetic_field.settle_time
        << ",\"magnetic_reference_time_s\":" << tuning.magnetic_field.reference_time
        << ",\"magnetic_new_field_time_s\":" << tuning.magnetic_field.new_field_time
        << ",\"magnetic_new_field_rotation_rad_s\":" << tuning.magnetic_field.new_field_rotation
        << "}\n";
}

/*
Run one estimator over all the samples. The first attitude comes from the
first sample alone; each later one includes the sample at its own time. The
time spent in the estimator is measured, excluding file input and output
(unless diagnostics are being written).
*/
template <typename AttitudeFilter, typename HeadingFilter>
static bool run_estimator(const std::vector<ImuSample>& samples, bool with_magnetometer,
        std::ostream* diagnostics, std::vector<UKF::Quaternion>& attitudes, double& elapsed_seconds,
        std::string& error) {
    Attitude_Estimator<AttitudeFilter, HeadingFilter> estimator;
    if (!estimator.initialise(samples[0], with_magnetometer)) {
        error = "The first accelerometer reading is zero, so the initial tilt is unknown";
        return false;
    }

    attitudes.clear();
    attitudes.reserve(samples.size());
    auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < samples.size(); i++) {
        if (i > 0 && !estimator.update(samples[i - 1], samples[i])) {
            error = "Filter step failed at sample " + std::to_string(i);
            return false;
        }
        attitudes.push_back(estimator.attitude());

        if (diagnostics) {
            UKF::Vector<3> bias = estimator.gyro_bias();
            UKF::Vector<3> velocity = estimator.velocity();
            *diagnostics << samples[i].time << ','
                << bias[0] << ',' << bias[1] << ',' << bias[2] << ','
                << velocity[0] << ',' << velocity[1] << ',' << velocity[2] << ','
                << estimator.resting() << ','
                << estimator.heading_offset() << ','
                << estimator.magnetometer_rejected() << '\n';
        }
    }
    elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    return true;
}

int main(int argc, char** argv) {
    /* Every option takes exactly one value. */
    std::map<std::string, std::string> options;
    for (int i = 1; i < argc; i += 2) {
        std::string option = argv[i];
        if (option == "--help") {
            std::cout <<
                "Usage: attitude_benchmark --input sensors.csv --output attitude.csv\n"
                "           --filter ukf|srukf --axes 6|9 [--diagnostics states.csv]\n"
                "       attitude_benchmark --describe-filter ukf|srukf\n"
                "Input:  time (s), gyroscope (rad/s), accelerometer (m/s^2, +g upwards at rest)\n"
                "        and magnetometer (any unit; ignored with --axes 6).\n"
                "Output: time,qw,qx,qy,qz. The quaternion rotates the body frame to a world\n"
                "        frame with +Z up and, with --axes 9, magnetic north along +X.\n";
            return 0;
        }
        if (option != "--input" && option != "--output" && option != "--filter" && option != "--axes" &&
                option != "--diagnostics" && option != "--describe-filter") {
            return fail("Unknown option " + option);
        }
        if (i + 1 >= argc) {
            return fail("Missing value for " + option);
        }
        options[option] = argv[i + 1];
    }

    if (options.count("--describe-filter")) {
        std::string filter = options["--describe-filter"];
        if (filter != "ukf" && filter != "srukf") {
            return fail("Unknown filter " + filter);
        }
        describe_filter(filter);
        return 0;
    }

    std::string input = options["--input"], output = options["--output"];
    std::string filter = options["--filter"], diagnostics = options["--diagnostics"];
    int axes = options["--axes"] == "6" ? 6 : options["--axes"] == "9" ? 9 : 0;
    if (input.empty() || output.empty() || axes == 0 || (filter != "ukf" && filter != "srukf")) {
        return fail("Required: --input PATH --output PATH --filter ukf|srukf --axes 6|9");
    }
    if (input == output || diagnostics == input || diagnostics == output) {
        return fail("Input, output and diagnostics paths must differ");
    }

    std::string error;
    std::vector<ImuSample> samples;
    if (!read_imu_csv(input, samples, error)) {
        return fail(error);
    }

    std::ofstream diagnostics_file;
    if (!diagnostics.empty()) {
        diagnostics_file.open(diagnostics);
        diagnostics_file << "time,bias_x,bias_y,bias_z,velocity_x,velocity_y,velocity_z,"
            "rest,heading_offset,magnetometer_rejected\n" << std::setprecision(17);
    }

    std::vector<UKF::Quaternion> attitudes;
    double elapsed_seconds = 0;
    bool succeeded = filter == "ukf" ?
        run_estimator<Attitude_UKF, Heading_UKF>(
            samples, axes == 9, diagnostics.empty() ? nullptr : &diagnostics_file,
            attitudes, elapsed_seconds, error) :
        run_estimator<Attitude_SRUKF, Heading_SRUKF>(
            samples, axes == 9, diagnostics.empty() ? nullptr : &diagnostics_file,
            attitudes, elapsed_seconds, error);
    if (!succeeded) {
        return fail(error);
    }

    if (!diagnostics.empty()) {
        diagnostics_file.close();
        if (!diagnostics_file) {
            return fail("Cannot write diagnostics: " + diagnostics);
        }
    }
    if (!write_attitude_csv(output, samples, attitudes, error)) {
        return fail(error);
    }

    /* Report the time spent in the estimator, as JSON. */
    std::cerr << std::setprecision(17)
        << "{\"filter\":\"" << filter << "\",\"axes\":" << axes
        << ",\"samples\":" << samples.size() << ",\"updates\":" << samples.size() - 1
        << ",\"elapsed_seconds\":" << elapsed_seconds << "}\n";

    return 0;
}
