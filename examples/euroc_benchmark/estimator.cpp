/*
Command-line program for the EuRoC MAV benchmark. It reads recorded IMU
samples from a CSV file, runs the multirotor attitude estimator over them,
and writes one attitude quaternion per sample. The filter itself is in
mav_estimator.h.

Usage:
    euroc_benchmark --input sensors.csv --output attitude.csv
        --filter mav-srukf|mav-srukf-gravity --axes 6 [--diagnostics states.csv]
    euroc_benchmark --describe-filter mav-srukf|mav-srukf-gravity

The mav-srukf-gravity filter is the same estimator with its rotor model
switched off, as a control. Failures are reported on standard error with a
non-zero exit code.
*/

#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "mav_estimator.h"
#include "imu_csv.h"

static int fail(const std::string& message) {
    std::cerr << "euroc_benchmark: " << message << '\n';
    return 1;
}

/* Print the filter configuration as JSON, for the benchmark reports. */
static void describe_filter(const std::string& filter, const MAV_Tuning& tuning) {
    std::cout << std::setprecision(17)
        << "{\"stage\":\"" << filter << "\""
        << ",\"core\":\"UKF::SquareRootCore\""
        << ",\"state\":\"attitude quaternion, IMU-frame velocity, gyroscope bias\""
        << ",\"measurements\":\""
        << (tuning.use_rotor_model ? "rotor drag in powered flight, zero velocity otherwise" : "zero velocity")
        << "; gyroscope at rest\""
        << ",\"integrator\":\"UKF::IntegratorRK4\""
        << ",\"rotor_model\":" << (tuning.use_rotor_model ? "true" : "false")
        << ",\"drag_rate_s_inv\":" << MAV_DRAG_RATE
        << ",\"thrust_tilt_rad\":" << MAV_THRUST_TILT
        << ",\"thrust_axis_imu\":[" << MAV_ThrustAxis[0] << ',' << MAV_ThrustAxis[1] << ',' << MAV_ThrustAxis[2] << "]"
        << ",\"thrust_axis_source\":\"calibrate_mount.py on the development recordings\""
        << ",\"initial_tilt_rad\":" << tuning.initial_tilt
        << ",\"initial_velocity_m_s\":" << tuning.initial_velocity
        << ",\"initial_gyro_bias_rad_s\":" << tuning.initial_gyro_bias
        << ",\"gyroscope_noise_rad_sqrt_s\":" << tuning.gyroscope_noise
        << ",\"gyroscope_scale_error\":" << tuning.gyroscope_scale_error
        << ",\"accelerometer_noise_m_s_sqrt_s\":" << tuning.accelerometer_noise
        << ",\"rotor_model_noise_m_s_sqrt_s\":" << tuning.rotor_model_noise
        << ",\"gyro_bias_walk_rad_s_sqrt_s\":" << tuning.gyro_bias_walk
        << ",\"rotor_drag_noise_m_s2_sqrt_s\":" << tuning.rotor_drag_noise
        << ",\"zero_velocity_noise_m_s_sqrt_s\":" << tuning.zero_velocity_noise
        << ",\"rest_gyroscope_noise_rad_s_sqrt_s\":" << tuning.rest_gyroscope_noise
        << ",\"heading_sigma_limit_rad\":" << tuning.heading_sigma_limit
        << ",\"rest_cutoff_frequency_hz\":" << tuning.rest_cutoff_frequency
        << ",\"rest_averaging_time_s\":" << tuning.rest.averaging_time
        << ",\"rest_gyroscope_change_rad_s\":" << tuning.rest.gyroscope_change
        << ",\"rest_accelerometer_change_m_s2\":" << tuning.rest.accelerometer_change
        << ",\"rest_accelerometer_minimum_m_s2\":" << tuning.rest.accelerometer_minimum
        << ",\"rest_accelerometer_maximum_m_s2\":" << tuning.rest.accelerometer_maximum
        << ",\"rest_gyroscope_bias_rad_s\":" << tuning.rest.gyroscope_bias
        << ",\"rest_minimum_time_s\":" << tuning.rest.minimum_time
        << ",\"rotor_vibration_cutoff_frequency_hz\":" << tuning.rotors.cutoff_frequency
        << ",\"rotor_vibration_averaging_time_s\":" << tuning.rotors.averaging_time
        << ",\"rotors_on_vibration_m_s2\":" << tuning.rotors.on_vibration
        << ",\"rotors_off_vibration_m_s2\":" << tuning.rotors.off_vibration
        << ",\"minimum_thrust_m_s2\":" << tuning.minimum_thrust
        << ",\"maximum_thrust_m_s2\":" << tuning.maximum_thrust
        << "}\n";
}

int main(int argc, char** argv) {
    /* Every option takes exactly one value. */
    std::map<std::string, std::string> options;
    for (int i = 1; i < argc; i += 2) {
        std::string option = argv[i];
        if (option == "--help") {
            std::cout <<
                "Usage: euroc_benchmark --input sensors.csv --output attitude.csv\n"
                "           --filter mav-srukf|mav-srukf-gravity --axes 6 [--diagnostics states.csv]\n"
                "       euroc_benchmark --describe-filter mav-srukf|mav-srukf-gravity\n"
                "Input:  time (s), gyroscope (rad/s), accelerometer (m/s^2, +g upwards at rest)\n"
                "        and three unused magnetometer columns.\n"
                "Output: time,qw,qx,qy,qz. The quaternion rotates the IMU frame to a world\n"
                "        frame with +Z up; the heading is arbitrary.\n";
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

    /* The control filter is the same estimator without its rotor model. */
    std::string filter = options.count("--describe-filter") ? options["--describe-filter"] : options["--filter"];
    if (filter != "mav-srukf" && filter != "mav-srukf-gravity") {
        return fail("Required: --filter mav-srukf|mav-srukf-gravity");
    }
    MAV_Tuning tuning;
    tuning.use_rotor_model = filter == "mav-srukf";

    if (options.count("--describe-filter")) {
        describe_filter(filter, tuning);
        return 0;
    }

    std::string input = options["--input"], output = options["--output"];
    std::string diagnostics = options["--diagnostics"];
    if (input.empty() || output.empty() || options["--axes"] != "6") {
        return fail("Required: --input PATH --output PATH --axes 6");
    }
    if (input == output || diagnostics == input || diagnostics == output) {
        return fail("Input, output and diagnostics paths must differ");
    }

    std::string error;
    std::vector<ImuSample> samples;
    if (!read_imu_csv(input, samples, error)) {
        return fail(error);
    }

    MAV_Estimator estimator(tuning);
    if (!estimator.initialise(samples[0])) {
        return fail("The first accelerometer reading is zero, so the initial tilt is unknown");
    }

    std::ofstream diagnostics_file;
    if (!diagnostics.empty()) {
        diagnostics_file.open(diagnostics);
        diagnostics_file << "time,bias_x,bias_y,bias_z,velocity_x,velocity_y,velocity_z,"
            "rest,rotor_model_active,rotor_vibration\n" << std::setprecision(17);
    }

    /*
    The first attitude comes from the first sample alone; each later one
    includes the sample at its own time. The time spent in the estimator is
    measured, excluding file input and output (unless diagnostics are being
    written).
    */
    std::vector<UKF::Quaternion> attitudes;
    attitudes.reserve(samples.size());
    auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < samples.size(); i++) {
        if (i > 0 && !estimator.update(samples[i - 1], samples[i])) {
            return fail("Filter step failed at sample " + std::to_string(i));
        }
        attitudes.push_back(estimator.attitude());

        if (!diagnostics.empty()) {
            UKF::Vector<3> bias = estimator.gyro_bias();
            UKF::Vector<3> velocity = estimator.velocity();
            diagnostics_file << samples[i].time << ','
                << bias[0] << ',' << bias[1] << ',' << bias[2] << ','
                << velocity[0] << ',' << velocity[1] << ',' << velocity[2] << ','
                << estimator.resting() << ','
                << estimator.rotor_model_active() << ','
                << estimator.rotor_vibration() << '\n';
        }
    }
    double elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

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
        << "{\"filter\":\"" << filter << "\",\"samples\":" << samples.size()
        << ",\"updates\":" << samples.size() - 1
        << ",\"elapsed_seconds\":" << elapsed_seconds << "}\n";

    return 0;
}
