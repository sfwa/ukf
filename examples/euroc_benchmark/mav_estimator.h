#pragma once
#include <array>
#include <string>
#include <vector>
struct MavSample {
    double time;
    std::array<double, 3> gyro, acc;
};
struct MavResult {
    std::vector<std::array<double, 4>> quaternions;
    double elapsed_seconds;
};
struct MavParameters {
    double drag_rate = .25;
    double thrust_tilt_deg = 20;
    double gyro_noise = .002;
    double velocity_noise = .2;
    double drag_noise = .3;
    bool planar = false;
    bool gravity_baseline = false;
};
// Both functions report failure by returning false with a message in error;
// no exceptions are used.
bool estimate_mav(const std::vector<MavSample> &samples, const MavParameters &parameters,
                  const std::string &diagnostics, MavResult &result, std::string &error);
bool mav_parameters(const MavParameters &parameters, std::string &description, std::string &error);
