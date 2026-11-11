#pragma once

#include <array>
#include <string>
#include <vector>

struct RobustSample {
    double time;
    std::array<double, 3> gyro, acc, mag;
};

struct RobustResult {
    std::vector<std::array<double, 4>> quaternions;
    double elapsed_seconds;
};

RobustResult estimate_robust(const std::vector<RobustSample>& samples, int axes,
                             const std::string& stage, const std::string& diagnostics);
std::string robust_parameters(const std::string& stage);
