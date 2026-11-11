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

// Both functions report failure by returning false with a message in error;
// no exceptions are used.
bool estimate_robust(const std::vector<RobustSample>& samples, int axes,
                     const std::string& stage, const std::string& diagnostics,
                     RobustResult& result, std::string& error);
bool robust_parameters(const std::string& stage, std::string& description);
