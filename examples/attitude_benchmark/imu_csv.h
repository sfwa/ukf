/*
Reading recorded IMU samples and writing attitude estimates as CSV. This is
plumbing shared by the benchmark executables; it has nothing to do with the
filter itself.

Input header:  time,gx,gy,gz,ax,ay,az,mx,my,mz
Output header: time,qw,qx,qy,qz

Errors are reported by returning false with a message in "error". The
examples, like the library, do not use exceptions.
*/

#ifndef IMU_CSV_H
#define IMU_CSV_H

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "imu_sample.h"

/* Parse a complete, finite decimal number. */
inline bool parse_number(const std::string& text, double& value) {
    const char* begin = text.c_str();
    char* end = nullptr;
    errno = 0;
    value = std::strtod(begin, &end);
    return end != begin && end == begin + text.size() && errno != ERANGE && std::isfinite(value);
}

inline bool read_imu_csv(const std::string& path, std::vector<ImuSample>& samples, std::string& error) {
    std::ifstream file(path);
    if (!file) {
        error = "Cannot open input: " + path;
        return false;
    }

    std::string line;
    if (!std::getline(file, line)) {
        error = "Empty input";
        return false;
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    if (line != "time,gx,gy,gz,ax,ay,az,mx,my,mz") {
        error = "Expected header time,gx,gy,gz,ax,ay,az,mx,my,mz";
        return false;
    }

    samples.clear();
    std::size_t line_number = 1;
    while (std::getline(file, line)) {
        line_number++;
        const std::string where = " at line " + std::to_string(line_number);
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        /* Exactly ten comma-separated numbers per row. */
        std::istringstream row(line);
        std::string token;
        double values[10];
        for (double& value : values) {
            if (!std::getline(row, token, ',') || !parse_number(token, value)) {
                error = "Expected ten finite numbers" + where;
                return false;
            }
        }
        if (std::getline(row, token, ',') || line.back() == ',') {
            error = "Extra CSV column" + where;
            return false;
        }
        if (!samples.empty() && values[0] <= samples.back().time) {
            error = "Time must increase" + where;
            return false;
        }

        ImuSample sample;
        sample.time = values[0];
        sample.gyroscope << values[1], values[2], values[3];
        sample.accelerometer << values[4], values[5], values[6];
        sample.magnetometer << values[7], values[8], values[9];
        samples.push_back(sample);
    }

    if (samples.size() < 2) {
        error = "At least two samples are required";
        return false;
    }
    return true;
}

/* One row per sample: time, then the body-to-world attitude quaternion. */
inline bool write_attitude_csv(const std::string& path, const std::vector<ImuSample>& samples,
        const std::vector<UKF::Quaternion>& attitudes, std::string& error) {
    std::ofstream file(path);
    file << "time,qw,qx,qy,qz\n" << std::setprecision(17);
    for (std::size_t i = 0; i < samples.size(); i++) {
        const UKF::Quaternion& q = attitudes[i];
        file << samples[i].time << ',' << q.w() << ',' << q.x() << ',' << q.y() << ',' << q.z() << '\n';
    }
    file.close();
    if (!file) {
        error = "Cannot write output: " + path;
        return false;
    }
    return true;
}

#endif
