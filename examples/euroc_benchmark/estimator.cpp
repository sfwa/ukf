// MIT license; see the repository LICENSE.
#include "mav_estimator.h"
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
// Errors are reported through return values; the benchmark uses no exceptions.
namespace {
// Parse a complete finite decimal number, as std::stod with a full-length and
// finiteness check did, but without exceptions. Out-of-range values fail.
bool parse_finite(const std::string &text, double &value) {
    const char *begin = text.c_str();
    char *end = nullptr;
    errno = 0;
    value = std::strtod(begin, &end);
    return end != begin && end == begin + text.size() && errno != ERANGE && std::isfinite(value);
}
int report_failure(const std::string &message) {
    std::cerr << "euroc_benchmark: " << message << '\n';
    return 1;
}
} // namespace
int main(int argc, char **argv) {
    std::string input, output, diagnostics, description, error;
    MavParameters p;
    bool axes_selected = false, filter_selected = false;
    for (int i = 1; i < argc; ++i) {
        std::string option = argv[i];
        if (option == "--help") {
            std::cout << "euroc_benchmark --input PATH --output PATH --filter mav-srukf --axes 6 "
                         "[--diagnostics PATH]\n";
            return 0;
        }
        if (i + 1 == argc)
            return report_failure("Missing option value");
        std::string value = argv[++i];
        double *parameter = nullptr;
        if (option == "--input")
            input = value;
        else if (option == "--output")
            output = value;
        else if (option == "--diagnostics")
            diagnostics = value;
        else if (option == "--filter" && (value == "mav-srukf" || value == "mav-srukf-gravity")) {
            p.gravity_baseline = value == "mav-srukf-gravity";
            filter_selected = true;
        } else if (option == "--axes" && value == "6") {
            axes_selected = true;
        } else if (option == "--describe-filter" &&
                   (value == "mav-srukf" || value == "mav-srukf-gravity")) {
            description = value;
            p.gravity_baseline = value == "mav-srukf-gravity";
        } else if (option == "--drag-rate")
            parameter = &p.drag_rate;
        else if (option == "--thrust-tilt")
            parameter = &p.thrust_tilt_deg;
        else if (option == "--gyro-noise")
            parameter = &p.gyro_noise;
        else if (option == "--velocity-noise")
            parameter = &p.velocity_noise;
        else if (option == "--drag-noise")
            parameter = &p.drag_noise;
        else if (option == "--planar" && value == "1")
            p.planar = true;
        else
            return report_failure("Unknown option/value: " + option + " " + value);
        if (parameter && !parse_finite(value, *parameter))
            return report_failure("Expected finite numeric parameter");
    }
    if (!description.empty()) {
        std::string parameters;
        if (!mav_parameters(p, parameters, error))
            return report_failure(error);
        std::cout << parameters << '\n';
        return 0;
    }
    if (!axes_selected || !filter_selected || input.empty() || output.empty() || input == output ||
        (!diagnostics.empty() && (diagnostics == input || diagnostics == output)))
        return report_failure("Required --filter, --axes 6, and distinct input/output/diagnostics paths");
    std::ifstream f(input);
    if (!f)
        return report_failure("Cannot open input");
    std::string line;
    if (!std::getline(f, line))
        return report_failure("Empty input");
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
    if (line != "time,gx,gy,gz,ax,ay,az,mx,my,mz")
        return report_failure("Expected time,gx,gy,gz,ax,ay,az,mx,my,mz header");
    std::vector<MavSample> samples;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        std::istringstream row(line);
        double v[10];
        for (int j = 0; j < 10; ++j) {
            std::string token;
            if (!std::getline(row, token, ','))
                return report_failure("Missing CSV column");
            if (!parse_finite(token, v[j]))
                return report_failure("Invalid CSV value");
        }
        std::string extra;
        if (std::getline(row, extra, ',') || line.back() == ',')
            return report_failure("Extra CSV column");
        if (!samples.empty() && v[0] <= samples.back().time)
            return report_failure("Non-increasing timestamp");
        samples.push_back({v[0], {v[1], v[2], v[3]}, {v[4], v[5], v[6]}});
    }
    MavResult r;
    if (!estimate_mav(samples, p, diagnostics, r, error))
        return report_failure(error);
    std::ofstream out(output);
    if (!out)
        return report_failure("Cannot open output");
    out << "time,qw,qx,qy,qz\n" << std::setprecision(17);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        out << samples[i].time;
        for (double q : r.quaternions[i])
            out << ',' << q;
        out << '\n';
    }
    out.close();
    if (!out)
        return report_failure("Failed writing output");
    std::cerr << std::setprecision(17) << "{\"filter\":\""
              << (p.gravity_baseline ? "mav-srukf-gravity" : "mav-srukf")
              << "\",\"samples\":" << samples.size() << ",\"elapsed_seconds\":" << r.elapsed_seconds
              << "}\n";
    return 0;
}
