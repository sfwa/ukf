// MIT license; see the repository LICENSE.
#include "mav_estimator.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
double number(const std::string &text) {
    std::size_t used = 0;
    double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value))
        throw std::runtime_error("Expected finite numeric parameter");
    return value;
}
int main(int argc, char **argv) {
    try {
        std::string input, output, diagnostics, description;
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
                throw std::runtime_error("Missing option value");
            std::string value = argv[++i];
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
                p.drag_rate = number(value);
            else if (option == "--thrust-tilt")
                p.thrust_tilt_deg = number(value);
            else if (option == "--gyro-noise")
                p.gyro_noise = number(value);
            else if (option == "--velocity-noise")
                p.velocity_noise = number(value);
            else if (option == "--drag-noise")
                p.drag_noise = number(value);
            else if (option == "--planar" && value == "1")
                p.planar = true;
            else
                throw std::runtime_error("Unknown option/value: " + option + " " + value);
        }
        if (!description.empty()) {
            std::cout << mav_parameters(p) << '\n';
            return 0;
        }
        if (!axes_selected || !filter_selected || input.empty() || output.empty() || input == output ||
            (!diagnostics.empty() && (diagnostics == input || diagnostics == output)))
            throw std::runtime_error(
                "Required --filter, --axes 6, and distinct input/output/diagnostics paths");
        std::ifstream f(input);
        if (!f)
            throw std::runtime_error("Cannot open input");
        std::string line;
        if (!std::getline(f, line))
            throw std::runtime_error("Empty input");
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line != "time,gx,gy,gz,ax,ay,az,mx,my,mz")
            throw std::runtime_error("Expected time,gx,gy,gz,ax,ay,az,mx,my,mz header");
        std::vector<MavSample> samples;
        while (std::getline(f, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            std::istringstream row(line);
            double v[10];
            for (int j = 0; j < 10; ++j) {
                std::string token;
                if (!std::getline(row, token, ','))
                    throw std::runtime_error("Missing CSV column");
                std::size_t used;
                v[j] = std::stod(token, &used);
                if (used != token.size() || !std::isfinite(v[j]))
                    throw std::runtime_error("Invalid CSV value");
            }
            std::string extra;
            if (std::getline(row, extra, ',') || line.back() == ',')
                throw std::runtime_error("Extra CSV column");
            if (!samples.empty() && v[0] <= samples.back().time)
                throw std::runtime_error("Non-increasing timestamp");
            samples.push_back({v[0], {v[1], v[2], v[3]}, {v[4], v[5], v[6]}});
        }
        auto r = estimate_mav(samples, p, diagnostics);
        std::ofstream out(output);
        if (!out)
            throw std::runtime_error("Cannot open output");
        out << "time,qw,qx,qy,qz\n" << std::setprecision(17);
        for (std::size_t i = 0; i < samples.size(); ++i) {
            out << samples[i].time;
            for (double q : r.quaternions[i])
                out << ',' << q;
            out << '\n';
        }
        out.close();
        if (!out)
            throw std::runtime_error("Failed writing output");
        std::cerr << std::setprecision(17) << "{\"filter\":\""
                  << (p.gravity_baseline ? "mav-srukf-gravity" : "mav-srukf")
                  << "\",\"samples\":" << samples.size() << ",\"elapsed_seconds\":" << r.elapsed_seconds
                  << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "euroc_benchmark: " << e.what() << '\n';
        return 1;
    }
}
