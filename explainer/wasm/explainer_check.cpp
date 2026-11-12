/*
Native check of the explainer's interface (explainer_api.h).

The explainer must show what the benchmark's estimator does, and nothing
else. This program reads the same sensor CSV as euroc_benchmark, runs the
recording through the explainer's interface, and writes the same attitude
CSV. The two output files must be identical. It also takes a number of
filter steps apart with ex_explain() and reports how closely the rebuilt
steps agree with the real ones.

Usage:
    explainer_check --input sensors.csv --output attitude.csv
        [--filter mav-srukf|mav-srukf-gravity] [--explain-every N]

A summary is printed as JSON on standard output. The exit code is non-zero
if anything fails.
*/

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "imu_csv.h"
#include "explainer_api.h"

static int fail(const std::string& message) {
    std::cerr << "explainer_check: " << message << '\n';
    return 1;
}

/* Find a named value in one of the interface's self-describing tables. */
static int find_name(int count, const char* (*name)(int), const std::string& wanted) {
    for (int i = 0; i < count; i++) {
        if (wanted == name(i)) {
            return i;
        }
    }
    return -1;
}

int main(int argc, char** argv) {
    std::map<std::string, std::string> options;
    for (int i = 1; i + 1 < argc; i += 2) {
        options[argv[i]] = argv[i + 1];
    }
    std::string input = options["--input"], output = options["--output"];
    std::string filter = options.count("--filter") ? options["--filter"] : "mav-srukf";
    int explain_every = options.count("--explain-every") ? std::atoi(options["--explain-every"].c_str()) : 997;
    if (input.empty() || output.empty() || input == output || explain_every < 1 ||
            (filter != "mav-srukf" && filter != "mav-srukf-gravity")) {
        return fail("Required: --input PATH --output PATH [--filter mav-srukf|mav-srukf-gravity]");
    }
    if (ex_abi_version() != EX_ABI_VERSION) {
        return fail("Interface version mismatch");
    }

    std::string error;
    std::vector<ImuSample> samples;
    if (!read_imu_csv(input, samples, error)) {
        return fail(error);
    }

    /* Hand the recording to the interface, as the page does. */
    const int count = static_cast<int>(samples.size());
    double* values = ex_samples(count);
    if (values == nullptr) {
        return fail("No room for the recording");
    }
    for (int i = 0; i < count; i++) {
        double* row = values + i * EX_SAMPLE_STRIDE;
        row[0] = samples[i].time;
        for (int axis = 0; axis < 3; axis++) {
            row[1 + axis] = samples[i].gyroscope[axis];
            row[4 + axis] = samples[i].accelerometer[axis];
        }
    }

    ex_tuning_reset();
    const int rotor_model = find_name(ex_tuning_count(), ex_tuning_name, "use_rotor_model");
    if (rotor_model < 0) {
        return fail("The tuning table has no use_rotor_model");
    }
    ex_tuning()[rotor_model] = filter == "mav-srukf" ? 1.0 : 0.0;

    if (ex_run() != count) {
        return fail("The estimator did not finish the recording");
    }

    /* Write the attitudes in the benchmark's format. */
    const int attitude = find_name(ex_trace_field_count(), ex_trace_field_name, "attitude");
    if (attitude < 0) {
        return fail("The trace has no attitude");
    }
    const double* trace = ex_trace();
    const int stride = ex_trace_stride(), offset = ex_trace_field_offset(attitude);
    std::vector<UKF::Quaternion> attitudes;
    attitudes.reserve(samples.size());
    for (int i = 0; i < count; i++) {
        const double* q = trace + static_cast<std::size_t>(i) * stride + offset;
        attitudes.push_back(UKF::Quaternion(q[0], q[1], q[2], q[3]));
    }
    if (!write_attitude_csv(output, samples, attitudes, error)) {
        return fail(error);
    }

    /* Take steps apart throughout the recording. */
    const int check = find_name(ex_step_field_count(), ex_step_field_name, "check");
    if (check < 0) {
        return fail("The step has no check");
    }
    double worst = 0.0;
    int steps = 0;
    for (int i = 1; i < count; i += explain_every) {
        if (ex_explain(i) != 0) {
            return fail("Cannot take step " + std::to_string(i) + " apart");
        }
        const double difference = ex_step()[ex_step_field_offset(check)];
        if (!std::isfinite(difference)) {
            return fail("Step " + std::to_string(i) + " has no check value");
        }
        worst = std::max(worst, difference);
        steps++;
    }

    std::cout << std::setprecision(17)
        << "{\"filter\":\"" << filter << "\",\"samples\":" << count
        << ",\"explained_steps\":" << steps << ",\"worst_step_difference\":" << worst << "}\n";

    return 0;
}
