/* Functional validation only: host timings are not ESP32-S31 estimates. */
#include <chrono>
#include <cstring>
#include "benchmark.h"

static uint64_t now_us() {
    return uint64_t(std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
}

static void yield_cpu() {}

static bool check_timing() {
    AttitudeBenchmark::Timing timing;
    const uint64_t samples[] = {1999, 2000, 2001};
    for (uint64_t sample : samples) {
        timing.add(sample);
    }
    return timing.count == 3 && timing.over_budget == 1 && timing.mean_us() == 2000.0 &&
        timing.min_us == 1999 && timing.max_us == 2001 &&
        AttitudeBenchmark::percentile(samples, 3, 50) == 2000 &&
        AttitudeBenchmark::percentile(samples, 3, 95) == 2001 &&
        AttitudeBenchmark::percentile(samples, 3, 99) == 2001 &&
        AttitudeBenchmark::percentile(samples, 0, 99) == 0;
}

int main(int argc, char** argv) {
    if (argc == 2 && std::strcmp(argv[1], "--check-timing") == 0) {
        return check_timing() ? 0 : 1;
    }
    if (argc != 1) {
        std::fprintf(stderr, "Usage: esp32s31_benchmark_host [--check-timing]\n");
        return 2;
    }
    std::fprintf(stderr, "Host functional check only; these are not ESP32-S31 performance measurements.\n");
    return AttitudeBenchmark::run_suite("host", now_us, yield_cpu) ? 0 : 1;
}
