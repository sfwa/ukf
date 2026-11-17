#ifndef ATTITUDE_BENCHMARK_TIMING_H
#define ATTITUDE_BENCHMARK_TIMING_H

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace AttitudeBenchmark {

inline constexpr uint64_t budget_us = 2000;

struct Timing {
    std::size_t count = 0;
    std::size_t over_budget = 0;
    uint64_t total_us = 0;
    uint64_t min_us = 0;
    uint64_t max_us = 0;

    void add(uint64_t elapsed_us) {
        min_us = count ? std::min(min_us, elapsed_us) : elapsed_us;
        max_us = std::max(max_us, elapsed_us);
        total_us += elapsed_us;
        ++count;
        over_budget += elapsed_us > budget_us;
    }

    double mean_us() const { return count ? double(total_us) / count : 0.0; }
};

/* Nearest-rank percentile of an already sorted, nonempty timing array. */
inline uint64_t percentile(const uint64_t* sorted, std::size_t count, unsigned percent) {
    return count ? sorted[(count * percent + 99) / 100 - 1] : 0;
}

}
#endif
