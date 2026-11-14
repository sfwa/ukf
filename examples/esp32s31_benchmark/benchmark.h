/* Shared target/host harness. Only update() and attitude() are timed. */
#ifndef ESP32S31_ATTITUDE_BENCHMARK_H
#define ESP32S31_ATTITUDE_BENCHMARK_H

#include <algorithm>
#include <atomic>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <type_traits>
#include "../attitude_benchmark/attitude_estimator.h"
#include "timing.h"

#ifdef UKF_BENCHMARK_RECORDING_HEADER
#include UKF_BENCHMARK_RECORDING_HEADER
#endif

namespace AttitudeBenchmark {

static_assert(std::is_same<real_t, float>::value, "Benchmark requires UKF_SINGLE_PRECISION");
static_assert(Attitude_StateVector::covariance_size() == 9, "Update benchmark metadata for a changed model");
static_assert(Attitude_StateVector::num_sigma() == 19, "Update benchmark metadata for a changed model");

using Clock = uint64_t (*)();
using Yield = void (*)();
using StandardEstimator = Attitude_Estimator<Attitude_UKF, Heading_UKF>;
using SquareRootEstimator = Attitude_Estimator<Attitude_SRUKF, Heading_SRUKF>;

/* Generation and loading are deliberately outside the measured region. */
struct SyntheticSource {
    static constexpr std::size_t size = 6001;
    static constexpr std::size_t start_row = 0;
    const char* name() const { return "synthetic_500hz_rest_motion_magnetic_disturbance"; }
    const char* sha256() const { return ""; }

    ImuSample sample(std::size_t index) const {
        ImuSample result;
        result.time = double(index) / 500.0;
        const real_t t = real_t(result.time);
        const UKF::Vector<3> omega(real_t(0.4), real_t(-0.25), real_t(0.65));
        const UKF::Quaternion initial =
            UKF::Quaternion(Eigen::AngleAxis<real_t>(real_t(0.35), UKF::Vector<3>::UnitX())) *
            UKF::Quaternion(Eigen::AngleAxis<real_t>(real_t(-0.25), UKF::Vector<3>::UnitY()));
        const real_t moving_time = std::clamp(t - real_t(3), real_t(0), real_t(5));
        const UKF::Quaternion attitude = initial * UKF::Quaternion(
            Eigen::AngleAxis<real_t>(omega.norm() * moving_time, omega.normalized()));
        result.gyroscope = UKF::Vector<3>(real_t(0.005), real_t(-0.008), real_t(0.004));
        if (t >= real_t(3) && t < real_t(8)) {
            result.gyroscope += omega;
        }
        result.accelerometer = attitude.conjugate() * UKF::Vector<3>(0, 0, real_t(9.80665));
        UKF::Vector<3> field(real_t(0.6), 0, real_t(-0.8));
        if (t >= real_t(5) && t < real_t(6)) {
            field = UKF::Vector<3>(real_t(1.5), real_t(0.8), real_t(0.2));
        } else if (t >= real_t(6) && t < real_t(6.5)) {
            field.setZero();
        }
        result.magnetometer = attitude.conjugate() * field;
        return result;
    }
};

#ifdef UKF_BENCHMARK_RECORDING_HEADER
struct RecordedSource {
    static constexpr std::size_t size = sizeof(Recording::samples) / sizeof(Recording::samples[0]);
    static constexpr std::size_t start_row = Recording::source_start_row;
    const char* name() const { return Recording::name; }
    const char* sha256() const { return Recording::sha256; }
    ImuSample sample(std::size_t index) const {
        const auto& row = Recording::samples[index];
        ImuSample result;
        result.time = row.time;
        result.gyroscope = Eigen::Map<const UKF::Vector<3>>(row.gyroscope);
        result.accelerometer = Eigen::Map<const UKF::Vector<3>>(row.accelerometer);
        result.magnetometer = Eigen::Map<const UKF::Vector<3>>(row.magnetometer);
        return result;
    }
};
using Source = RecordedSource;
#else
using Source = SyntheticSource;
#endif

static_assert(Source::size >= 2 && Source::size <= 8193, "Use a recording of 2 to 8193 samples");
/* One reusable buffer in internal .bss, never on the estimator task stack. */
inline uint64_t durations[Source::size - 1];

inline void json_string(const char* value) {
    std::putchar('"');
    for (const unsigned char* p = reinterpret_cast<const unsigned char*>(value); *p; ++p) {
        if (*p == '"' || *p == '\\') {
            std::putchar('\\');
            std::putchar(*p);
        } else if (*p < 32) {
            std::printf("\\u%04x", unsigned(*p));
        } else {
            std::putchar(*p);
        }
    }
    std::putchar('"');
}

inline void print_timing(const Timing& timing) {
    std::printf("{\"count\":%zu,\"mean_us\":%.3f,\"min_us\":%" PRIu64
        ",\"max_us\":%" PRIu64 ",\"over_budget\":%zu}",
        timing.count, timing.mean_us(), timing.min_us, timing.max_us, timing.over_budget);
}

inline bool valid_attitude(const UKF::Quaternion& q) {
    return q.coeffs().allFinite() && std::abs(q.squaredNorm() - real_t(1)) < real_t(1e-3);
}

/* Force the output calculation to finish before reading the end timestamp. */
inline void timing_barrier(const UKF::Quaternion& q, bool updated) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(q.coeffs().data()), "g"(updated) : "memory");
#else
    volatile real_t consumed = q.w() + q.x() + q.y() + q.z() + real_t(updated);
    (void)consumed;
    std::atomic_signal_fence(std::memory_order_seq_cst);
#endif
}

struct Warmup {
    bool valid = false;
    std::size_t updates = 0;
    std::size_t failed_sample = 0;
};

template <typename Estimator>
Warmup warm_up(const Source& source, Yield yield_cpu) {
    Warmup result;
    Estimator estimator;
    ImuSample previous = source.sample(0);
    if (!estimator.initialise(previous, true)) {
        return result;
    }
    const std::size_t end = std::min(Source::size, std::size_t(257));
    for (std::size_t i = 1; i < end; ++i) {
        const ImuSample current = source.sample(i);
        ++result.updates;
        if (!estimator.update(previous, current) || !valid_attitude(estimator.attitude())) {
            result.failed_sample = i;
            return result;
        }
        previous = current;
        if (i % 32 == 0) {
            yield_cpu();
        }
    }
    result.valid = true;
    return result;
}

template <typename Estimator>
bool run_one(const Source& source, const char* platform, const char* filter,
        unsigned pass, Clock now_us, Yield yield_cpu) {
    const Warmup warmup = warm_up<Estimator>(source, yield_cpu);
    /* A fresh state preserves initialization, rest entry and recovery in each replay. */
    Estimator estimator;
    ImuSample previous = source.sample(0);
    const uint64_t init_begin = now_us();
    const bool initialised = estimator.initialise(previous, true);
    const uint64_t init_us = now_us() - init_begin;
    bool valid = initialised && warmup.valid;
    Timing all, rest, moving, magnetic_accepted, magnetic_rejected;
    UKF::Quaternion attitude = UKF::Quaternion::Identity();
    std::size_t failed_sample = 0;
    double checksum = 0;

    for (std::size_t i = 1; valid && i < Source::size; ++i) {
        const ImuSample current = source.sample(i);
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const uint64_t begin = now_us();
        const bool updated = estimator.update(previous, current);
        attitude = estimator.attitude();
        timing_barrier(attitude, updated);
        const uint64_t elapsed = now_us() - begin;
        /* Checking/consuming output, accounting and yielding are not timed. */
        durations[all.count] = elapsed;
        all.add(elapsed);
        (estimator.resting() ? rest : moving).add(elapsed);
        (estimator.magnetometer_rejected() ? magnetic_rejected : magnetic_accepted).add(elapsed);
        valid = updated && valid_attitude(attitude);
        if (!valid) {
            failed_sample = i;
        } else {
            checksum += double(attitude.w()) + double(attitude.x()) + double(attitude.y()) + double(attitude.z());
        }
        previous = current;
        if (i % 32 == 0) {
            yield_cpu();
        }
    }

    std::sort(durations, durations + all.count);
    std::printf("{\"type\":\"attitude_benchmark\",\"platform\":");
    json_string(platform);
    std::printf(",\"filter\":\"%s\",\"pass\":%u,\"precision\":\"float32\",\"axes\":9,"
        "\"timed_scope\":\"update_and_attitude\",\"budget_us\":%" PRIu64
        ",\"expected_updates\":%zu,\"warmup_updates\":%zu,\"warmup_valid\":%s,\"warmup_failed_sample\":%zu,"
        "\"initialised\":%s,\"initialise_us\":%" PRIu64
        ",\"replay_valid\":%s,\"failed_sample\":%zu,\"timing\":",
        filter, pass, budget_us, Source::size - 1, warmup.updates,
        warmup.valid ? "true" : "false", warmup.failed_sample, initialised ? "true" : "false",
        init_us, valid ? "true" : "false", failed_sample);
    print_timing(all);
    std::printf(",\"p50_us\":%" PRIu64 ",\"p95_us\":%" PRIu64 ",\"p99_us\":%" PRIu64
        ",\"mean_budget_percent\":%.3f,\"max_budget_percent\":%.3f,\"min_headroom_us\":%.3f,\"observed_compute_within_budget\":%s,\"rest\":",
        percentile(durations, all.count, 50), percentile(durations, all.count, 95),
        percentile(durations, all.count, 99), all.mean_us() * 100.0 / budget_us,
        double(all.max_us) * 100.0 / budget_us,
        double(budget_us) - double(all.max_us),
        valid && all.count == Source::size - 1 && all.over_budget == 0 ? "true" : "false");
    print_timing(rest);
    std::printf(",\"moving\":"); print_timing(moving);
    std::printf(",\"magnetic_accepted\":"); print_timing(magnetic_accepted);
    std::printf(",\"magnetic_rejected\":"); print_timing(magnetic_rejected);
    /* Avoid non-JSON NaN/Inf if a broken estimator produced invalid output. */
    if (all.count && valid_attitude(attitude)) {
        std::printf(",\"final_attitude_wxyz\":[%.9g,%.9g,%.9g,%.9g]",
            double(attitude.w()), double(attitude.x()), double(attitude.y()), double(attitude.z()));
    } else {
        std::printf(",\"final_attitude_wxyz\":null");
    }
    std::printf(",\"checksum\":%.12g}\n", checksum);
    yield_cpu();
    return valid;
}

inline bool run_suite(const char* platform, Clock now_us, Yield yield_cpu) {
    const Source source;
    double dt_min = 0;
    double dt_max = 0;
    double previous_time = source.sample(0).time;
    for (std::size_t i = 0; i < Source::size; ++i) {
        const ImuSample sample = source.sample(i);
        const double dt = sample.time - previous_time;
        if (!std::isfinite(sample.time) || (i && (!std::isfinite(real_t(dt)) || real_t(dt) <= 0)) ||
                !sample.gyroscope.allFinite() || !sample.accelerometer.allFinite() ||
                !sample.magnetometer.allFinite()) {
            std::printf("{\"type\":\"benchmark_error\",\"reason\":\"invalid_input\",\"sample\":%zu}\n", i);
            return false;
        }
        if (i) {
            dt_min = i == 1 ? dt : std::min(dt_min, dt);
            dt_max = i == 1 ? dt : std::max(dt_max, dt);
        }
        previous_time = sample.time;
        if (i % 32 == 0) {
            yield_cpu();
        }
    }
    std::printf("{\"type\":\"benchmark_input\",\"platform\":"); json_string(platform);
    std::printf(",\"source\":"); json_string(source.name());
    std::printf(",\"source_sha256\":"); json_string(source.sha256());
    std::printf(",\"samples\":%zu,\"source_start_row\":%zu,\"dt_min_s\":%.12g,\"dt_max_s\":%.12g,\"dt_mean_s\":%.12g,"
        "\"covariance_dimensions\":9,\"sigma_points\":19,\"integrator\":\"RK4\","
        "\"heading_filter_dimensions\":1,\"passes\":2,\"timing_buffer_bytes\":%zu}\n",
        Source::size, Source::start_row, dt_min, dt_max,
        (source.sample(Source::size - 1).time - source.sample(0).time) / (Source::size - 1), sizeof(durations));

    /* Reverse order on pass two to expose simple order/cache effects. */
    bool valid = run_one<StandardEstimator>(source, platform, "ukf", 1, now_us, yield_cpu);
    valid = run_one<SquareRootEstimator>(source, platform, "srukf", 1, now_us, yield_cpu) && valid;
    valid = run_one<SquareRootEstimator>(source, platform, "srukf", 2, now_us, yield_cpu) && valid;
    valid = run_one<StandardEstimator>(source, platform, "ukf", 2, now_us, yield_cpu) && valid;
    return valid;
}

}
#endif
