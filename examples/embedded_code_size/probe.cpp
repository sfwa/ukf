/*
Link-only size probe for the complete recorded-data attitude estimator.

The three C entry points are retained by the measurement script. Inputs and
state storage come from the caller, so the compiler cannot specialize a
synthetic trajectory or replace the estimator with constant output. Nothing
here supplies a startup routine, sensor driver, timer, logging, or benchmark
harness. The resulting ELF is a size artifact, not runnable firmware.
*/

#include <cstddef>
#include <cstdint>
#include <new>
#include "../attitude_benchmark/attitude_estimator.h"

#if !defined(UKF_SINGLE_PRECISION)
#error "The embedded size probe requires UKF_SINGLE_PRECISION"
#endif
#if !defined(UKF_SIZE_VARIANT) || UKF_SIZE_VARIANT < 0 || UKF_SIZE_VARIANT > 2
#error "Set UKF_SIZE_VARIANT to 0 (empty), 1 (UKF), or 2 (SR-UKF)"
#endif

static_assert(sizeof(real_t) == sizeof(float), "The probe must use float32");

#if UKF_SIZE_VARIANT == 1
using SizeEstimator = Attitude_Estimator<Attitude_UKF, Heading_UKF>;
#elif UKF_SIZE_VARIANT == 2
using SizeEstimator = Attitude_Estimator<Attitude_SRUKF, Heading_SRUKF>;
#endif

/*
Fixed-width metadata lets the host read the required caller-owned state size
and alignment directly from either 32-bit little-endian target ELF. It is
identical in length for all variants and is not estimator executable code.
*/
extern "C" const std::uint32_t ukf_size_state_layout[2]
    __attribute__((section(".probe_metadata"), used)) = {
#if UKF_SIZE_VARIANT == 0
        0, 1
#else
        static_cast<std::uint32_t>(sizeof(SizeEstimator)),
        static_cast<std::uint32_t>(alignof(SizeEstimator))
#endif
    };

/* storage must have the size/alignment reported by ukf_size_state_layout. */
extern "C" bool ukf_size_initialise(void* storage, const ImuSample* sample) {
#if UKF_SIZE_VARIANT == 0
    (void)storage;
    (void)sample;
    return true;
#else
    auto* estimator = new (storage) SizeEstimator;
    return estimator->initialise(*sample, true);
#endif
}

extern "C" bool ukf_size_update(void* storage,
        const ImuSample* previous, const ImuSample* sample) {
#if UKF_SIZE_VARIANT == 0
    (void)storage;
    (void)previous;
    (void)sample;
    return true;
#else
    return static_cast<SizeEstimator*>(storage)->update(*previous, *sample);
#endif
}

/* The output order matches the benchmark JSON: w, x, y, z. */
extern "C" void ukf_size_attitude(const void* storage, float* output_wxyz) {
#if UKF_SIZE_VARIANT == 0
    (void)storage;
    output_wxyz[0] = 1.0f;
    output_wxyz[1] = output_wxyz[2] = output_wxyz[3] = 0.0f;
#else
    const UKF::Quaternion attitude = static_cast<const SizeEstimator*>(storage)->attitude();
    output_wxyz[0] = attitude.w();
    output_wxyz[1] = attitude.x();
    output_wxyz[2] = attitude.y();
    output_wxyz[3] = attitude.z();
#endif
}
