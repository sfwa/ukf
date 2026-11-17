#include <cstdint>
#include <cstdio>
#include "sdkconfig.h"
#include "esp_chip_info.h"
#include "esp_cpu.h"
#include "esp_heap_caps.h"
#include "esp_system.h"
#include "esp_timer.h"
#include "freertos/FreeRTOS.h"
#include "freertos/idf_additions.h"
#include "freertos/task.h"
#include "benchmark.h"

#if !defined(CONFIG_IDF_TARGET_ESP32S31)
#error "Build this application for ESP32-S31"
#endif
#if CONFIG_ESP_DEFAULT_CPU_FREQ_MHZ != 320
#error "Use 320 MHz for the ESP32-S31 benchmark"
#endif
#if !defined(CONFIG_COMPILER_OPTIMIZATION_PERF)
#error "Use CONFIG_COMPILER_OPTIMIZATION_PERF (-O2) for this benchmark"
#endif
#if defined(CONFIG_PM_ENABLE) || defined(CONFIG_SPIRAM) || defined(CONFIG_FREERTOS_UNICORE)
#error "Use the benchmark defaults: dual core, fixed CPU frequency, internal SRAM"
#endif

namespace {

constexpr unsigned kStackBytes = 64 * 1024;
constexpr int kBenchmarkCore = 1;

uint64_t now_us() {
    return static_cast<uint64_t>(esp_timer_get_time());
}

void yield_cpu() {
    // Called outside the measured update. Allow the idle task/watchdog to run.
    vTaskDelay(1);
}

void benchmark_task(void*) {
    esp_chip_info_t chip;
    esp_chip_info(&chip);
    std::printf("{\"type\":\"platform\",\"platform\":\"esp32s31\","
        "\"idf\":\"%s\",\"compiler\":\"%s\",\"chip_revision\":%u,"
        "\"chip_cores\":%u,\"task_core\":%d,\"cpu_mhz\":%d,"
        "\"optimization\":\"O2\",\"stack_bytes\":%u,"
        "\"internal_heap_free_bytes\":%zu,\"timer\":\"esp_timer_get_time\","
        "\"timer_resolution_us\":1,\"wifi_started\":false,\"bluetooth_started\":false}\n",
        esp_get_idf_version(), __VERSION__, static_cast<unsigned>(chip.revision),
        static_cast<unsigned>(chip.cores), esp_cpu_get_core_id(),
        CONFIG_ESP_DEFAULT_CPU_FREQ_MHZ, kStackBytes,
        heap_caps_get_free_size(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT));

    const bool ok = AttitudeBenchmark::run_suite("esp32s31", now_us, yield_cpu);
    std::printf("{\"type\":\"platform_complete\",\"suite_ok\":%s,"
        "\"stack_min_free_bytes\":%u,\"internal_heap_min_free_bytes\":%zu}\n",
        ok ? "true" : "false",
        static_cast<unsigned>(uxTaskGetStackHighWaterMark(nullptr)),
        heap_caps_get_minimum_free_size(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT));
    std::fflush(stdout);
    vTaskDelete(nullptr);
}

}

extern "C" void app_main() {
    const BaseType_t created = xTaskCreatePinnedToCore(benchmark_task, "ukf_benchmark",
        kStackBytes, nullptr, 5, nullptr, kBenchmarkCore);
    if (created != pdPASS) {
        std::printf("{\"type\":\"error\",\"message\":\"Cannot allocate benchmark task\"}\n");
        std::fflush(stdout);
    }
}
