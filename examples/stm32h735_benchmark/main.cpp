#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <sys/stat.h>
#include "stm32h7xx_hal.h"
#include "../esp32s31_benchmark/benchmark.h"

/* DTCM is not cached, so these symbols can be read directly over SWD. */
extern "C" {
volatile char benchmark_log[16384] __attribute__((section(".capture")));
volatile uint32_t benchmark_log_length __attribute__((section(".capture")));
volatile uint32_t benchmark_log_overflow __attribute__((section(".capture")));
/* 0 = running, 1 = complete/valid, 2 = replay failure, 3 = setup/fault. */
volatile uint32_t benchmark_complete __attribute__((section(".capture")));
volatile uint32_t benchmark_fault_cfsr __attribute__((section(".capture")));
volatile uint32_t benchmark_fault_hfsr __attribute__((section(".capture")));
extern char __heap_start[], __heap_end[];
extern uint32_t __stack_limit[];
}

namespace {
constexpr uint32_t cpu_hz = 550000000;
constexpr unsigned stack_bytes = 64 * 1024;
UART_HandleTypeDef uart{};
bool uart_ready = false;
uint32_t previous_cycles = 0;
uint64_t total_cycles = 0;

[[noreturn]] void stop() {
    __DSB();
    while (true) { __NOP(); }
}

[[noreturn]] void fail(const char* reason) {
    std::printf("{\"type\":\"benchmark_error\",\"reason\":\"%s\"}\n", reason);
    std::fflush(stdout);
    benchmark_complete = 3;
    stop();
}

void clock_init() {
    /* STM32H735G-DK's solder-bridge default is direct SMPS, also selected
       before RAM initialization by the vendor ExitRun0Mode startup routine. */
    if (HAL_PWREx_ConfigSupply(PWR_DIRECT_SMPS_SUPPLY) != HAL_OK) {
        fail("direct_smps_configuration");
    }
    __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE0);
    const uint32_t start = HAL_GetTick();
    while (!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {
        if (HAL_GetTick() - start > 1000) { fail("voltage_scale_timeout"); }
    }
    RCC_OscInitTypeDef oscillator{};
    oscillator.OscillatorType = RCC_OSCILLATORTYPE_HSE;
    oscillator.HSEState = RCC_HSE_ON;
    oscillator.PLL.PLLState = RCC_PLL_ON;
    oscillator.PLL.PLLSource = RCC_PLLSOURCE_HSE;
    /* 25 MHz / 5 * 110 / 1 = 550 MHz, VOS0, as permitted by STM32H735. */
    oscillator.PLL.PLLM = 5;
    oscillator.PLL.PLLN = 110;
    oscillator.PLL.PLLP = 1;
    oscillator.PLL.PLLQ = 4;
    oscillator.PLL.PLLR = 2;
    oscillator.PLL.PLLRGE = RCC_PLL1VCIRANGE_2;
    oscillator.PLL.PLLVCOSEL = RCC_PLL1VCOWIDE;
    oscillator.PLL.PLLFRACN = 0;
    if (HAL_RCC_OscConfig(&oscillator) != HAL_OK) { fail("pll_configuration"); }

    RCC_ClkInitTypeDef clocks{};
    clocks.ClockType = RCC_CLOCKTYPE_SYSCLK | RCC_CLOCKTYPE_HCLK |
        RCC_CLOCKTYPE_D1PCLK1 | RCC_CLOCKTYPE_PCLK1 | RCC_CLOCKTYPE_PCLK2 | RCC_CLOCKTYPE_D3PCLK1;
    clocks.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
    clocks.SYSCLKDivider = RCC_SYSCLK_DIV1;
    clocks.AHBCLKDivider = RCC_HCLK_DIV2;
    clocks.APB3CLKDivider = RCC_APB3_DIV2;
    clocks.APB1CLKDivider = RCC_APB1_DIV2;
    clocks.APB2CLKDivider = RCC_APB2_DIV2;
    clocks.APB4CLKDivider = RCC_APB4_DIV2;
    if (HAL_RCC_ClockConfig(&clocks, FLASH_LATENCY_3) != HAL_OK) { fail("bus_clock_configuration"); }
    SystemCoreClockUpdate();
    if (SystemCoreClock != cpu_hz || HAL_RCC_GetHCLKFreq() != cpu_hz / 2) {
        fail("unexpected_cpu_clock");
    }
}

void uart_init() {
    __HAL_RCC_GPIOD_CLK_ENABLE();
    __HAL_RCC_USART3_CLK_ENABLE();
    GPIO_InitTypeDef gpio{};
    gpio.Pin = GPIO_PIN_8 | GPIO_PIN_9;
    gpio.Mode = GPIO_MODE_AF_PP;
    gpio.Pull = GPIO_PULLUP;
    gpio.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
    gpio.Alternate = GPIO_AF7_USART3;
    HAL_GPIO_Init(GPIOD, &gpio);
    uart.Instance = USART3;
    uart.Init.BaudRate = 115200;
    uart.Init.WordLength = UART_WORDLENGTH_8B;
    uart.Init.StopBits = UART_STOPBITS_1;
    uart.Init.Parity = UART_PARITY_NONE;
    uart.Init.Mode = UART_MODE_TX_RX;
    uart.Init.HwFlowCtl = UART_HWCONTROL_NONE;
    uart.Init.OverSampling = UART_OVERSAMPLING_16;
    uart.Init.OneBitSampling = UART_ONE_BIT_SAMPLE_DISABLE;
    uart.Init.ClockPrescaler = UART_PRESCALER_DIV1;
    if (HAL_UART_Init(&uart) != HAL_OK) { fail("uart_configuration"); }
    uart_ready = true;
}

void timer_init() {
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->LAR = 0xC5ACCE55;
    if (DWT->CTRL & DWT_CTRL_NOCYCCNT_Msk) { fail("cycle_counter_unavailable"); }
    DWT->CYCCNT = 0;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
    __DSB();
    __ISB();
    const uint32_t first = DWT->CYCCNT;
    for (unsigned i = 0; i < 100; ++i) { __NOP(); }
    if (DWT->CYCCNT == first) { fail("cycle_counter_stopped"); }
}

uint64_t now_us() {
    __DSB();
    __ISB();
    const uint32_t cycles = DWT->CYCCNT;
    total_cycles += static_cast<uint32_t>(cycles - previous_cycles);
    previous_cycles = cycles;
    /* Extend wraps before converting to integer microseconds. Calls occur at
       every update and yield, well inside the 7.8-second wrap interval. */
    return total_cycles / (cpu_hz / 1000000);
}

void yield_cpu() {
    (void)now_us();
    /* No RTOS or sensor workload; the 1 kHz SysTick interrupt remains enabled. */
    __NOP();
}

void paint_stack() {
    volatile uint32_t* begin = __stack_limit;
    const auto limit = (__get_MSP() - 1024U) & ~uint32_t(3);
    while (reinterpret_cast<uintptr_t>(begin) < limit) { *begin++ = 0xA5A5A5A5; }
}

uint32_t stack_free_bytes() {
    const volatile uint32_t* current = __stack_limit;
    unsigned count = 0;
    while (count < stack_bytes / sizeof(uint32_t) && current[count] == 0xA5A5A5A5) { ++count; }
    return count * sizeof(uint32_t);
}
}

extern "C" void SysTick_Handler() { HAL_IncTick(); }
extern "C" void HardFault_Handler() {
    benchmark_fault_cfsr = SCB->CFSR;
    benchmark_fault_hfsr = SCB->HFSR;
    benchmark_complete = 3;
    stop();
}
extern "C" void MemManage_Handler() { HardFault_Handler(); }
extern "C" void BusFault_Handler() { HardFault_Handler(); }
extern "C" void UsageFault_Handler() { HardFault_Handler(); }

extern "C" int _write(int, char* data, int size) {
    for (int i = 0; i < size; ++i) {
        const uint32_t at = benchmark_log_length;
        if (at < sizeof(benchmark_log) - 1) {
            benchmark_log[at] = data[i];
            benchmark_log[at + 1] = '\0';
            benchmark_log_length = at + 1;
        } else {
            benchmark_log_overflow = 1;
        }
    }
    if (uart_ready && size > 0) {
        (void)HAL_UART_Transmit(&uart, reinterpret_cast<uint8_t*>(data), static_cast<uint16_t>(size), 1000);
    }
    return size;
}

extern "C" void* _sbrk(ptrdiff_t increment) {
    static char* next = __heap_start;
    if (increment < 0 || static_cast<size_t>(increment) >
            reinterpret_cast<uintptr_t>(__heap_end) - reinterpret_cast<uintptr_t>(next)) {
        errno = ENOMEM;
        return reinterpret_cast<void*>(-1);
    }
    char* result = next;
    next += increment;
    return result;
}
extern "C" int _close(int) { return -1; }
extern "C" int _fstat(int, struct stat* status) { status->st_mode = S_IFCHR; return 0; }
extern "C" int _isatty(int) { return 1; }
extern "C" int _lseek(int, int, int) { return 0; }
extern "C" int _read(int, char*, int) { return 0; }
extern "C" int _getpid() { return 1; }
extern "C" int _kill(int, int) { errno = EINVAL; return -1; }
extern "C" void _init() {}
extern "C" void _fini() {}
extern "C" [[noreturn]] void _exit(int) { benchmark_complete = 3; stop(); }

int main() {
    benchmark_log_length = 0;
    benchmark_log_overflow = 0;
    benchmark_complete = 0;
    benchmark_fault_cfsr = 0;
    benchmark_fault_hfsr = 0;
    benchmark_log[0] = '\0';
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    paint_stack();
    SCB_EnableICache();
    SCB_EnableDCache();
    HAL_Init();
    clock_init();
    uart_init();
    timer_init();
    std::printf("{\"type\":\"platform\",\"platform\":\"stm32h735\","
        "\"board\":\"STM32H735G-DK\",\"compiler\":\"%s\",\"stm32cube\":\"%s\","
        "\"device_id\":%lu,\"chip_revision\":%lu,\"cpu_mhz\":%lu,\"hclk_mhz\":%.1f,"
        "\"optimization\":\"O2\",\"float_abi\":\"hard\",\"fpu\":\"fpv5-d16\","
        "\"icache\":true,\"dcache\":true,\"supply\":\"direct_smps\","
        "\"stack_bytes\":%u,\"stack_memory\":\"DTCM\",\"data_memory\":\"AXI_SRAM\","
        "\"timer\":\"DWT_CYCCNT_extended\",\"timer_resolution_us\":1,"
        "\"systick_hz\":1000,\"rtos\":false,\"sensor_io\":false}\n",
        __VERSION__, STM32CUBE_VERSION, static_cast<unsigned long>(DBGMCU->IDCODE & 0xFFF),
        static_cast<unsigned long>(DBGMCU->IDCODE >> 16), static_cast<unsigned long>(SystemCoreClock / 1000000),
        double(HAL_RCC_GetHCLKFreq()) / 1e6, stack_bytes);
    const bool replay_ok = AttitudeBenchmark::run_suite("stm32h735", now_us, yield_cpu);
    const uint32_t free_stack = stack_free_bytes();
    const bool ok = replay_ok && free_stack != 0 && benchmark_log_overflow == 0;
    std::printf("{\"type\":\"platform_complete\",\"suite_ok\":%s,\"stack_min_free_bytes\":%lu,"
        "\"log_overflow\":%s}\n", ok ? "true" : "false", static_cast<unsigned long>(free_stack),
        benchmark_log_overflow ? "true" : "false");
    std::fflush(stdout);
    benchmark_complete = ok ? 1 : 2;
    stop();
}
