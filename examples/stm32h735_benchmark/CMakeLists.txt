cmake_minimum_required(VERSION 3.16)
# Repository host builds discover every example. This target needs its own toolchain.
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    return()
endif()
project(ukf_stm32h735_benchmark LANGUAGES C CXX ASM)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
get_filename_component(UKF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
get_filename_component(SFWA_ROOT "${UKF_ROOT}/.." ABSOLUTE)
set(STM32CUBE_ROOT "${SFWA_ROOT}/toolchains/STM32CubeH7" CACHE PATH "Pinned STM32CubeH7 checkout")
set(UKF_EIGEN_DIR "${UKF_ROOT}/build/src/eigen3" CACHE PATH "Repository-patched Eigen 3.2.8 headers")
set(UKF_BENCHMARK_RECORDING_HEADER "" CACHE FILEPATH "Optional generated IMU recording header")
set(UKF_BENCHMARK_OPTIMIZATION "O2" CACHE STRING "Benchmark compiler optimization: O2 or Os")
set_property(CACHE UKF_BENCHMARK_OPTIMIZATION PROPERTY STRINGS O2 Os)
if(NOT UKF_BENCHMARK_OPTIMIZATION MATCHES "^(O2|Os)$")
    message(FATAL_ERROR "UKF_BENCHMARK_OPTIMIZATION must be O2 or Os")
endif()
set(CMSIS "${STM32CUBE_ROOT}/Drivers/CMSIS")
set(DEVICE "${CMSIS}/Device/ST/STM32H7xx")
set(HAL "${STM32CUBE_ROOT}/Drivers/STM32H7xx_HAL_Driver")
if(NOT EXISTS "${UKF_EIGEN_DIR}/Eigen/Core")
    message(FATAL_ERROR "Build the repository's eigen3 target first, or set UKF_EIGEN_DIR to its patched headers")
endif()
add_executable(ukf_stm32h735_benchmark
    main.cpp
    "${DEVICE}/Source/Templates/gcc/startup_stm32h735xx.s"
    "${DEVICE}/Source/Templates/system_stm32h7xx.c"
    "${HAL}/Src/stm32h7xx_hal.c"
    "${HAL}/Src/stm32h7xx_hal_cortex.c"
    "${HAL}/Src/stm32h7xx_hal_dma.c"
    "${HAL}/Src/stm32h7xx_hal_flash.c"
    "${HAL}/Src/stm32h7xx_hal_gpio.c"
    "${HAL}/Src/stm32h7xx_hal_pwr.c"
    "${HAL}/Src/stm32h7xx_hal_pwr_ex.c"
    "${HAL}/Src/stm32h7xx_hal_rcc.c"
    "${HAL}/Src/stm32h7xx_hal_rcc_ex.c"
    "${HAL}/Src/stm32h7xx_hal_uart.c"
    "${HAL}/Src/stm32h7xx_hal_uart_ex.c")
set_target_properties(ukf_stm32h735_benchmark PROPERTIES SUFFIX ".elf" CXX_STANDARD 17 CXX_EXTENSIONS OFF)
target_include_directories(ukf_stm32h735_benchmark PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${UKF_ROOT}/include")
target_include_directories(ukf_stm32h735_benchmark SYSTEM PRIVATE
    "${UKF_EIGEN_DIR}" "${CMSIS}/Include" "${DEVICE}/Include" "${HAL}/Inc")
target_compile_definitions(ukf_stm32h735_benchmark PRIVATE
    STM32H735xx USE_HAL_DRIVER USE_PWR_DIRECT_SMPS_SUPPLY USER_VECT_TAB_ADDRESS
    HSE_VALUE=25000000U UKF_SINGLE_PRECISION EIGEN_NO_MALLOC
    STM32CUBE_VERSION="1.13.0"
    UKF_BENCHMARK_OPTIMIZATION="${UKF_BENCHMARK_OPTIMIZATION}")
if(UKF_BENCHMARK_RECORDING_HEADER)
    get_filename_component(recording "${UKF_BENCHMARK_RECORDING_HEADER}" ABSOLUTE)
    if(NOT EXISTS "${recording}")
        message(FATAL_ERROR "Recording header does not exist: ${recording}")
    endif()
    target_compile_definitions(ukf_stm32h735_benchmark PRIVATE UKF_BENCHMARK_RECORDING_HEADER="${recording}")
endif()
set(ARCH_FLAGS -mcpu=cortex-m7 -mthumb -mfpu=fpv5-d16 -mfloat-abi=hard)
target_compile_options(ukf_stm32h735_benchmark PRIVATE ${ARCH_FLAGS} -${UKF_BENCHMARK_OPTIMIZATION} -g3
    -ffunction-sections -fdata-sections -Wall -Wextra
    $<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-rtti;-fno-threadsafe-statics;-fno-use-cxa-atexit;-include${UKF_ROOT}/cmake/EigenGccWarnings.h>)
target_link_options(ukf_stm32h735_benchmark PRIVATE ${ARCH_FLAGS} -nostartfiles
    "-T${CMAKE_CURRENT_SOURCE_DIR}/stm32h735.ld"
    "-Wl,-Map=${CMAKE_CURRENT_BINARY_DIR}/ukf_stm32h735_benchmark.map"
    -Wl,--gc-sections -Wl,--print-memory-usage)
# Full newlib retains %zu and 64-bit integer formatting used by the shared JSON harness.
target_link_libraries(ukf_stm32h735_benchmark PRIVATE m c gcc nosys)
set_property(TARGET ukf_stm32h735_benchmark APPEND PROPERTY LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/stm32h735.ld")
add_custom_command(TARGET ukf_stm32h735_benchmark POST_BUILD
    COMMAND "${CMAKE_OBJCOPY}" -O binary $<TARGET_FILE:ukf_stm32h735_benchmark> "${CMAKE_CURRENT_BINARY_DIR}/ukf_stm32h735_benchmark.bin"
    COMMAND "${CMAKE_SIZE}" $<TARGET_FILE:ukf_stm32h735_benchmark>)
