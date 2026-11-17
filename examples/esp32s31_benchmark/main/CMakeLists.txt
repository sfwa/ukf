idf_component_register(SRCS "main.cpp" INCLUDE_DIRS ".."
    REQUIRES esp_timer esp_system esp_hw_support freertos heap)

# Consume the same header library and patched Eigen as the desktop examples.
# Disable its programs here to avoid recursively configuring this application.
get_filename_component(UKF_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
set(UKF_BUILD_TESTS OFF CACHE BOOL "Configure UKF unit tests" FORCE)
set(UKF_BUILD_BENCHMARKS OFF CACHE BOOL "Configure UKF microbenchmarks" FORCE)
set(UKF_BUILD_EXAMPLES OFF CACHE BOOL "Configure UKF desktop examples" FORCE)
add_subdirectory("${UKF_ROOT}" "${CMAKE_BINARY_DIR}/ukf")
target_link_libraries(${COMPONENT_LIB} PRIVATE ukf::single)
target_compile_definitions(${COMPONENT_LIB} PRIVATE EIGEN_NO_MALLOC)

# The warning helper uses this path to include its GCC/Eigen compatibility file.
set(ukf_project_source_dir "${UKF_ROOT}")
ukf_project_warnings(${COMPONENT_LIB})

set(UKF_BENCHMARK_RECORDING_HEADER "" CACHE FILEPATH "Optional generated IMU recording header")
if(UKF_BENCHMARK_RECORDING_HEADER)
    get_filename_component(recording_header "${UKF_BENCHMARK_RECORDING_HEADER}" ABSOLUTE
        BASE_DIR "${PROJECT_DIR}")
    if(NOT EXISTS "${recording_header}")
        message(FATAL_ERROR "Recording header does not exist: ${recording_header}")
    endif()
    target_compile_definitions(${COMPONENT_LIB} PRIVATE
        UKF_BENCHMARK_RECORDING_HEADER="${recording_header}")
endif()
