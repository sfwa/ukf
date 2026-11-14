cmake_minimum_required(VERSION 3.16)

# This directory is both an ESP-IDF application and a normal host example.
# Configuring the repository itself must not require an installed ESP-IDF.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(IDF_TARGET esp32s31 CACHE STRING "ESP-IDF target")
    if(NOT IDF_TARGET STREQUAL "esp32s31")
        message(FATAL_ERROR "This benchmark targets ESP32-S31; use idf.py set-target esp32s31")
    endif()
    if(NOT DEFINED ENV{IDF_PATH})
        message(FATAL_ERROR "Source the ESP-IDF export script before configuring this application")
    endif()
    include($ENV{IDF_PATH}/tools/cmake/project.cmake)
    project(ukf_esp32s31_benchmark)
else()
    add_executable(esp32s31_benchmark_host host.cpp)
    target_link_libraries(esp32s31_benchmark_host PRIVATE ukf::single)
    target_compile_definitions(esp32s31_benchmark_host PRIVATE EIGEN_NO_MALLOC)
    ukf_project_warnings(esp32s31_benchmark_host)

    set(UKF_BENCHMARK_RECORDING_HEADER "" CACHE FILEPATH "Optional generated IMU recording header")
    if(UKF_BENCHMARK_RECORDING_HEADER)
        get_filename_component(recording_header "${UKF_BENCHMARK_RECORDING_HEADER}" ABSOLUTE
            BASE_DIR "${CMAKE_SOURCE_DIR}")
        if(NOT EXISTS "${recording_header}")
            message(FATAL_ERROR "Recording header does not exist: ${recording_header}")
        endif()
        target_compile_definitions(esp32s31_benchmark_host PRIVATE
            UKF_BENCHMARK_RECORDING_HEADER="${recording_header}")
    endif()
endif()
