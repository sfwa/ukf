cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(ukf LANGUAGES CXX)

option(UKF_USE_SYSTEM_EIGEN "Use an installed Eigen3 package instead of patched Eigen 3.2.8" OFF)
option(UKF_WARNINGS_AS_ERRORS "Treat warnings in UKF's own targets as errors" ON)
set(ukf_project_source_dir ${CMAKE_CURRENT_SOURCE_DIR})

# A consumer needs only the header targets unless it opts into the project tools.
set(ukf_standalone OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(ukf_standalone ON)
endif()
option(UKF_BUILD_TESTS "Configure the UKF unit tests" ${ukf_standalone})
option(UKF_BUILD_BENCHMARKS "Configure the UKF benchmarks" ${ukf_standalone})
option(UKF_BUILD_EXAMPLES "Configure the UKF examples" ${ukf_standalone})

add_library(ukf INTERFACE)
add_library(ukf::ukf ALIAS ukf)
target_compile_features(ukf INTERFACE cxx_std_17)
target_include_directories(ukf INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(UKF_USE_SYSTEM_EIGEN)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(ukf INTERFACE Eigen3::Eigen)
else()
    include(ExternalProject)
    set_directory_properties(PROPERTIES EP_PREFIX .)
    ExternalProject_Add(eigen3
        URL https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.gz
        TIMEOUT 30
        # Preserve the allocation and C++17/libc++ compatibility fixes.
        PATCH_COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-no-malloc-in-llt.patch
        COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-cxx17.patch
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON)
    ExternalProject_Get_Property(eigen3 source_dir)
    target_include_directories(ukf SYSTEM INTERFACE ${source_dir})
    # INTERFACE dependencies propagate to consumers, including parallel builds.
    add_dependencies(ukf eigen3)
endif()

# Keep ukf precision-neutral for applications that already define a precision.
add_library(ukf_single INTERFACE)
add_library(ukf::single ALIAS ukf_single)
target_link_libraries(ukf_single INTERFACE ukf)
target_compile_definitions(ukf_single INTERFACE UKF_SINGLE_PRECISION)
add_library(ukf_double INTERFACE)
add_library(ukf::double ALIAS ukf_double)
target_link_libraries(ukf_double INTERFACE ukf)
target_compile_definitions(ukf_double INTERFACE UKF_DOUBLE_PRECISION)

# Warning policy is private to project executables/libraries, never consumers.
# Project targets are also built without C++ exceptions: the library reports
# failures through return values, and this keeps throw/try out of the tree.
function(ukf_project_warnings target)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(UKF_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -fno-exceptions)
        if(UKF_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
                CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13 AND
                (NOT UKF_USE_SYSTEM_EIGEN OR Eigen3_VERSION VERSION_EQUAL 3.4.0))
            # GCC diagnoses unreachable packet/QR/LLT paths in these Eigen
            # versions even with SYSTEM includes. Limit this to dependency headers.
            target_compile_options(${target} PRIVATE
                "-include${ukf_project_source_dir}/cmake/EigenGccWarnings.h")
        endif()
    endif()
endfunction()

if(UKF_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test EXCLUDE_FROM_ALL)
endif()
if(UKF_BUILD_BENCHMARKS)
    add_subdirectory(benchmark EXCLUDE_FROM_ALL)
endif()
if(UKF_BUILD_EXAMPLES)
    add_subdirectory(examples EXCLUDE_FROM_ALL)
endif()
