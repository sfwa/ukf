cmake_minimum_required(VERSION 3.15)
project(ukf)

include(CheckCXXCompilerFlag)

# Options for consumers / packagers
option(UKF_USE_SYSTEM_EIGEN "Use system Eigen3 (via find_package) instead of downloading it" ON)
option(UKF_DEFAULT_DOUBLE "Default to double precision for UKF (when not overridden by consumer)" ON)

# Notes for maintainers:
# - We prefer system Eigen by default to make packaging and editor tooling easier.
# - When UKF_USE_SYSTEM_EIGEN=OFF we fall back to downloading Eigen 3.2.8 via
#   ExternalProject (this preserves original upstream behavior while offering
#   a reproducible fallback for environments without system Eigen).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Werror)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register -Wno-ignored-attributes")

# If requested prefer system Eigen (easier for packaging and editor tooling)
if(UKF_USE_SYSTEM_EIGEN)
    find_package(Eigen3 3.2 REQUIRED CONFIG)
    set(eigen_dir ${EIGEN3_INCLUDE_DIR})
else()
    include(ExternalProject)
    # Set default ExternalProject root directory
    set_directory_properties(PROPERTIES EP_PREFIX .)

    # Add eigen3.2.8 as a fallback ExternalProject
    ExternalProject_Add(
        eigen3
        URL https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.gz
        TIMEOUT 30
        # Apply the LLT malloc fix patch, necessary to avoid dynamic memory allocation on embedded systems.
        PATCH_COMMAND patch -p1 -t -N < ${PROJECT_SOURCE_DIR}/patch/eigen-3.2.8-no-malloc-in-llt.patch
        # Disable install step
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON)

    ExternalProject_Get_Property(eigen3 source_dir)
    set(eigen_dir ${source_dir})
endif()

enable_testing()

add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(benchmark EXCLUDE_FROM_ALL)
add_subdirectory(examples EXCLUDE_FROM_ALL)

# Create the INTERFACE target and provide include dirs and default compile defs
add_library(ukf INTERFACE)
target_include_directories(ukf INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${eigen_dir}
)

if(UKF_DEFAULT_DOUBLE)
    target_compile_definitions(ukf INTERFACE UKF_DOUBLE_PRECISION)
else()
    target_compile_definitions(ukf INTERFACE UKF_SINGLE_PRECISION)
endif()

# If system Eigen was used and the imported target exists, propagate it
if(TARGET Eigen3::Eigen)
    target_link_libraries(ukf INTERFACE Eigen3::Eigen)
endif()
