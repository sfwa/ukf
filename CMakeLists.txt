cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# The version is written once, in include/UKF/Version.h.
include(${CMAKE_CURRENT_LIST_DIR}/cmake/ReadVersion.cmake)
read_version(${CMAKE_CURRENT_LIST_DIR}/include/UKF/Version.h UKF ukf_version ukf_version_prerelease)
project(ukf VERSION ${ukf_version} LANGUAGES CXX)

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

option(UKF_INSTALL "Generate install rules and a ukf CMake package" ${ukf_standalone})
set(UKF_EIGEN_VERSION "3.2.8" CACHE STRING
    "Eigen release to download when UKF_USE_SYSTEM_EIGEN is OFF (3.2.8 is patched; 3.3 or newer is used as released)")
set(UKF_EIGEN_URL "" CACHE STRING
    "Eigen source archive (URL or local file) to use instead of the GitLab release of UKF_EIGEN_VERSION")

include(GNUInstallDirs)
add_library(ukf INTERFACE)
add_library(ukf::ukf ALIAS ukf)
target_compile_features(ukf INTERFACE cxx_std_17)
target_include_directories(ukf INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(UKF_USE_SYSTEM_EIGEN)
    # Any installed Eigen 3.3 or newer, including Eigen 5. Point CMake at a
    # particular installation with Eigen3_DIR or CMAKE_PREFIX_PATH. The version
    # is checked here rather than in find_package because Eigen 5's package
    # accepts a single requested version only within the same minor release.
    find_package(Eigen3 REQUIRED NO_MODULE)
    if(Eigen3_VERSION VERSION_LESS 3.3)
        message(FATAL_ERROR "UKF needs Eigen 3.3 or newer when UKF_USE_SYSTEM_EIGEN is ON; found ${Eigen3_VERSION}")
    endif()
    message(STATUS "UKF: using installed Eigen ${Eigen3_VERSION}")
    target_link_libraries(ukf INTERFACE Eigen3::Eigen)
else()
    if(NOT UKF_EIGEN_VERSION STREQUAL "3.2.8" AND
            (NOT UKF_EIGEN_VERSION MATCHES "^[0-9]+\\.[0-9]+(\\.[0-9]+)?$" OR UKF_EIGEN_VERSION VERSION_LESS 3.3))
        message(FATAL_ERROR "UKF_EIGEN_VERSION must be 3.2.8 or a release of 3.3 or newer, such as 3.4.0 or 5.0.0")
    endif()
    message(STATUS "UKF: using downloaded Eigen ${UKF_EIGEN_VERSION}")
    if(UKF_EIGEN_URL STREQUAL "")
        set(ukf_eigen_url
            https://gitlab.com/libeigen/eigen/-/archive/${UKF_EIGEN_VERSION}/eigen-${UKF_EIGEN_VERSION}.tar.gz)
    else()
        set(ukf_eigen_url ${UKF_EIGEN_URL})
    endif()
    if(UKF_EIGEN_VERSION STREQUAL "3.2.8")
        # Preserve the allocation and C++17/libc++ compatibility fixes.
        set(ukf_eigen_patch
            PATCH_COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-no-malloc-in-llt.patch
            COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-cxx17.patch)
    else()
        # Newer releases need no patches: square-root rank updates use
        # UKF-owned workspace, so EIGEN_NO_MALLOC checks still apply.
        set(ukf_eigen_patch "")
    endif()
    include(ExternalProject)
    set_directory_properties(PROPERTIES EP_PREFIX .)
    ExternalProject_Add(eigen3
        URL ${ukf_eigen_url}
        TIMEOUT 30
        ${ukf_eigen_patch}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON)
    ExternalProject_Get_Property(eigen3 source_dir)
    # An installed package carries the downloaded headers with it. CMake
    # does not prefix a relative SYSTEM install path, so do it here unless
    # the include directory is already absolute.
    if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
        set(ukf_installed_eigen_dir "${CMAKE_INSTALL_INCLUDEDIR}/ukf/eigen3")
    else()
        set(ukf_installed_eigen_dir "$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/ukf/eigen3")
    endif()
    target_include_directories(ukf SYSTEM INTERFACE
        $<BUILD_INTERFACE:${source_dir}>
        $<INSTALL_INTERFACE:${ukf_installed_eigen_dir}>)
    # INTERFACE dependencies propagate to consumers, including parallel builds.
    add_dependencies(ukf eigen3)
endif()

# Eigen unrolls fixed-size expressions up to EIGEN_UNROLLING_LIMIT (100 by
# default). A lower limit trades speed for flash: for the float attitude
# estimator on a Cortex-M7, 0 saves about 13 KB but takes 13 to 15 % more
# instructions an update. Leave it empty for Eigen's default. It is defined
# for every target that links ukf, which must agree with any other code in the
# program that includes Eigen.
set(UKF_EIGEN_UNROLLING_LIMIT "" CACHE STRING
    "EIGEN_UNROLLING_LIMIT for targets that use UKF (empty: Eigen's default of 100; 0 disables unrolling)")
if(NOT UKF_EIGEN_UNROLLING_LIMIT STREQUAL "")
    if(NOT UKF_EIGEN_UNROLLING_LIMIT MATCHES "^[0-9]+$")
        message(FATAL_ERROR "UKF_EIGEN_UNROLLING_LIMIT must be empty or a nonnegative integer")
    endif()
    target_compile_definitions(ukf INTERFACE EIGEN_UNROLLING_LIMIT=${UKF_EIGEN_UNROLLING_LIMIT})
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

# An installed tree is used with find_package(ukf), giving the same ukf::ukf,
# ukf::single and ukf::double targets as add_subdirectory.
if(UKF_INSTALL)
    include(CMakePackageConfigHelpers)
    set_target_properties(ukf PROPERTIES EXPORT_NAME ukf)
    set_target_properties(ukf_single PROPERTIES EXPORT_NAME single)
    set_target_properties(ukf_double PROPERTIES EXPORT_NAME double)
    install(TARGETS ukf ukf_single ukf_double EXPORT ukfTargets)
    install(DIRECTORY include/UKF DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    if(NOT UKF_USE_SYSTEM_EIGEN)
        # The downloaded (and, for 3.2.8, patched) Eigen headers. Build the
        # project before installing so that they have been downloaded.
        install(DIRECTORY ${source_dir}/Eigen ${source_dir}/unsupported
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ukf/eigen3)
    endif()
    set(ukf_cmake_dir ${CMAKE_INSTALL_LIBDIR}/cmake/ukf)
    install(EXPORT ukfTargets NAMESPACE ukf:: DESTINATION ${ukf_cmake_dir})
    configure_package_config_file(cmake/ukfConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/ukfConfig.cmake
        INSTALL_DESTINATION ${ukf_cmake_dir})
    # Within a major version, later releases are source compatible
    # (RELEASING.md). The package is header-only, so any architecture will do.
    set(ukf_version_file_options COMPATIBILITY SameMajorVersion)
    if(NOT CMAKE_VERSION VERSION_LESS 3.14)
        list(APPEND ukf_version_file_options ARCH_INDEPENDENT)
    endif()
    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ukfConfigVersion.cmake
        VERSION ${PROJECT_VERSION} ${ukf_version_file_options})
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/ukfConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/ukfConfigVersion.cmake
        DESTINATION ${ukf_cmake_dir})
endif()

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
                CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
            # GCC diagnoses unreachable packet/QR/LLT paths and possibly
            # uninitialized values inside Eigen even with SYSTEM includes.
            # Limit this to dependency headers.
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
