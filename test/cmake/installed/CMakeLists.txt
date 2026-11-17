cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(ukf_installed_consumer LANGUAGES CXX)

# Consume an installed UKF through its CMake package, from a parent that has
# not enabled C++17. Configure with CMAKE_PREFIX_PATH (or ukf_DIR) naming the
# install prefix.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()
# A request for this major version must be met; the package is source
# compatible within a major version.
find_package(ukf 1.0 REQUIRED CONFIG)
set(ukf_version_string ${ukf_VERSION})
if(NOT ukf_VERSION_PRERELEASE STREQUAL "")
    set(ukf_version_string ${ukf_VERSION}-${ukf_VERSION_PRERELEASE})
endif()
message(STATUS "Found ukf ${ukf_version_string} at ${ukf_DIR}")

foreach(target ukf::ukf ukf::single ukf::double)
    if(NOT TARGET ${target})
        message(FATAL_ERROR "The installed ukf package must provide ${target}")
    endif()
endforeach()

add_executable(installed_double ../main.cpp)
target_link_libraries(installed_double PRIVATE ukf::double)
target_compile_definitions(installed_double PRIVATE EXPECT_DOUBLE)
add_executable(installed_single ../main.cpp)
target_link_libraries(installed_single PRIVATE ukf::single)
add_executable(installed_manual ../main.cpp)
target_link_libraries(installed_manual PRIVATE ukf::ukf)
target_compile_definitions(installed_manual PRIVATE UKF_SINGLE_PRECISION)
add_executable(installed_no_exceptions ../no_exceptions.cpp)
target_link_libraries(installed_no_exceptions PRIVATE ukf::double)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(installed_no_exceptions PRIVATE -fno-exceptions)
elseif(MSVC)
    target_compile_definitions(installed_no_exceptions PRIVATE _HAS_EXCEPTIONS=0)
endif()

foreach(target installed_double installed_single installed_manual)
    target_compile_definitions(${target} PRIVATE "EXPECT_VERSION_STRING=\"${ukf_version_string}\"")
endforeach()

foreach(target installed_double installed_single installed_manual installed_no_exceptions)
    add_test(NAME ${target} COMMAND ${target})
endforeach()
