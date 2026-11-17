cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(ukf_installed_consumer LANGUAGES CXX)

# Consume an installed UKF through its CMake package, from a parent that has
# not enabled C++17. Configure with CMAKE_PREFIX_PATH (or ukf_DIR) naming the
# install prefix.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()
# A release meets a request for its major version, since the package is
# source compatible within a major version. A pre-release, such as master's
# 1.0.0-dev, comes before 1.0.0 and is found without a version.
find_package(ukf REQUIRED CONFIG)
if(ukf_VERSION_PRERELEASE STREQUAL "")
    find_package(ukf ${ukf_VERSION_MAJOR}.0 REQUIRED CONFIG)
endif()

# Ask the installed version file directly, as find_package does, whether it
# accepts a request (PACKAGE_VERSION_COMPATIBLE) and whether it is an exact
# match (PACKAGE_VERSION_EXACT, which find_package ... EXACT requires).
function(expect_request request expected_COMPATIBLE expected_EXACT)
    set(PACKAGE_FIND_VERSION ${request})
    string(REGEX MATCH "^[0-9]+" PACKAGE_FIND_VERSION_MAJOR ${request})
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
    set(PACKAGE_VERSION_EXACT FALSE)
    include(${ukf_DIR}/ukfConfigVersion.cmake)
    foreach(result COMPATIBLE EXACT)
        set(got FALSE)
        if(PACKAGE_VERSION_${result})
            set(got TRUE)
        endif()
        if(NOT got STREQUAL expected_${result})
            message(FATAL_ERROR "Asking for ukf ${request}: PACKAGE_VERSION_${result} is ${got}, expected ${expected_${result}}")
        endif()
    endforeach()
endfunction()
math(EXPR next_major "${ukf_VERSION_MAJOR} + 1")
expect_request(${next_major}.0 FALSE FALSE)
if(ukf_VERSION_PRERELEASE STREQUAL "")
    expect_request(${ukf_VERSION} TRUE TRUE)
else()
    # 1.0.0-rc.1 must not pass for 1.0.0, exactly or otherwise.
    expect_request(${ukf_VERSION} FALSE FALSE)
    expect_request(${ukf_VERSION_MAJOR} FALSE FALSE)
endif()
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
