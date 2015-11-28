CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(unittest)

# Enable ExternalProject CMake module
INCLUDE(ExternalProject)

set(CMAKE_CXX_FLAGS "-O3 -mssse3 -oaddmul")
set(CMAKE_C_FLAGS "-O3 -mssse3 -oaddmul")

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# Add gtest
ExternalProject_Add(
googletest
SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
TIMEOUT 30
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
SET(googletest_dir ${source_dir})

INCLUDE_DIRECTORIES(${googletest_dir}/include ../include ${eigen_dir})

# Add test executable target
ADD_EXECUTABLE(unittest
    test_state.cpp
    test_dynamics.cpp
    test_integrator.cpp
    test_sensors.cpp
    test_ukf.cpp)

# Create dependency of test on googletest
ADD_DEPENDENCIES(unittest googletest ukf)

# Specify test's link libraries
ExternalProject_Get_Property(googletest binary_dir)
TARGET_LINK_LIBRARIES(unittest
    ukf
    ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
    ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a
    pthread)

ADD_TEST(unittest unittest)
ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS unittest)
