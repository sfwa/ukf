CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(unittest)

# Enable ExternalProject CMake module
INCLUDE(ExternalProject)

set(CMAKE_CXX_FLAGS "-isystem /usr/local/include/eigen3 -O3")

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# Add gtest
ExternalProject_Add(
googletest
SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
TIMEOUT 10
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
INCLUDE_DIRECTORIES(${source_dir}/include ../include /usr/local/include/eigen3)

# Add test executable target
ADD_EXECUTABLE(unittest
    src/test_state.cpp
    src/test_dynamics.cpp
    src/test_integrator.cpp
    src/test_sensors.cpp
    src/test_ukf.cpp)

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
