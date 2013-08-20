CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(ukf)
INCLUDE(ExternalProject)

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# Add eigen3.2
ExternalProject_Add(
eigen3
URL https://bitbucket.org/eigen/eigen/get/3.2.0.tar.gz
TIMEOUT 30
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

ExternalProject_Get_Property(eigen3 source_dir)
SET(eigen_dir ${source_dir})

INCLUDE_DIRECTORIES(include ${eigen_dir})

SET(CMAKE_CXX_FLAGS "-isystem ${eigen_dir} -O3 -Weverything -Werror -Wno-documentation -Wno-system-headers -Wno-padded")

ADD_LIBRARY(ukf STATIC
    src/state.cpp
    src/dynamics.cpp
    src/sensors.cpp
    src/ukf.cpp
    src/ukf-estimates.cpp
    src/ukf-gain.cpp)

ADD_SUBDIRECTORY(c EXCLUDE_FROM_ALL)

ADD_DEPENDENCIES(ukf eigen3)

ENABLE_TESTING()

ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)
