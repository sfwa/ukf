CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(ukf)

INCLUDE_DIRECTORIES(include /usr/local/include/eigen3)

set(CMAKE_CXX_FLAGS "-O3")

ADD_LIBRARY(ukf STATIC
    src/state.cpp
    src/dynamics.cpp
    src/sensors.cpp
    src/ukf.cpp
    src/ukf-estimates.cpp)

ADD_SUBDIRECTORY(c EXCLUDE_FROM_ALL)

ENABLE_TESTING()

ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)
