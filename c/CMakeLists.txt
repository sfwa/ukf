CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(cukf)

INCLUDE_DIRECTORIES(../include /usr/local/include/eigen3)

set(CMAKE_CXX_FLAGS "-isystem /usr/local/include/eigen3 -O3 -Weverything -Werror -Wno-documentation -Wno-system-headers -Wno-error=padded -Wno-global-constructors -Wno-exit-time-destructors")

ADD_LIBRARY(cukf SHARED interface.cpp)

ADD_DEPENDENCIES(cukf ukf)

TARGET_LINK_LIBRARIES(cukf ukf)
