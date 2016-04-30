CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4 FATAL_ERROR)
PROJECT(sfwa_ukf)

# Disable dynamic memory allocation in Eigen
ADD_DEFINITIONS(-DEIGEN_NO_MALLOC)

INCLUDE_DIRECTORIES(../../include ${eigen_dir})

# Add executable target
ADD_EXECUTABLE(sfwa_ukf
    sfwa_ukf.cpp)
