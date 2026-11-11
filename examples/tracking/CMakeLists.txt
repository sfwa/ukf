CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)

ADD_EXECUTABLE(ukf_tracking_example tracking.cpp)
TARGET_INCLUDE_DIRECTORIES(ukf_tracking_example PRIVATE ../../include ${eigen_dir})
TARGET_COMPILE_DEFINITIONS(ukf_tracking_example PRIVATE UKF_DOUBLE_PRECISION EIGEN_NO_MALLOC)
ADD_DEPENDENCIES(ukf_tracking_example eigen3)
