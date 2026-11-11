cmake_minimum_required(VERSION 3.10)

add_executable(ukf_attitude_example attitude.cpp)
target_compile_definitions(ukf_attitude_example PRIVATE UKF_DOUBLE_PRECISION EIGEN_NO_MALLOC)
target_include_directories(ukf_attitude_example PRIVATE ../../include ${eigen_dir})
add_dependencies(ukf_attitude_example eigen3)
