add_executable(attitude_benchmark estimator.cpp)
target_compile_definitions(attitude_benchmark PRIVATE UKF_DOUBLE_PRECISION)
target_include_directories(attitude_benchmark PRIVATE ../../include ${eigen_dir})
add_dependencies(attitude_benchmark eigen3)
