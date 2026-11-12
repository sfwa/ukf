add_executable(euroc_benchmark estimator.cpp)
target_link_libraries(euroc_benchmark PRIVATE ukf::double)
# The CSV and sensor-check helpers are shared with the attitude benchmark.
target_include_directories(euroc_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../attitude_benchmark)
target_compile_definitions(euroc_benchmark PRIVATE EIGEN_NO_MALLOC)
ukf_project_warnings(euroc_benchmark)
