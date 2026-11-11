add_executable(euroc_benchmark estimator.cpp mav_estimator.cpp)
target_link_libraries(euroc_benchmark PRIVATE ukf::double)
ukf_project_warnings(euroc_benchmark)

add_executable(mav_geometry_test test_mav_geometry.cpp)
target_link_libraries(mav_geometry_test PRIVATE ukf::double)
# Must precede every Eigen include, including the GCC warning header that
# ukf_project_warnings force-includes, to enable Eigen's allocation check.
target_compile_definitions(mav_geometry_test PRIVATE EIGEN_RUNTIME_NO_MALLOC)
ukf_project_warnings(mav_geometry_test)
add_dependencies(euroc_benchmark mav_geometry_test)
