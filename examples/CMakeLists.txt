CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
PROJECT(examples)

# List examples explicitly so generated directories such as CMakeFiles are
# never mistaken for source projects (including after an in-source build).
ADD_SUBDIRECTORY(ahrs)
ADD_SUBDIRECTORY(sfwa_ukf)
ADD_SUBDIRECTORY(tracking)
ADD_SUBDIRECTORY(smoothing)
ADD_SUBDIRECTORY(attitude)
ADD_SUBDIRECTORY(attitude_benchmark)
ADD_SUBDIRECTORY(euroc_benchmark)
ADD_SUBDIRECTORY(esp32s31_benchmark)
ADD_SUBDIRECTORY(stm32h735_benchmark)
ADD_SUBDIRECTORY(embedded_code_size)

ADD_CUSTOM_TARGET(ukf_examples DEPENDS
    ahrs sfwa_ukf ukf_tracking_example ukf_attitude_example ukf_smoothing_example)

ADD_TEST(NAME tracking_example COMMAND ukf_tracking_example)
ADD_TEST(NAME attitude_example COMMAND ukf_attitude_example)
ADD_TEST(NAME smoothing_example COMMAND ukf_smoothing_example)
SET_TESTS_PROPERTIES(tracking_example attitude_example smoothing_example PROPERTIES
    LABELS examples TIMEOUT 30)

# The check target, when tests are configured, runs all registered CTest tests.
IF(TARGET check)
    ADD_DEPENDENCIES(check ukf_tracking_example ukf_attitude_example ukf_smoothing_example)
ENDIF()

IF(UKF_BUILD_TESTS)
    ADD_EXECUTABLE(ahrs_api_test ahrs/test_api.cpp)
    TARGET_LINK_LIBRARIES(ahrs_api_test PRIVATE ahrs ukf::single)
    ukf_project_warnings(ahrs_api_test)
    ADD_EXECUTABLE(sfwa_api_test sfwa_ukf/test_api.cpp)
    TARGET_LINK_LIBRARIES(sfwa_api_test PRIVATE sfwa_ukf ukf::double)
    ukf_project_warnings(sfwa_api_test)
    ADD_TEST(NAME ahrs_api COMMAND ahrs_api_test)
    ADD_TEST(NAME sfwa_api COMMAND sfwa_api_test)
    SET_TESTS_PROPERTIES(ahrs_api sfwa_api PROPERTIES LABELS examples TIMEOUT 30)
    ADD_DEPENDENCIES(ukf_examples ahrs_api_test sfwa_api_test)
    IF(TARGET check)
        ADD_DEPENDENCIES(check ahrs_api_test sfwa_api_test)
    ENDIF()
ENDIF()
