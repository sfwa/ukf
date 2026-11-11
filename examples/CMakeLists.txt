CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
PROJECT(examples)

# List examples explicitly so generated directories such as CMakeFiles are
# never mistaken for source projects (including after an in-source build).
ADD_SUBDIRECTORY(ahrs)
ADD_SUBDIRECTORY(sfwa_ukf)
ADD_SUBDIRECTORY(tracking)
ADD_SUBDIRECTORY(attitude)

ADD_CUSTOM_TARGET(ukf_examples DEPENDS
    ahrs sfwa_ukf ukf_tracking_example ukf_attitude_example)

ADD_TEST(NAME tracking_example COMMAND ukf_tracking_example)
ADD_TEST(NAME attitude_example COMMAND ukf_attitude_example)
SET_TESTS_PROPERTIES(tracking_example attitude_example PROPERTIES
    LABELS examples TIMEOUT 30)

# The existing check target runs all registered CTest tests.
ADD_DEPENDENCIES(check ukf_tracking_example ukf_attitude_example)
