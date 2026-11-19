CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
PROJECT(unittest)
INCLUDE(ExternalProject)

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# The static libraries that building googletest produces. They are named as
# its byproducts so that a generator which checks every input before it
# builds anything, as Ninja does, knows where they come from.
SET(googletest_binary_dir ${CMAKE_CURRENT_BINARY_DIR}/src/googletest-build)
SET(googletest_libraries
    ${googletest_binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
    ${googletest_binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)

# Add gtest
ExternalProject_Add(
googletest
URL https://github.com/google/googletest/archive/release-1.7.0.zip
TIMEOUT 30
CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DCMAKE_POLICY_VERSION_MINIMUM:STRING=3.10
BINARY_DIR ${googletest_binary_dir}
BUILD_BYPRODUCTS ${googletest_libraries}
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
SET(googletest_dir ${source_dir})

# The same tests build in double precision (unittest) and in single
# precision (unittest_single), with tolerances for each precision; single
# precision is what most embedded targets run.
SET(unittest_sources
    TestAngle.cpp
    TestStateVector.cpp
    TestMeasurementStep.cpp
    TestFixedMeasurementVector.cpp
    TestDynamicMeasurementVector.cpp
    TestDenseMeasurementCovariance.cpp
    TestCore.cpp
    TestSmoother.cpp
    TestSquareRootCore.cpp
    TestInnovationGating.cpp
    TestInnovationScaling.cpp
    TestQuaternionMeasurement.cpp)
FIND_PACKAGE(Threads REQUIRED)

FOREACH(precision double single)
    IF(precision STREQUAL "double")
        SET(target unittest)
    ELSE()
        SET(target unittest_${precision})
    ENDIF()
    ADD_EXECUTABLE(${target} ${unittest_sources})
    TARGET_COMPILE_DEFINITIONS(${target} PRIVATE EIGEN_NO_MALLOC)
    TARGET_INCLUDE_DIRECTORIES(${target} SYSTEM PRIVATE ${googletest_dir}/include)
    ukf_project_warnings(${target})

    # Create dependency of test on googletest
    ADD_DEPENDENCIES(${target} googletest)

    # Specify test's link libraries
    TARGET_LINK_LIBRARIES(${target} PRIVATE
        ukf::${precision}
        ${googletest_libraries}
        Threads::Threads)

    ADD_TEST(${target} ${target})
ENDFOREACH()
ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS unittest unittest_single)
