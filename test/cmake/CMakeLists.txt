cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(ukf_consumer LANGUAGES CXX)

# Exercise consumption from a parent that has not already enabled C++17.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()
add_subdirectory(../.. ukf)

if(TARGET unittest OR TARGET benchmark OR TARGET ahrs OR TARGET sfwa_ukf)
    message(FATAL_ERROR "UKF project tools must be opt-in when used as a subdirectory")
endif()
if(UKF_USE_SYSTEM_EIGEN AND TARGET eigen3)
    message(FATAL_ERROR "System Eigen must not configure the bundled download")
endif()

add_executable(consumer_double main.cpp)
target_link_libraries(consumer_double PRIVATE ukf::double)
target_compile_definitions(consumer_double PRIVATE EXPECT_DOUBLE)
add_executable(consumer_single main.cpp)
target_link_libraries(consumer_single PRIVATE ukf::single)
add_executable(consumer_manual main.cpp)
target_link_libraries(consumer_manual PRIVATE ukf::ukf)
target_compile_definitions(consumer_manual PRIVATE UKF_SINGLE_PRECISION)

# A sibling that does not use UKF must retain the parent's language and flags.
add_executable(consumer_sibling sibling.cpp)
foreach(target consumer_double consumer_single consumer_manual consumer_sibling)
    add_test(NAME ${target} COMMAND ${target})
endforeach()
