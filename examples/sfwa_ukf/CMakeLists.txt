CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4 FATAL_ERROR)
PROJECT(sfwa_ukf)

# Disable dynamic memory allocation in Eigen
ADD_DEFINITIONS(-DEIGEN_NO_MALLOC)

INCLUDE_DIRECTORIES(../../ukf ${eigen_dir} ${PROJECT_BINARY_DIR})

# Add benchmark executable target
#ADD_EXECUTABLE(sfwa_ukf
#    sfwa_ukf.cpp)
