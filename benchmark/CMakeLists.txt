CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
PROJECT(benchmark)
INCLUDE(ExternalProject)

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# Add benchmark
ExternalProject_Add(
googlebenchmark
URL https://github.com/google/benchmark/archive/v1.3.0.zip
TIMEOUT 30
CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DCMAKE_POLICY_VERSION_MINIMUM:STRING=3.10
    -DBENCHMARK_ENABLE_TESTING:BOOL=OFF
    # This old dependency enables -Werror in Release independently of UKF.
    # Keep strict warnings on our benchmarks, not the vendored library.
    -DHAVE_CXX_FLAG_WERROR:BOOL=OFF
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

# Specify include dir
ExternalProject_Get_Property(googlebenchmark source_dir)
SET(googlebenchmark_dir ${source_dir})

# Add benchmark executable target
ADD_EXECUTABLE(benchmark
    StateVectorBenchmark.cpp
    MeasurementVectorBenchmark.cpp
    CoreBenchmark.cpp
    SquareRootCoreBenchmark.cpp
    BenchmarkMain.cpp)

TARGET_COMPILE_DEFINITIONS(benchmark PRIVATE EIGEN_NO_MALLOC)
TARGET_INCLUDE_DIRECTORIES(benchmark SYSTEM PRIVATE ${googlebenchmark_dir}/include)
ukf_project_warnings(benchmark)

# Create dependency of benchmark on googlebenchmark
ADD_DEPENDENCIES(benchmark googlebenchmark)

# Specify benchmark's link libraries
ExternalProject_Get_Property(googlebenchmark binary_dir)
FIND_PACKAGE(Threads REQUIRED)
TARGET_LINK_LIBRARIES(benchmark PRIVATE
    ukf::double
    ${binary_dir}/src/${CMAKE_FIND_LIBRARY_PREFIXES}benchmark.a
    Threads::Threads)
