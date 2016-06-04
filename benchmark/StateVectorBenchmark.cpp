#include <benchmark/benchmark.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"

using SV16_FourQuaternions = UKF::StateVector<
    UKF::Field<1, UKF::Quaternion>,
    UKF::Field<2, UKF::Quaternion>,
    UKF::Field<3, UKF::Quaternion>,
    UKF::Field<4, UKF::Quaternion>
>;

using SV16_FourVectors = UKF::StateVector<
    UKF::Field<1, UKF::Vector<4>>,
    UKF::Field<2, UKF::Vector<4>>,
    UKF::Field<3, UKF::Vector<4>>,
    UKF::Field<4, UKF::Vector<4>>
>;

using SV16_OneVector = UKF::StateVector<
    UKF::Field<1, UKF::Vector<16>>
>;

/*
Tests to ensure setField() and getField() don't hurt performance compared to
segment().
*/
void StateVector_SetGetUsingSegment(benchmark::State& state) {
    SV16_FourVectors test_state;
    test_state.segment<4>(0) << UKF::Vector<4>(1, 2, 3, 4);
    while(state.KeepRunning()) {
        test_state.segment<4>(0) << UKF::Vector<4>(1, 2, 3, 4);
        benchmark::DoNotOptimize(test_state.segment<4>(0));
    }
}

void StateVector_SetGetUsingSetField(benchmark::State& state) {
    SV16_FourVectors test_state;
    test_state.segment<4>(0) << UKF::Vector<4>(1, 2, 3, 4);
    while(state.KeepRunning()) {
        test_state.set_field<1>(UKF::Vector<4>(1, 2, 3, 4));
        benchmark::DoNotOptimize(test_state.get_field<1>());
    }
}

BENCHMARK(StateVector_SetGetUsingSegment);
BENCHMARK(StateVector_SetGetUsingSetField);

template <typename T>
void StateVector_SigmaPointGeneration(benchmark::State& state) {
    T test_state;
    test_state << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
    typename T::CovarianceMatrix covariance = T::CovarianceMatrix::Zero();
    covariance.diagonal() << UKF::Vector<T::covariance_size()>::Ones();
    while(state.KeepRunning()) {
        benchmark::DoNotOptimize(test_state.calculate_sigma_point_distribution(covariance));
    }
}

BENCHMARK_TEMPLATE(StateVector_SigmaPointGeneration, SV16_OneVector);
BENCHMARK_TEMPLATE(StateVector_SigmaPointGeneration, SV16_FourVectors);
BENCHMARK_TEMPLATE(StateVector_SigmaPointGeneration, SV16_FourQuaternions);

template <typename T>
void StateVector_SigmaPointMean(benchmark::State& state) {
    T test_state;
    test_state << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
    typename T::CovarianceMatrix covariance = T::CovarianceMatrix::Zero();
    covariance.diagonal() << UKF::Vector<T::covariance_size()>::Ones();
    typename T::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);
    while(state.KeepRunning()) {
        benchmark::DoNotOptimize(T::calculate_sigma_point_mean(sigma_points));
    }
}

BENCHMARK_TEMPLATE(StateVector_SigmaPointMean, SV16_OneVector);
BENCHMARK_TEMPLATE(StateVector_SigmaPointMean, SV16_FourVectors);
BENCHMARK_TEMPLATE(StateVector_SigmaPointMean, SV16_FourQuaternions);

template <typename T>
void StateVector_SigmaPointDeltas(benchmark::State& state) {
    T test_state;
    test_state << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
    typename T::CovarianceMatrix covariance = T::CovarianceMatrix::Zero();
    covariance.diagonal() << UKF::Vector<T::covariance_size()>::Ones();
    typename T::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);
    T test_mean = T::calculate_sigma_point_mean(sigma_points);
    while(state.KeepRunning()) {
        benchmark::DoNotOptimize(test_mean.calculate_sigma_point_deltas(sigma_points));
    }
}

BENCHMARK_TEMPLATE(StateVector_SigmaPointDeltas, SV16_OneVector);
BENCHMARK_TEMPLATE(StateVector_SigmaPointDeltas, SV16_FourVectors);
BENCHMARK_TEMPLATE(StateVector_SigmaPointDeltas, SV16_FourQuaternions);

template <typename T>
void StateVector_SigmaPointCovariance(benchmark::State& state) {
    T test_state;
    test_state << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
    typename T::CovarianceMatrix covariance = T::CovarianceMatrix::Zero();
    covariance.diagonal() << UKF::Vector<T::covariance_size()>::Ones();
    typename T::SigmaPointDistribution sigma_points = test_state.calculate_sigma_point_distribution(covariance);
    T test_mean = T::calculate_sigma_point_mean(sigma_points);
    typename T::SigmaPointDeltas sigma_point_deltas = test_mean.calculate_sigma_point_deltas(sigma_points);
    while(state.KeepRunning()) {
        benchmark::DoNotOptimize(T::calculate_sigma_point_covariance(sigma_point_deltas));
    }
}

BENCHMARK_TEMPLATE(StateVector_SigmaPointCovariance, SV16_OneVector);
BENCHMARK_TEMPLATE(StateVector_SigmaPointCovariance, SV16_FourVectors);
BENCHMARK_TEMPLATE(StateVector_SigmaPointCovariance, SV16_FourQuaternions);

template <typename T>
void StateVector_FullAPrioriCalculation(benchmark::State& state) {
    T test_state;
    test_state << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
    typename T::CovarianceMatrix covariance = T::CovarianceMatrix::Zero();
    covariance.diagonal() << UKF::Vector<T::covariance_size()>::Ones();
    typename T::SigmaPointDistribution sigma_points;
    typename T::SigmaPointDeltas sigma_point_deltas;
    T test_mean;
    typename T::CovarianceMatrix apriori_covariance;
    while(state.KeepRunning()) {
        sigma_points = test_state.calculate_sigma_point_distribution(covariance);
        test_mean = T::calculate_sigma_point_mean(sigma_points);
        sigma_point_deltas = test_mean.calculate_sigma_point_deltas(sigma_points);
        apriori_covariance = T::calculate_sigma_point_covariance(sigma_point_deltas);
    }
}

BENCHMARK_TEMPLATE(StateVector_FullAPrioriCalculation, SV16_OneVector);
BENCHMARK_TEMPLATE(StateVector_FullAPrioriCalculation, SV16_FourVectors);
BENCHMARK_TEMPLATE(StateVector_FullAPrioriCalculation, SV16_FourQuaternions);

template <typename T>
void StateVector_UpdateDelta(benchmark::State& state) {
    T test_state;
    test_state << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
    typename T::StateVectorDelta test_delta = T::StateVectorDelta::Ones();
    while(state.KeepRunning()) {
        test_state.apply_delta(test_delta);
    }
}

BENCHMARK_TEMPLATE(StateVector_UpdateDelta, SV16_OneVector);
BENCHMARK_TEMPLATE(StateVector_UpdateDelta, SV16_FourVectors);
BENCHMARK_TEMPLATE(StateVector_UpdateDelta, SV16_FourQuaternions);
