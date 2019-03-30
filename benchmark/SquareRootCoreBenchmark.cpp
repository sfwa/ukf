#include <benchmark/benchmark.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"

/*
Set up state vector class. The order of these is changed to prevent
linker collisions with the ones in CoreBenchmark.cpp.
*/
enum MyStateFields {
    Attitude,
    AngularVelocity,
    Position,
    Velocity
};

using MyStateVector = UKF::StateVector<
    UKF::Field<Position, UKF::Vector<3>>,
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>,
    UKF::Field<AngularVelocity, UKF::Vector<3>>
>;

namespace UKF {
namespace Parameters {
template <> constexpr real_t AlphaSquared<MyStateVector> = 1e-6;
}


/*
State vector process model. One version takes body frame kinematic
acceleration and angular acceleration as inputs, the other doesn't (assumes
zero accelerations).
*/
template <> template <>
MyStateVector MyStateVector::derivative<UKF::Vector<3>, UKF::Vector<3>>(
        const UKF::Vector<3>& acceleration, const UKF::Vector<3>& angular_acceleration) const {
    MyStateVector temp;
    
    /* Position derivative. */
    temp.set_field<Position>(get_field<Velocity>());

    /* Velocity derivative. */
    temp.set_field<Velocity>(get_field<Attitude>().conjugate() * acceleration);

    /* Attitude derivative. */
    UKF::Quaternion temp_q;
    temp_q.vec() = get_field<AngularVelocity>();
    temp_q.w() = 0;
    temp.set_field<Attitude>(temp_q);

    /* Angular velocity derivative. */
    temp.set_field<AngularVelocity>(angular_acceleration);

    return temp;
}

template <> template <>
MyStateVector MyStateVector::derivative<>() const {
    return derivative(UKF::Vector<3>(0, 0, 0), UKF::Vector<3>(0, 0, 0));
}

}

/*
Set up measurement vector class. The order of these is changed to prevent
linker collisions with the ones in CoreBenchmark.cpp.
*/
enum MyMeasurementFields {
    Accelerometer,
    Magnetometer,
    Gyroscope,
    GPS_Position,
    GPS_Velocity
};

using MyMeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<GPS_Position, UKF::Vector<3>>,
    UKF::Field<GPS_Velocity, UKF::Vector<3>>,
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Magnetometer, UKF::FieldVector>,
    UKF::Field<Gyroscope, UKF::Vector<3>>
>;

using MyCore = UKF::SquareRootCore<
    MyStateVector,
    MyMeasurementVector,
    UKF::IntegratorRK4
>;

namespace UKF {
/*
Define measurement model to be used in tests. NOTE: These are just for
testing, don't expect them to make any physical sense whatsoever.
*/
template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, GPS_Position>(const MyStateVector& state) {
    return state.get_field<Position>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, GPS_Velocity>(const MyStateVector& state) {
    return state.get_field<Velocity>();
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Accelerometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -9.8);
}

template <> template <>
UKF::FieldVector MyMeasurementVector::expected_measurement
<MyStateVector, Magnetometer>(const MyStateVector& state) {
    return state.get_field<Attitude>() * UKF::FieldVector(1, 0, 0);
}

template <> template <>
UKF::Vector<3> MyMeasurementVector::expected_measurement
<MyStateVector, Gyroscope>(const MyStateVector& state) {
    return state.get_field<AngularVelocity>();
}

}

/*
Initialise covariances as square roots to be comparable with CoreBenchmark.cpp.
*/
MyCore create_initialised_sr_test_filter() {
    MyCore test_filter;
    test_filter.state.set_field<Position>(UKF::Vector<3>(0, 0, 0));
    test_filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
    test_filter.state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    test_filter.state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));
    test_filter.root_covariance = MyStateVector::CovarianceMatrix::Zero();
    test_filter.root_covariance.diagonal() <<
        10000, 10000, 10000, 100, 100, 100, 1, 1, 5, 10, 10, 10;
    test_filter.root_covariance = test_filter.root_covariance.llt().matrixU();
    test_filter.measurement_root_covariance << 
        10, 10, 10, 1, 1, 1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0.05, 0.05, 0.05;
    test_filter.measurement_root_covariance = test_filter.measurement_root_covariance.array().sqrt();

    real_t a, b;
    real_t dt = 0.01;
    a = std::sqrt(0.1*dt*dt);
    b = std::sqrt(0.1*dt);
    test_filter.process_noise_root_covariance << a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, b, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, b;
    test_filter.process_noise_root_covariance = test_filter.process_noise_root_covariance.llt().matrixU();

    return test_filter;
}

void SquareRootCore_APrioriStep(benchmark::State& state) {
    MyCore test_filter = create_initialised_sr_test_filter();

    while(state.KeepRunning()) {
        test_filter.a_priori_step(0.01);
    }
}

BENCHMARK(SquareRootCore_APrioriStep);

void SquareRootCore_InnovationStep(benchmark::State& state) {
    MyCore test_filter = create_initialised_sr_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01);

    while(state.KeepRunning()) {
        test_filter.innovation_step(m);
    }
}

BENCHMARK(SquareRootCore_InnovationStep);

void SquareRootCore_APosterioriStep(benchmark::State& state) {
    MyCore test_filter = create_initialised_sr_test_filter();
    MyMeasurementVector m;

    m.set_field<GPS_Position>(UKF::Vector<3>(100, 10, -50));
    m.set_field<GPS_Velocity>(UKF::Vector<3>(20, 0, 0));
    m.set_field<Accelerometer>(UKF::Vector<3>(0, 0, -9.8));
    m.set_field<Magnetometer>(UKF::FieldVector(0, -1, 0));
    m.set_field<Gyroscope>(UKF::Vector<3>(0.5, 0, 0));

    test_filter.a_priori_step(0.01);
    test_filter.innovation_step(m);

    MyStateVector::CovarianceMatrix initial_cov = test_filter.root_covariance;
    MyStateVector initial_state = test_filter.state;

    while(state.KeepRunning()) {
        test_filter.root_covariance = initial_cov;
        test_filter.state = initial_state;
        test_filter.a_posteriori_step();
    }
}

BENCHMARK(SquareRootCore_APosterioriStep);
