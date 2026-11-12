// Instantiate every filter core with exceptions disabled by the build.
#include "UKF/Core.h"
#include "UKF/MeasurementVector.h"

#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
#if !defined(_MSC_VER)
#error "This consumer must be compiled with exceptions disabled"
#endif
#endif

using State = UKF::StateVector<UKF::Field<0, UKF::Vector<2>>>;
using Measurement = UKF::DynamicMeasurementVector<UKF::Field<1, real_t>>;

namespace UKF {
template <> template <>
State State::process_model<>(real_t dt) const {
    State next;
    next << (*this)(0) + dt * (*this)(1), (*this)(1);
    return next;
}
template <> template <>
State State::derivative<>() const {
    State derivative;
    derivative << (*this)(1), 0;
    return derivative;
}
template <> template <>
real_t Measurement::expected_measurement<State, 1>(const State& state) {
    return state(0);
}
}

template <typename Filter>
bool run(Filter& filter, Measurement& z) {
    filter.state << 0, 1;
    filter.process_noise_root_covariance = State::CovarianceMatrix::Identity() * real_t(0.01);
    filter.measurement_root_covariance << real_t(0.1);
    return filter.step(real_t(0.1), z);
}

int main() {
    Measurement z;
    z.set_field<1>(real_t(0.1));

    UKF::Core<State, Measurement> standard;
    standard.state << 0, 1;
    standard.covariance = State::CovarianceMatrix::Identity();
    standard.process_noise_covariance = State::CovarianceMatrix::Identity() * real_t(1e-4);
    standard.measurement_covariance << real_t(0.01);

    UKF::SquareRootCore<State, Measurement> direct;
    direct.root_covariance = State::CovarianceMatrix::Identity();
    UKF::SquareRootCore<State, Measurement, UKF::IntegratorRK4> integrated;
    integrated.root_covariance = State::CovarianceMatrix::Identity();
    UKF::SquareRootParameterEstimationCore<State, Measurement> parameter;
    parameter.state << 0, 1;
    parameter.root_covariance = State::CovarianceMatrix::Identity();
    parameter.process_noise_root_covariance = State::CovarianceMatrix::Identity() * real_t(0.01);
    parameter.measurement_root_covariance << real_t(0.1);

    const bool ok = standard.step(real_t(0.1), z) && run(direct, z) && run(integrated, z) &&
        parameter.step(z);
    return ok ? 0 : 1;
}
