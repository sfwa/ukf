#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <stdexcept>
#include <type_traits>
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"

using GatingState = UKF::StateVector<UKF::Field<230, UKF::Vector<2>>>;
using GatingMeasurement = UKF::FixedMeasurementVector<UKF::Field<231, UKF::Vector<2>>>;
using GatingDynamicMeasurement = UKF::DynamicMeasurementVector<UKF::Field<231, UKF::Vector<2>>>;

namespace UKF {
template <> template <>
GatingState GatingState::derivative<>() const {
    return GatingState::Zero();
}
template <> template <>
GatingState GatingState::derivative<UKF::Vector<2>>(const UKF::Vector<2>& velocity) const {
    return velocity;
}
template <> template <>
UKF::Vector<2> GatingMeasurement::expected_measurement<GatingState, 231>(const GatingState& state) {
    return state;
}
template <> template <>
UKF::Vector<2> GatingMeasurement::expected_measurement<GatingState, 231, UKF::Vector<2>>(
        const GatingState& state, const UKF::Vector<2>& input) {
    return state + real_t(0.2) * input;
}
template <> template <>
UKF::Vector<2> GatingDynamicMeasurement::expected_measurement<GatingState, 231>(const GatingState& state) {
    return state;
}
template <> template <>
UKF::Vector<2> GatingDynamicMeasurement::expected_measurement<GatingState, 231, UKF::Vector<2>>(
        const GatingState& state, const UKF::Vector<2>& input) {
    return state + real_t(0.2) * input;
}
}

namespace {
template <typename M> using Standard = UKF::Core<GatingState, M, UKF::IntegratorEuler>;
template <typename M> using SquareRoot = UKF::SquareRootCore<GatingState, M, UKF::IntegratorEuler>;
template <typename M> using Parameter = UKF::SquareRootParameterEstimationCore<GatingState, M>;

template <typename F, typename M>
F make_filter() {
    F filter;
    filter.state << 2, -1;
    GatingState::CovarianceMatrix root;
    root << 1, 0, real_t(0.25), real_t(1.5);
    if constexpr (std::is_same_v<F, Standard<M>>) {
        filter.covariance = root * root.transpose();
        filter.process_noise_covariance.setZero();
        filter.measurement_covariance << real_t(1.96), real_t(1.21);
    } else {
        filter.root_covariance = root;
        filter.process_noise_root_covariance.setZero();
        filter.measurement_root_covariance << real_t(1.4), real_t(1.1);
    }
    return filter;
}

template <typename F, typename M>
GatingState::CovarianceMatrix covariance(const F& filter) {
    if constexpr (std::is_same_v<F, Standard<M>>) {
        return filter.covariance;
    } else {
        return filter.root_covariance * filter.root_covariance.transpose();
    }
}

template <typename F, typename M, typename... U>
void predict(F& filter, const U&... input) {
    if constexpr (std::is_same_v<F, Parameter<M>>) {
        filter.a_priori_step();
    } else {
        filter.a_priori_step(real_t(0.5), input...);
    }
}

template <typename F, typename M, typename... U>
void step(F& filter, const M& measurement, const U&... input) {
    if constexpr (std::is_same_v<F, Parameter<M>>) {
        filter.step(measurement, input...);
    } else {
        filter.step(real_t(0.5), measurement, input...);
    }
}

template <typename F, typename M, typename Gate, typename... U>
bool gated_step(F& filter, const M& measurement, Gate&& gate, const U&... input) {
    if constexpr (std::is_same_v<F, Parameter<M>>) {
        return filter.step_with_gate(measurement, std::forward<Gate>(gate), input...);
    } else {
        return filter.step_with_gate(real_t(0.5), measurement, std::forward<Gate>(gate), input...);
    }
}

template <typename F, typename M>
real_t normalized_innovation_squared(const F& filter) {
    UKF::Vector<2> whitened;
    if constexpr (std::is_same_v<F, Standard<M>>) {
        whitened = filter.innovation_covariance.llt().solve(filter.innovation);
        return filter.innovation.dot(whitened);
    } else {
        whitened = filter.innovation_root_covariance.template triangularView<Eigen::Lower>().solve(filter.innovation);
        return whitened.squaredNorm();
    }
}

struct NoncopyableGate {
    int calls = 0;
    NoncopyableGate() = default;
    NoncopyableGate(const NoncopyableGate&) = delete;
    template <typename F, typename M, typename... U>
    bool operator()(const F&, const M&, const U&...) {
        ++calls;
        return true;
    }
};

struct RvalueGate {
    template <typename F>
    bool operator()(const F&, const bool& accept) && { return accept; }
};

template <typename F, typename M>
void check_full_step() {
    F filter = make_filter<F, M>();
    F predicted = filter;
    UKF::Vector<2> input(real_t(0.5), real_t(-0.25));
    predict<F, M>(predicted, input);
    M measurement;
    measurement.template set_field<231>(UKF::Vector<2>(50, -50));
    int calls = 0;
    auto reject_outlier = [&](const auto& current, const M& sample, const UKF::Vector<2>& context) {
        static_assert(std::is_const_v<std::remove_reference_t<decltype(current)>>, "Gate must receive const filter");
        ++calls;
        EXPECT_EQ(&filter, &current);
        EXPECT_EQ(&measurement, &sample);
        EXPECT_EQ(&input, &context);
        EXPECT_EQ(real_t(0), (current.state - predicted.state).norm());
        EXPECT_EQ(real_t(0), (covariance<F, M>(current) - covariance<F, M>(predicted)).norm());
        UKF::Vector<2> residual = measurement - predicted.state - real_t(0.2) * input;
        EXPECT_LT((current.innovation - residual).norm(), real_t(1e-5));
        return normalized_innovation_squared<F, M>(current) < real_t(9);
    };
    EXPECT_FALSE(gated_step(filter, measurement, reject_outlier, input));
    EXPECT_EQ(1, calls);
    EXPECT_EQ(real_t(0), (filter.state - predicted.state).norm());
    EXPECT_EQ(real_t(0), (covariance<F, M>(filter) - covariance<F, M>(predicted)).norm());
    EXPECT_GT((normalized_innovation_squared<F, M>(filter)), real_t(9));

    F reference = filter;
    measurement.template set_field<231>(UKF::Vector<2>(3, -2));
    step(reference, measurement, input);
    NoncopyableGate accept;
    EXPECT_TRUE(gated_step(filter, measurement, accept, input));
    EXPECT_EQ(1, accept.calls);
    EXPECT_EQ(real_t(0), (filter.state - reference.state).norm());
    EXPECT_EQ(real_t(0), (covariance<F, M>(filter) - covariance<F, M>(reference)).norm());
}

template <typename F, typename M>
void check_split_step_and_exception() {
    F filter = make_filter<F, M>();
    predict<F, M>(filter);
    M measurement;
    measurement.template set_field<231>(UKF::Vector<2>(4, -3));
    filter.innovation_step(measurement);
    const F predicted = filter;
    real_t threshold = 0;
    auto gate = [&](const F& current, const real_t& limit) {
        EXPECT_EQ(&threshold, &limit);
        return normalized_innovation_squared<F, M>(current) < limit;
    };
    EXPECT_FALSE(filter.a_posteriori_step(gate, threshold));
    EXPECT_EQ(real_t(0), (filter.state - predicted.state).norm());
    EXPECT_EQ(real_t(0), (covariance<F, M>(filter) - covariance<F, M>(predicted)).norm());
    EXPECT_THROW(filter.a_posteriori_step([](const F&) -> bool {
        throw std::runtime_error("gate failed");
    }), std::runtime_error);
    EXPECT_EQ(real_t(0), (filter.state - predicted.state).norm());
    EXPECT_EQ(real_t(0), (covariance<F, M>(filter) - covariance<F, M>(predicted)).norm());

    F reference = filter;
    reference.a_posteriori_step();
    EXPECT_TRUE(filter.a_posteriori_step(RvalueGate{}, true));
    EXPECT_EQ(real_t(0), (filter.state - reference.state).norm());
    EXPECT_EQ(real_t(0), (covariance<F, M>(filter) - covariance<F, M>(reference)).norm());
}

template <typename F>
void check_empty() {
    using M = GatingDynamicMeasurement;
    F filter = make_filter<F, M>();
    M measurement;
    measurement.template set_field<231>(UKF::Vector<2>(3, -2));
    NoncopyableGate accept;
    ASSERT_TRUE(gated_step(filter, measurement, accept));
    ASSERT_EQ(1, accept.calls);
    F predicted = filter;
    predict<F, M>(predicted);
    M empty;
    EXPECT_FALSE(gated_step(filter, empty, accept));
    EXPECT_EQ(1, accept.calls);
    EXPECT_EQ(0, filter.innovation.size());
    if constexpr (std::is_same_v<F, Standard<M>>) {
        EXPECT_EQ(0, filter.innovation_covariance.size());
    } else {
        EXPECT_EQ(0, filter.innovation_root_covariance.size());
    }
    EXPECT_EQ(real_t(0), (filter.state - predicted.state).norm());
    EXPECT_EQ(real_t(0), (covariance<F, M>(filter) - covariance<F, M>(predicted)).norm());
    EXPECT_FALSE(filter.a_posteriori_step(RvalueGate{}, true));
    F reference = filter;
    step(reference, measurement);
    EXPECT_TRUE(gated_step(filter, measurement, accept));
    EXPECT_EQ(2, accept.calls);
    EXPECT_EQ(real_t(0), (filter.state - reference.state).norm());
    EXPECT_EQ(real_t(0), (covariance<F, M>(filter) - covariance<F, M>(reference)).norm());
}
}

TEST(InnovationGatingTest, FullStepAcrossFilterVariants) {
    check_full_step<Standard<GatingMeasurement>, GatingMeasurement>();
    check_full_step<SquareRoot<GatingMeasurement>, GatingMeasurement>();
    check_full_step<Parameter<GatingMeasurement>, GatingMeasurement>();
}

TEST(InnovationGatingTest, SplitStepAcrossFilterVariants) {
    check_split_step_and_exception<Standard<GatingMeasurement>, GatingMeasurement>();
    check_split_step_and_exception<SquareRoot<GatingMeasurement>, GatingMeasurement>();
    check_split_step_and_exception<Parameter<GatingMeasurement>, GatingMeasurement>();
}

TEST(InnovationGatingTest, EmptyDynamicMeasurementsAcrossFilterVariants) {
    check_empty<Standard<GatingDynamicMeasurement>>();
    check_empty<SquareRoot<GatingDynamicMeasurement>>();
    check_empty<Parameter<GatingDynamicMeasurement>>();
}
