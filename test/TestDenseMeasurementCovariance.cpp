#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include "UKF/Types.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"

namespace {

using MixedDenseMeasurement = UKF::FixedMeasurementVector<
    UKF::Field<200, UKF::Vector<2>>,
    UKF::Field<201, UKF::FieldVector>,
    UKF::Field<202, real_t>>;
using DynamicDenseMeasurement = UKF::DynamicMeasurementVector<
    UKF::Field<200, UKF::Vector<2>>,
    UKF::Field<201, UKF::FieldVector>,
    UKF::Field<202, real_t>>;

template <typename Actual, typename Expected>
void expect_matrix_near(const Actual& actual, const Expected& expected) {
    ASSERT_EQ(expected.rows(), actual.rows());
    ASSERT_EQ(expected.cols(), actual.cols());
    EXPECT_LT((actual - expected).norm(), real_t(1e-4) * (real_t(1) + expected.norm()));
}

UKF::Matrix<6, 6> correlated_noise_root() {
    UKF::Matrix<6, 6> root = UKF::Matrix<6, 6>::Identity();
    root.diagonal() << 1, 2, 3, 4, 5, 6;
    root(1, 0) = real_t(0.4);
    root(2, 0) = real_t(-0.2);
    root(3, 1) = real_t(0.3);
    root(4, 2) = real_t(-0.5);
    root(5, 2) = real_t(0.7);
    root(5, 3) = real_t(0.6);
    return root;
}

}

TEST(DenseMeasurementCovarianceTest, TransformsCorrelationsAcrossFieldBoundaries) {
    MixedDenseMeasurement z, predicted;
    z.set_field<200>(UKF::Vector<2>(2, -1));
    z.set_field<201>(UKF::FieldVector(1, 2, 3));
    z.set_field<202>(4);
    predicted = z;
    predicted.set_field<201>(UKF::FieldVector(2, 1, 3));

    const MixedDenseMeasurement::DenseCovarianceMatrix root = correlated_noise_root();
    const MixedDenseMeasurement::DenseCovarianceMatrix covariance = root * root.transpose();
    UKF::Matrix<6, 6> transform = UKF::Matrix<6, 6>::Identity();
    transform.block<3, 3>(2, 2) = UKF::Detail::calculate_rotation_vector_jacobian<MixedDenseMeasurement>(
        z.get_field<201>(), predicted.get_field<201>());
    const UKF::Matrix<6, 6> expected_root = transform * root;
    const UKF::Matrix<6, 6> expected_covariance = expected_root * expected_root.transpose();

    expect_matrix_near(z.calculate_measurement_covariance(covariance, predicted), expected_covariance);
    expect_matrix_near(z.calculate_measurement_root_covariance(root, predicted), expected_root);
    EXPECT_GT((expected_covariance.block<2, 3>(0, 2).norm()), real_t(0.01));
}

TEST(DenseMeasurementCovarianceTest, DiagonalInputMatchesExistingInterface) {
    MixedDenseMeasurement z, predicted;
    z.set_field<200>(UKF::Vector<2>(2, -1));
    z.set_field<201>(UKF::FieldVector(1, 2, 3));
    z.set_field<202>(4);
    predicted = z;
    predicted.set_field<201>(UKF::FieldVector(2, 1, 3));
    MixedDenseMeasurement::CovarianceVector diagonal;
    diagonal << 1, 2, 3, 4, 5, 6;
    MixedDenseMeasurement::DenseCovarianceMatrix dense = MixedDenseMeasurement::DenseCovarianceMatrix::Zero();
    dense.diagonal() = diagonal;

    expect_matrix_near(z.calculate_measurement_covariance(dense, predicted),
        z.calculate_measurement_covariance(diagonal, predicted));
    expect_matrix_near(z.calculate_measurement_root_covariance(dense, predicted),
        z.calculate_measurement_root_covariance(diagonal, predicted));
}

TEST(DenseMeasurementCovarianceTest, DynamicSubsetReordersRowsAndRetainsLatentColumns) {
    // The declared scalar is last, but it is inserted first. The FieldVector
    // is absent, while its latent noise columns still contribute to the scalar.
    DynamicDenseMeasurement z;
    z.set_field<202>(4);
    z.set_field<200>(UKF::Vector<2>(2, -1));
    const DynamicDenseMeasurement::DenseCovarianceMatrix root = correlated_noise_root();
    const DynamicDenseMeasurement::DenseCovarianceMatrix covariance = root * root.transpose();
    UKF::Matrix<3, 6> expected_root;
    expected_root.row(0) = root.row(5);
    expected_root.row(1) = root.row(0);
    expected_root.row(2) = root.row(1);
    const UKF::Matrix<3, 3> expected_covariance = expected_root * expected_root.transpose();

    expect_matrix_near(z.calculate_measurement_root_covariance(root, z), expected_root);
    expect_matrix_near(z.calculate_measurement_covariance(covariance, z), expected_covariance);
    EXPECT_GT(expected_covariance(0, 0), root(5, 5) * root(5, 5));
}

TEST(DenseMeasurementCovarianceTest, EmptyDynamicTransformsHaveNoRows) {
    DynamicDenseMeasurement z;
    const DynamicDenseMeasurement::DenseCovarianceMatrix noise =
        DynamicDenseMeasurement::DenseCovarianceMatrix::Identity();
    const auto covariance = z.calculate_measurement_covariance(noise, z);
    const auto root = z.calculate_measurement_root_covariance(noise, z);
    EXPECT_EQ(0, covariance.rows());
    EXPECT_EQ(0, covariance.cols());
    EXPECT_EQ(0, root.rows());
    EXPECT_EQ(6, root.cols());
}

TEST(DenseMeasurementCovarianceTest, DynamicFieldVectorTransformsCrossFieldCorrelations) {
    DynamicDenseMeasurement z, predicted;
    z.set_field<202>(4);
    z.set_field<201>(UKF::FieldVector(1, 2, 3));
    // Prediction field access is independent of its insertion order as well.
    predicted.set_field<201>(UKF::FieldVector(2, 1, 3));
    predicted.set_field<202>(4);
    const DynamicDenseMeasurement::DenseCovarianceMatrix root = correlated_noise_root();
    const DynamicDenseMeasurement::DenseCovarianceMatrix covariance = root * root.transpose();
    UKF::Matrix<4, 6> transform = UKF::Matrix<4, 6>::Zero();
    transform(0, 5) = real_t(1);
    transform.block<3, 3>(1, 2) = UKF::Detail::calculate_rotation_vector_jacobian<DynamicDenseMeasurement>(
        z.get_field<201>(), predicted.get_field<201>());
    const UKF::Matrix<4, 6> expected_root = transform * root;
    const UKF::Matrix<4, 4> expected_covariance = expected_root * expected_root.transpose();

    expect_matrix_near(z.calculate_measurement_covariance(covariance, predicted), expected_covariance);
    expect_matrix_near(z.calculate_measurement_root_covariance(root, predicted), expected_root);
    EXPECT_GT((expected_covariance.block<1, 3>(0, 1).norm()), real_t(0.01));
}

TEST(DenseMeasurementCovarianceTest, ScalarDenseNoiseIsDistinctFromDiagonalStorage) {
    using Fixed = UKF::FixedMeasurementVector<UKF::Field<220, real_t>>;
    using Dynamic = UKF::DynamicMeasurementVector<UKF::Field<220, real_t>>;
    Fixed fixed;
    fixed.set_field<220>(real_t(3));
    Fixed::DenseCovarianceMatrix dense;
    dense(0, 0) = real_t(4);
    Fixed::CovarianceVector diagonal;
    diagonal.set_field<220>(real_t(4));
    EXPECT_EQ(real_t(4), fixed.calculate_measurement_covariance(dense, fixed)(0, 0));
    expect_matrix_near(fixed.calculate_measurement_covariance(dense, fixed),
        fixed.calculate_measurement_covariance(diagonal, fixed));
    dense(0, 0) = real_t(2);
    EXPECT_EQ(real_t(2), fixed.calculate_measurement_root_covariance(dense, fixed)(0, 0));

    Dynamic dynamic;
    dynamic.set_field<220>(real_t(3));
    EXPECT_EQ(real_t(2), dynamic.calculate_measurement_root_covariance(dense, dynamic)(0, 0));
    dense(0, 0) = real_t(4);
    EXPECT_EQ(real_t(4), dynamic.calculate_measurement_covariance(dense, dynamic)(0, 0));
}

using DenseNoiseTestState = UKF::StateVector<UKF::Field<210, UKF::Vector<2>>>;
using DenseNoiseTestMeasurement = UKF::FixedMeasurementVector<UKF::Field<211, UKF::Vector<2>>>;
using DenseNoiseSubsetMeasurement = UKF::DynamicMeasurementVector<
    UKF::Field<212, real_t>, UKF::Field<211, UKF::Vector<2>>>;

namespace UKF {

template <> template <>
DenseNoiseTestState DenseNoiseTestState::derivative<>() const {
    return DenseNoiseTestState::Zero();
}

template <> template <>
UKF::Vector<2> DenseNoiseTestMeasurement::expected_measurement<DenseNoiseTestState, 211>(
        const DenseNoiseTestState& state) {
    return state;
}

template <> template <>
UKF::Vector<2> DenseNoiseSubsetMeasurement::expected_measurement<DenseNoiseTestState, 211>(
        const DenseNoiseTestState& state) {
    return state;
}

template <> template <>
real_t DenseNoiseSubsetMeasurement::expected_measurement<DenseNoiseTestState, 212>(
        const DenseNoiseTestState& state) {
    return state(0);
}

}

namespace {

template <typename Measurement>
void check_linear_filter(const Measurement& z, const typename Measurement::DenseCovarianceMatrix& noise_root,
        const UKF::Matrix<2, 2>& active_noise) {
    using Dense = typename Measurement::DenseCovarianceMatrix;
    using Standard = UKF::Core<DenseNoiseTestState, Measurement, UKF::IntegratorEuler, Dense>;
    using SquareRoot = UKF::SquareRootCore<DenseNoiseTestState, Measurement, UKF::IntegratorEuler, Dense>;
    using Parameter = UKF::SquareRootParameterEstimationCore<DenseNoiseTestState, Measurement, Dense>;
    const UKF::Vector<2> prior_mean(1, -2);
    const UKF::Matrix<2, 2> prior_covariance = UKF::Matrix<2, 2>::Identity();
    const UKF::Matrix<2, 2> innovation_covariance = prior_covariance + active_noise;
    const UKF::Matrix<2, 2> gain = innovation_covariance.inverse();
    const UKF::Vector<2> expected_mean = prior_mean + gain * (UKF::Vector<2>(z) - prior_mean);
    const UKF::Matrix<2, 2> expected_covariance = prior_covariance - gain;

    Standard standard;
    standard.state = prior_mean;
    standard.covariance = prior_covariance;
    standard.process_noise_covariance.setZero();
    standard.measurement_covariance = noise_root * noise_root.transpose();
    standard.step(real_t(1), z);
    expect_matrix_near(standard.innovation_covariance, innovation_covariance);
    expect_matrix_near(standard.state, expected_mean);
    expect_matrix_near(standard.covariance, expected_covariance);

    SquareRoot square_root;
    square_root.state = prior_mean;
    square_root.root_covariance.setIdentity();
    square_root.process_noise_root_covariance.setZero();
    square_root.measurement_root_covariance = noise_root;
    square_root.step(real_t(1), z);
    expect_matrix_near(square_root.innovation_root_covariance * square_root.innovation_root_covariance.transpose(),
        innovation_covariance);
    expect_matrix_near(square_root.state, expected_mean);
    expect_matrix_near(square_root.root_covariance * square_root.root_covariance.transpose(), expected_covariance);

    Parameter parameter;
    parameter.state = prior_mean;
    parameter.root_covariance.setIdentity();
    parameter.process_noise_root_covariance.setZero();
    parameter.measurement_root_covariance = noise_root;
    parameter.step(z);
    expect_matrix_near(parameter.innovation_root_covariance * parameter.innovation_root_covariance.transpose(),
        innovation_covariance);
    expect_matrix_near(parameter.state, expected_mean);
    expect_matrix_near(parameter.root_covariance * parameter.root_covariance.transpose(), expected_covariance);

    // Correlation must affect the answer, rather than silently becoming diagonal.
    UKF::Matrix<2, 2> diagonal_innovation = UKF::Matrix<2, 2>::Zero();
    diagonal_innovation.diagonal() = innovation_covariance.diagonal();
    const UKF::Vector<2> diagonal_mean = prior_mean + diagonal_innovation.inverse() *
        (UKF::Vector<2>(z) - prior_mean);
    EXPECT_GT((expected_mean - diagonal_mean).norm(), real_t(0.01));
}

}

TEST(DenseMeasurementCovarianceTest, AllFilterVariantsMatchCorrelatedLinearKalmanUpdate) {
    DenseNoiseTestMeasurement z;
    z << 2, 1;
    DenseNoiseTestMeasurement::DenseCovarianceMatrix root;
    root << 1, 0, real_t(0.4), real_t(1.2);
    const UKF::Matrix<2, 2> covariance = root * root.transpose();
    check_linear_filter(z, root, covariance);
}

TEST(DenseMeasurementCovarianceTest, AllFilterVariantsUseEveryColumnOfDynamicNoiseRoot) {
    DenseNoiseSubsetMeasurement z;
    z.set_field<211>(UKF::Vector<2>(2, 1));
    DenseNoiseSubsetMeasurement::DenseCovarianceMatrix root;
    root << 1, 0, 0,
            real_t(0.5), 1, 0,
            real_t(-0.3), real_t(0.4), real_t(1.2);
    const UKF::Matrix<3, 3> covariance = root * root.transpose();
    const UKF::Matrix<2, 2> active_noise = covariance.block<2, 2>(1, 1);
    check_linear_filter(z, root, active_noise);
}
