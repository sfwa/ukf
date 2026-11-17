#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "UKF/Core.h"
#include "UKF/Integrator.h"
#include "UKF/MeasurementVector.h"
#include "UKF/StateVector.h"
#include "UKF/Types.h"
#include "comparisons.h"

enum MyFields {
  LatLon,
  Altitude,
  Velocity,
  Acceleration,
  Attitude,
  AngularVelocity,
  Position
};

using MyStateVector = UKF::StateVector<
    UKF::Field<LatLon, UKF::Vector<2>>, UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<Attitude, UKF::Quaternion>, UKF::Field<Altitude, real_t>>;

TEST(StateVectorTest, Instantiation) {
  MyStateVector test_state;

  EXPECT_EQ(10, MyStateVector::MaxRowsAtCompileTime);
  EXPECT_EQ(10, test_state.size());
  EXPECT_EQ(9, MyStateVector::CovarianceMatrix::MaxRowsAtCompileTime);
  EXPECT_EQ(9, MyStateVector::CovarianceMatrix::MaxColsAtCompileTime);
  EXPECT_EQ(10, MyStateVector::SigmaPointDistribution::MaxRowsAtCompileTime);
  EXPECT_EQ(19, MyStateVector::SigmaPointDistribution::MaxColsAtCompileTime);
}

TEST(StateVectorTest, Assignment) {
  MyStateVector test_state;

  test_state.set_field<LatLon>(UKF::Vector<2>(-37.8136, 144.9631));
  test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
  test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(10);

  EXPECT_VECTOR_EQ(UKF::Vector<2>(-37.8136, 144.9631),
                   test_state.get_field<LatLon>());
  EXPECT_VECTOR_EQ(UKF::Vector<3>(1, 2, 3), test_state.get_field<Velocity>());
  EXPECT_QUATERNION_EQ(UKF::Quaternion(1, 0, 0, 0),
                       test_state.get_field<Attitude>());
  EXPECT_EQ(10, test_state.get_field<Altitude>());
}

TEST(StateVectorTest, Arithmetic) {}

TEST(StateVectorTest, DefaultParameters) {
  EXPECT_FLOAT_EQ(1.0, UKF::Parameters::AlphaSquared<MyStateVector>);
  EXPECT_FLOAT_EQ(0.0, UKF::Parameters::Beta<MyStateVector>);
  EXPECT_FLOAT_EQ(3.0, UKF::Parameters::Kappa<MyStateVector>);
  EXPECT_FLOAT_EQ(3.0, UKF::Parameters::Lambda<MyStateVector>);
  EXPECT_FLOAT_EQ(0.0, UKF::Parameters::MRP_A<MyStateVector>);
  EXPECT_FLOAT_EQ(2.0, UKF::Parameters::MRP_F<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WM0<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WC0<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WMI<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WCI<MyStateVector>);
}

using AlternateStateVector =
    UKF::StateVector<UKF::Field<Attitude, UKF::Quaternion>,
                     UKF::Field<AngularVelocity, UKF::Vector<3>>,
                     UKF::Field<Acceleration, UKF::Vector<3>>>;

namespace UKF {
namespace Parameters {
template <>
constexpr real_t AlphaSquared<AlternateStateVector> = 1e-6;
template <>
constexpr real_t Beta<AlternateStateVector> = 2.0;
template <>
constexpr real_t Kappa<AlternateStateVector> = 3.0;
}  // namespace Parameters
}  // namespace UKF

TEST(StateVectorTest, CustomParameters) {
  EXPECT_FLOAT_EQ(1e-6, UKF::Parameters::AlphaSquared<AlternateStateVector>);
  EXPECT_FLOAT_EQ(2.0, UKF::Parameters::Beta<AlternateStateVector>);
  EXPECT_FLOAT_EQ(3.0, UKF::Parameters::Kappa<AlternateStateVector>);
  EXPECT_FLOAT_EQ(12e-6 - 9, UKF::Parameters::Lambda<AlternateStateVector>);
  EXPECT_FLOAT_EQ(0.0, UKF::Parameters::MRP_A<AlternateStateVector>);
  EXPECT_FLOAT_EQ(2.0, UKF::Parameters::MRP_F<AlternateStateVector>);
  EXPECT_FLOAT_EQ(-749999, UKF::Parameters::Sigma_WM0<AlternateStateVector>);
  EXPECT_FLOAT_EQ(-749996, UKF::Parameters::Sigma_WC0<AlternateStateVector>);
  EXPECT_FLOAT_EQ(41666.667, UKF::Parameters::Sigma_WMI<AlternateStateVector>);
  EXPECT_FLOAT_EQ(41666.667, UKF::Parameters::Sigma_WCI<AlternateStateVector>);

  EXPECT_FLOAT_EQ(1.0, UKF::Parameters::AlphaSquared<MyStateVector>);
  EXPECT_FLOAT_EQ(0.0, UKF::Parameters::Beta<MyStateVector>);
  EXPECT_FLOAT_EQ(3.0, UKF::Parameters::Kappa<MyStateVector>);
  EXPECT_FLOAT_EQ(3.0, UKF::Parameters::Lambda<MyStateVector>);
  EXPECT_FLOAT_EQ(0.0, UKF::Parameters::MRP_A<MyStateVector>);
  EXPECT_FLOAT_EQ(2.0, UKF::Parameters::MRP_F<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WM0<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 4.0, UKF::Parameters::Sigma_WC0<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WMI<MyStateVector>);
  EXPECT_FLOAT_EQ(1.0 / 24.0, UKF::Parameters::Sigma_WCI<MyStateVector>);
}

TEST(StateVectorTest, SigmaPointGeneration) {
  MyStateVector test_state;

  test_state.set_field<LatLon>(UKF::Vector<2>(45, 135));
  test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
  test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(10);

  MyStateVector::CovarianceMatrix covariance =
      MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points, target_sigma_points;

  target_sigma_points << 45, 48.464, 45, 45, 45, 45, 45, 45, 45, 45, 41.536, 45,
      45, 45, 45, 45, 45, 45, 45, 135, 135, 138.46, 135, 135, 135, 135, 135,
      135, 135, 135, 131.54, 135, 135, 135, 135, 135, 135, 135, 1, 1, 1, 4.464,
      1, 1, 1, 1, 1, 1, 1, 1, -2.464, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 5.464, 2, 2,
      2, 2, 2, 2, 2, 2, -1.464, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 6.464, 3, 3, 3, 3,
      3, 3, 3, 3, -0.464, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0.866, 0, 0, 0, 0, 0, 0,
      0, 0, -0.866, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.866, 0, 0, 0, 0, 0, 0, 0, 0,
      -0.866, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.866, 0, 0, 0, 0, 0, 0, 0, 0,
      -0.866, 0, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 0.5, 0.5,
      0.5, 1, 10, 10, 10, 10, 10, 10, 10, 10, 10, 13.464, 10, 10, 10, 10, 10,
      10, 10, 10, 6.536;
  sigma_points = test_state.calculate_sigma_point_distribution(
      (covariance * (MyStateVector::covariance_size() +
                     UKF::Parameters::Lambda<MyStateVector>))
          .llt()
          .matrixL());

  EXPECT_VECTOR_EQ(target_sigma_points.col(0), sigma_points.col(0));
  EXPECT_VECTOR_EQ(target_sigma_points.col(1), sigma_points.col(1));
  EXPECT_VECTOR_EQ(target_sigma_points.col(2), sigma_points.col(2));
  EXPECT_VECTOR_EQ(target_sigma_points.col(3), sigma_points.col(3));
  EXPECT_VECTOR_EQ(target_sigma_points.col(4), sigma_points.col(4));
  EXPECT_VECTOR_EQ(target_sigma_points.col(5), sigma_points.col(5));
  EXPECT_VECTOR_EQ(target_sigma_points.col(6), sigma_points.col(6));
  EXPECT_VECTOR_EQ(target_sigma_points.col(7), sigma_points.col(7));
  EXPECT_VECTOR_EQ(target_sigma_points.col(8), sigma_points.col(8));
  EXPECT_VECTOR_EQ(target_sigma_points.col(9), sigma_points.col(9));
  EXPECT_VECTOR_EQ(target_sigma_points.col(10), sigma_points.col(10));
  EXPECT_VECTOR_EQ(target_sigma_points.col(11), sigma_points.col(11));
  EXPECT_VECTOR_EQ(target_sigma_points.col(12), sigma_points.col(12));
  EXPECT_VECTOR_EQ(target_sigma_points.col(13), sigma_points.col(13));
  EXPECT_VECTOR_EQ(target_sigma_points.col(14), sigma_points.col(14));
  EXPECT_VECTOR_EQ(target_sigma_points.col(15), sigma_points.col(15));
  EXPECT_VECTOR_EQ(target_sigma_points.col(16), sigma_points.col(16));
  EXPECT_VECTOR_EQ(target_sigma_points.col(17), sigma_points.col(17));
  EXPECT_VECTOR_EQ(target_sigma_points.col(18), sigma_points.col(18));
}

TEST(StateVectorTest, SigmaPointMean) {
  MyStateVector test_state;

  test_state.set_field<LatLon>(UKF::Vector<2>(-37.8136, 144.9631));
  test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
  test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(10);

  MyStateVector::CovarianceMatrix covariance =
      MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points =
      test_state.calculate_sigma_point_distribution(
          (covariance * (MyStateVector::covariance_size() +
                         UKF::Parameters::Lambda<MyStateVector>))
              .llt()
              .matrixL());

  EXPECT_VECTOR_EQ(test_state,
                   MyStateVector::calculate_sigma_point_mean(sigma_points));
}

TEST(StateVectorTest, SigmaPointDeltas) {
  MyStateVector test_state;

  test_state.set_field<LatLon>(UKF::Vector<2>(-37.8136, 144.9631));
  test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
  test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(10);

  MyStateVector::CovarianceMatrix covariance =
      MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points =
      test_state.calculate_sigma_point_distribution(
          (covariance * (MyStateVector::covariance_size() +
                         UKF::Parameters::Lambda<MyStateVector>))
              .llt()
              .matrixL());
  MyStateVector test_mean =
      MyStateVector::calculate_sigma_point_mean(sigma_points);
  MyStateVector::SigmaPointDeltas sigma_point_deltas, target_sigma_point_deltas;

  target_sigma_point_deltas << 0, 3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      3.464, 0, 0, 0, 0, 0, 0, 0, 0, -3.464;
  sigma_point_deltas = test_mean.calculate_sigma_point_deltas(sigma_points);

  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(0), sigma_point_deltas.col(0));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(1), sigma_point_deltas.col(1));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(2), sigma_point_deltas.col(2));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(3), sigma_point_deltas.col(3));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(4), sigma_point_deltas.col(4));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(5), sigma_point_deltas.col(5));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(6), sigma_point_deltas.col(6));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(7), sigma_point_deltas.col(7));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(8), sigma_point_deltas.col(8));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(9), sigma_point_deltas.col(9));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(10),
                   sigma_point_deltas.col(10));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(11),
                   sigma_point_deltas.col(11));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(12),
                   sigma_point_deltas.col(12));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(13),
                   sigma_point_deltas.col(13));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(14),
                   sigma_point_deltas.col(14));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(15),
                   sigma_point_deltas.col(15));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(16),
                   sigma_point_deltas.col(16));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(17),
                   sigma_point_deltas.col(17));
  EXPECT_VECTOR_EQ(target_sigma_point_deltas.col(18),
                   sigma_point_deltas.col(18));
}

TEST(StateVectorTest, SigmaPointCovariance) {
  MyStateVector test_state;

  test_state.set_field<LatLon>(UKF::Vector<2>(-37.8136, 144.9631));
  test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
  test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(10);

  MyStateVector::CovarianceMatrix covariance =
      MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  MyStateVector::SigmaPointDistribution sigma_points =
      test_state.calculate_sigma_point_distribution(
          (covariance * (MyStateVector::covariance_size() +
                         UKF::Parameters::Lambda<MyStateVector>))
              .llt()
              .matrixL());
  MyStateVector test_mean =
      MyStateVector::calculate_sigma_point_mean(sigma_points);
  MyStateVector::SigmaPointDeltas sigma_point_deltas =
      test_mean.calculate_sigma_point_deltas(sigma_points);
  MyStateVector::CovarianceMatrix calculated_covariance =
      MyStateVector::calculate_sigma_point_covariance(sigma_point_deltas);

  EXPECT_VECTOR_EQ(covariance.col(0), calculated_covariance.col(0));
  EXPECT_VECTOR_EQ(covariance.col(1), calculated_covariance.col(1));
  EXPECT_VECTOR_EQ(covariance.col(2), calculated_covariance.col(2));
  EXPECT_VECTOR_EQ(covariance.col(3), calculated_covariance.col(3));
  EXPECT_VECTOR_EQ(covariance.col(4), calculated_covariance.col(4));
  EXPECT_VECTOR_EQ(covariance.col(5), calculated_covariance.col(5));
  EXPECT_VECTOR_EQ(covariance.col(6), calculated_covariance.col(6));
  EXPECT_VECTOR_EQ(covariance.col(7), calculated_covariance.col(7));
  EXPECT_VECTOR_EQ(covariance.col(8), calculated_covariance.col(8));
}

TEST(StateVectorTest, SmallAlphaSigmaPoints) {
  AlternateStateVector test_state;

  test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
  test_state.set_field<AngularVelocity>(UKF::Vector<3>(1, 2, 3));
  test_state.set_field<Acceleration>(UKF::Vector<3>(1, 2, 3));

  AlternateStateVector::CovarianceMatrix covariance =
      MyStateVector::CovarianceMatrix::Zero();
  covariance.diagonal() << 1000.0, 1000.0, 1000.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  /*
  Scale by ScaledDimension, as the cores do: n + Lambda cancels most of
  alpha^2 when alpha is small, by 3 % in single precision.
  */
  AlternateStateVector::SigmaPointDistribution sigma_points =
      test_state.calculate_sigma_point_distribution(
          (covariance * UKF::Parameters::ScaledDimension<AlternateStateVector>)
              .llt()
              .matrixL());
  AlternateStateVector test_mean =
      AlternateStateVector::calculate_sigma_point_mean(sigma_points);
  AlternateStateVector::SigmaPointDeltas sigma_point_deltas =
      test_mean.calculate_sigma_point_deltas(sigma_points);
  AlternateStateVector::CovarianceMatrix calculated_covariance =
      AlternateStateVector::calculate_sigma_point_covariance(
          sigma_point_deltas);

  EXPECT_VECTOR_EQ(covariance.col(0), calculated_covariance.col(0));
  EXPECT_VECTOR_EQ(covariance.col(1), calculated_covariance.col(1));
  EXPECT_VECTOR_EQ(covariance.col(2), calculated_covariance.col(2));
  EXPECT_VECTOR_EQ(covariance.col(3), calculated_covariance.col(3));
  EXPECT_VECTOR_EQ(covariance.col(4), calculated_covariance.col(4));
  EXPECT_VECTOR_EQ(covariance.col(5), calculated_covariance.col(5));
  EXPECT_VECTOR_EQ(covariance.col(6), calculated_covariance.col(6));
  EXPECT_VECTOR_EQ(covariance.col(7), calculated_covariance.col(7));
  EXPECT_VECTOR_EQ(covariance.col(8), calculated_covariance.col(8));

  /*
  Each sigma point is stored to within epsilon of its value, 2 here, and
  the mean multiplies those roundings by the outer weight, 1 / (2 alpha^2
  (n + kappa)) = 4e4: in single precision the mean is good to only a few
  parts in a thousand.
  */
  const real_t mean_tolerance = precision_tolerance(1e-4, 1e-2);
  EXPECT_QUATERNION_EQ(test_state.get_field<Attitude>(),
                       test_mean.get_field<Attitude>());
  EXPECT_LT((test_state.get_field<AngularVelocity>() -
             test_mean.get_field<AngularVelocity>())
                .norm(),
            mean_tolerance * test_state.get_field<AngularVelocity>().norm());
  EXPECT_LT((test_state.get_field<Acceleration>() -
             test_mean.get_field<Acceleration>())
                .norm(),
            mean_tolerance * test_state.get_field<Acceleration>().norm());
}

/*
These are linear, but they don't have to be – just makes it easier to
calculate the predicted output.
*/

using ProcessModelTestStateVector =
    UKF::StateVector<UKF::Field<Position, UKF::Vector<3>>,
                     UKF::Field<Velocity, UKF::Vector<3>>>;

namespace UKF {
template <>
template <>
ProcessModelTestStateVector ProcessModelTestStateVector::derivative<>() const {
  ProcessModelTestStateVector temp;
  /* Position derivative. */
  temp.set_field<Position>(get_field<Velocity>());

  /* Velocity derivative. */
  temp.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));

  return temp;
}

template <>
template <>
ProcessModelTestStateVector
ProcessModelTestStateVector::derivative<UKF::Vector<3>>(
    const UKF::Vector<3>& acceleration) const {
  ProcessModelTestStateVector temp;
  /* Position derivative. */
  temp.set_field<Position>(get_field<Velocity>());

  /* Velocity derivative. */
  temp.set_field<Velocity>(acceleration);

  return temp;
}

}  // namespace UKF

TEST(StateVectorTest, ProcessModel) {
  ProcessModelTestStateVector test_state;
  UKF::Vector<6> expected_state;

  test_state.set_field<Position>(UKF::Vector<3>(0, 0, 0));
  test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));

  expected_state << 0.1, 0.2, 0.3, 1, 2, 3;
  EXPECT_VECTOR_EQ(expected_state,
                   test_state.process_model<UKF::IntegratorRK4>(0.1));
  EXPECT_VECTOR_EQ(expected_state,
                   test_state.process_model<UKF::IntegratorHeun>(0.1));
  EXPECT_VECTOR_EQ(expected_state,
                   test_state.process_model<UKF::IntegratorEuler>(0.1));

  expected_state << 0.1 + (0.5 * 3 * 0.1 * 0.1), 0.2 + (0.5 * 2 * 0.1 * 0.1),
      0.3 + (0.5 * 1 * 0.1 * 0.1), 1.3, 2.2, 3.1;
  EXPECT_VECTOR_EQ(expected_state, test_state.process_model<UKF::IntegratorRK4>(
                                       0.1, UKF::Vector<3>(3, 2, 1)));
  EXPECT_VECTOR_EQ(expected_state,
                   test_state.process_model<UKF::IntegratorHeun>(
                       0.1, UKF::Vector<3>(3, 2, 1)));
  expected_state << 0.1, 0.2, 0.3, 1.3, 2.2, 3.1;
  EXPECT_VECTOR_EQ(expected_state,
                   test_state.process_model<UKF::IntegratorEuler>(
                       0.1, UKF::Vector<3>(3, 2, 1)));
}

TEST(StateVectorTest, UpdateDelta) {
  MyStateVector test_state;

  test_state.set_field<LatLon>(UKF::Vector<2>(45, 135));
  test_state.set_field<Velocity>(UKF::Vector<3>(1, 2, 3));
  test_state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
  test_state.set_field<Altitude>(10);

  MyStateVector::StateVectorDelta test_delta;
  test_delta << 10, -20, 3, 2, 1, 0, 0, 0.5, 1;

  test_state.apply_delta(test_delta);

  MyStateVector expected_state;
  expected_state.set_field<LatLon>(UKF::Vector<2>(55, 115));
  expected_state.set_field<Velocity>(UKF::Vector<3>(4, 4, 4));
  expected_state.set_field<Attitude>(UKF::Quaternion(0.9692, 0, 0, 0.2462));
  expected_state.set_field<Altitude>(11);

  EXPECT_VECTOR_EQ(expected_state, test_state);
}

/* A direct model deliberately has no derivative specialization. */
using DirectModelTestState = UKF::StateVector<UKF::Field<101, UKF::Vector<2>>>;
using DirectModelTestMeasurement =
    UKF::FixedMeasurementVector<UKF::Field<102, UKF::Vector<2>>>;

namespace UKF {
template <>
template <>
DirectModelTestState DirectModelTestState::process_model<>(real_t dt) const {
  return *this + dt * UKF::Vector<2>(1, -2);
}

template <>
template <>
DirectModelTestState DirectModelTestState::process_model<>(
    real_t dt, const real_t& speed) const {
  return *this + dt * speed * UKF::Vector<2>(1, -2);
}

template <>
template <>
UKF::Vector<2>
DirectModelTestMeasurement::expected_measurement<DirectModelTestState, 102>(
    const DirectModelTestState& state) {
  return state;
}

template <>
template <>
UKF::Vector<2> DirectModelTestMeasurement::expected_measurement<
    DirectModelTestState, 102, real_t>(const DirectModelTestState& state,
                                       const real_t&) {
  return state;
}
}  // namespace UKF

TEST(StateVectorTest, DirectProcessModel) {
  DirectModelTestState state;
  state << 2, 3;
  EXPECT_VECTOR_EQ(UKF::Vector<2>(2.5, 2), state.process_model(real_t(0.5)));
  EXPECT_VECTOR_EQ(UKF::Vector<2>(3.5, 0),
                   state.process_model(real_t(0.5), real_t(3)));
}

TEST(StateVectorTest, CoreDispatchesDirectProcessModel) {
  UKF::Core<DirectModelTestState, DirectModelTestMeasurement> filter;
  filter.state << 2, 3;
  filter.covariance.setIdentity();
  filter.process_noise_covariance.setZero();
  filter.measurement_covariance.setOnes();
  DirectModelTestMeasurement z;
  z << 2.5, 2;
  filter.step(real_t(0.5), z);
  EXPECT_VECTOR_EQ(UKF::Vector<2>(2.5, 2), filter.state);
  EXPECT_VECTOR_EQ(UKF::Vector<2>(0.5, 0.5), filter.covariance.diagonal());

  z << 4, -1;
  filter.step(real_t(0.5), z, real_t(3));
  EXPECT_VECTOR_EQ(UKF::Vector<2>(4, -1), filter.state);
  EXPECT_VECTOR_EQ(UKF::Vector<2>::Constant(real_t(1.0 / 3.0)),
                   filter.covariance.diagonal());
}

/*
q and -q are the same rotation. A process model which keeps w >= 0 flips the
sigma points that cross 180 degrees, which must change neither the mean nor
the covariance. Gibbs vectors (the default MRP_A = 0) do not depend on the
sign of a delta; with MRP_A = 1 the shorter rotation must be taken.
*/
using SignTestState = UKF::StateVector<UKF::Field<210, UKF::Quaternion>>;
using SignTestMrpState = UKF::StateVector<UKF::Field<211, UKF::Quaternion>>;
namespace UKF {
namespace Parameters {
template <>
constexpr real_t MRP_A<SignTestMrpState> = 1.0;
}
}  // namespace UKF

template <typename S, int Key>
void check_quaternion_sign_invariance() {
  const real_t epsilon = std::numeric_limits<real_t>::epsilon();
  S state;
  state.template set_field<Key>(UKF::Quaternion(
      Eigen::AngleAxis<real_t>(real_t(3.13), UKF::Vector<3>::UnitZ())));
  const typename S::CovarianceMatrix root =
      S::CovarianceMatrix::Identity() * real_t(0.3);
  const typename S::SigmaPointDistribution sigma_points =
      state.calculate_sigma_point_distribution(
          root * std::sqrt(UKF::Parameters::ScaledDimension<S>));
  const S mean = S::calculate_sigma_point_mean(sigma_points);
  const typename S::CovarianceMatrix covariance =
      S::calculate_sigma_point_covariance(
          mean.calculate_sigma_point_deltas(sigma_points));

  /* Canonicalise each sigma point, then also flip the whole set. */
  typename S::SigmaPointDistribution canonical = sigma_points;
  int flipped = 0;
  for (int i = 0; i < canonical.cols(); ++i) {
    if (canonical(3, i) < 0) {
      canonical.col(i) = -canonical.col(i);
      ++flipped;
    }
  }
  ASSERT_GT(flipped, 0);
  for (const typename S::SigmaPointDistribution& points :
       {canonical, typename S::SigmaPointDistribution(-canonical)}) {
    const S other_mean = S::calculate_sigma_point_mean(points);
    const UKF::Quaternion a = mean.template get_field<Key>(),
                          b = other_mean.template get_field<Key>();
    EXPECT_LT(real_t(2) * (a * b.conjugate()).vec().norm(),
              real_t(8) * epsilon);
    const typename S::CovarianceMatrix other_covariance =
        S::calculate_sigma_point_covariance(
            other_mean.calculate_sigma_point_deltas(points));
    EXPECT_LT((other_covariance - covariance).norm(),
              real_t(16) * epsilon * covariance.norm());
  }

  /* The spread is that of the root: the deltas are the short rotations. */
  EXPECT_LT((covariance - root * root.transpose()).norm(),
            real_t(0.1) * covariance.norm());
}

TEST(StateVectorTest, QuaternionMeanIsSignInvariant) {
  check_quaternion_sign_invariance<SignTestState, 210>();
}

TEST(StateVectorTest, QuaternionMeanIsSignInvariantWithMrp) {
  check_quaternion_sign_invariance<SignTestMrpState, 211>();
}
