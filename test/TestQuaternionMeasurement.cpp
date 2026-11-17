#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <type_traits>

#include "UKF/Core.h"
#include "UKF/MeasurementVector.h"
#include "UKF/StateVector.h"
#include "UKF/Types.h"

/*
Quaternion-valued measurements: an orientation sensor reports a quaternion,
and its innovation and noise are three-dimensional rotations in the same
convention as quaternion states.
*/
enum QuaternionMeasurementKeys {
  QAttitude = 300,
  QBias,
  QAttitudeSensor,
  QBiasSensor
};

using QState = UKF::StateVector<UKF::Field<QAttitude, UKF::Quaternion>,
                                UKF::Field<QBias, real_t>>;
using QFixed =
    UKF::FixedMeasurementVector<UKF::Field<QAttitudeSensor, UKF::Quaternion>,
                                UKF::Field<QBiasSensor, real_t>>;
using QDynamic =
    UKF::DynamicMeasurementVector<UKF::Field<QBiasSensor, real_t>,
                                  UKF::Field<QAttitudeSensor, UKF::Quaternion>>;
using QAttitudeOnly =
    UKF::FixedMeasurementVector<UKF::Field<QAttitudeSensor, UKF::Quaternion>>;

namespace UKF {
template <>
template <>
QState QState::process_model<>(real_t) const {
  return *this;
}
template <>
template <>
UKF::Quaternion QFixed::expected_measurement<QState, QAttitudeSensor>(
    const QState& state) {
  return state.get_field<QAttitude>();
}
template <>
template <>
real_t QFixed::expected_measurement<QState, QBiasSensor>(const QState& state) {
  return state.get_field<QBias>();
}
template <>
template <>
UKF::Quaternion QDynamic::expected_measurement<QState, QAttitudeSensor>(
    const QState& state) {
  return state.get_field<QAttitude>();
}
template <>
template <>
real_t QDynamic::expected_measurement<QState, QBiasSensor>(
    const QState& state) {
  return state.get_field<QBias>();
}
template <>
template <>
UKF::Quaternion QAttitudeOnly::expected_measurement<QState, QAttitudeSensor>(
    const QState& state) {
  return state.get_field<QAttitude>();
}
}  // namespace UKF

namespace {

UKF::Quaternion rotation(real_t angle, const UKF::Vector<3>& axis) {
  return UKF::Quaternion(Eigen::AngleAxis<real_t>(angle, axis.normalized()));
}

/* The rotation taking b to a, as an angle-axis vector. */
UKF::Vector<3> rotation_between(const UKF::Quaternion& a,
                                const UKF::Quaternion& b) {
  UKF::Quaternion d = a * b.conjugate();
  if (d.w() < real_t(0)) {
    d.coeffs() = -d.coeffs();
  }
  const Eigen::AngleAxis<real_t> angle_axis(d);
  return angle_axis.angle() * angle_axis.axis();
}

constexpr real_t prior_variance = real_t(0.04);
constexpr real_t attitude_variance = real_t(0.01);
constexpr real_t bias_variance = real_t(0.25);
const UKF::Vector<3> attitude_noise =
    UKF::Vector<3>::Constant(attitude_variance);

}  // namespace

TEST(QuaternionMeasurementTest, FieldsAndTypes) {
  static_assert(QFixed::size() == 5,
                "A quaternion is stored as four coefficients");
  static_assert(QFixed::covariance_size() == 4,
                "A quaternion's uncertainty is a rotation");
  static_assert(std::is_same<QFixed::CovarianceVector,
                             UKF::FixedMeasurementVector<
                                 UKF::Field<QAttitudeSensor, UKF::Vector<3>>,
                                 UKF::Field<QBiasSensor, real_t>>>::value,
                "Quaternion noise is a three-vector");
  static_assert(
      std::is_same<QFixed::InnovationVector, QFixed::CovarianceVector>::value,
      "The innovation has the uncertainty's layout");
  using Plain = UKF::FixedMeasurementVector<UKF::Field<1, UKF::Vector<2>>,
                                            UKF::Field<2, real_t>>;
  static_assert(
      std::is_same<Plain::InnovationVector, Plain>::value,
      "Without quaternion fields the innovation is the measurement type");
  static_assert(
      std::is_same<UKF::Core<QState, Plain>::InnovationVectorType,
                   Plain>::value,
      "Without quaternion fields the innovation is the measurement type");

  const UKF::Quaternion q = rotation(real_t(0.7), UKF::Vector<3>(1, -2, 0.5));
  QFixed fixed;
  fixed.set_field<QAttitudeSensor>(q);
  fixed.set_field<QBiasSensor>(real_t(3));
  EXPECT_LT((fixed.get_field<QAttitudeSensor>().coeffs() - q.coeffs()).norm(),
            real_t(1e-6));
  EXPECT_EQ(real_t(3), fixed.get_field<QBiasSensor>());

  QDynamic dynamic;
  dynamic.set_field<QAttitudeSensor>(q);
  dynamic.set_field<QBiasSensor>(real_t(3));
  EXPECT_EQ(5, dynamic.size());
  EXPECT_EQ(4u, dynamic.covariance_size());
  EXPECT_LT((dynamic.get_field<QAttitudeSensor>().coeffs() - q.coeffs()).norm(),
            real_t(1e-6));
  EXPECT_EQ(real_t(3), dynamic.get_field<QBiasSensor>());
}

/*
The innovation of a small rotation is that rotation's vector, whichever sign
either quaternion has, and fields after a quaternion keep their values.
*/
TEST(QuaternionMeasurementTest, InnovationIsTheRotation) {
  const UKF::Quaternion predicted =
      rotation(real_t(1.2), UKF::Vector<3>(0, 1, 1));
  const UKF::Vector<3> error(real_t(0.01), real_t(-0.02), real_t(0.015));
  const UKF::Quaternion measured = rotation(error.norm(), error) * predicted;

  QFixed z_pred, z;
  z_pred.set_field<QAttitudeSensor>(predicted);
  z_pred.set_field<QBiasSensor>(real_t(1));
  z.set_field<QAttitudeSensor>(measured);
  z.set_field<QBiasSensor>(real_t(1.5));
  const QFixed::InnovationVector innovation = z_pred.calculate_innovation(z);
  EXPECT_LT((innovation.get_field<QAttitudeSensor>() - error).norm(),
            real_t(1e-5));
  EXPECT_NEAR(real_t(0.5), innovation.get_field<QBiasSensor>(), real_t(1e-6));

  QFixed flipped = z;
  flipped.set_field<QAttitudeSensor>(UKF::Quaternion(-measured.coeffs()));
  EXPECT_LT((z_pred.calculate_innovation(flipped).get_field<QAttitudeSensor>() -
             error)
                .norm(),
            real_t(1e-5));

  /* A dynamic vector keeps the order in which its fields were set. */
  for (int order = 0; order < 2; ++order) {
    QDynamic dz_pred, dz;
    if (order == 0) {
      dz_pred.set_field<QAttitudeSensor>(predicted);
      dz_pred.set_field<QBiasSensor>(real_t(1));
      dz.set_field<QAttitudeSensor>(measured);
      dz.set_field<QBiasSensor>(real_t(1.5));
    } else {
      dz_pred.set_field<QBiasSensor>(real_t(1));
      dz_pred.set_field<QAttitudeSensor>(predicted);
      dz.set_field<QBiasSensor>(real_t(1.5));
      dz.set_field<QAttitudeSensor>(measured);
    }
    const QDynamic::InnovationVector dynamic_innovation =
        dz_pred.calculate_innovation(dz);
    ASSERT_EQ(4, dynamic_innovation.size());
    EXPECT_LT((dynamic_innovation.get_field<QAttitudeSensor>() - error).norm(),
              real_t(1e-5));
    EXPECT_NEAR(real_t(0.5), dynamic_innovation.get_field<QBiasSensor>(),
                real_t(1e-6));
    EXPECT_NEAR(real_t(0.5), dynamic_innovation(order == 0 ? 3 : 0),
                real_t(1e-6));
  }
}

template <typename Measurement>
Measurement direct_observation(bool attitude_first) {
  Measurement z;
  const UKF::Quaternion measured =
      rotation(real_t(0.05), UKF::Vector<3>(0, 0, 1));
  if (attitude_first) {
    z.template set_field<QAttitudeSensor>(measured);
    z.template set_field<QBiasSensor>(real_t(1.5));
  } else {
    z.template set_field<QBiasSensor>(real_t(1.5));
    z.template set_field<QAttitudeSensor>(measured);
  }
  return z;
}

/*
At exactly half a turn the error quaternion's scalar part is +0 or -0
depending on the measurement's sign; both signs give the same innovation.
*/
TEST(QuaternionMeasurementTest, HalfTurnInnovationIgnoresSign) {
  QFixed z_pred, z, flipped;
  z_pred.set_field<QAttitudeSensor>(UKF::Quaternion::Identity());
  z_pred.set_field<QBiasSensor>(real_t(1));
  for (const UKF::Quaternion& half_turn :
       {UKF::Quaternion(0, 1, 0, 0), UKF::Quaternion(0, 0, -1, 0),
        UKF::Quaternion(0, 0, 0, 1)}) {
    z.set_field<QAttitudeSensor>(half_turn);
    z.set_field<QBiasSensor>(real_t(1));
    flipped = z;
    flipped.set_field<QAttitudeSensor>(UKF::Quaternion(-half_turn.coeffs()));
    const UKF::Vector<3> a =
        z_pred.calculate_innovation(z).get_field<QAttitudeSensor>();
    const UKF::Vector<3> b =
        z_pred.calculate_innovation(flipped).get_field<QAttitudeSensor>();
    EXPECT_TRUE(a == b);
    EXPECT_TRUE(UKF::Detail::quaternion_to_rotation_vector<QFixed>(half_turn) ==
                UKF::Detail::quaternion_to_rotation_vector<QFixed>(
                    UKF::Quaternion(-half_turn.coeffs())));
  }
}

/*
A direct observation of a small rotation behaves as a linear one: the
correction and the posterior uncertainty are those of a Kalman filter with
the prior and noise variances, in each core and with either vector type.
*/
template <template <typename...> class Core, typename Measurement>
void check_core_with(bool attitude_first) {
  Core<QState, Measurement> filter;
  filter.state.template set_field<QAttitude>(UKF::Quaternion::Identity());
  filter.state.template set_field<QBias>(real_t(1));
  const QState::CovarianceMatrix prior =
      QState::CovarianceMatrix::Identity() * prior_variance;
  typename Measurement::CovarianceVector noise;
  noise.template set_field<QAttitudeSensor>(attitude_noise);
  noise.template set_field<QBiasSensor>(bias_variance);
  if constexpr (std::is_same<Core<QState, Measurement>,
                             UKF::Core<QState, Measurement>>::value) {
    filter.covariance = prior;
    filter.process_noise_covariance.setZero();
    filter.measurement_covariance = noise;
    ASSERT_TRUE(filter.a_priori_step(real_t(1)));
  } else {
    filter.root_covariance = prior.llt().matrixL();
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance = noise.array().sqrt().matrix();
    if constexpr (std::is_same<Core<QState, Measurement>,
                               UKF::SquareRootParameterEstimationCore<
                                   QState, Measurement>>::value) {
      ASSERT_TRUE(filter.a_priori_step());
    } else {
      ASSERT_TRUE(filter.a_priori_step(real_t(1)));
    }
  }
  const Measurement z = direct_observation<Measurement>(attitude_first);
  ASSERT_TRUE(filter.innovation_step(z));
  EXPECT_EQ(4, filter.innovation.size());
  ASSERT_TRUE(filter.a_posteriori_step());

  const real_t angle = real_t(0.05);
  const real_t gain = prior_variance / (prior_variance + attitude_variance);
  const real_t posterior =
      prior_variance * attitude_variance / (prior_variance + attitude_variance);
  const real_t bias_gain = prior_variance / (prior_variance + bias_variance);
  const UKF::Vector<3> corrected =
      rotation_between(filter.state.template get_field<QAttitude>(),
                       UKF::Quaternion::Identity());
  EXPECT_LT((corrected - UKF::Vector<3>(0, 0, gain * angle)).norm(),
            real_t(2e-4));
  EXPECT_NEAR(real_t(1) + bias_gain * real_t(0.5),
              filter.state.template get_field<QBias>(), real_t(1e-4));
  EXPECT_NEAR(real_t(1), filter.state.template get_field<QAttitude>().norm(),
              real_t(1e-5));

  QState::CovarianceMatrix covariance;
  if constexpr (std::is_same<Core<QState, Measurement>,
                             UKF::Core<QState, Measurement>>::value) {
    covariance = filter.covariance;
  } else {
    covariance = filter.root_covariance * filter.root_covariance.transpose();
  }
  QState::CovarianceMatrix expected = QState::CovarianceMatrix::Zero();
  expected.diagonal() << posterior, posterior, posterior,
      prior_variance * bias_variance / (prior_variance + bias_variance);
  EXPECT_LT((covariance - expected).cwiseAbs().maxCoeff(), real_t(2e-5));
}

template <typename S, typename M>
using StandardCore = UKF::Core<S, M>;
template <typename S, typename M>
using RootCore = UKF::SquareRootCore<S, M>;
template <typename S, typename M>
using ParameterCore = UKF::SquareRootParameterEstimationCore<S, M>;

TEST(QuaternionMeasurementTest, DirectObservationMatchesLinearFilter) {
  for (bool attitude_first : {true, false}) {
    SCOPED_TRACE(attitude_first ? "attitude set first" : "bias set first");
    check_core_with<StandardCore, QFixed>(attitude_first);
    check_core_with<RootCore, QFixed>(attitude_first);
    check_core_with<ParameterCore, QFixed>(attitude_first);
    check_core_with<StandardCore, QDynamic>(attitude_first);
    check_core_with<RootCore, QDynamic>(attitude_first);
    check_core_with<ParameterCore, QDynamic>(attitude_first);
  }
}

/*
Repeated observations of an orientation far from the prior converge on it,
in the standard and square-root cores, and the innovation's rotation stays
the shorter one.
*/
TEST(QuaternionMeasurementTest, ConvergesFromAFarPrior) {
  const UKF::Quaternion truth = rotation(real_t(2.5), UKF::Vector<3>(1, 2, -1));
  QAttitudeOnly z;
  z.set_field<QAttitudeSensor>(truth);
  QAttitudeOnly::CovarianceVector noise;
  noise.set_field<QAttitudeSensor>(attitude_noise);

  UKF::Core<QState, QAttitudeOnly> filter;
  UKF::SquareRootCore<QState, QAttitudeOnly> root_filter;
  filter.state.set_field<QAttitude>(UKF::Quaternion::Identity());
  filter.state.set_field<QBias>(real_t(0));
  filter.covariance = QState::CovarianceMatrix::Identity();
  filter.process_noise_covariance =
      QState::CovarianceMatrix::Identity() * real_t(1e-6);
  filter.measurement_covariance = noise;
  root_filter.state = filter.state;
  root_filter.root_covariance = filter.covariance.llt().matrixL();
  root_filter.process_noise_root_covariance =
      filter.process_noise_covariance.llt().matrixL();
  root_filter.measurement_root_covariance = noise.array().sqrt().matrix();

  for (int i = 0; i < 30; ++i) {
    ASSERT_TRUE(filter.step(real_t(1), z));
    ASSERT_TRUE(root_filter.step(real_t(1), z));
  }
  EXPECT_LT(rotation_between(filter.state.get_field<QAttitude>(), truth).norm(),
            real_t(1e-3));
  EXPECT_LT(
      rotation_between(root_filter.state.get_field<QAttitude>(), truth).norm(),
      real_t(1e-3));
  EXPECT_LT(filter.innovation.norm(), real_t(1e-3));
  EXPECT_LT(root_filter.innovation.norm(), real_t(1e-3));
}

/*
Dense noise for a quaternion field is given in its three rotation
coordinates, and gives the same correction as the diagonal noise it equals.
*/
template <typename Measurement>
void check_dense_noise_matches_diagonal() {
  using Dense = typename Measurement::DenseCovarianceMatrix;
  UKF::Core<QState, Measurement> diagonal;
  UKF::Core<QState, Measurement, UKF::NotUsedIntegrator, Dense> dense;
  UKF::SquareRootCore<QState, Measurement, UKF::NotUsedIntegrator, Dense>
      dense_root;
  diagonal.state.template set_field<QAttitude>(
      rotation(real_t(0.3), UKF::Vector<3>(1, 0, 0)));
  diagonal.state.template set_field<QBias>(real_t(1));
  diagonal.covariance = QState::CovarianceMatrix::Identity() * prior_variance;
  diagonal.process_noise_covariance.setZero();
  diagonal.measurement_covariance.template set_field<QAttitudeSensor>(
      UKF::Vector<3>(0.01, 0.02, 0.03));
  diagonal.measurement_covariance.template set_field<QBiasSensor>(
      bias_variance);
  dense.state = diagonal.state;
  dense.covariance = diagonal.covariance;
  dense.process_noise_covariance = diagonal.process_noise_covariance;
  dense.measurement_covariance =
      Dense(diagonal.measurement_covariance.asDiagonal());
  dense_root.state = diagonal.state;
  dense_root.root_covariance = diagonal.covariance.llt().matrixL();
  dense_root.process_noise_root_covariance.setZero();
  dense_root.measurement_root_covariance = Dense(
      diagonal.measurement_covariance.array().sqrt().matrix().asDiagonal());

  const Measurement z = direct_observation<Measurement>(false);
  ASSERT_TRUE(diagonal.step(real_t(1), z));
  ASSERT_TRUE(dense.step(real_t(1), z));
  ASSERT_TRUE(dense_root.step(real_t(1), z));
  EXPECT_LT((diagonal.state - dense.state).norm(), real_t(1e-5));
  EXPECT_LT((diagonal.state - dense_root.state).norm(), real_t(1e-5));
  EXPECT_LT((diagonal.covariance - dense.covariance).norm(), real_t(1e-5));
  EXPECT_LT((diagonal.covariance - dense_root.root_covariance *
                                       dense_root.root_covariance.transpose())
                .norm(),
            real_t(1e-5));
}

TEST(QuaternionMeasurementTest, DenseNoiseMatchesDiagonal) {
  check_dense_noise_matches_diagonal<QFixed>();
  check_dense_noise_matches_diagonal<QDynamic>();
}
