#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <limits>
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
With the attitude known, a quaternion sensor's innovation is its own noise,
so the innovation covariance is the sensor's variance for any reading, and a
normalised innovation gate at the 95 % point of three degrees of freedom
rejects readings far from the prediction.
*/
template <template <typename...> class Core>
void check_known_attitude_gate() {
  const real_t threshold = real_t(7.815);
  const real_t pi = real_t(3.14159265358979323846);
  for (real_t degrees : {real_t(1), real_t(90), real_t(170), real_t(179)}) {
    SCOPED_TRACE(degrees);
    const real_t angle = degrees * pi / real_t(180);
    QAttitudeOnly z;
    z.set_field<QAttitudeSensor>(rotation(angle, UKF::Vector<3>(1, 2, -1)));
    Core<QState, QAttitudeOnly> filter;
    filter.state.template set_field<QAttitude>(UKF::Quaternion::Identity());
    filter.state.template set_field<QBias>(real_t(0));
    real_t nis = std::numeric_limits<real_t>::quiet_NaN();
    bool applied;
    UKF::Matrix<3, 3> innovation_covariance;
    if constexpr (std::is_same_v<Core<QState, QAttitudeOnly>,
                                 StandardCore<QState, QAttitudeOnly>>) {
      filter.covariance.setZero();
      filter.process_noise_covariance.setZero();
      filter.measurement_covariance.template set_field<QAttitudeSensor>(
          attitude_noise);
      auto gate = [&](const auto& current, const auto&, const auto&...) {
        const auto weighted = current.innovation_covariance.llt()
                                  .solve(current.innovation)
                                  .eval();
        nis = current.innovation.dot(weighted);
        return std::isfinite(nis) && nis <= threshold;
      };
      applied = filter.step_with_gate(real_t(0), z, gate);
      innovation_covariance = filter.innovation_covariance;
    } else {
      filter.root_covariance.setZero();
      filter.process_noise_root_covariance.setZero();
      filter.measurement_root_covariance.template set_field<QAttitudeSensor>(
          attitude_noise.array().sqrt().matrix());
      auto gate = [&](const auto& current, const auto&, const auto&...) {
        const auto whitened = current.innovation_root_covariance
                                  .template triangularView<Eigen::Lower>()
                                  .solve(current.innovation)
                                  .eval();
        nis = whitened.squaredNorm();
        return std::isfinite(nis) && nis <= threshold;
      };
      applied = filter.step_with_gate(real_t(0), z, gate);
      innovation_covariance = filter.innovation_root_covariance *
                              filter.innovation_root_covariance.transpose();
    }

    /* In the default coordinates the innovation is 2 tan(angle / 2). */
    const real_t length = real_t(2) * std::tan(angle / real_t(2));
    EXPECT_NEAR(length * length / attitude_variance, nis,
                real_t(1e-3) * length * length / attitude_variance);
    EXPECT_LT((innovation_covariance -
               UKF::Matrix<3, 3>::Identity() * attitude_variance)
                  .cwiseAbs()
                  .maxCoeff(),
              real_t(1e-5) * attitude_variance);
    EXPECT_EQ(degrees == real_t(1), applied);
    EXPECT_NEAR(real_t(1),
                std::abs(filter.state.template get_field<QAttitude>().w()),
                real_t(1e-6));
  }
}

TEST(QuaternionMeasurementTest, KnownAttitudeGateKeepsThePredictiveNoise) {
  check_known_attitude_gate<StandardCore>();
  check_known_attitude_gate<RootCore>();
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

/*
Sensor noise that is exact in one direction leaves the innovation only the
state's small variance there, and a correction keeps that variance between
zero and its prior value, in single precision too.
*/
enum { QGibbsAngle = 320, QGibbsSensor };
using QGibbsState = UKF::StateVector<UKF::Field<QGibbsAngle, real_t>>;
using QGibbsMeasurement =
    UKF::FixedMeasurementVector<UKF::Field<QGibbsSensor, UKF::Quaternion>>;

namespace UKF {
template <>
template <>
QGibbsState QGibbsState::process_model<>(real_t) const {
  return *this;
}
template <>
template <>
UKF::Quaternion
QGibbsMeasurement::expected_measurement<QGibbsState, QGibbsSensor>(
    const QGibbsState& state) {
  return UKF::Detail::rotation_vector_to_quaternion<QGibbsState>(
      UKF::Vector<3>(state.get_field<QGibbsAngle>(), 0, 0));
}
}  // namespace UKF

template <template <typename...> class Core>
void check_exact_noise_direction() {
  using Dense = QGibbsMeasurement::DenseCovarianceMatrix;
  using Filter =
      Core<QGibbsState, QGibbsMeasurement, UKF::NotUsedIntegrator, Dense>;
  constexpr bool standard =
      std::is_same_v<Filter, UKF::Core<QGibbsState, QGibbsMeasurement,
                                       UKF::NotUsedIntegrator, Dense>>;
  UKF::Matrix<3, 3> factor;
  factor << real_t(0.125), 0, 0, real_t(-0.125), 0, 0, 0, 0, 2;
  Filter filter;
  filter.state.template set_field<QGibbsAngle>(real_t(0));
  if constexpr (standard) {
    filter.covariance << real_t(1e-8);
    filter.process_noise_covariance << 0;
    filter.measurement_covariance = factor * factor.transpose();
  } else {
    filter.root_covariance << real_t(1e-4);
    filter.process_noise_root_covariance << 0;
    filter.measurement_root_covariance = factor;
  }
  const real_t pi = real_t(3.14159265358979323846);
  QGibbsMeasurement z;
  z.set_field<QGibbsSensor>(rotation(pi / real_t(2), UKF::Vector<3>::UnitZ()));
  ASSERT_TRUE(filter.step(real_t(0), z));
  real_t variance;
  if constexpr (standard) {
    variance = filter.covariance(0, 0);
  } else {
    variance = filter.root_covariance(0, 0) * filter.root_covariance(0, 0);
  }
  EXPECT_GE(variance, real_t(0));
  EXPECT_LE(variance, real_t(1e-8));
}

TEST(QuaternionMeasurementTest, ExactNoiseDirectionKeepsTheStateVariance) {
  check_exact_noise_direction<UKF::Core>();
  check_exact_noise_direction<UKF::SquareRootCore>();
}

/*
A split-step caller's innovation covariance is checked before the correction
uses it, here behind a 90-degree quaternion reading with large noise.
*/
TEST(QuaternionMeasurementTest, SplitPosteriorChecksThePublicCovariance) {
  UKF::Core<QState, QAttitudeOnly> filter;
  filter.state.set_field<QAttitude>(UKF::Quaternion::Identity());
  filter.state.set_field<QBias>(real_t(0));
  filter.covariance = QState::CovarianceMatrix::Identity() * prior_variance;
  filter.process_noise_covariance.setZero();
  filter.measurement_covariance.set_field<QAttitudeSensor>(
      UKF::Vector<3>::Constant(real_t(4)));
  const real_t pi = real_t(3.14159265358979323846);
  QAttitudeOnly z;
  z.set_field<QAttitudeSensor>(
      rotation(pi / real_t(2), UKF::Vector<3>::UnitZ()));
  ASSERT_TRUE(filter.a_priori_step(real_t(0)));
  ASSERT_TRUE(filter.innovation_step(z));
  const QState predicted = filter.state;
  const QState::CovarianceMatrix predicted_covariance = filter.covariance;
  const QAttitudeOnly::CovarianceMatrix valid = filter.innovation_covariance;
  filter.innovation_covariance(0, 0) = real_t(-1);
  EXPECT_FALSE(filter.a_posteriori_step());
  EXPECT_TRUE(predicted == filter.state);
  EXPECT_TRUE(predicted_covariance == filter.covariance);
  filter.innovation_covariance = valid;
  EXPECT_TRUE(filter.a_posteriori_step());
}

/*
A split-step caller's edit to the public innovation covariance or root takes
effect in the correction, whether or not the reading sets the quaternion
field that its type declares. The bias reading is read first and is
independent of the attitude, so its gain is the edited one either way.
*/
template <typename Filter>
void check_split_edit(bool attitude) {
  Filter filter;
  filter.state.template set_field<QAttitude>(UKF::Quaternion::Identity());
  filter.state.template set_field<QBias>(real_t(0));
  typename QDynamic::CovarianceVector noise;
  noise.template set_field<QAttitudeSensor>(UKF::Vector<3>::Ones());
  noise.template set_field<QBiasSensor>(real_t(1));
  QDynamic z;
  z.template set_field<QBiasSensor>(real_t(1));
  if (attitude) {
    z.template set_field<QAttitudeSensor>(rotation(
        real_t(3.14159265358979323846) / real_t(2), UKF::Vector<3>::UnitZ()));
  }
  real_t variance;
  if constexpr (std::is_same_v<Filter, UKF::Core<QState, QDynamic>>) {
    filter.covariance.setIdentity();
    filter.process_noise_covariance.setZero();
    filter.measurement_covariance = noise;
    ASSERT_TRUE(filter.a_priori_step(real_t(0)));
    ASSERT_TRUE(filter.innovation_step(z));
    filter.innovation_covariance(0, 0) = real_t(8);
    ASSERT_TRUE(filter.a_posteriori_step());
    variance = filter.covariance(3, 3);
  } else {
    filter.root_covariance.setIdentity();
    filter.process_noise_root_covariance.setZero();
    filter.measurement_root_covariance = noise;
    if constexpr (std::is_same_v<Filter, UKF::SquareRootParameterEstimationCore<
                                             QState, QDynamic>>) {
      ASSERT_TRUE(filter.a_priori_step());
    } else {
      ASSERT_TRUE(filter.a_priori_step(real_t(0)));
    }
    ASSERT_TRUE(filter.innovation_step(z));
    filter.innovation_root_covariance(0, 0) = std::sqrt(real_t(8));
    ASSERT_TRUE(filter.a_posteriori_step());
    variance =
        (filter.root_covariance * filter.root_covariance.transpose())(3, 3);
  }
  EXPECT_NEAR(real_t(0.125), filter.state.template get_field<QBias>(),
              real_t(1e-5));
  EXPECT_NEAR(real_t(0.875), variance, real_t(1e-5));
}

TEST(QuaternionMeasurementTest, SplitEditsTakeEffect) {
  for (bool attitude : {false, true}) {
    SCOPED_TRACE(attitude);
    check_split_edit<UKF::Core<QState, QDynamic>>(attitude);
    check_split_edit<UKF::SquareRootCore<QState, QDynamic>>(attitude);
    check_split_edit<UKF::SquareRootParameterEstimationCore<QState, QDynamic>>(
        attitude);
  }
}
