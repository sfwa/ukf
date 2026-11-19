/*
Copyright (C) 2016 Thiemar Pty Ltd

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef MEASUREMENTVECTOR_H
#define MEASUREMENTVECTOR_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

#include "UKF/StateVector.h"
#include "UKF/Types.h"

namespace UKF {

namespace Detail {

/*
Store a field value in its segment of a vector: a quaternion as its
coefficients x, y, z, w, as in a state vector.
*/
template <typename Segment, typename T>
inline void assign_field(Segment&& segment, const T& in) {
  if constexpr (std::is_same<T, Quaternion>::value) {
    segment << in.vec(), in.w();
  } else {
    segment << in;
  }
}

/* A deterministic unit axis perpendicular to a nonzero field vector. */
inline Vector<3> antiparallel_axis(const Vector<3>& v) {
  if (std::abs(v(0)) > std::abs(v(1))) {
    return Vector<3>(-v(2), 0, v(0)).normalized();
  }
  return Vector<3>(0, -v(2), v(1)).normalized();
}

/*
Scale a field vector to unit length, returning its length, or zero for a
zero vector. Dividing by the largest component first keeps the squared
norm clear of underflow and overflow at any physical scale.
*/
inline real_t unit_direction(const Vector<3>& v, Vector<3>& direction) {
  const real_t scale = v.cwiseAbs().maxCoeff();
  if (scale == real_t(0)) {
    direction = v;
    return real_t(0);
  }
  direction = v / scale;
  const real_t length = direction.norm();
  direction /= length;
  return scale * length;
}

/*
For unit directions the shortest rotation has quaternion components
(u1 cross u2, 1 + u1 dot u2), up to a common scale. Near a half turn,
|u1 + u2|^2 / 2 gives the same scalar without subtracting nearly equal
numbers. Only directions indistinguishable from opposites need a chosen axis.
*/
inline void direction_rotation(const Vector<3>& u1, const Vector<3>& u2,
                               Vector<3>& axis, real_t& w) {
  const real_t dot = u1.dot(u2);
  if (dot < real_t(0)) {
    const Vector<3> sum = u1 + u2;
    if (sum.cwiseAbs().maxCoeff() <=
        real_t(4) * std::numeric_limits<real_t>::epsilon()) {
      axis = antiparallel_axis(u1);
      w = real_t(0);
    } else {
      axis = u1.cross(sum);
      w = real_t(0.5) * sum.squaredNorm();
    }
  } else {
    axis = u1.cross(u2);
    w = real_t(1) + dot;
  }
}

/*
Calculate the smallest rotation vector which transforms v1 into v2. The
rotation depends only on the directions of the two vectors, so the
calculation is done on unit vectors: its tolerances are then relative,
and a field given in teslas (about 5e-5) is handled like one in gauss. A
zero vector has no direction, and gives no rotation.
*/
template <typename T>
inline Vector<3> calculate_rotation_vector(const Vector<3>& v2,
                                           const Vector<3>& v1) {
  Vector<3> u1, u2;
  if (unit_direction(v1, u1) == real_t(0) ||
      unit_direction(v2, u2) == real_t(0)) {
    return Vector<3>::Zero();
  }
  Vector<3> axis;
  real_t w;
  direction_rotation(u1, u2, axis, w);
  if constexpr (Parameters::MRP_A<T> == real_t(0)) {
    // The common quaternion scale cancels for the default Gibbs coordinates.
    return Parameters::MRP_F<T> * axis /
           (w == real_t(0) ? std::numeric_limits<real_t>::epsilon() : w);
  }
  const real_t norm = std::sqrt(axis.squaredNorm() + w * w);
  const real_t denominator = Parameters::MRP_A<T> * norm + w;
  const real_t minimum = std::numeric_limits<real_t>::epsilon() * norm;
  return Parameters::MRP_F<T> * axis /
         (std::abs(denominator) > minimum ? denominator : minimum);
}

/*
Calculate the transformation matrix to transform a covariance matrix of
field vector noise for the field vector v2 into a covariance matrix for the
rotation vector corresponding to the transformation between the field
vector v1 and v2.

Differentiate r = f*c/(a*n + w), where c = v1 cross v2,
w = |v1|*|v2| + v1 dot v2, and n = sqrt(c dot c + w*w).
At an exactly antiparallel pair the rotation axis is not unique; use
the same deterministic axis convention as calculate_rotation_vector.
*/
template <typename T>
inline Matrix<3, 3> rotation_vector_jacobian_of(const Vector<3>& v2,
                                                const Vector<3>& v1) {
  Vector<3> c = Vector<3>::Zero();
  real_t w = real_t(0);
  Vector<3> dw = v1;
  if (!v1.isZero(real_t(0)) && !v2.isZero(real_t(0))) {
    // Both arguments here are unit vectors when both directions exist.
    direction_rotation(v1, v2, c, w);
    dw += v2;
  }
  Matrix<3, 3> dc;
  dc << 0, -v1(2), v1(1), v1(2), 0, -v1(0), -v1(1), v1(0), 0;

  const real_t quaternion_norm = std::sqrt(c.squaredNorm() + w * w);
  real_t denominator = Parameters::MRP_A<T> * quaternion_norm + w;
  Vector<3> denominator_gradient = dw;
  if (quaternion_norm > real_t(0)) {
    denominator_gradient +=
        Parameters::MRP_A<T> * (dc.transpose() * c + w * dw) / quaternion_norm;
  }
  const real_t minimum =
      std::numeric_limits<real_t>::epsilon() *
      (quaternion_norm > real_t(0) ? quaternion_norm : real_t(1));
  if (std::abs(denominator) < minimum) {
    denominator = minimum;
  }
  return Parameters::MRP_F<T> *
         (dc / denominator -
          c * denominator_gradient.transpose() / (denominator * denominator));
}

/*
The rotation vector does not change when either field vector is scaled,
so its Jacobian at v2 is the Jacobian at the unit vectors divided by
|v2|: the scale of the field noise carries through to the rotation. The
unit vectors keep the tolerances above relative. A zero vector, which
has no direction, is passed through as it is.
*/
template <typename T>
inline Matrix<3, 3> calculate_rotation_vector_jacobian(const Vector<3>& v2,
                                                       const Vector<3>& v1) {
  Vector<3> u1, u2;
  const real_t length = unit_direction(v2, u2);
  if (unit_direction(v1, u1) == real_t(0) || length == real_t(0)) {
    return rotation_vector_jacobian_of<T>(v2, v1);
  }
  return rotation_vector_jacobian_of<T>(u2, u1) / length;
}

/*
This helper class links a MeasurementVector type and a StateVector type
in order to allow function overloading without ambiguity. These functions
are also shared between the FixedMeasurementVector and
DynamicMeasurementVector to reduce duplicated code.
*/
template <typename M, typename S = StateVector<UKF::Field<0, real_t>>>
class MeasurementStateHelper {
 public:
  /* Map sensor noise into the corresponding innovation coordinates. */
  template <typename T>
  static Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>
  noise_jacobian([[maybe_unused]] const T& z_pred,
                 [[maybe_unused]] const T& z) {
    return Matrix<Detail::CovarianceDimension<T>,
                  Detail::CovarianceDimension<T>>::Identity();
  }

  static Matrix<3, 3> noise_jacobian(const FieldVector& z_pred,
                                     const FieldVector& z) {
    return Detail::calculate_rotation_vector_jacobian<M>(z, z_pred);
  }

  /*
  A quaternion sensor's noise n is a left perturbation of the measured
  rotation, z = q(n) h(x), and its variances are given in that perturbation.
  The innovation is u = r(z z_pred^-1). With the prediction error
  e = r(h(x) z_pred^-1), q(u) = q(n) q(e), so n = r(q(u) q(e)^-1), whose
  derivative in u at u = e is the attitude reset Jacobian J(e) of a state
  correction.

  Before the observation, e is zero on average and the noise enters the
  innovation unmapped, so the innovation covariance that gates and
  diagnostics read keeps the sensor's own variances, through the generic
  noise_jacobian above. Given the observation, the correction maps the noise
  through J(u)^-1 at the observed innovation, which widens it as the
  innovation grows; see calculate_correction_noise_transform. Judging an
  observation by noise widened by that same observation would let an outlier
  pass a gate, so the widening is kept out of the innovation covariance.
  */
  static Matrix<3, 3> correction_noise_jacobian(const Quaternion& z_pred,
                                                const Quaternion& z) {
    return Detail::quaternion_reset_jacobian<M>(
               Detail::quaternion_to_rotation_vector<M>(z * z_pred.conjugate()))
        .inverse();
  }

  /*
  Functions for calculating the covariance of a field in the
  measurement vector. The noise p has the field's covariance type: for a
  quaternion, the variances of the rotation's three components.
  */
  template <typename P, typename T>
  static Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>
  field_covariance(const P& p, [[maybe_unused]] const T& z_pred,
                   [[maybe_unused]] const T& z) {
    Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>
        temp = Matrix<Detail::CovarianceDimension<T>,
                      Detail::CovarianceDimension<T>>::Zero();
    if constexpr (std::is_arithmetic_v<P>) {
      temp(0, 0) = p;
    } else {
      // Avoid the nested diagonal block created by Eigen 3.2's comma
      // initializer, which GCC 13 can miscompile in optimized builds.
      temp.diagonal() = p;
    }
    return temp;
  }

  static Matrix<3, 3> field_covariance(const FieldVector& p,
                                       const FieldVector& z_pred,
                                       const FieldVector& z) {
    Matrix<3, 3> T = Detail::calculate_rotation_vector_jacobian<M>(z, z_pred);

    /*
    To calculate the covariance, pre-multiply by the transformation
    matrix and then post-multiply by the transformation matrix
    transpose.
    */
    return T * Eigen::DiagonalMatrix<real_t, 3>(p) * T.transpose();
  }

  template <typename P, typename T>
  static Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>
  field_root_covariance(const P& p, [[maybe_unused]] const T& z_pred,
                        [[maybe_unused]] const T& z) {
    Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>
        temp = Matrix<Detail::CovarianceDimension<T>,
                      Detail::CovarianceDimension<T>>::Zero();
    if constexpr (std::is_arithmetic_v<P>) {
      temp(0, 0) = p;
    } else {
      // Use direct assignment here for the same reason as above.
      temp.diagonal() = p;
    }
    return temp;
  }

  static Matrix<3, 3> field_root_covariance(const FieldVector& p,
                                            const FieldVector& z_pred,
                                            const FieldVector& z) {
    /*
    To calculate the root covariance, we simply pre-multiply it by
    the transformation matrix. This will yield a positive-indefinite
    and possibly non-triangular matrix, but it can be shown that
    since it multiplied by its transpose it equal to the covariance
    matrix, it gives the correct result in QR decomposition used by
    the square-root filter.

    A proof of this is left as an exercise to the reader.
    */
    return Detail::calculate_rotation_vector_jacobian<M>(z, z_pred) *
           Eigen::DiagonalMatrix<real_t, 3>(p);
  }

  /*
  Functions for calculating the mean of each field in a sigma point
  distribution. The field type is an explicit template argument and
  the rows of the sigma points are an Eigen expression, so no overload
  has to be chosen from a converted block or a dummy field value; every
  field mean is a vector of the field's stored dimension.

  Means are accumulated as offsets from the central point; see
  StateVector. The field vector mean is found by first calculating the
  set of rotation vectors which transforms the central sigma point into
  the set of sigma points. Then, calculate the mean of these rotations
  and apply it to the central sigma point. A quaternion's mean is the
  renormalised barycentric mean, with each point taken in the
  hemisphere of the central one, as for a quaternion state.
  */
  template <typename T, typename Derived>
  static Vector<Detail::StateVectorDimension<T>> sigma_point_mean(
      const Eigen::MatrixBase<Derived>& sigma) {
    if constexpr (std::is_same<T, Quaternion>::value) {
      Array<1, S::num_sigma() - 1> signs;
      for (std::size_t i = 0; i < S::num_sigma() - 1; i++) {
        signs(i) = sigma.col(i + 1).dot(sigma.col(0)) < real_t(0) ? real_t(-1)
                                                                  : real_t(1);
      }
      const Vector<4> temp =
          sigma.col(0) + Parameters::Sigma_WMI<S> *
                             ((sigma.template rightCols<S::num_sigma() - 1>()
                                   .array()
                                   .rowwise() *
                               signs)
                                  .matrix()
                                  .colwise() -
                              sigma.col(0))
                                 .rowwise()
                                 .sum();
      const Quaternion temp_q = Quaternion(temp).normalized();
      return Vector<4>(temp_q.x(), temp_q.y(), temp_q.z(), temp_q.w());
    } else if constexpr (std::is_same<T, FieldVector>::value) {
      Vector<3> temp = Vector<3>::Zero();

      for (std::size_t i = 0; i < S::num_sigma() - 1; i++) {
        temp += Parameters::Sigma_WMI<S> * Detail::calculate_rotation_vector<M>(
                                               sigma.col(i + 1), sigma.col(0));
      }

      return Detail::rotation_vector_to_quaternion<M>(temp) *
             Vector<3>(sigma.col(0));
    } else {
      return sigma.col(0) +
             Parameters::Sigma_WMI<S> *
                 (sigma.template rightCols<S::num_sigma() - 1>().colwise() -
                  sigma.col(0))
                     .rowwise()
                     .sum();
    }
  }

  /*
  Functions for calculating the difference between two fields in a
  measurement vector, used when calculating the innovation.
  */
  template <typename T>
  static Detail::CovarianceFieldType<T> measurement_delta(const T& z,
                                                          const T& z_pred) {
    if constexpr (std::is_same<T, Quaternion>::value) {
      /*
      A quaternion's innovation is the rotation vector of the
      error quaternion z * z_pred^-1, in the convention used for
      quaternion states.
      */
      return Detail::quaternion_to_rotation_vector<M>(z * z_pred.conjugate());
    } else if constexpr (std::is_same<T, FieldVector>::value) {
      return Detail::calculate_rotation_vector<M>(z, z_pred);
    } else {
      return z - z_pred;
    }
  }

  /*
  Functions for calculating the difference between each point in a
  sigma point distribution and the mean.
  */
  template <typename T, typename Derived>
  static Matrix<Detail::CovarianceDimension<T>, S::num_sigma()>
  sigma_point_deltas(const T& mean, const Eigen::MatrixBase<Derived>& Z) {
    if constexpr (std::is_same<T, Quaternion>::value) {
      Matrix<3, S::num_sigma()> temp;

      for (std::size_t i = 0; i < S::num_sigma(); i++) {
        temp.col(i) = Detail::quaternion_to_rotation_vector<M>(
            Quaternion(Z.col(i)) * mean.conjugate());
      }

      return temp;
    } else if constexpr (std::is_same<T, FieldVector>::value) {
      Matrix<3, S::num_sigma()> temp;

      for (std::size_t i = 0; i < S::num_sigma(); i++) {
        temp.col(i) = Detail::calculate_rotation_vector<M>(Z.col(i), mean);
      }

      return temp;
    } else if constexpr (std::is_arithmetic<T>::value) {
      return Z.array() - mean;
    } else {
      return Z.colwise() - mean;
    }
  }
};

}  // namespace Detail

/* Alias for the Eigen type from which FixedMeasurementVector inherits. */
template <typename... Fields>
using MeasurementVectorFixedBaseType =
    Vector<Detail::get_composite_vector_dimension<Fields...>()>;

/* Alias for the Eigen type from which DynamicMeasurementVector inherits. */
template <typename... Fields>
using MeasurementVectorDynamicBaseType =
    VectorDynamic<Detail::get_composite_vector_dimension<Fields...>()>;

/*
This class provides a fixed measurement vector, to be used when the same
measurement are available every time step.
*/
template <typename... Fields>
class FixedMeasurementVector
    : public MeasurementVectorFixedBaseType<typename Fields::type...> {
 public:
  using Base = MeasurementVectorFixedBaseType<typename Fields::type...>;
  using Base::Base;
  using Base::operator=;

  /* Get size of measurement vector. */
  static constexpr std::size_t size() {
    return Detail::get_composite_vector_dimension<typename Fields::type...>();
  }

  /* Get size of measurement vector covariance. */
  static constexpr std::size_t covariance_size() {
    return Detail::get_covariance_dimension<typename Fields::type...>();
  }

  /* Aliases for types needed during filter iteration. */
  template <typename S>
  using SigmaPointDistribution =
      Matrix<FixedMeasurementVector::size(), S::num_sigma()>;
  template <typename S>
  using SigmaPointDeltas =
      Matrix<FixedMeasurementVector::covariance_size(), S::num_sigma()>;
  using CovarianceMatrix = Matrix<covariance_size(), covariance_size()>;
  /*
  Diagonal sensor noise, and the innovation, have one entry for each
  component of a field's uncertainty: three, the rotation, for a
  quaternion. Without quaternion fields these are this vector type.
  Experimental: quaternion measurement fields, and these types where they
  differ from the measurement type, may change in a minor release (see
  RELEASING.md).
  */
  using CovarianceVector = FixedMeasurementVector<Field<
      Fields::key, Detail::CovarianceFieldType<typename Fields::type>>...>;
  using InnovationVector = CovarianceVector;
  // Dense sensor noise always follows the declared field order.
  using DenseCovarianceMatrix = CovarianceMatrix;

  /* Functions for accessing individual fields. */
  template <int Key>
  typename Detail::FieldTypes<Key, Fields...>::type get_field() const {
    static_assert(Detail::get_field_offset<0, Fields...>(Key) !=
                      std::numeric_limits<std::size_t>::max(),
                  "Specified key not present in state vector");
    return Detail::convert_from_segment<
        typename Detail::FieldTypes<Key, Fields...>::type>(
        Base::template segment<Detail::get_field_size<Fields...>(Key)>(
            Detail::get_field_offset<0, Fields...>(Key)));
  }

  template <int Key, typename T>
  void set_field(T in) {
    static_assert(Detail::get_field_offset<0, Fields...>(Key) !=
                      std::numeric_limits<std::size_t>::max(),
                  "Specified key not present in state vector");
    Detail::assign_field(
        Base::template segment<Detail::get_field_size<Fields...>(Key)>(
            Detail::get_field_offset<0, Fields...>(Key)),
        in);
  }

  /* Calculate the mean from a measurement sigma point distribution. */
  template <typename S>
  FixedMeasurementVector calculate_sigma_point_mean(
      const SigmaPointDistribution<S>& Z) const {
    FixedMeasurementVector mean;
    calculate_field_mean<S, Fields...>(Z, mean);

    return mean;
  }

  /*
  Calculate the set of sigma point delta vectors; these are used for
  calculating the measurement covariance and the Kalman gain.
  */
  template <typename S>
  SigmaPointDeltas<S> calculate_sigma_point_deltas(
      const SigmaPointDistribution<S>& Z) const {
    SigmaPointDeltas<S> z_prime;

    /* Calculate the delta vectors. */
    calculate_field_deltas<S, Fields...>(Z, z_prime);

    return z_prime;
  }

  /*
  Calculate the expected measurement covariance covariance as described in
  equation 68 of the Kraft papers.
  */
  template <typename S>
  CovarianceMatrix calculate_sigma_point_covariance(
      const SigmaPointDeltas<S>& z_prime) const {
    CovarianceMatrix cov;

    /* Calculate the covariance using equation 64 from the Kraft paper. */
    cov.noalias() =
        (Parameters::Sigma_WCI<S> *
         z_prime.template rightCols<S::num_sigma() - 1>()) *
        z_prime.template rightCols<S::num_sigma() - 1>().transpose();
    cov.noalias() += (Parameters::Sigma_WC0<S> * z_prime.col(0)) *
                     z_prime.col(0).transpose();

    return cov;
  }

  /*
  Create a measurement sigma point distribution using the sigma points.
  */
  template <typename S, typename... U>
  SigmaPointDistribution<S> calculate_sigma_point_distribution(
      const typename S::SigmaPointDistribution& X, const U&... input) const {
    SigmaPointDistribution<S> Z(size(), S::num_sigma());

    /* Reference the inputs once rather than copying them per sigma point. */
    const std::tuple<const U&...> inputs{input...};
    for (std::size_t i = 0; i < S::num_sigma(); i++) {
      FixedMeasurementVector temp;
      calculate_field_measurements<S, const std::tuple<const U&...>&,
                                   Fields...>(temp, X.col(i), inputs);
      Z.col(i) = temp;
    }

    return Z;
  }

  /*
  Returns a matrix containing the measurement covariance. Used for the
  standard UKF.
  */
  CovarianceMatrix calculate_measurement_covariance(
      const CovarianceVector& measurement_covariance,
      const FixedMeasurementVector& z_pred) const {
    CovarianceMatrix temp = CovarianceMatrix::Zero();

    calculate_field_covariance<Fields...>(temp, measurement_covariance, z_pred);
    return temp;
  }

  /*
  Returns a matrix containing the measurement root covariance. Used for the
  square-root UKF.
  */
  CovarianceMatrix calculate_measurement_root_covariance(
      const CovarianceVector& measurement_root_covariance,
      const FixedMeasurementVector& z_pred) const {
    CovarianceMatrix temp = CovarianceMatrix::Zero();

    calculate_field_root_covariance<Fields...>(
        temp, measurement_root_covariance, z_pred);
    return temp;
  }

  /* Dense noise is expressed in raw sensor coordinates and declared field
     order. */
  CovarianceMatrix calculate_measurement_covariance(
      const DenseCovarianceMatrix& measurement_covariance,
      const FixedMeasurementVector& z_pred) const {
    if constexpr (has_field_vector) {
      DenseCovarianceMatrix transform = DenseCovarianceMatrix::Zero();
      (calculate_field_noise_transform<Fields>(transform, z_pred), ...);
      return transform * measurement_covariance * transform.transpose();
    } else {
      return measurement_covariance;
    }
  }

  /* Keep every factor column: selecting columns would discard noise variance.
   */
  DenseCovarianceMatrix calculate_measurement_root_covariance(
      const DenseCovarianceMatrix& measurement_root_covariance,
      const FixedMeasurementVector& z_pred) const {
    if constexpr (has_field_vector) {
      DenseCovarianceMatrix transform = DenseCovarianceMatrix::Zero();
      (calculate_field_noise_transform<Fields>(transform, z_pred), ...);
      return transform * measurement_root_covariance;
    } else {
      return measurement_root_covariance;
    }
  }

  /* Resolve Eigen expressions before conversion can match both noise types. */
  template <typename Derived>
  CovarianceMatrix calculate_measurement_covariance(
      const Eigen::MatrixBase<Derived>& measurement_covariance,
      const FixedMeasurementVector& z_pred) const {
    if constexpr (Derived::ColsAtCompileTime == 1) {
      return calculate_measurement_covariance(
          CovarianceVector(measurement_covariance), z_pred);
    } else {
      if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (measurement_covariance.cols() == 1) {
          return calculate_measurement_covariance(
              CovarianceVector(measurement_covariance), z_pred);
        }
      }
      return calculate_measurement_covariance(
          DenseCovarianceMatrix(measurement_covariance), z_pred);
    }
  }

  template <typename Derived>
  CovarianceMatrix calculate_measurement_root_covariance(
      const Eigen::MatrixBase<Derived>& measurement_root_covariance,
      const FixedMeasurementVector& z_pred) const {
    if constexpr (Derived::ColsAtCompileTime == 1) {
      return calculate_measurement_root_covariance(
          CovarianceVector(measurement_root_covariance), z_pred);
    } else {
      if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (measurement_root_covariance.cols() == 1) {
          return calculate_measurement_root_covariance(
              CovarianceVector(measurement_root_covariance), z_pred);
        }
      }
      return calculate_measurement_root_covariance(
          DenseCovarianceMatrix(measurement_root_covariance), z_pred);
    }
  }

  /*
  Whether a field is a quaternion, whose noise the correction maps to the
  observed innovation.
  */
  static constexpr bool has_quaternion_field =
      (std::is_same_v<typename Fields::type, Quaternion> || ...);

  /* Whether this reading holds a quaternion field. */
  constexpr bool contains_quaternion_field() const {
    return has_quaternion_field;
  }

  /*
  The transform from the noise in the innovation covariance to the noise
  that the correction uses: J(u)^-1 at the observed innovation u for each
  quaternion field and the identity elsewhere, so that a noise covariance R
  becomes T R T^T and a noise factor L becomes T L. Called on the observed
  vector with the predicted one.
  */
  CovarianceMatrix calculate_correction_noise_transform(
      const FixedMeasurementVector& z_pred) const {
    CovarianceMatrix transform = CovarianceMatrix::Identity();
    (calculate_field_correction_transform<Fields>(transform, z_pred), ...);
    return transform;
  }

  /* Return the innovation using the supplied measurement vector. */
  InnovationVector calculate_innovation(const FixedMeasurementVector& z) const {
    InnovationVector temp;

    (calculate_field_innovation<Fields>(z, temp), ...);
    return temp;
  }

 private:
  // Only FieldVector noise needs conversion from sensor to innovation
  // coordinates.
  static constexpr bool has_field_vector =
      (std::is_same_v<typename Fields::type, FieldVector> || ...);

  /*
  This function is intended to be specialised by the user for each field in
  the measurement vector, and allows the user to specify how a particular
  state vector is transformed into a measurement vector.

  Template parameters are a StateVector type, a field key and a parameter
  pack of input vectors.
  */
  template <typename S, int Key, typename... U>
  static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement(
      const S& state, const U&... input);

  template <typename S, int Key, typename U, size_t... I>
  static typename Detail::FieldTypes<Key, Fields...>::type
  expected_measurement_helper(const S& state, U&& input,
                              std::index_sequence<I...>) {
    return expected_measurement<S, Key>(state, std::get<I>(input)...);
  }

  /*
  These functions build the measurement estimate from the expected
  measurement of each individual field.
  */
  template <typename S, typename U, typename T>
  static void calculate_field_measurements(FixedMeasurementVector& expected,
                                           const S& state, U&& input) {
    constexpr std::size_t len =
        std::tuple_size<typename std::remove_reference<U>::type>::value;

    expected.template set_field<T::key>(
        expected_measurement_helper<S, T::key, U>(
            state, std::forward<U>(input), std::make_index_sequence<len>()));
  }

  template <typename S, typename U, typename T1, typename T2, typename... Tail>
  static void calculate_field_measurements(FixedMeasurementVector& expected,
                                           const S& state, U&& input) {
    calculate_field_measurements<S, U, T1>(expected, state,
                                           std::forward<U>(input));
    calculate_field_measurements<S, U, T2, Tail...>(expected, state,
                                                    std::forward<U>(input));
  }

  /*
  These functions build the measurement covariance from all fields in the
  measurement vector.
  */
  template <typename T>
  void calculate_field_correction_transform(
      CovarianceMatrix& transform, const FixedMeasurementVector& z_pred) const {
    if constexpr (std::is_same_v<typename T::type, Quaternion>) {
      constexpr std::size_t offset =
          Detail::get_field_covariance_offset<0, Fields...>(T::key);
      transform.template block<3, 3>(offset, offset) =
          Detail::MeasurementStateHelper<FixedMeasurementVector>::
              correction_noise_jacobian(z_pred.template get_field<T::key>(),
                                        get_field<T::key>());
    }
  }

  template <typename T>
  void calculate_field_noise_transform(
      DenseCovarianceMatrix& transform,
      const FixedMeasurementVector& z_pred) const {
    constexpr std::size_t offset =
        Detail::get_field_covariance_offset<0, Fields...>(T::key);
    transform.template block<Detail::CovarianceDimension<typename T::type>,
                             Detail::CovarianceDimension<typename T::type>>(
        offset, offset) =
        Detail::MeasurementStateHelper<FixedMeasurementVector>::noise_jacobian(
            z_pred.template get_field<T::key>(), get_field<T::key>());
  }

  template <typename T>
  void calculate_field_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_covariance,
      const FixedMeasurementVector& z_pred) const {
    P.template block<Detail::CovarianceDimension<typename T::type>,
                     Detail::CovarianceDimension<typename T::type>>(
        Detail::get_field_covariance_offset<0, Fields...>(T::key),
        Detail::get_field_covariance_offset<0, Fields...>(T::key)) = Detail::
        MeasurementStateHelper<FixedMeasurementVector>::field_covariance(
            measurement_covariance.template get_field<T::key>(),
            z_pred.get_field<T::key>(), get_field<T::key>());
  }

  template <typename T1, typename T2, typename... Tail>
  void calculate_field_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_covariance,
      const FixedMeasurementVector& z_pred) const {
    calculate_field_covariance<T1>(P, measurement_covariance, z_pred);
    calculate_field_covariance<T2, Tail...>(P, measurement_covariance, z_pred);
  }

  /*
  These functions build the measurement root covariance from all fields in
  the measurement vector.
  */
  template <typename T>
  void calculate_field_root_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_root_covariance,
      const FixedMeasurementVector& z_pred) const {
    P.template block<Detail::CovarianceDimension<typename T::type>,
                     Detail::CovarianceDimension<typename T::type>>(
        Detail::get_field_covariance_offset<0, Fields...>(T::key),
        Detail::get_field_covariance_offset<0, Fields...>(T::key)) = Detail::
        MeasurementStateHelper<FixedMeasurementVector>::field_root_covariance(
            measurement_root_covariance.template get_field<T::key>(),
            z_pred.get_field<T::key>(), get_field<T::key>());
  }

  template <typename T1, typename T2, typename... Tail>
  void calculate_field_root_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_root_covariance,
      const FixedMeasurementVector& z_pred) const {
    calculate_field_root_covariance<T1>(P, measurement_root_covariance, z_pred);
    calculate_field_root_covariance<T2, Tail...>(P, measurement_root_covariance,
                                                 z_pred);
  }

  /*
  Private functions for calculating the mean of each field in a sigma point
  distribution.
  */
  template <typename S, typename T>
  void calculate_field_mean(const SigmaPointDistribution<S>& Z,
                            FixedMeasurementVector& mean) const {
    mean.template segment<Detail::StateVectorDimension<typename T::type>>(
        Detail::get_field_offset<0, Fields...>(T::key)) =
        Detail::MeasurementStateHelper<FixedMeasurementVector, S>::
            template sigma_point_mean<typename T::type>(
                Z.template block<Detail::StateVectorDimension<typename T::type>,
                                 S::num_sigma()>(
                    Detail::get_field_offset<0, Fields...>(T::key), 0));
  }

  template <typename S, typename T1, typename T2, typename... Tail>
  void calculate_field_mean(const SigmaPointDistribution<S>& Z,
                            FixedMeasurementVector& mean) const {
    calculate_field_mean<S, T1>(Z, mean);
    calculate_field_mean<S, T2, Tail...>(Z, mean);
  }

  /* These functions build the innovation from all populated fields. */
  template <typename T>
  void calculate_field_innovation(const FixedMeasurementVector& z,
                                  InnovationVector& innovation) const {
    innovation.template set_field<T::key>(
        Detail::MeasurementStateHelper<FixedMeasurementVector>::
            measurement_delta(z.template get_field<T::key>(),
                              get_field<T::key>()));
  }

  /*
  Private functions for calculating the delta vectors between each sigma
  point and the mean.
  */
  template <typename S, typename T>
  void calculate_field_deltas(const SigmaPointDistribution<S>& Z,
                              SigmaPointDeltas<S>& z_prime) const {
    z_prime.template block<Detail::CovarianceDimension<typename T::type>,
                           S::num_sigma()>(
        Detail::get_field_covariance_offset<0, Fields...>(T::key), 0) = Detail::
        MeasurementStateHelper<FixedMeasurementVector, S>::sigma_point_deltas(
            get_field<T::key>(),
            Z.template block<Detail::StateVectorDimension<typename T::type>,
                             S::num_sigma()>(
                Detail::get_field_offset<0, Fields...>(T::key), 0));
  }

  template <typename S, typename T1, typename T2, typename... Tail>
  void calculate_field_deltas(const SigmaPointDistribution<S>& Z,
                              SigmaPointDeltas<S>& z_prime) const {
    calculate_field_deltas<S, T1>(Z, z_prime);
    calculate_field_deltas<S, T2, Tail...>(Z, z_prime);
  }
};

/*
This class provides a dynamic measurement vector, to be used when not all
measurements are available every time step.
*/
template <typename... Fields>
class DynamicMeasurementVector
    : public MeasurementVectorDynamicBaseType<typename Fields::type...> {
 public:
  using Base = MeasurementVectorDynamicBaseType<typename Fields::type...>;
  using Base::Base;
  using Base::operator=;

  /* Get maximum size of dynamic measurement vector. */
  static constexpr std::size_t max_size() {
    return Detail::get_composite_vector_dimension<typename Fields::type...>();
  }

  /* Get maximum size of dynamic measurement vector covariance. */
  static constexpr std::size_t max_covariance_size() {
    return Detail::get_covariance_dimension<typename Fields::type...>();
  }

  /* Aliases for types needed during filter iteration. */
  template <typename S>
  using SigmaPointDistribution =
      MatrixDynamic<DynamicMeasurementVector::max_size(), S::num_sigma()>;
  template <typename S>
  using SigmaPointDeltas =
      MatrixDynamic<DynamicMeasurementVector::max_covariance_size(),
                    S::num_sigma()>;
  using CovarianceMatrix =
      MatrixDynamic<max_covariance_size(), max_covariance_size()>;
  /*
  Diagonal sensor noise covers every declared field, and the innovation
  the fields that are set, with one entry for each component of a field's
  uncertainty: three, the rotation, for a quaternion. Without quaternion
  fields the innovation is this vector type. Experimental: quaternion
  measurement fields, and these types where they differ from the
  measurement type, may change in a minor release (see RELEASING.md).
  */
  using CovarianceVector = FixedMeasurementVector<Field<
      Fields::key, Detail::CovarianceFieldType<typename Fields::type>>...>;
  using InnovationVector = DynamicMeasurementVector<Field<
      Fields::key, Detail::CovarianceFieldType<typename Fields::type>>...>;
  // Dense sensor noise includes all declared fields, even when some are absent.
  using DenseCovarianceMatrix =
      Matrix<max_covariance_size(), max_covariance_size()>;
  using DenseNoiseRootMatrix =
      Eigen::Matrix<real_t, Eigen::Dynamic, max_covariance_size(),
                    Eigen::ColMajor, max_covariance_size(),
                    max_covariance_size()>;

  /* Test if an individual field is set. */
  template <int Key>
  bool is_field_set() const {
    static_assert(Detail::get_field_size<Fields...>(Key) !=
                      std::numeric_limits<std::size_t>::max(),
                  "Specified key not present in measurement vector");

    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(Key)>(field_offsets);

    return (offset != std::numeric_limits<std::size_t>::max());
  }

  /*
  The size of the set fields' uncertainty, and so of the innovation: the
  size, less one for each quaternion field set.
  */
  std::size_t covariance_size() const {
    return covariance_offset(static_cast<std::size_t>(Base::size()));
  }

  /* Functions for accessing individual fields. */
  template <int Key>
  typename Detail::FieldTypes<Key, Fields...>::type get_field() const {
    static_assert(Detail::get_field_size<Fields...>(Key) !=
                      std::numeric_limits<std::size_t>::max(),
                  "Specified key not present in measurement vector");

    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(Key)>(field_offsets);

    assert(offset != std::numeric_limits<std::size_t>::max() &&
           "Specified key not present in measurement vector");

    return Detail::convert_from_segment<
        typename Detail::FieldTypes<Key, Fields...>::type>(
        Base::template segment<Detail::get_field_size<Fields...>(Key)>(offset));
  }

  template <int Key, typename T>
  void set_field(T in) {
    static_assert(Detail::get_field_size<Fields...>(Key) !=
                      std::numeric_limits<std::size_t>::max(),
                  "Specified key not present in measurement vector");

    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(Key)>(field_offsets);

    /* Check if this field has already been set. If so, replace it. */
    if (offset < static_cast<std::size_t>(Base::size())) {
      Detail::assign_field(
          Base::template segment<Detail::get_field_size<Fields...>(Key)>(
              offset),
          in);
    } else {
      /*
      Otherwise, resize the measurement vector to fit it and store the
      order in which fields have been set.
      */
      std::size_t previous_size = Base::size();
      Base::conservativeResize(previous_size +
                               Detail::get_field_size<Fields...>(Key));

      /* Assign the value to the field. */
      Detail::assign_field(
          Base::template segment<Detail::get_field_size<Fields...>(Key)>(
              previous_size),
          in);

      /* Store the offset in field_offsets. */
      std::get<Detail::get_field_order<0, Fields...>(Key)>(field_offsets) =
          previous_size;
    }
  }

  /*
  Calculate the mean from a measurement sigma point distribution. Ensure
  that the returned object has the same field_offsets array so that its
  field accessors work.
  */
  template <typename S>
  DynamicMeasurementVector calculate_sigma_point_mean(
      const SigmaPointDistribution<S>& Z) const {
    DynamicMeasurementVector mean(Base::size());
    calculate_field_mean<S, Fields...>(Z, mean);

    mean.field_offsets = field_offsets;

    return mean;
  }

  /*
  Calculate the set of sigma point delta vectors; these are used for
  calculating the measurement covariance and the Kalman gain.
  The function isn't static; it uses the current measurement vector as the
  mean.
  */
  template <typename S>
  SigmaPointDeltas<S> calculate_sigma_point_deltas(
      const SigmaPointDistribution<S>& Z) const {
    SigmaPointDeltas<S> z_prime(covariance_size(), S::num_sigma());

    /* Calculate the delta vectors. */
    calculate_field_deltas<S, Fields...>(Z, z_prime);

    return z_prime;
  }

  /*
  Calculate the expected measurement covariance covariance as described in
  equation 68 of the Kraft papers.
  */
  template <typename S>
  CovarianceMatrix calculate_sigma_point_covariance(
      const SigmaPointDeltas<S>& z_prime) const {
    CovarianceMatrix cov(covariance_size(), covariance_size());

    /* Calculate the covariance using equation 64 from the Kraft paper. */
    cov.noalias() =
        (Parameters::Sigma_WCI<S> *
         z_prime.template rightCols<S::num_sigma() - 1>()) *
        z_prime.template rightCols<S::num_sigma() - 1>().transpose();
    cov.noalias() += (Parameters::Sigma_WC0<S> * z_prime.col(0)) *
                     z_prime.col(0).transpose();

    return cov;
  }

  /*
  Create a measurement sigma point distribution using the sigma points.
  */
  template <typename S, typename... U>
  SigmaPointDistribution<S> calculate_sigma_point_distribution(
      const typename S::SigmaPointDistribution& X, const U&... input) const {
    SigmaPointDistribution<S> Z(Base::size(), S::num_sigma());

    /* Reference the inputs once rather than copying them per sigma point. */
    const std::tuple<const U&...> inputs{input...};
    for (std::size_t i = 0; i < S::num_sigma(); i++) {
      DynamicMeasurementVector temp(Base::size());
      calculate_field_measurements<S, const std::tuple<const U&...>&,
                                   Fields...>(temp, X.col(i), inputs);
      Z.col(i) = temp;
    }

    return Z;
  }

  /*
  Returns a diagonal matrix containing the measurement covariance. This is
  used for the standard UKF.
  */
  CovarianceMatrix calculate_measurement_covariance(
      const CovarianceVector& measurement_covariance,
      const DynamicMeasurementVector& z_pred) const {
    CovarianceMatrix temp =
        CovarianceMatrix::Zero(covariance_size(), covariance_size());

    calculate_field_covariance<Fields...>(temp, measurement_covariance, z_pred);
    return temp;
  }

  /*
  Returns a diagonal matrix containing the measurement root covariance.
  This is used for the square-root UKF.
  */
  CovarianceMatrix calculate_measurement_root_covariance(
      const CovarianceVector& measurement_root_covariance,
      const DynamicMeasurementVector& z_pred) const {
    CovarianceMatrix temp =
        CovarianceMatrix::Zero(covariance_size(), covariance_size());

    calculate_field_root_covariance<Fields...>(
        temp, measurement_root_covariance, z_pred);
    return temp;
  }

  /* Dense noise is expressed in raw sensor coordinates and declared field
     order. */
  CovarianceMatrix calculate_measurement_covariance(
      const DenseCovarianceMatrix& measurement_covariance,
      const DynamicMeasurementVector& z_pred) const {
    DenseNoiseRootMatrix transform =
        DenseNoiseRootMatrix::Zero(covariance_size(), max_covariance_size());
    (calculate_field_noise_transform<Fields>(transform, z_pred), ...);
    return transform * measurement_covariance * transform.transpose();
  }

  /* Keep every factor column: selecting columns would discard noise variance.
   */
  DenseNoiseRootMatrix calculate_measurement_root_covariance(
      const DenseCovarianceMatrix& measurement_root_covariance,
      const DynamicMeasurementVector& z_pred) const {
    DenseNoiseRootMatrix transform =
        DenseNoiseRootMatrix::Zero(covariance_size(), max_covariance_size());
    (calculate_field_noise_transform<Fields>(transform, z_pred), ...);
    return transform * measurement_root_covariance;
  }

  /* Resolve Eigen expressions before conversion can match both noise types. */
  template <typename Derived>
  CovarianceMatrix calculate_measurement_covariance(
      const Eigen::MatrixBase<Derived>& measurement_covariance,
      const DynamicMeasurementVector& z_pred) const {
    if constexpr (Derived::ColsAtCompileTime == 1) {
      return calculate_measurement_covariance(
          CovarianceVector(measurement_covariance), z_pred);
    } else {
      if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (measurement_covariance.cols() == 1) {
          return calculate_measurement_covariance(
              CovarianceVector(measurement_covariance), z_pred);
        }
      }
      return calculate_measurement_covariance(
          DenseCovarianceMatrix(measurement_covariance), z_pred);
    }
  }

  /* Bounded dynamic columns accommodate both diagonal and dense root shapes. */
  template <typename Derived>
  CovarianceMatrix calculate_measurement_root_covariance(
      const Eigen::MatrixBase<Derived>& measurement_root_covariance,
      const DynamicMeasurementVector& z_pred) const {
    if constexpr (Derived::ColsAtCompileTime == 1) {
      return calculate_measurement_root_covariance(
          CovarianceVector(measurement_root_covariance), z_pred);
    } else {
      if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (measurement_root_covariance.cols() == 1) {
          return calculate_measurement_root_covariance(
              CovarianceVector(measurement_root_covariance), z_pred);
        }
      }
      return calculate_measurement_root_covariance(
          DenseCovarianceMatrix(measurement_root_covariance), z_pred);
    }
  }

  /*
  Whether a field is a quaternion, whose noise the correction maps to the
  observed innovation.
  */
  static constexpr bool has_quaternion_field =
      (std::is_same_v<typename Fields::type, Quaternion> || ...);

  /* Whether this reading holds a quaternion field; one may be left unset. */
  bool contains_quaternion_field() const {
    return (contains_field_of<Fields, Quaternion>() || ...);
  }

  /*
  The transform from the noise in the innovation covariance to the noise
  that the correction uses: J(u)^-1 at the observed innovation u for each
  quaternion field and the identity elsewhere, so that a noise covariance R
  becomes T R T^T and a noise factor L becomes T L. Called on the observed
  vector with the predicted one.
  */
  CovarianceMatrix calculate_correction_noise_transform(
      const DynamicMeasurementVector& z_pred) const {
    CovarianceMatrix transform =
        CovarianceMatrix::Identity(covariance_size(), covariance_size());
    (calculate_field_correction_transform<Fields>(transform, z_pred), ...);
    return transform;
  }

  /*
  Return the innovation using the supplied measurement vector. Its fields
  are in the order in which they were set, as the uncertainty is.
  */
  InnovationVector calculate_innovation(
      const DynamicMeasurementVector& z) const {
    InnovationVector temp(covariance_size());
    (calculate_field_innovation<Fields>(z, temp), ...);
    return temp;
  }

 private:
  /*
  This vector keeps track of which fields have been set in the measurement
  vector, and the offset within the measurement vector of each field. This
  allows any combination of measurements to be supplied in any order and
  they will be handled correctly.
  */
  std::array<std::size_t, sizeof...(Fields)> field_offsets =
      Detail::create_array<sizeof...(Fields)>(
          std::numeric_limits<std::size_t>::max());

  /* The innovation vector type sets its field offsets from this one. */
  template <typename... OtherFields>
  friend class DynamicMeasurementVector;

  /*
  The offset in the uncertainty of the field stored at the given offset:
  the offset less one for each quaternion field stored before it.
  */
  std::size_t covariance_offset(std::size_t offset) const {
    if constexpr ((std::is_same<typename Fields::type, Quaternion>::value ||
                   ...)) {
      std::size_t quaternions_before = 0;
      ((quaternions_before +=
        std::is_same<typename Fields::type, Quaternion>::value &&
        std::get<Detail::get_field_order<0, Fields...>(Fields::key)>(
            field_offsets) < offset),
       ...);
      return offset - quaternions_before;
    } else {
      return offset;
    }
  }

  /*
  This function is intended to be specialised by the user for each field in
  the measurement vector, and allows the user to specify how a particular
  state vector is transformed into a measurement vector.

  Template parameters are a StateVector type, a field key and a parameter
  pack of input vectors.
  */
  template <typename S, int Key, typename... U>
  static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement(
      const S& state, const U&... input);

  template <typename S, int Key, typename U, size_t... I>
  static typename Detail::FieldTypes<Key, Fields...>::type
  expected_measurement_helper(const S& state, U&& input,
                              std::index_sequence<I...>) {
    return expected_measurement<S, Key>(state, std::get<I>(input)...);
  }

  /*
  These functions build the measurement estimate from the expected
  measurement of each individual field.
  */
  template <typename S, typename U, typename T>
  void calculate_field_measurements(DynamicMeasurementVector& expected,
                                    const S& state, U&& input) const {
    constexpr std::size_t len =
        std::tuple_size<typename std::remove_reference<U>::type>::value;

    /*
    If this field has been set, then generate an expected measurement for
    it. Otherwise, do nothing.
    */
    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
    if (offset != std::numeric_limits<std::size_t>::max()) {
      Detail::assign_field(
          expected
              .template segment<Detail::StateVectorDimension<typename T::type>>(
                  offset),
          expected_measurement_helper<S, T::key, U>(
              state, std::forward<U>(input), std::make_index_sequence<len>()));
    } else {
      return;
    }
  }

  template <typename S, typename U, typename T1, typename T2, typename... Tail>
  void calculate_field_measurements(DynamicMeasurementVector& expected,
                                    const S& state, U&& input) const {
    calculate_field_measurements<S, U, T1>(expected, state,
                                           std::forward<U>(input));
    calculate_field_measurements<S, U, T2, Tail...>(expected, state,
                                                    std::forward<U>(input));
  }

  /*
  These functions build the measurement covariance from all populated
  fields.
  */
  /* Whether field T is set and of type Type. */
  template <typename T, typename Type>
  bool contains_field_of() const {
    if constexpr (std::is_same_v<typename T::type, Type>) {
      return std::get<Detail::get_field_order<0, Fields...>(T::key)>(
                 field_offsets) != std::numeric_limits<std::size_t>::max();
    } else {
      return false;
    }
  }

  template <typename T>
  void calculate_field_correction_transform(
      CovarianceMatrix& transform,
      const DynamicMeasurementVector& z_pred) const {
    if constexpr (std::is_same_v<typename T::type, Quaternion>) {
      const std::size_t offset =
          std::get<Detail::get_field_order<0, Fields...>(T::key)>(
              field_offsets);
      if (offset != std::numeric_limits<std::size_t>::max()) {
        transform.template block<3, 3>(covariance_offset(offset),
                                       covariance_offset(offset)) =
            Detail::MeasurementStateHelper<DynamicMeasurementVector>::
                correction_noise_jacobian(z_pred.template get_field<T::key>(),
                                          get_field<T::key>());
      }
    }
  }

  template <typename T>
  void calculate_field_noise_transform(
      DenseNoiseRootMatrix& transform,
      const DynamicMeasurementVector& z_pred) const {
    const std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
    if (offset != std::numeric_limits<std::size_t>::max()) {
      constexpr std::size_t noise_offset =
          Detail::get_field_covariance_offset<0, Fields...>(T::key);
      transform.template block<Detail::CovarianceDimension<typename T::type>,
                               Detail::CovarianceDimension<typename T::type>>(
          covariance_offset(offset), noise_offset) = Detail::
          MeasurementStateHelper<DynamicMeasurementVector>::noise_jacobian(
              z_pred.template get_field<T::key>(), get_field<T::key>());
    }
  }

  template <typename T>
  void calculate_field_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_covariance,
      const DynamicMeasurementVector& z_pred) const {
    /*
    If this field has been set, then fill the measurement covariance for
    it. Otherwise, do nothing.
    */
    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
    if (offset != std::numeric_limits<std::size_t>::max()) {
      P.template block<Detail::CovarianceDimension<typename T::type>,
                       Detail::CovarianceDimension<typename T::type>>(
          covariance_offset(offset), covariance_offset(offset)) = Detail::
          MeasurementStateHelper<DynamicMeasurementVector>::field_covariance(
              measurement_covariance.template get_field<T::key>(),
              z_pred.get_field<T::key>(), get_field<T::key>());
    } else {
      return;
    }
  }

  template <typename T1, typename T2, typename... Tail>
  void calculate_field_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_covariance,
      const DynamicMeasurementVector& z_pred) const {
    calculate_field_covariance<T1>(P, measurement_covariance, z_pred);
    calculate_field_covariance<T2, Tail...>(P, measurement_covariance, z_pred);
  }

  template <typename T>
  void calculate_field_root_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_root_covariance,
      const DynamicMeasurementVector& z_pred) const {
    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
    if (offset != std::numeric_limits<std::size_t>::max()) {
      P.template block<Detail::CovarianceDimension<typename T::type>,
                       Detail::CovarianceDimension<typename T::type>>(
          covariance_offset(offset), covariance_offset(offset)) =
          Detail::MeasurementStateHelper<DynamicMeasurementVector>::
              field_root_covariance(
                  measurement_root_covariance.template get_field<T::key>(),
                  z_pred.get_field<T::key>(), get_field<T::key>());
    } else {
      return;
    }
  }

  template <typename T1, typename T2, typename... Tail>
  void calculate_field_root_covariance(
      CovarianceMatrix& P, const CovarianceVector& measurement_root_covariance,
      const DynamicMeasurementVector& z_pred) const {
    calculate_field_root_covariance<T1>(P, measurement_root_covariance, z_pred);
    calculate_field_root_covariance<T2, Tail...>(P, measurement_root_covariance,
                                                 z_pred);
  }

  /*
  Private functions for calculating the mean of each field in a sigma point
  distribution.
  */
  template <typename S, typename T>
  void calculate_field_mean(const SigmaPointDistribution<S>& Z,
                            DynamicMeasurementVector& mean) const {
    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
    if (offset != std::numeric_limits<std::size_t>::max()) {
      mean.template segment<Detail::StateVectorDimension<typename T::type>>(
          offset) = Detail::MeasurementStateHelper<DynamicMeasurementVector,
                                                   S>::
          template sigma_point_mean<typename T::type>(
              Z.template block<Detail::StateVectorDimension<typename T::type>,
                               S::num_sigma()>(offset, 0));
    } else {
      return;
    }
  }

  template <typename S, typename T1, typename T2, typename... Tail>
  void calculate_field_mean(const SigmaPointDistribution<S>& Z,
                            DynamicMeasurementVector& mean) const {
    calculate_field_mean<S, T1>(Z, mean);
    calculate_field_mean<S, T2, Tail...>(Z, mean);
  }

  /* These functions build the innovation from all populated fields. */
  template <typename T>
  void calculate_field_innovation(const DynamicMeasurementVector& z,
                                  InnovationVector& innovation) const {
    constexpr std::size_t order = Detail::get_field_order<0, Fields...>(T::key);
    const std::size_t offset = std::get<order>(field_offsets);
    /* A field contributes only when both vectors hold it. */
    if (offset != std::numeric_limits<std::size_t>::max() &&
        std::get<order>(z.field_offsets) !=
            std::numeric_limits<std::size_t>::max()) {
      const std::size_t innovation_offset = covariance_offset(offset);
      std::get<order>(innovation.field_offsets) = innovation_offset;
      Detail::assign_field(
          innovation
              .template segment<Detail::CovarianceDimension<typename T::type>>(
                  innovation_offset),
          Detail::MeasurementStateHelper<DynamicMeasurementVector>::
              measurement_delta(z.template get_field<T::key>(),
                                get_field<T::key>()));
    }
  }

  /*
  Private functions for calculating the delta vectors between each sigma
  point and the mean.
  */
  template <typename S, typename T>
  void calculate_field_deltas(const SigmaPointDistribution<S>& Z,
                              SigmaPointDeltas<S>& z_prime) const {
    std::size_t offset =
        std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
    if (offset != std::numeric_limits<std::size_t>::max()) {
      z_prime.template block<Detail::CovarianceDimension<typename T::type>,
                             S::num_sigma()>(
          covariance_offset(offset),
          0) = Detail::MeasurementStateHelper<DynamicMeasurementVector, S>::
          sigma_point_deltas(
              get_field<T::key>(),
              Z.template block<Detail::StateVectorDimension<typename T::type>,
                               S::num_sigma()>(offset, 0));
    } else {
      return;
    }
  }

  template <typename S, typename T1, typename T2, typename... Tail>
  void calculate_field_deltas(const SigmaPointDistribution<S>& Z,
                              SigmaPointDeltas<S>& z_prime) const {
    calculate_field_deltas<S, T1>(Z, z_prime);
    calculate_field_deltas<S, T2, Tail...>(Z, z_prime);
  }
};

}  // namespace UKF

#endif
