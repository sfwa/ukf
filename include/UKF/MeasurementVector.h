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

#include <array>
#include <limits>
#include <tuple>
#include <cassert>
#include <cstddef>
#include <utility>
#include <Eigen/Core>
#include "UKF/Types.h"
#include "UKF/StateVector.h"

namespace UKF {

    namespace Detail {

    /* Calculate the smallest rotation vector which transforms v1 into v2. */
    template <typename T>
    inline Vector<3> calculate_rotation_vector(const Vector<3>& v2, const Vector<3>& v1) {
        Vector<3> axis = Vector<3>::Zero();
        real_t norm = std::sqrt(v1.squaredNorm() * v2.squaredNorm());
        real_t q_w = norm + v1.dot(v2);

        /* Check whether the vectors are antiparallel or too small. */
        if(q_w > std::numeric_limits<real_t>::epsilon()) {
            axis = v1.cross(v2);
        } else if(norm > std::numeric_limits<real_t>::epsilon()) {
            axis = Vector<3>(-v1(2), v1(1), v1(0));
        }

        return Parameters::MRP_F<T> * axis /
            (std::abs(Parameters::MRP_A<T> + q_w) > std::numeric_limits<real_t>::epsilon() ?
            Parameters::MRP_A<T> + q_w : std::numeric_limits<real_t>::epsilon());
    }

    /*
    Calculate the transformation matrix to transform a covariance matrix of
    field vector noise for the field vector v2 into a covariance matrix for the
    rotation vector corresponding to the transformation between the field
    vector v1 and v2.

    This is done by using an expression for the exact Jacobian of the rotation
    vector, as a function of the field vector.

    The expression for the exact Jacobian was derived using the following
    MATLAB commands (with the Symbolic Toolbox):

    v1_ = sym('v1_', [3 1], 'real');
    v2_ = sym('v2_', [3 1], 'real');
    f = sym('f', 1, 'real');
    a = sym('a', 1, 'real');
    Ra = cross(v1_, v2_) * (f/(a + sqrt(dot(v1_, v1_) * dot(v2_, v2_)) + dot(v1_, v2_)));
    J_Ra = jacobian(Ra, v2_);
    pretty(J_Ra);
    */
    template <typename T>
    inline Matrix<3, 3> calculate_rotation_vector_jacobian(const Vector<3>& v2, const Vector<3>& v1) {
        Matrix<3, 3> j;
        real_t c8 = std::sqrt(v1.squaredNorm() * v2.squaredNorm());
        real_t q_w = c8 + v1.dot(v2);
        real_t c9_8 = real_t(0.0);

        /* Check whether the vectors are antiparallel or too small. */
        Vector<3> c123 = Vector<3>::Zero();
        if(c8 > std::numeric_limits<real_t>::epsilon()) {
            c9_8 = v1.squaredNorm() / c8;

            if(q_w > std::numeric_limits<real_t>::epsilon()) {
                c123 = v1.cross(v2);
            } else {
                c123 = Vector<3>(-v1(2), v1(1), v1(0));
            }
            c123(1) = -c123(1);
        }

        real_t c4 = Parameters::MRP_A<T> + q_w;
        if(std::abs(c4) < std::numeric_limits<real_t>::epsilon()) {
            c4 = std::numeric_limits<real_t>::epsilon();
        }
        real_t c4_2 = c4*c4;

        Vector<3> c765 = v1 + v2*c9_8;

        j <<            -c123(0)*c765(0)/c4_2, -v1(2)/c4 - c123(0)*c765(1)/c4_2,  v1(1)/c4 - c123(0)*c765(2)/c4_2,
              v1(2)/c4 + c123(1)*c765(0)/c4_2,             c123(1)*c765(1)/c4_2, -v1(0)/c4 + c123(1)*c765(2)/c4_2,
             -v1(1)/c4 - c123(2)*c765(0)/c4_2,  v1(0)/c4 - c123(2)*c765(1)/c4_2,            -c123(2)*c765(2)/c4_2;
        j *= Parameters::MRP_F<T>;

        return j;
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
        /*
        Functions for calculating the covariance of a field in the
        measurement vector.
        */
        template <typename T>
        static Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>> field_covariance(
                const T& p, const T& z_pred, const T& z) {
            Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>> temp =
                Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>::Zero();
            temp.diagonal() << p;
            return temp;
        }

        static Matrix<3, 3> field_covariance(const FieldVector& p, const FieldVector& z_pred, const FieldVector& z) {
            Matrix<3, 3> T = Detail::calculate_rotation_vector_jacobian<M>(z, z_pred);

            /*
            To calculate the covariance, pre-multiply by the transformation
            matrix and then post-multiply by the transformation matrix
            transpose.
            */
            return T * Eigen::DiagonalMatrix<real_t, 3>(p) * T.transpose();
        }

        template <typename T>
        static Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>> field_root_covariance(
                const T& p, const T& z_pred, const T& z) {
            Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>> temp =
                Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>::Zero();
            temp.diagonal() << p;
            return temp;
        }

        static Matrix<3, 3> field_root_covariance(
                const FieldVector& p, const FieldVector& z_pred, const FieldVector& z) {
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
        distribution. Note that sigma_point_mean takes a dummy argument so
        that the overrides work properly.
        */
        template <typename T>
        static T sigma_point_mean(
                const Matrix<Detail::StateVectorDimension<T>, S::num_sigma()>& sigma, const T& field) {
            return Parameters::Sigma_WMI<S>*sigma.template block<Detail::StateVectorDimension<T>, S::num_sigma()-1>(
                0, 1).rowwise().sum() + Parameters::Sigma_WM0<S>*sigma.col(0);
        }

        static real_t sigma_point_mean(const Matrix<1, S::num_sigma()>& sigma, const real_t& field) {
            return Parameters::Sigma_WMI<S>*sigma.template segment<S::num_sigma()-1>(1).sum()
                + Parameters::Sigma_WM0<S>*sigma(0);
        }

        /*
        Calculate the field vector mean by first calculating the set of
        rotation vectors which transforms the central sigma point into the
        set of sigma points. Then, calculate the mean of these rotations and
        apply it to the central sigma point.
        */
        static FieldVector sigma_point_mean(const Matrix<3, S::num_sigma()>& sigma, const FieldVector& field) {
            Vector<3> temp = Vector<3>::Zero();

            for(std::size_t i = 0; i < S::num_sigma()-1; i++) {
                temp += Parameters::Sigma_WMI<S>*Detail::calculate_rotation_vector<M>(sigma.col(i+1), sigma.col(0));
            }

            return Detail::rotation_vector_to_quaternion<M>(temp) * sigma.col(0);
        }

        /*
        Functions for calculating the difference between two fields in a
        measurement vector, used when calculating the innovation.
        */
        template <typename T>
        static T measurement_delta(const T& z, const T& z_pred) {
            return z - z_pred;
        }

        static real_t measurement_delta(const real_t& z, const real_t& z_pred) {
            return z - z_pred;
        }

        static FieldVector measurement_delta(const FieldVector& z, const FieldVector& z_pred) {
            return Detail::calculate_rotation_vector<M>(z, z_pred);
        }

        /*
        Functions for calculating the difference between each point in a
        sigma point distribution and the mean.
        */
        template <typename T>
        static Matrix<Detail::CovarianceDimension<T>, S::num_sigma()> sigma_point_deltas(
                const T& mean, const Matrix<Detail::StateVectorDimension<T>, S::num_sigma()>& Z) {
            return Z.colwise() - mean;
        }

        static Matrix<1, S::num_sigma()> sigma_point_deltas(real_t mean, const Matrix<1, S::num_sigma()>& Z) {
            return Z.array() - mean;
        }

        static Matrix<3, S::num_sigma()> sigma_point_deltas(
                const FieldVector& mean, const Matrix<3, S::num_sigma()>& Z) {
            Matrix<3, S::num_sigma()> temp;

            for(std::size_t i = 0; i < S::num_sigma(); i++) {
                temp.col(i) = Detail::calculate_rotation_vector<M>(Z.col(i), mean);
            }

            return temp;
        }
    };

    }

/* Alias for the Eigen type from which FixedMeasurementVector inherits. */
template <typename... Fields>
using MeasurementVectorFixedBaseType = Vector<Detail::get_composite_vector_dimension<Fields...>()>;

/* Alias for the Eigen type from which DynamicMeasurementVector inherits. */
template <typename... Fields>
using MeasurementVectorDynamicBaseType = VectorDynamic<Detail::get_composite_vector_dimension<Fields...>()>;

/*
This class provides a fixed measurement vector, to be used when the same
measurement are available every time step.
*/
template <typename... Fields>
class FixedMeasurementVector : public MeasurementVectorFixedBaseType<typename Fields::type...> {
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
    using SigmaPointDistribution = Matrix<FixedMeasurementVector::size(), S::num_sigma()>;
    template <typename S>
    using SigmaPointDeltas = Matrix<FixedMeasurementVector::covariance_size(), S::num_sigma()>;
    using CovarianceMatrix = Matrix<covariance_size(), covariance_size()>;
    using CovarianceVector = FixedMeasurementVector<Fields...>;

    /*
    Measurement noise covariance. This is defined by the user and can be
    adjusted between iterations. This is used for the standard UKF.
    */
    static CovarianceVector measurement_covariance;

    /*
    Measurement noise root covariance. This is used for the square-root UKF.
    */
    static CovarianceVector measurement_root_covariance;

    /* Functions for accessing individual fields. */
    template <int Key>
    typename Detail::FieldTypes<Key, Fields...>::type get_field() const {
        static_assert(Detail::get_field_offset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in state vector");
        return Detail::convert_from_segment<typename Detail::FieldTypes<Key, Fields...>::type>(
            Base::template segment<Detail::get_field_size<Fields...>(Key)>(
                Detail::get_field_offset<0, Fields...>(Key)));
    }

    template <int Key, typename T>
    void set_field(T in) {
        static_assert(Detail::get_field_offset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in state vector");
        Base::template segment<Detail::get_field_size<Fields...>(Key)>(
            Detail::get_field_offset<0, Fields...>(Key)) << in;
    }

    /* Calculate the mean from a measurement sigma point distribution. */
    template <typename S>
    FixedMeasurementVector calculate_sigma_point_mean(const SigmaPointDistribution<S>& Z) const {
        FixedMeasurementVector mean;
        calculate_field_mean<S, Fields...>(Z, mean);

        return mean;
    }

    /*
    Calculate the set of sigma point delta vectors; these are used for
    calculating the measurement covariance and the Kalman gain.
    */
    template <typename S>
    SigmaPointDeltas<S> calculate_sigma_point_deltas(const SigmaPointDistribution<S>& Z) const {
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
    CovarianceMatrix calculate_sigma_point_covariance(const SigmaPointDeltas<S>& z_prime) const {
        CovarianceMatrix cov;

        /* Calculate the covariance using equation 64 from the Kraft paper. */
        cov = CovarianceMatrix::Zero();
        for(std::size_t i = 1; i < S::num_sigma(); i++) {
            cov.noalias() += Parameters::Sigma_WCI<S> * (z_prime.col(i) * z_prime.col(i).transpose());
        }
        cov.noalias() += Parameters::Sigma_WC0<S> * (z_prime.col(0) * z_prime.col(0).transpose());

        return cov;
    }

    /*
    Create a measurement sigma point distribution using the sigma points.
    */
    template <typename S, typename... U>
    SigmaPointDistribution<S> calculate_sigma_point_distribution(
            const typename S::SigmaPointDistribution& X, const U&... input) const {
        SigmaPointDistribution<S> Z(size(), S::num_sigma());

        for(std::size_t i = 0; i < S::num_sigma(); i++) {
            FixedMeasurementVector temp;
            calculate_field_measurements<S, std::tuple<U...>, Fields...>(temp, X.col(i), std::make_tuple(input...));
            Z.col(i) = temp;
        }

        return Z;
    }

    /*
    Returns a matrix containing the measurement covariance. Used for the
    standard UKF.
    */
    CovarianceMatrix calculate_measurement_covariance(const FixedMeasurementVector& z_pred) const {
        CovarianceMatrix temp = CovarianceMatrix::Zero();

        calculate_field_covariance<Fields...>(temp, z_pred);
        return temp;
    }

    /*
    Returns a matrix containing the measurement root covariance. Used for the
    square-root UKF.
    */
    CovarianceMatrix calculate_measurement_root_covariance(const FixedMeasurementVector& z_pred) const {
        CovarianceMatrix temp = CovarianceMatrix::Zero();

        calculate_field_root_covariance<Fields...>(temp, z_pred);
        return temp;
    }

    /* Return the innovation using the supplied measurement vector. */
    FixedMeasurementVector calculate_innovation(const FixedMeasurementVector& z) const {
        FixedMeasurementVector temp;

        temp = z;
        calculate_field_innovation<Fields...>(temp);
        return temp;
    }

private:
    /*
    This function is intended to be specialised by the user for each field in
    the measurement vector, and allows the user to specify how a particular
    state vector is transformed into a measurement vector.

    Template parameters are a StateVector type, a field key and a parameter
    pack of input vectors.
    */
    template <typename S, int Key, typename... U>
    static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement(const S& state, const U&... input);

    template <typename S, int Key, typename U, size_t... I>
    static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement_helper(
            const S& state, U&& input, std::index_sequence<I...>) {
        return expected_measurement<S, Key>(state, std::get<I>(input)...);
    }

    /*
    These functions build the measurement estimate from the expected
    measurement of each individual field.
    */
    template <typename S, typename U, typename T>
    static void calculate_field_measurements(FixedMeasurementVector& expected, const S& state, U&& input) {
        constexpr std::size_t len = std::tuple_size<typename std::remove_reference<U>::type>::value;

        expected.template segment<Detail::StateVectorDimension<typename T::type>>(
            Detail::get_field_offset<0, Fields...>(T::key)) << expected_measurement_helper<S, T::key, U>(
                state, std::forward<U>(input), std::make_index_sequence<len>());
    }

    template <typename S, typename U, typename T1, typename T2, typename... Tail>
    static void calculate_field_measurements(FixedMeasurementVector& expected, const S& state, U&& input) {
        calculate_field_measurements<S, U, T1>(expected, state, std::forward<U>(input));
        calculate_field_measurements<S, U, T2, Tail...>(expected, state, std::forward<U>(input));
    }

    /*
    These functions build the measurement covariance from all fields in the
    measurement vector.
    */
    template <typename T>
    void calculate_field_covariance(CovarianceMatrix& P, const FixedMeasurementVector& z_pred) const {
        P.template block<Detail::CovarianceDimension<typename T::type>, Detail::CovarianceDimension<typename T::type>>(
            Detail::get_field_covariance_offset<0, Fields...>(T::key),
            Detail::get_field_covariance_offset<0, Fields...>(T::key)) =
                Detail::MeasurementStateHelper<FixedMeasurementVector>::field_covariance(
                    measurement_covariance.template get_field<T::key>(),
                    z_pred.get_field<T::key>(), get_field<T::key>());
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_covariance(CovarianceMatrix& P, const FixedMeasurementVector& z_pred) const {
        calculate_field_covariance<T1>(P, z_pred);
        calculate_field_covariance<T2, Tail...>(P, z_pred);
    }

    /*
    These functions build the measurement root covariance from all fields in
    the measurement vector.
    */
    template <typename T>
    void calculate_field_root_covariance(CovarianceMatrix& P, const FixedMeasurementVector& z_pred) const {
        P.template block<Detail::CovarianceDimension<typename T::type>, Detail::CovarianceDimension<typename T::type>>(
            Detail::get_field_covariance_offset<0, Fields...>(T::key),
            Detail::get_field_covariance_offset<0, Fields...>(T::key)) = 
                Detail::MeasurementStateHelper<FixedMeasurementVector>::field_root_covariance(
                    measurement_root_covariance.template get_field<T::key>(),
                    z_pred.get_field<T::key>(), get_field<T::key>());
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_root_covariance(CovarianceMatrix& P, const FixedMeasurementVector& z_pred) const {
        calculate_field_root_covariance<T1>(P, z_pred);
        calculate_field_root_covariance<T2, Tail...>(P, z_pred);
    }

    /*
    Private functions for calculating the mean of each field in a sigma point
    distribution.
    */
    template <typename S, typename T>
    void calculate_field_mean(const SigmaPointDistribution<S>& Z, FixedMeasurementVector& mean) const {
        mean.template segment<Detail::StateVectorDimension<typename T::type>>(
            Detail::get_field_offset<0, Fields...>(T::key)) <<
                Detail::MeasurementStateHelper<FixedMeasurementVector, S>::sigma_point_mean(
                    Z.template block<Detail::StateVectorDimension<typename T::type>, S::num_sigma()>(
                    Detail::get_field_offset<0, Fields...>(T::key), 0), typename T::type());
    }

    template <typename S, typename T1, typename T2, typename... Tail>
    void calculate_field_mean(const SigmaPointDistribution<S>& Z, FixedMeasurementVector& mean) const {
        calculate_field_mean<S, T1>(Z, mean);
        calculate_field_mean<S, T2, Tail...>(Z, mean);
    }

    /* These functions build the innovation from all populated fields. */
    template <typename T>
    void calculate_field_innovation(FixedMeasurementVector& z) const {
        z.set_field<T::key>(Detail::MeasurementStateHelper<FixedMeasurementVector>::measurement_delta(
            z.get_field<T::key>(), get_field<T::key>()));
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_innovation(FixedMeasurementVector& z) const {
        calculate_field_innovation<T1>(z);
        calculate_field_innovation<T2, Tail...>(z);
    }

    /*
    Private functions for calculating the delta vectors between each sigma
    point and the mean.
    */
    template <typename S, typename T>
    void calculate_field_deltas(const SigmaPointDistribution<S>& Z, SigmaPointDeltas<S>& z_prime) const {
        z_prime.template block<Detail::CovarianceDimension<typename T::type>, S::num_sigma()>(
            Detail::get_field_covariance_offset<0, Fields...>(T::key), 0) =
                Detail::MeasurementStateHelper<FixedMeasurementVector, S>::sigma_point_deltas(
                    get_field<T::key>(), Z.template block<Detail::StateVectorDimension<typename T::type>,
                    S::num_sigma()>(Detail::get_field_offset<0, Fields...>(T::key), 0));
    }

    template <typename S, typename T1, typename T2, typename... Tail>
    void calculate_field_deltas(const SigmaPointDistribution<S>& Z, SigmaPointDeltas<S>& z_prime) const {
        calculate_field_deltas<S, T1>(Z, z_prime);
        calculate_field_deltas<S, T2, Tail...>(Z, z_prime);
    }
};

/*
This class provides a dynamic measurement vector, to be used when not all
measurements are available every time step.
*/
template <typename... Fields>
class DynamicMeasurementVector : public MeasurementVectorDynamicBaseType<typename Fields::type...> {
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
    using SigmaPointDistribution = MatrixDynamic<DynamicMeasurementVector::max_size(), S::num_sigma()>;
    template <typename S>
    using SigmaPointDeltas = MatrixDynamic<DynamicMeasurementVector::max_covariance_size(), S::num_sigma()>;
    using CovarianceMatrix = MatrixDynamic<max_covariance_size(), max_covariance_size()>;
    using CovarianceVector = FixedMeasurementVector<Fields...>;

    /*
    Measurement noise covariance. This is used for the standard UKF.
    */
    static CovarianceVector measurement_covariance;

    /*
    Measurement noise root covariance. This is used for the square-root UKF.
    */
    static CovarianceVector measurement_root_covariance;

    /* Functions for accessing individual fields. */
    template <int Key>
    typename Detail::FieldTypes<Key, Fields...>::type get_field() const {
        static_assert(Detail::get_field_size<Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");

        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(Key)>(field_offsets);

        assert(offset != std::numeric_limits<std::size_t>::max() &&
            "Specified key not present in measurement vector");

        return Detail::convert_from_segment<typename Detail::FieldTypes<Key, Fields...>::type>(
            Base::template segment<Detail::get_field_size<Fields...>(Key)>(offset));
    }

    template <int Key, typename T>
    void set_field(T in) {
        static_assert(Detail::get_field_size<Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");

        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(Key)>(field_offsets);

        /* Check if this field has already been set. If so, replace it. */
        if(offset < Base::size()) {
            Base::template segment<Detail::get_field_size<Fields...>(Key)>(offset) << in;
        } else {
            /*
            Otherwise, resize the measurement vector to fit it and store the
            order in which fields have been set.
            */
            std::size_t previous_size = Base::size();
            Base::conservativeResize(previous_size + Detail::get_field_size<Fields...>(Key));

            /* Assign the value to the field. */
            Base::template segment<Detail::get_field_size<Fields...>(Key)>(previous_size) << in;

            /* Store the offset in field_offsets. */
            std::get<Detail::get_field_order<0, Fields...>(Key)>(field_offsets) = previous_size;
        }
    }

    /*
    Calculate the mean from a measurement sigma point distribution. Ensure
    that the returned object has the same field_offsets array so that its
    field accessors work.
    */
    template <typename S>
    DynamicMeasurementVector calculate_sigma_point_mean(const SigmaPointDistribution<S>& Z) const {
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
    SigmaPointDeltas<S> calculate_sigma_point_deltas(const SigmaPointDistribution<S>& Z) const {
        SigmaPointDeltas<S> z_prime(Base::size(), S::num_sigma());

        /* Calculate the delta vectors. */
        calculate_field_deltas<S, Fields...>(Z, z_prime);

        return z_prime;
    }

    /*
    Calculate the expected measurement covariance covariance as described in
    equation 68 of the Kraft papers.
    */
    template <typename S>
    CovarianceMatrix calculate_sigma_point_covariance(const SigmaPointDeltas<S>& z_prime) const {
        CovarianceMatrix cov(Base::size(), Base::size());

        /* Calculate the covariance using equation 64 from the Kraft paper. */
        cov = CovarianceMatrix::Zero(Base::size(), Base::size());
        for(std::size_t i = 1; i < S::num_sigma(); i++) {
            cov.noalias() += Parameters::Sigma_WCI<S> * (z_prime.col(i) * z_prime.col(i).transpose());
        }
        cov.noalias() += Parameters::Sigma_WC0<S> * (z_prime.col(0) * z_prime.col(0).transpose());

        return cov;
    }

    /*
    Create a measurement sigma point distribution using the sigma points.
    */
    template <typename S, typename... U>
    SigmaPointDistribution<S> calculate_sigma_point_distribution(
            const typename S::SigmaPointDistribution& X, const U&... input) const {
        SigmaPointDistribution<S> Z(Base::size(), S::num_sigma());

        for(std::size_t i = 0; i < S::num_sigma(); i++) {
            DynamicMeasurementVector temp(Base::size());
            calculate_field_measurements<S, std::tuple<U...>, Fields...>(temp, X.col(i), std::make_tuple(input...));
            Z.col(i) = temp;
        }

        return Z;
    }

    /*
    Returns a diagonal matrix containing the measurement covariance. This is
    used for the standard UKF.
    */
    CovarianceMatrix calculate_measurement_covariance(const DynamicMeasurementVector& z_pred) const {
        CovarianceMatrix temp = CovarianceMatrix::Zero(Base::size(), Base::size());

        calculate_field_covariance<Fields...>(temp, z_pred);
        return temp;
    }

    /*
    Returns a diagonal matrix containing the measurement root covariance.
    This is used for the square-root UKF.
    */
    CovarianceMatrix calculate_measurement_root_covariance(const DynamicMeasurementVector& z_pred) const {
        CovarianceMatrix temp = CovarianceMatrix::Zero(Base::size(), Base::size());

        calculate_field_root_covariance<Fields...>(temp, z_pred);
        return temp;
    }

    /* Return the innovation using the supplied measurement vector. */
    DynamicMeasurementVector calculate_innovation(const DynamicMeasurementVector& z) const {
        DynamicMeasurementVector temp(Base::size());

        temp = z;
        calculate_field_innovation<Fields...>(temp);
        return temp;
    }

private:
    /*
    This vector keeps track of which fields have been set in the measurement
    vector, and the offset within the measurement vector of each field. This
    allows any combination of measurements to be supplied in any order and
    they will be handled correctly.
    */
    std::array<std::size_t, sizeof...(Fields)> field_offsets = Detail::create_array<sizeof...(Fields)>(
        std::numeric_limits<std::size_t>::max());

    /*
    This function is intended to be specialised by the user for each field in
    the measurement vector, and allows the user to specify how a particular
    state vector is transformed into a measurement vector.

    Template parameters are a StateVector type, a field key and a parameter
    pack of input vectors.
    */
    template <typename S, int Key, typename... U>
    static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement(const S& state, const U&... input);

    template <typename S, int Key, typename U, size_t... I>
    static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement_helper(
            const S& state, U&& input, std::index_sequence<I...>) {
        return expected_measurement<S, Key>(state, std::get<I>(input)...);
    }

    /*
    These functions build the measurement estimate from the expected
    measurement of each individual field.
    */
    template <typename S, typename U, typename T>
    void calculate_field_measurements(DynamicMeasurementVector& expected, const S& state, U&& input) const {
        constexpr std::size_t len = std::tuple_size<typename std::remove_reference<U>::type>::value;

        /*
        If this field has been set, then generate an expected measurement for
        it. Otherwise, do nothing.
        */
        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
        if(offset != std::numeric_limits<std::size_t>::max()) {
            expected.template segment<Detail::StateVectorDimension<typename T::type>>(offset) <<
                expected_measurement_helper<S, T::key, U>(state, std::forward<U>(input),
                    std::make_index_sequence<len>());
        } else {
            return;
        }
    }

    template <typename S, typename U, typename T1, typename T2, typename... Tail>
    void calculate_field_measurements(DynamicMeasurementVector& expected, const S& state, U&& input) const {
        calculate_field_measurements<S, U, T1>(expected, state, std::forward<U>(input));
        calculate_field_measurements<S, U, T2, Tail...>(expected, state, std::forward<U>(input));
    }

    /*
    These functions build the measurement covariance from all populated
    fields.
    */
    template <typename T>
    void calculate_field_covariance(CovarianceMatrix& P, const DynamicMeasurementVector& z_pred) const {
        /*
        If this field has been set, then fill the measurement covariance for
        it. Otherwise, do nothing.
        */
        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
        if(offset != std::numeric_limits<std::size_t>::max()) {
            P.template block<Detail::CovarianceDimension<typename T::type>,
                Detail::CovarianceDimension<typename T::type>>(offset, offset) =
                    Detail::MeasurementStateHelper<DynamicMeasurementVector>::field_covariance(
                        measurement_covariance.template get_field<T::key>(),
                        z_pred.get_field<T::key>(), get_field<T::key>());
        } else {
            return;
        }
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_covariance(CovarianceMatrix& P, const DynamicMeasurementVector& z_pred) const {
        calculate_field_covariance<T1>(P, z_pred);
        calculate_field_covariance<T2, Tail...>(P, z_pred);
    }

    /*
    These functions build the measurement root covariance from all populated
    fields.
    */
    template <typename T>
    static Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>> field_root_covariance(
            const T& p, const T& z_pred, const T& z) {
        Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>> temp =
            Matrix<Detail::CovarianceDimension<T>, Detail::CovarianceDimension<T>>::Zero();
        temp.diagonal() << p;
        return temp;
    }

    static Matrix<3, 3> field_root_covariance(const FieldVector& p, const FieldVector& z_pred, const FieldVector& z) {
        return Detail::calculate_rotation_vector_jacobian<DynamicMeasurementVector>(z, z_pred) *
            Eigen::DiagonalMatrix<real_t, 3>(p);
    }

    template <typename T>
    void calculate_field_root_covariance(CovarianceMatrix& P, const DynamicMeasurementVector& z_pred) const {
        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
        if(offset != std::numeric_limits<std::size_t>::max()) {
            P.template block<Detail::CovarianceDimension<typename T::type>,
                Detail::CovarianceDimension<typename T::type>>(offset, offset) =
                Detail::MeasurementStateHelper<DynamicMeasurementVector>::field_root_covariance(
                    measurement_root_covariance.template get_field<T::key>(),
                    z_pred.get_field<T::key>(), get_field<T::key>());
        } else {
            return;
        }
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_root_covariance(CovarianceMatrix& P, const DynamicMeasurementVector& z_pred) const {
        calculate_field_root_covariance<T1>(P, z_pred);
        calculate_field_root_covariance<T2, Tail...>(P, z_pred);
    }

    /*
    Private functions for calculating the mean of each field in a sigma point
    distribution.
    */
    template <typename S, typename T>
    void calculate_field_mean(const SigmaPointDistribution<S>& Z, DynamicMeasurementVector& mean) const {
        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
        if(offset != std::numeric_limits<std::size_t>::max()) {
            mean.template segment<Detail::StateVectorDimension<typename T::type>>(offset) <<
                Detail::MeasurementStateHelper<DynamicMeasurementVector, S>::sigma_point_mean(
                    Z.template block<Detail::StateVectorDimension<typename T::type>, S::num_sigma()>(
                    offset, 0), typename T::type());
        } else {
            return;
        }
    }

    template <typename S, typename T1, typename T2, typename... Tail>
    void calculate_field_mean(const SigmaPointDistribution<S>& Z, DynamicMeasurementVector& mean) const {
        calculate_field_mean<S, T1>(Z, mean);
        calculate_field_mean<S, T2, Tail...>(Z, mean);
    }

    /* These functions build the innovation from all populated fields. */
    template <typename T>
    void calculate_field_innovation(DynamicMeasurementVector& z) const {
        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
        if(offset != std::numeric_limits<std::size_t>::max()) {
            z.set_field<T::key>(Detail::MeasurementStateHelper<DynamicMeasurementVector>::measurement_delta(
                z.get_field<T::key>(), get_field<T::key>()));
        } else {
            return;
        }
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_innovation(DynamicMeasurementVector& z) const {
        calculate_field_innovation<T1>(z);
        calculate_field_innovation<T2, Tail...>(z);
    }

    /*
    Private functions for calculating the delta vectors between each sigma
    point and the mean.
    */
    template <typename S, typename T>
    void calculate_field_deltas(const SigmaPointDistribution<S>& Z, SigmaPointDeltas<S>& z_prime) const {
        std::size_t offset = std::get<Detail::get_field_order<0, Fields...>(T::key)>(field_offsets);
        if(offset != std::numeric_limits<std::size_t>::max()) {
            z_prime.template block<Detail::CovarianceDimension<typename T::type>, S::num_sigma()>(offset, 0) =
                Detail::MeasurementStateHelper<DynamicMeasurementVector, S>::sigma_point_deltas(get_field<T::key>(), 
                    Z.template block<Detail::StateVectorDimension<typename T::type>, S::num_sigma()>(offset, 0));
        } else {
            return;
        }
    }

    template <typename S, typename T1, typename T2, typename... Tail>
    void calculate_field_deltas(const SigmaPointDistribution<S>& Z, SigmaPointDeltas<S>& z_prime) const {
        calculate_field_deltas<S, T1>(Z, z_prime);
        calculate_field_deltas<S, T2, Tail...>(Z, z_prime);
    }
};

}

#endif
