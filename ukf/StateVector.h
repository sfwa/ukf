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

#ifndef STATEVECTOR_H
#define STATEVECTOR_H

#include <limits>
#include <tuple>
#include <cstddef>
#include <utility>
#include <Eigen/Core>
#include "Config.h"
#include "Types.h"

namespace UKF {

    namespace Detail {

    /*
    This variable template defines the dimensionality of a particular
    StateVector field. For a normal Eigen column vector, the default
    implementation returns the number of rows at compile time.

    Template specialisations can be used to define the size of fields which
    don't have a MaxRowsAtCompileTime member, or need a custom size for some
    reason (eg. the Eigen Quaternion classes).
    */
    template <typename T>
    constexpr std::size_t StateVectorDimension = T::MaxRowsAtCompileTime;

    template <>
    constexpr std::size_t StateVectorDimension<Quaternion> = 4;

    template <>
    constexpr std::size_t StateVectorDimension<real_t> = 1;

    /*
    This variable template defines the dimensionality of a particular
    StateVector field as it applies to the covariance matrix. For most types,
    this is just the same as the StateVectorDimension.

    An example of a type which needs a special case is a quaternion – it
    needs more special treatment than a simple additive update, so its
    representation in the covariance matrix is different to its
    representation in the state vector.
    */
    template <typename T>
    constexpr std::size_t CovarianceDimension = StateVectorDimension<T>;

    template <>
    constexpr std::size_t CovarianceDimension<Quaternion> = 3;

    template <typename T>
    constexpr T Adder(T v) {
        return v;
    }

    template <typename T, typename... Args>
    constexpr T Adder(T first, Args... args) {
        return first + Adder(args...);
    }

    /*
    Get the dimension of the state vector by summing the dimension of all
    fields.
    */
    template <typename... Fields>
    constexpr std::size_t GetCompositeVectorDimension() {
        return Adder(StateVectorDimension<Fields>...);
    }

    /*
    Get the dimension of the covariance matrix by summing the dimension of
    all fields.
    */
    template <typename... Fields>
    constexpr std::size_t GetCovarianceDimension() {
        return Adder(CovarianceDimension<Fields>...);
    }

    /*
    Get the offset of a particular field in a parameter pack of Field objects
    by matching the provided key parameter.
    */
    template <std::size_t Offset, typename T>
    constexpr std::size_t GetFieldOffset(int Key) {
        return Key == T::key ? Offset : std::numeric_limits<std::size_t>::max();
    }

    template <std::size_t Offset, typename T1, typename T2, typename... Fields>
    constexpr std::size_t GetFieldOffset(int Key) {
        return Key == T1::key ? Offset : GetFieldOffset<
            Offset + StateVectorDimension<typename T1::type>, T2, Fields...>(Key);
    }

    /* Do the same as above, but for the covariance matrix. */
    template <std::size_t Offset, typename T>
    constexpr std::size_t GetFieldCovarianceOffset(int Key) {
        return Key == T::key ? Offset : std::numeric_limits<std::size_t>::max();
    }

    template <std::size_t Offset, typename T1, typename T2, typename... Fields>
    constexpr std::size_t GetFieldCovarianceOffset(int Key) {
        return Key == T1::key ? Offset : GetFieldCovarianceOffset<
            Offset + CovarianceDimension<typename T1::type>, T2, Fields...>(Key);
    }

    template <typename T>
    constexpr T ConvertFromSegment(const Vector<StateVectorDimension<T>> &state) {
        return static_cast<T>(state);
    }

    template <>
    constexpr real_t ConvertFromSegment<real_t>(const Vector<1> &state) {
        return static_cast<real_t>(state(0));
    }

    /*
    Get the size of a particular field in a parameter pack of Field objects
    by matching the provided key parameter.
    */
    template <typename T>
    constexpr std::size_t GetFieldSize(int Key) {
        return Key == T::key ? StateVectorDimension<typename T::type> : std::numeric_limits<std::size_t>::max();
    }

    template <typename T1, typename T2, typename... Fields>
    constexpr std::size_t GetFieldSize(int Key) {
        return Key == T1::key ? StateVectorDimension<typename T1::type> : GetFieldSize<T2, Fields...>(Key);
    }

    /* Do the same as above, but for the covariance matrix. */
    template <typename T>
    constexpr std::size_t GetFieldCovarianceSize(int Key) {
        return Key == T::key ? CovarianceDimension<typename T::type> : std::numeric_limits<std::size_t>::max();
    }

    template <typename T1, typename T2, typename... Fields>
    constexpr std::size_t GetFieldCovarianceSize(int Key) {
        return Key == T1::key ? CovarianceDimension<typename T1::type> : GetFieldCovarianceSize<T2, Fields...>(Key);
    }

    }

    namespace Parameters {

    /*
    This namespace contains the compile-time adjustable parameters used by
    various routines, such as sigma point weights and MRP scaling factors.

    They are implemented using variable templates, so to change them, the
    user simply provides a specialisation for their own version of the state
    vector class.

    These are the default parameters used for the scaled unscented transform.
    See "The Unscented Kalman Filter for Nonlinear Estimation" by Eric A. Wan
    and Rudolph van der Merwe for details.

    Note that alpha^2 here should be small for correct operation of the
    filter. Most literature seems to quote about 1e-3 for alpha (1e-6 for
    alpha^2), but the parameters suggested in "Gaussian Processes for State
    Space Models and Change Point Detection" by Ryan Tuner (2011) provide a
    more stable filter.
    */
    template <typename T> constexpr real_t AlphaSquared = 1.0;
    template <typename T> constexpr real_t Beta = 0.0;
    template <typename T> constexpr real_t Kappa = 3.0;
    template <typename T> constexpr real_t Lambda =
        AlphaSquared<T> * (T::covariance_size() + Kappa<T>) - T::covariance_size();

    /*
    Definitions for parameters used to calculated MRP vectors.
    See the Markley paper for further details.
    */
    template <typename T> constexpr real_t MRP_A = 1.0;
    template <typename T> constexpr real_t MRP_F = 2.0 * (MRP_A<T> + 1.0);

    /*
    Definitions for sigma point weights. The naming convention follows that used
    in in the paper given above.
    */
    template <typename T> constexpr real_t Sigma_WM0 = Lambda<T>/(T::covariance_size() + Lambda<T>);
    template <typename T> constexpr real_t Sigma_WC0 = Sigma_WM0<T> + (1.0 - AlphaSquared<T> + Beta<T>);
    template <typename T> constexpr real_t Sigma_WMI = 1.0 / (2.0 * (T::covariance_size() + Lambda<T>));
    template <typename T> constexpr real_t Sigma_WCI = Sigma_WMI<T>;

    }

/*
A StateVector object takes a variable number of Fields as template parameters.
Each field has a type (from which the size of the field is inferred) and a
key (used to access the field).
*/
template <int Key, typename T>
class Field {
public:
    using type = T;
    static const int key = Key;
};

/* Alias for the Eigen type from which StateVector inherits. */
template <typename... Fields>
using StateVectorBaseType = Vector<Detail::GetCompositeVectorDimension<Fields...>()>;

/*
Templated state vector class. A particular UKF implementation should
specialise this class with a list of fields that make up the state vector.

By default, fields can be Eigen vectors (including Quaternions) or scalar
floating point types (real_t). Support for other types can be added by
specialising the StateVectorDimension variable for the desired class.
*/
template <typename IntegratorType, typename... Fields>
class StateVector : public StateVectorBaseType<typename Fields::type...> {
public:
    /* Inherit Eigen::Matrix constructors and assignment operators. */
    using Base = StateVectorBaseType<typename Fields::type...>;
    using Self = StateVector<IntegratorType, Fields...>;
    using Base::Base;
    using Base::operator=;

    /* Get size of state vector. */
    static constexpr std::size_t size() {
        return Detail::GetCompositeVectorDimension<typename Fields::type...>();
    }

    /* Get size of state vector delta. */
    static constexpr std::size_t covariance_size() {
        return Detail::GetCovarianceDimension<typename Fields::type...>();
    }

    static constexpr std::size_t num_sigma = 2*covariance_size() + 1;

    /* Aliases for types needed during filter iteration. */
    using CovarianceMatrix = Matrix<covariance_size(), covariance_size()>;
    using SigmaPointDistribution = Matrix<size(), num_sigma>;

    /* Methods for accessing individual fields. */
    template <int Key>
    auto field() {
        static_assert(Detail::GetFieldOffset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in state vector");
        return Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(Detail::GetFieldOffset<0, Fields...>(Key));
    }

    template <int Key>
    auto field() const {
        static_assert(Detail::GetFieldOffset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in state vector");
        return Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(Detail::GetFieldOffset<0, Fields...>(Key));
    }

    /* Calculate the mean from a sigma point distribution. */
    static Self calculate_sigma_point_mean(const SigmaPointDistribution &X) {
        Self mean;
        calculate_field_mean<Fields...>(X, mean);

        return mean;
    }

    /* Calculate the covariance matrix from a sigma point distribution. */
    static CovarianceMatrix calculate_sigma_point_covariance(const SigmaPointDistribution &X) {

    }

    /*
    Create a sigma point distribution using the provided covariance matrix.
    Return value optimisation will ensure this does not involve a copy.
    */
    SigmaPointDistribution calculate_sigma_point_distribution(const CovarianceMatrix &P) const {
        /* Calculate the LLT decomposition of the scaled covariance matrix. */
        assert((P * (covariance_size() + Parameters::Lambda<Self>)).llt().info() == Eigen::Success &&
            "Covariance matrix is not positive definite");
        CovarianceMatrix S = (P * (covariance_size() + Parameters::Lambda<Self>)).llt().matrixL();

        /* Calculate the sigma point distribution from all the fields. */
        SigmaPointDistribution X;
        calculate_field_sigmas<Fields...>(S, X);

        return X;
    }

private:
    static IntegratorType integrator;

    /* Private functions for creating a sigma point distribution. */
    template <typename T>
    Matrix<Detail::CovarianceDimension<T>, num_sigma> perturb_state(
            const T &state, const Matrix<Detail::CovarianceDimension<T>, covariance_size()> &cov) const {
        Matrix<Detail::CovarianceDimension<T>, num_sigma> temp;
        temp.col(0) = state;
        temp.block(0, 1, Detail::CovarianceDimension<T>, covariance_size()) = cov.colwise() + state;
        temp.block(0, covariance_size()+1, Detail::CovarianceDimension<T>, covariance_size()) = -(cov.colwise() - state);

        return temp;
    }

    Matrix<1, num_sigma> perturb_state(real_t state, const Matrix<1, covariance_size()> &cov) const {
        Matrix<1, num_sigma> temp;
        temp(0) = state;
        temp.segment(1, covariance_size()) = cov.array() + state;
        temp.segment(covariance_size()+1, covariance_size()) = -(cov.array() - state);

        return temp;
    }

    /*
    Construct error quaternions using the MRP method, equation 34 from the
    Markley paper.
    */
    Matrix<4, num_sigma> perturb_state(
            const Quaternion &state,
            const Matrix<Detail::CovarianceDimension<Quaternion>, covariance_size()> &cov) const {
        Matrix<4, num_sigma> temp;
        temp.col(0) << state.vec(), state.w();

        Array<1, covariance_size()> x_2 = cov.colwise().squaredNorm();
        Array<1, covariance_size()> err_w =
            (-Parameters::MRP_A<Self> * x_2 + Parameters::MRP_F<Self> * (
                x_2 * (1.0 - Parameters::MRP_A<Self>*Parameters::MRP_A<Self>)
                + Parameters::MRP_F<Self> * Parameters::MRP_F<Self>).sqrt())
            / (Parameters::MRP_F<Self> * Parameters::MRP_F<Self> + x_2);
        Array<3, covariance_size()> err_xyz =
            cov.array().rowwise() * (err_w + Parameters::MRP_A<Self>) * (1.0 / Parameters::MRP_F<Self>);

        Quaternion temp_q;
        for(int i = 0; i < covariance_size(); i++) {
            temp_q = Quaternion(err_w(i), err_xyz(0, i), err_xyz(1, i), err_xyz(2, i)) * state;
            temp.col(i+1) << temp_q.vec(), temp_q.w();
        }

        for(int i = 0; i < covariance_size(); i++) {
            temp_q = Quaternion(err_w(i), err_xyz(0, i), err_xyz(1, i), err_xyz(2, i)).conjugate() * state;
            temp.col(i+covariance_size()+1) << temp_q.vec(), temp_q.w();
        }

        return temp;
    }

    template <typename T>
    void calculate_field_sigmas(const CovarianceMatrix &S, SigmaPointDistribution &X) const {
        X.block(Detail::GetFieldOffset<0, Fields...>(T::key), 0,
            Detail::StateVectorDimension<typename T::type>, num_sigma) = perturb_state(
            Detail::ConvertFromSegment<typename T::type>(field<T::key>()),
            S.block(Detail::GetFieldCovarianceOffset<0, Fields...>(T::key), 0,
                Detail::CovarianceDimension<typename T::type>, covariance_size()));
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_sigmas(const CovarianceMatrix &S, SigmaPointDistribution &X) const {
        calculate_field_sigmas<T1>(S, X);
        calculate_field_sigmas<T2, Tail...>(S, X);
    }

    /*
    Private functions for calculating the mean of each field in a sigma point
    distribution. Note that sigma_point_mean takes a dummy argument so that
    the overrides work properly.
    */
    template <typename T>
    static T sigma_point_mean(const Matrix<Detail::StateVectorDimension<T>, num_sigma> &sigma, const T &field) {
        return Parameters::Sigma_WMI<Self>*sigma.block(0, 1, Detail::StateVectorDimension<T>, num_sigma-1).rowwise().sum()
            + Parameters::Sigma_WM0<Self>*sigma.col(0);
    }

    static real_t sigma_point_mean(const Matrix<1, num_sigma> &sigma, const real_t &field) {
        return Parameters::Sigma_WMI<Self>*sigma.segment(1, num_sigma-1).sum()
            + Parameters::Sigma_WM0<Self>*sigma(0);
    }

    static Vector<4> sigma_point_mean(const Matrix<4, num_sigma> &sigma, const Quaternion &field) {
        Vector<4> temp = Parameters::Sigma_WMI<Self>*sigma.block(0, 1, 4, num_sigma-1).rowwise().sum()
            + Parameters::Sigma_WM0<Self>*sigma.col(0);
        Quaternion temp_q = Quaternion(temp(3), temp(0), temp(1), temp(2)).normalized();
        return Vector<4>(temp_q.x(), temp_q.y(), temp_q.z(), temp_q.w());
    }

    template <typename T>
    static void calculate_field_mean(const SigmaPointDistribution &X, StateVector &mean) {
        mean.segment(
            Detail::GetFieldOffset<0, Fields...>(T::key),
            Detail::StateVectorDimension<typename T::type>) << sigma_point_mean(
                X.block(Detail::GetFieldOffset<0, Fields...>(T::key), 0,
                Detail::StateVectorDimension<typename T::type>, num_sigma), typename T::type());
    }

    template <typename T1, typename T2, typename... Tail>
    static void calculate_field_mean(const SigmaPointDistribution &X, StateVector &mean) {
        calculate_field_mean<T1>(X, mean);
        calculate_field_mean<T2, Tail...>(X, mean);
    }
};

}

#endif