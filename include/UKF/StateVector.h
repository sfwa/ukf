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
#include "UKF/Types.h"
#include "UKF/Integrator.h"

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
    inline constexpr std::size_t StateVectorDimension<Quaternion> = 4;

    template <>
    inline constexpr std::size_t StateVectorDimension<real_t> = 1;

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
    inline constexpr std::size_t CovarianceDimension<Quaternion> = 3;

    template <typename T>
    constexpr T adder(T v) {
        return v;
    }

    template <typename T, typename... Args>
    constexpr T adder(T first, Args... args) {
        return first + adder(args...);
    }

    /*
    Get the dimension of the state vector by summing the dimension of all
    fields.
    */
    template <typename... Fields>
    constexpr std::size_t get_composite_vector_dimension() {
        return adder(StateVectorDimension<Fields>...);
    }

    /*
    Get the dimension of the covariance matrix by summing the dimension of
    all fields.
    */
    template <typename... Fields>
    constexpr std::size_t get_covariance_dimension() {
        return adder(CovarianceDimension<Fields>...);
    }

    /*
    Get the offset of a particular field in a parameter pack of Field objects
    by matching the provided key parameter.
    */
    template <std::size_t Offset, typename T>
    constexpr std::size_t get_field_offset(int Key) {
        return Key == T::key ? Offset : std::numeric_limits<std::size_t>::max();
    }

    template <std::size_t Offset, typename T1, typename T2, typename... Fields>
    constexpr std::size_t get_field_offset(int Key) {
        return Key == T1::key ? Offset : get_field_offset<
            Offset + StateVectorDimension<typename T1::type>, T2, Fields...>(Key);
    }

    /* Do the same as above, but for the covariance matrix. */
    template <std::size_t Offset, typename T>
    constexpr std::size_t get_field_covariance_offset(int Key) {
        return Key == T::key ? Offset : std::numeric_limits<std::size_t>::max();
    }

    template <std::size_t Offset, typename T1, typename T2, typename... Fields>
    constexpr std::size_t get_field_covariance_offset(int Key) {
        return Key == T1::key ? Offset : get_field_covariance_offset<
            Offset + CovarianceDimension<typename T1::type>, T2, Fields...>(Key);
    }

    /* Get the order of the specified key within the field pack. */
    template <std::size_t Offset, typename T>
    constexpr std::size_t get_field_order(int Key) {
        return Key == T::key ? Offset : std::numeric_limits<std::size_t>::max();
    }

    template <std::size_t Offset, typename T1, typename T2, typename... Fields>
    constexpr std::size_t get_field_order(int Key) {
        return Key == T1::key ? Offset : get_field_order<Offset + 1, T2, Fields...>(Key);
    }

    /*
    These helper structs allow various functions to determine field types at
    compile time based only on the integer key parameter.
    */
    template <bool condition, class T, class U>
    struct IfHelper {
      using type = U;
    };

    template <class T, class U>
    struct IfHelper<true, T, U> {
      using type = T;
    };

    template <int Key, typename... Fields>
    struct FieldTypesBase {
        using type = void;
    };

    template <int Key, typename T, typename... Fields>
    struct FieldTypesBase<Key, T, Fields...> {
        using type = typename IfHelper<
            Key == T::key,
            typename T::type,
            typename FieldTypesBase<Key, Fields...>::type>::type;
    };

    /*
    A `boost::mpl::identity`-like wrapper to avoid inference errors when
    specializing function templates.
    In some situations when using FieldTypesBase directly only the primary
    template gets resolved, resulting in failed substitutions.
    This is solved in C++17, e.g. when compiling by GCC8
    */
    template <int Key, typename... Fields>
    struct FieldTypes : FieldTypesBase<Key, Fields...> {};

    template <typename T>
    inline T convert_from_segment(const Vector<StateVectorDimension<T>>& state) {
        return static_cast<T>(state);
    }

    template <>
    inline real_t convert_from_segment<real_t>(const Vector<1>& state) {
        return static_cast<real_t>(state(0));
    }

    /*
    Get the size of a particular field in a parameter pack of Field objects
    by matching the provided key parameter.
    */
    template <typename T>
    constexpr std::size_t get_field_size(int Key) {
        return Key == T::key ? StateVectorDimension<typename T::type> : std::numeric_limits<std::size_t>::max();
    }

    template <typename T1, typename T2, typename... Fields>
    constexpr std::size_t get_field_size(int Key) {
        return Key == T1::key ? StateVectorDimension<typename T1::type> : get_field_size<T2, Fields...>(Key);
    }

    /* Do the same as above, but for the covariance matrix. */
    template <typename T>
    constexpr std::size_t get_field_covariance_size(int Key) {
        return Key == T::key ? CovarianceDimension<typename T::type> : std::numeric_limits<std::size_t>::max();
    }

    template <typename T1, typename T2, typename... Fields>
    constexpr std::size_t get_field_covariance_size(int Key) {
        return Key == T1::key ? CovarianceDimension<typename T1::type> : get_field_covariance_size<T2, Fields...>(Key);
    }

    /* Function for creating an array initialised to a specific value. */
    template <typename T, std::size_t... Indices>
    constexpr std::array<T, sizeof...(Indices)> create_array(T value, std::index_sequence<Indices...>) {
        // Cast Indices to void to remove the unused value warning.
        return {{(static_cast<void>(Indices), value)...}};
    }

    template <std::size_t N, typename T>
    constexpr std::array<T, N> create_array(const T& value) {
        return create_array(value, std::make_index_sequence<N>());
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

    Note: By default, we use Gibbs vectors, which have a singularity at 180
    degrees. This is to avoid numerical issues calculating the covariance
    matrix when the quaternion covariance vector is large and the SUT scaling
    parameter alpha is set very small.

    The singularity being 180 degrees instead of 360 is not a problem unless
    the attitude is expected to change by more than 180 degrees in a single
    filter iteration; if it is, setting the MRP_A parameter to 1.0 moves the
    singularity to 360 degrees.
    */
    template <typename T> constexpr real_t MRP_A = 0.0;
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

    /* Functions requiring things from the "Parameters" namespace. */
    namespace Detail {
        /* Convert a rotation vector into a quaternion. */
        template <typename T>
        inline Quaternion rotation_vector_to_quaternion(const Vector<3>& r) {
            real_t x_2 = r.squaredNorm();
            real_t d_q_w = (-UKF::Parameters::MRP_A<T> * x_2 + UKF::Parameters::MRP_F<T>
                * std::sqrt(x_2 * (real_t(1.0) - UKF::Parameters::MRP_A<T>*UKF::Parameters::MRP_A<T>)
                + UKF::Parameters::MRP_F<T> * UKF::Parameters::MRP_F<T>))
                / (UKF::Parameters::MRP_F<T> * UKF::Parameters::MRP_F<T> + x_2);
            Vector<3> d_q_xyz = r * (d_q_w + UKF::Parameters::MRP_A<T>) * (real_t(1.0) / UKF::Parameters::MRP_F<T>);
            Quaternion d_q;
            d_q.vec() = d_q_xyz;
            d_q.w() = d_q_w;

            return d_q;
        }
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
using StateVectorBaseType = Vector<Detail::get_composite_vector_dimension<Fields...>()>;

/*
Templated state vector class. A particular UKF implementation should
specialise this class with a list of fields that make up the state vector.

By default, fields can be Eigen vectors (including Quaternions) or scalar
floating point types (real_t). Support for other types can be added by
specialising the StateVectorDimension variable for the desired class.
*/
template <typename... Fields>
class StateVector : public StateVectorBaseType<typename Fields::type...> {
public:
    /* Inherit Eigen::Matrix constructors and assignment operators. */
    using Base = StateVectorBaseType<typename Fields::type...>;
    using Base::Base;
    using Base::operator=;

    /* Get size of state vector. */
    static constexpr std::size_t size() {
        return Detail::get_composite_vector_dimension<typename Fields::type...>();
    }

    /* Get size of state vector delta. */
    static constexpr std::size_t covariance_size() {
        return Detail::get_covariance_dimension<typename Fields::type...>();
    }

    static constexpr std::size_t num_sigma() {
        return 2*covariance_size() + 1;
    }

    /* Aliases for types needed during filter iteration. */
    using CovarianceMatrix = Matrix<covariance_size(), covariance_size()>;
    using SigmaPointDistribution = Matrix<size(), num_sigma()>;
    using SigmaPointDeltas = Matrix<covariance_size(), num_sigma()>;
    using StateVectorDelta = Vector<covariance_size()>;

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

    template <int Key>
    void set_field(Quaternion in) {
        static_assert(Detail::get_field_offset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in state vector");
        Base::template segment<Detail::get_field_size<Fields...>(Key)>(
            Detail::get_field_offset<0, Fields...>(Key)) << in.vec(), in.w();
    }

    /* Calculate the mean from a sigma point distribution. */
    static StateVector calculate_sigma_point_mean(const SigmaPointDistribution& X) {
        StateVector mean;
        calculate_field_mean<Fields...>(X, mean);

        return mean;
    }

    /*
    Calculate the set of sigma point delta vectors; these are used for
    calculating the a priori covariance and the Kalman gain.
    */
    SigmaPointDeltas calculate_sigma_point_deltas(const SigmaPointDistribution& X) const {
        SigmaPointDeltas w_prime;

        /* Calculate the delta vectors. */
        calculate_field_deltas<Fields...>(X, w_prime);

        return w_prime;
    }

    /*
    Calculate the covariance as described in section 3.5.1 of the Kraft
    paper. Note that we operate on the transformed sigma points; there
    appears to be a typographical error in equation 63, but the explanatory
    text makes it clear that this is what is intended.
    The function isn't static; it uses the current state vector as the mean.
    */
    static CovarianceMatrix calculate_sigma_point_covariance(const SigmaPointDeltas& w_prime) {
        CovarianceMatrix cov;

        /* Calculate the covariance using equation 64 from the Kraft paper. */
        cov = CovarianceMatrix::Zero();
        for(std::size_t i = 1; i < num_sigma(); i++) {
            cov.noalias() += Parameters::Sigma_WCI<StateVector> * (w_prime.col(i) * w_prime.col(i).transpose());
        }
        cov.noalias() += Parameters::Sigma_WC0<StateVector> * (w_prime.col(0) * w_prime.col(0).transpose());

        return cov;
    }

    /*
    Create a sigma point distribution using the provided covariance matrix.
    */
    SigmaPointDistribution calculate_sigma_point_distribution(const CovarianceMatrix& S) const {
        /* Calculate the sigma point distribution from all the fields. */
        SigmaPointDistribution X;
        calculate_field_sigmas<Fields...>(S, X);

        return X;
    }

    /*
    This function calculates the derivative of the state vector. The
    derivative can be a function of the current state and any number of user-
    supplied non-state inputs.
    */
    template <typename... U>
    StateVector derivative(const U&... input) const;

    /*
    Apply the process model to the state vector to calculate the predicted
    state based on the supplied time delta. This is achieved by using a
    numerical integrator and the derivative function.
    */
    template <typename IntegratorType = IntegratorRK4, typename... U>
    StateVector process_model(real_t delta, const U&... input) const {
        return IntegratorType::integrate(delta, *this, input...);
    }

    /* Update the state vector using a delta vector. */
    void apply_delta(const StateVectorDelta& delta) {
        apply_field_deltas<Fields...>(delta);
    }

private:
    /* Private functions for creating a sigma point distribution. */
    template <typename T>
    static Matrix<Detail::StateVectorDimension<T>, num_sigma()> perturb_state(
            const T& state, const Matrix<Detail::CovarianceDimension<T>, covariance_size()>& cov) {
        Matrix<Detail::StateVectorDimension<T>, num_sigma()> temp;
        temp.col(0) = state;
        temp.template block<Detail::StateVectorDimension<T>, covariance_size()>(0, 1) = cov.colwise() + state;
        temp.template block<Detail::StateVectorDimension<T>, covariance_size()>(0, covariance_size()+1) =
            -(cov.colwise() - state);

        return temp;
    }

    static Matrix<1, num_sigma()> perturb_state(real_t state, const Matrix<1, covariance_size()>& cov) {
        Matrix<1, num_sigma()> temp;
        temp(0) = state;
        temp.template segment<covariance_size()>(1) = cov.array() + state;
        temp.template segment<covariance_size()>(covariance_size()+1) = -(cov.array() - state);

        return temp;
    }

    /*
    Construct error quaternions using the MRP method, equation 34 from the
    Markley paper.
    */
    static Matrix<4, num_sigma()> perturb_state(const Quaternion& state, const Matrix<3, covariance_size()>& cov) {
        Matrix<4, num_sigma()> temp;
        temp.col(0) << state.vec(), state.w();

        Array<1, covariance_size()> x_2 = cov.colwise().squaredNorm();
        Array<1, covariance_size()> err_w =
            (-Parameters::MRP_A<StateVector> * x_2 + Parameters::MRP_F<StateVector>
                * (x_2 * (real_t(1.0) - Parameters::MRP_A<StateVector>*Parameters::MRP_A<StateVector>)
                + Parameters::MRP_F<StateVector> * Parameters::MRP_F<StateVector>).sqrt())
            / (Parameters::MRP_F<StateVector> * Parameters::MRP_F<StateVector> + x_2);
        Array<3, covariance_size()> err_xyz = cov.array().rowwise() * (err_w + Parameters::MRP_A<StateVector>)
            * (real_t(1.0) / Parameters::MRP_F<StateVector>);

        Quaternion temp_q;
        for(std::size_t i = 0; i < covariance_size(); i++) {
            temp_q = Quaternion(err_w(i), err_xyz(0, i), err_xyz(1, i), err_xyz(2, i)) * state;
            temp.col(i+1) << temp_q.vec(), temp_q.w();
        }

        for(std::size_t i = 0; i < covariance_size(); i++) {
            temp_q = Quaternion(err_w(i), err_xyz(0, i), err_xyz(1, i), err_xyz(2, i)).conjugate() * state;
            temp.col(i+covariance_size()+1) << temp_q.vec(), temp_q.w();
        }

        return temp;
    }

    template <typename T>
    void calculate_field_sigmas(const CovarianceMatrix& S, SigmaPointDistribution& X) const {
        X.template block<Detail::StateVectorDimension<typename T::type>, num_sigma()>(
            Detail::get_field_offset<0, Fields...>(T::key), 0) = perturb_state(get_field<T::key>(),
                S.template block<Detail::CovarianceDimension<typename T::type>, covariance_size()>(
                    Detail::get_field_covariance_offset<0, Fields...>(T::key), 0));
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_sigmas(const CovarianceMatrix& S, SigmaPointDistribution& X) const {
        calculate_field_sigmas<T1>(S, X);
        calculate_field_sigmas<T2, Tail...>(S, X);
    }

    /*
    Private functions for calculating the mean of each field in a sigma point
    distribution. Note that sigma_point_mean takes a dummy argument so that
    the overrides work properly.

    Calculate the mean as described in section 3.4 of the Kraft paper, using
    the weights from the scaled unscented transform.
    For calculating a priori attitude estimate, we use the algorithm for
    computing the quaternion mean found in "Unscented Filtering in a Unit
    Quaternion Space for Spacecraft Attitude Estimation" by Yee-Jin Cheon.
    The following algorithm implements equation (41d) from that paper.
    */
    template <typename T>
    static T sigma_point_mean(const Matrix<Detail::StateVectorDimension<T>, num_sigma()>& sigma, const T& field) {
        return Parameters::Sigma_WMI<StateVector>*sigma.template block<Detail::StateVectorDimension<T>, num_sigma()-1>(
            0, 1).rowwise().sum() + Parameters::Sigma_WM0<StateVector>*sigma.col(0);
    }

    static real_t sigma_point_mean(const Matrix<1, num_sigma()>& sigma, const real_t& field) {
        return Parameters::Sigma_WMI<StateVector>*sigma.template segment<num_sigma()-1>(1).sum()
            + Parameters::Sigma_WM0<StateVector>*sigma(0);
    }

    /*
    Calculate the quaternion barycentric mean with renormalisation. Note that
    this is not an ad-hoc renormalisation of an appoximation; see the paper
    mentioned above for details.
    */
    static Vector<4> sigma_point_mean(const Matrix<4, num_sigma()>& sigma, const Quaternion& field) {
        Vector<4> temp = Parameters::Sigma_WMI<StateVector>*sigma.template block<4, num_sigma()-1>(
            0, 1).rowwise().sum() + Parameters::Sigma_WM0<StateVector>*sigma.col(0);
        Quaternion temp_q = Quaternion(temp).normalized();
        return Vector<4>(temp_q.x(), temp_q.y(), temp_q.z(), temp_q.w());
    }

    template <typename T>
    static void calculate_field_mean(const SigmaPointDistribution& X, StateVector& mean) {
        mean.template segment<Detail::StateVectorDimension<typename T::type>>(
            Detail::get_field_offset<0, Fields...>(T::key)) << sigma_point_mean(
                X.template block<Detail::StateVectorDimension<typename T::type>, num_sigma()>(
                    Detail::get_field_offset<0, Fields...>(T::key), 0), typename T::type());
    }

    template <typename T1, typename T2, typename... Tail>
    static void calculate_field_mean(const SigmaPointDistribution& X, StateVector& mean) {
        calculate_field_mean<T1>(X, mean);
        calculate_field_mean<T2, Tail...>(X, mean);
    }

    /*
    Private functions for calculating the delta vectors between each sigma
    point and the mean.
    */
    template <typename T>
    static Matrix<Detail::CovarianceDimension<T>, num_sigma()> sigma_point_deltas(
            const T& mean, const Matrix<Detail::StateVectorDimension<T>, num_sigma()>& X) {
        return X.colwise() - mean;
    }

    static Matrix<1, num_sigma()> sigma_point_deltas(real_t mean, const Matrix<1, num_sigma()>& X) {
        return X.array() - mean;
    }

    static Matrix<3, num_sigma()> sigma_point_deltas(const Quaternion& mean, const Matrix<4, num_sigma()>& X) {
        Matrix<3, num_sigma()> temp;

        /*
        The attitude part of this set of vectors is calculated using equation
        45 from the Kraft paper.
        */
        for(std::size_t i = 0; i < num_sigma(); i++) {
            Quaternion delta_q = Quaternion(X.col(i)) * mean.conjugate();
            temp.col(i) = Parameters::MRP_F<StateVector> * delta_q.vec() /
                (std::abs(Parameters::MRP_A<StateVector> + delta_q.w()) > std::numeric_limits<real_t>::epsilon() ?
                    Parameters::MRP_A<StateVector> + delta_q.w() : std::numeric_limits<real_t>::epsilon());
        }

        return temp;
    }

    template <typename T>
    void calculate_field_deltas(const SigmaPointDistribution& X, SigmaPointDeltas& w_prime) const {
        w_prime.template block<Detail::CovarianceDimension<typename T::type>, num_sigma()>(
            Detail::get_field_covariance_offset<0, Fields...>(T::key), 0) = sigma_point_deltas(
                get_field<T::key>(), X.template block<Detail::StateVectorDimension<typename T::type>, num_sigma()>(
                    Detail::get_field_offset<0, Fields...>(T::key), 0));
    }

    template <typename T1, typename T2, typename... Tail>
    void calculate_field_deltas(const SigmaPointDistribution& X, SigmaPointDeltas& w_prime) const {
        calculate_field_deltas<T1>(X, w_prime);
        calculate_field_deltas<T2, Tail...>(X, w_prime);
    }

    /*
    Private functions to apply a delta vector to the state vector field by
    field.
    */
    template <typename T>
    static T update_field(const T& state, const Vector<Detail::CovarianceDimension<T>>& delta) {
        return state + delta;
    }

    static real_t update_field(const real_t& state, const Vector<1>& delta) {
        return state + delta(0);
    }

    static Quaternion update_field(const Quaternion& state, const Vector<3>& delta) {
        /*
        Update the attitude quaternion from the MRP vector, equation 45 from
        the Markley paper.
        */
        return Detail::rotation_vector_to_quaternion<StateVector>(delta) * state;
    }

    template <typename T>
    void apply_field_deltas(const StateVectorDelta& delta) {
        set_field<T::key>(update_field(get_field<T::key>(),
            delta.template segment<Detail::get_field_covariance_size<Fields...>(T::key)>(
                Detail::get_field_covariance_offset<0, Fields...>(T::key))));
    }

    template <typename T1, typename T2, typename... Tail>
    void apply_field_deltas(const StateVectorDelta& delta) {
        apply_field_deltas<T1>(delta);
        apply_field_deltas<T2, Tail...>(delta);
    }
};

}

#endif
