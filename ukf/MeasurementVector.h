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
#include "Config.h"
#include "Types.h"
#include "StateVector.h"

namespace UKF {

/* Alias for the Eigen type from which FixedMeasurementVector inherits. */
template <typename... Fields>
using MeasurementVectorFixedBaseType = Vector<Detail::GetCompositeVectorDimension<Fields...>()>;

/* Alias for the Eigen type from which DynamicMeasurementVector inherits. */
template <typename... Fields>
using MeasurementVectorDynamicBaseType = VectorDynamic<Detail::GetCompositeVectorDimension<Fields...>()>;

/* Templated measurement vector abstract base class. */
template <template<typename...> class B, typename... Fields>
class MeasurementVector : public B<typename Fields::type...> {
public:
    /* Inherit Eigen::Matrix constructors and assignment operators. */
    using Base = B<typename Fields::type...>;
    using Base::Base;
    using Base::operator=;

private:
    /*
    Measurement covariance is represented as a vector the same length as the
    measurement vector.
    */
    static Base measurement_covariance;
};

/*
This class provides a fixed measurement vector, to be used when the same
measurement are available every time step.
*/
template <typename... Fields>
class FixedMeasurementVector : public MeasurementVector<MeasurementVectorFixedBaseType, Fields...> {
public:
    using Base = MeasurementVector<MeasurementVectorFixedBaseType, Fields...>;
    using Base::Base;
    using Base::operator=;

    /* Get size of measurement vector. */
    static constexpr std::size_t size() {
        return Detail::GetCompositeVectorDimension<typename Fields::type...>();
    }

    /* Get size of measurement vector covariance. */
    static constexpr std::size_t covariance_size() {
        return Detail::GetCovarianceDimension<typename Fields::type...>();
    }

    /* Aliases for types needed during filter iteration. */
    template <typename S>
    using SigmaPointDistribution = Matrix<size(), S::num_sigma()>;
    template <typename S>
    using SigmaPointDeltas = Matrix<size(), S::num_sigma()>;
    using CovarianceMatrix = Matrix<covariance_size(), covariance_size()>;

    /* Functions for accessing individual fields. */
    template <int Key>
    typename Detail::FieldTypes<Key, Fields...>::type get_field() const {
        static_assert(Detail::GetFieldOffset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in state vector");
        return Detail::ConvertFromSegment<typename Detail::FieldTypes<Key, Fields...>::type>(
            Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(Detail::GetFieldOffset<0, Fields...>(Key)));
    }

    template <int Key, typename T>
    void set_field(T in) {
        static_assert(Detail::GetFieldOffset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in state vector");
        Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(Detail::GetFieldOffset<0, Fields...>(Key)) << in;
    }

    /* Calculate the mean from a measurement sigma point distribution. */
    template <typename S>
    FixedMeasurementVector calculate_sigma_point_mean(const SigmaPointDistribution<S> &Z) const {
        return Parameters::Sigma_WMI<S>*Z.block(0, 1, size(), S::num_sigma()-1).rowwise().sum()
            + Parameters::Sigma_WM0<S>*Z.col(0);
    }

    /*
    Calculate the set of sigma point delta vectors; these are used for
    calculating the measurement covariance and the Kalman gain.
    The function isn't static; it uses the current measurement vector as the mean.
    */
    template <typename S>
    SigmaPointDeltas<S> calculate_sigma_point_deltas(const SigmaPointDistribution<S> &Z) const {
        SigmaPointDeltas<S> z_prime;

        /* Calculate the delta vectors. */
        z_prime = Z.colwise() - *this;

        return z_prime;
    }

    /*
    Calculate the expected measurement covariance covariance as described in
    equation 68 of the Kraft papers.
    */
    template <typename S>
    CovarianceMatrix calculate_sigma_point_covariance(const SigmaPointDeltas<S> &z_prime) const {
        CovarianceMatrix cov;

        /* Calculate the covariance using equation 64 from the Kraft paper. */
        cov = Parameters::Sigma_WC0<S> * (z_prime.col(0) * z_prime.col(0).transpose());
        for(int i = 1; i < S::num_sigma(); i++) {
            cov += Parameters::Sigma_WCI<S> * (z_prime.col(i) * z_prime.col(i).transpose());
        }

        return cov;
    }

    /*
    Create a measurement sigma point distribution using the sigma points.
    */
    template <typename S>
    SigmaPointDistribution<S> calculate_sigma_point_distribution(const typename S::SigmaPointDistribution &X) const {
        SigmaPointDistribution<S> Z;

        for(int i = 0; i < S::num_sigma(); i++) {
            FixedMeasurementVector temp;
            calculate_field_measurements<S, Fields...>(X.col(i), temp);
            Z.col(i) = temp;
        }

        return Z;
    }

private:
    /*
    This function is intended to be specialised by the user for each field in
    the measurement vector, and allows the user to specify how a particular
    state vector is transformed into a measurement vector.

    Template parameters are a StateVector type and a field key.
    */
    template <typename S, int Key>
    static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement(const S &state);

    /*
    These functions build the measurement estimate from the expected
    measurement of each individual field.
    */
    template <typename S, typename T>
    static void calculate_field_measurements(const S &state, FixedMeasurementVector &expected) {
        expected.segment(Detail::GetFieldOffset<0, Fields...>(T::key),
            Detail::StateVectorDimension<typename T::type>) << expected_measurement<S, T::key>(state);
    }

    template <typename S, typename T1, typename T2, typename... Tail>
    static void calculate_field_measurements(const S &state, FixedMeasurementVector &expected) {
        calculate_field_measurements<S, T1>(state, expected);
        calculate_field_measurements<S, T2, Tail...>(state, expected);
    }
};

/*
This class provides a dynamic measurement vector, to be used when not all
measurements are available every time step.
*/
template <typename... Fields>
class DynamicMeasurementVector : public MeasurementVector<MeasurementVectorDynamicBaseType, Fields...> {
public:
    using Base = MeasurementVector<MeasurementVectorDynamicBaseType, Fields...>;
    using Base::Base;
    using Base::operator=;

    /* Get maximum size of dynamic measurement vector. */
    static constexpr std::size_t max_size() {
        return Detail::GetCompositeVectorDimension<typename Fields::type...>();
    }

    /* Get maximum size of dynamic measurement vector covariance. */
    static constexpr std::size_t max_covariance_size() {
        return Detail::GetCovarianceDimension<typename Fields::type...>();
    }

    /* Aliases for types needed during filter iteration. */
    template <typename S>
    using SigmaPointDistribution = MatrixDynamic<max_size(), S::num_sigma()>;
    template <typename S>
    using SigmaPointDeltas = MatrixDynamic<max_size(), S::num_sigma()>;
    using CovarianceMatrix = MatrixDynamic<max_covariance_size(), max_covariance_size()>;

    /* Functions for accessing individual fields. */
    template <int Key>
    typename Detail::FieldTypes<Key, Fields...>::type get_field() const {
        static_assert(Detail::GetFieldSize<Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");

        std::size_t offset = std::get<Detail::GetFieldOrder<0, Fields...>(Key)>(field_offsets);

        assert(offset != std::numeric_limits<std::size_t>::max() &&
            "Specified key not present in measurement vector");

        return Detail::ConvertFromSegment<typename Detail::FieldTypes<Key, Fields...>::type>(
            Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(offset));
    }

    template <int Key, typename T>
    void set_field(T in) {
        static_assert(Detail::GetFieldSize<Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");

        std::size_t offset = std::get<Detail::GetFieldOrder<0, Fields...>(Key)>(field_offsets);

        /* Check if this field has already been set. If so, replace it. */
        if(offset < Base::template size()) {
            Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(offset) << in;
        } else {
            /*
            Otherwise, resize the measurement vector to fit it and store the
            order in which fields have been set.
            */
            std::size_t previous_size = Base::template size();
            Base::template conservativeResize(previous_size + Detail::GetFieldSize<Fields...>(Key));

            /* Assign the value to the field. */
            Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(previous_size) << in;

            /* Store the offset in field_offsets. */
            std::get<Detail::GetFieldOrder<0, Fields...>(Key)>(field_offsets) = previous_size;
        }
    }

    /*
    Calculate the mean from a measurement sigma point distribution. Ensure
    that the returned object has the same field_offsets array so that its
    field accessors work.
    */
    template <typename S>
    DynamicMeasurementVector calculate_sigma_point_mean(const SigmaPointDistribution<S> &Z) const {
        DynamicMeasurementVector temp = DynamicMeasurementVector(
            Parameters::Sigma_WMI<S>*Z.block(0, 1, Base::template size(), S::num_sigma()-1).rowwise().sum()
            + Parameters::Sigma_WM0<S>*Z.col(0));

        temp.field_offsets = field_offsets;

        return temp;
    }

    /*
    Calculate the set of sigma point delta vectors; these are used for
    calculating the measurement covariance and the Kalman gain.
    The function isn't static; it uses the current measurement vector as the mean.
    */
    template <typename S>
    SigmaPointDeltas<S> calculate_sigma_point_deltas(const SigmaPointDistribution<S> &Z) const {
        SigmaPointDeltas<S> z_prime(Base::template size(), S::num_sigma());

        /* Calculate the delta vectors. */
        z_prime = Z.colwise() - *this;

        return z_prime;
    }

    /*
    Calculate the expected measurement covariance covariance as described in
    equation 68 of the Kraft papers.
    */
    template <typename S>
    CovarianceMatrix calculate_sigma_point_covariance(const SigmaPointDeltas<S> &z_prime) const {
        CovarianceMatrix cov(Base::template size(), Base::template size());
        
        /* Calculate the covariance using equation 64 from the Kraft paper. */
        cov = Parameters::Sigma_WC0<S> * (z_prime.col(0) * z_prime.col(0).transpose());
        for(int i = 1; i < S::num_sigma(); i++) {
            cov += Parameters::Sigma_WCI<S> * (z_prime.col(i) * z_prime.col(i).transpose());
        }

        return cov;
    }

    /*
    Create a measurement sigma point distribution using the sigma points.
    */
    template <typename S>
    SigmaPointDistribution<S> calculate_sigma_point_distribution(const typename S::SigmaPointDistribution &X) const {
        SigmaPointDistribution<S> Z(Base::template size(), S::num_sigma());

        for(int i = 0; i < S::num_sigma(); i++) {
            DynamicMeasurementVector temp(Base::template size());
            calculate_field_measurements<S, Fields...>(X.col(i), temp);
            Z.col(i) = temp;
        }

        return Z;
    }

private:
    /*
    This vector keeps track of which fields have been set in the measurement
    vector, and the offset within the measurement vector of each field. This
    allows any combination of measurements to be supplied in any order and
    they will be handled correctly.
    */
    std::array<std::size_t, sizeof...(Fields)> field_offsets = Detail::CreateArray<sizeof...(Fields)>(
        std::numeric_limits<std::size_t>::max());

    /*
    This function is intended to be specialised by the user for each field in
    the measurement vector, and allows the user to specify how a particular
    state vector is transformed into a measurement vector.

    Template parameters are a StateVector type and a field key.
    */
    template <typename S, int Key>
    static typename Detail::FieldTypes<Key, Fields...>::type expected_measurement(const S &state);

    /*
    These functions build the measurement estimate from the expected
    measurement of each individual field.
    */
    template <typename S, typename T>
    void calculate_field_measurements(const S &state, DynamicMeasurementVector &expected) const {
        /*
        If this field has been set, then generate an expected measurement for
        it. Otherwise, do nothing.
        */
        std::size_t offset = std::get<Detail::GetFieldOrder<0, Fields...>(T::key)>(field_offsets);
        if(offset != std::numeric_limits<std::size_t>::max()) {
            expected.segment(offset, Detail::StateVectorDimension<typename T::type>) <<
                expected_measurement<S, T::key>(state);
        } else {
            return;
        }
    }

    template <typename S, typename T1, typename T2, typename... Tail>
    void calculate_field_measurements(const S &state, DynamicMeasurementVector &expected) const {
        calculate_field_measurements<S, T1>(state, expected);
        calculate_field_measurements<S, T2, Tail...>(state, expected);
    }
};

}

#endif