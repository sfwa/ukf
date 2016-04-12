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
    using SigmaPointDistribution = Matrix<size(), S::num_sigma>;
    using CovarianceMatrix = Matrix<covariance_size(), size()>;

    template <int Key>
    auto field() {
        static_assert(Detail::GetFieldOffset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");
        return Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(Detail::GetFieldOffset<0, Fields...>(Key));
    }

    template <int Key>
    auto field() const {
        static_assert(Detail::GetFieldOffset<0, Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");
        return Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(Detail::GetFieldOffset<0, Fields...>(Key));
    }

    /* Calculate the mean from a measurement sigma point distribution. */
    template <typename S>
    static FixedMeasurementVector calculate_sigma_point_mean(const SigmaPointDistribution<S> &Z) {
        return Parameters::Sigma_WMI<S>*Z.block(0, 1, size(), S::num_sigma-1).rowwise().sum()
            + Parameters::Sigma_WM0<S>*Z.col(0);
    }

    /*
    Calculate the expected measurement covariance covariance as described in
    equation 68 of the Kraft papers.
    The function isn't static; it uses the current measurement vector as the mean.
    */
    template <typename S>
    CovarianceMatrix calculate_sigma_point_covariance(const SigmaPointDistribution<S> &Z) const {
        CovarianceMatrix cov;
        SigmaPointDistribution<S> z_prime;

        /* Calculate the delta vectors. */
        z_prime = Z.colwise() - *this;

        /* Calculate the covariance using equation 64 from the Kraft paper. */
        cov = Parameters::Sigma_WC0<S> * (z_prime.col(0) * z_prime.col(0).transpose());
        for(int i = 1; i < S::num_sigma; i++) {
            cov += Parameters::Sigma_WCI<S> * (z_prime.col(i) * z_prime.col(i).transpose());
        }

        return cov;
    }

    /*
    Create a measurement sigma point distribution using the sigma points.
    Return value optimisation will ensure this does not involve a copy.
    */
    template <typename S>
    static SigmaPointDistribution<S> calculate_sigma_point_distribution(const typename S::SigmaPointDistribution &X) {
        SigmaPointDistribution<S> Z;

        for(int i = 0; i < S::num_sigma; i++) {
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

    Template parameters are a StateVector type and a Field type.
    */
    template <typename S, typename T>
    static typename T::type expected_measurement(const S &state);

    /*
    These functions build the measurement estimate from the expected
    measurement of each individual field.
    */
    template <typename S, typename T>
    static void calculate_field_measurements(const S &state, FixedMeasurementVector &expected) {
        expected.segment(Detail::GetFieldOffset<0, Fields...>(T::key),
            Detail::StateVectorDimension<typename T::type>) << expected_measurement<S, T>(state);
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
    template <int N>
    using SigmaPointDistribution = MatrixDynamic<max_size(), N>;
    using CovarianceMatrix = MatrixDynamic<max_covariance_size(), max_size()>;

    template <int Key>
    auto field() {
        std::size_t offset = get_offset(Key);

        static_assert(Detail::GetFieldSize<Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");

        /* Check if this field has already been set. If so, replace it. */
        if(offset < Base::template size()) {
            return Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(offset);
        } else {
            /*
            Otherwise, resize the measurement vector to fit it and store the
            order in which fields have been set.
            */
            std::size_t previous_size = Base::template size();
            Base::template conservativeResize(previous_size + Detail::GetFieldSize<Fields...>(Key));

            /*
            Resize the current_measurements matrix and store the new key at
            the end.
            */
            std::size_t num_measurements = current_measurements.size();
            current_measurements.conservativeResize(num_measurements + 1);
            current_measurements(num_measurements) = Key;

            /* Assign the value to the field. */
            return Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(previous_size);
        }
    }

    /* Read-only version of the field accessor method. */
    template <int Key>
    auto field() const {
        std::size_t offset = get_offset(Key);

        static_assert(Detail::GetFieldSize<Fields...>(Key) != std::numeric_limits<std::size_t>::max(),
            "Specified key not present in measurement vector");

        assert(offset != std::numeric_limits<std::size_t>::max() &&
            "Specified key not present in measurement vector");

        return Base::template segment<Detail::GetFieldSize<Fields...>(Key)>(offset);
    }

private:
    /*
    This vector keeps track of which fields have been set in the measurement
    vector, and the order they were supplied in. This allows any combination
    of measurements to be supplied in any order and they will be handled
    correctly.
    */
    Eigen::Matrix<int, Eigen::Dynamic, 1, 0, sizeof...(Fields), 1> current_measurements;

    /*
    This method gets the offset of the specified key in the measurement
    vector, or returns std::numeric_limits<std::size_t>::max() if it's not
    present.
    */
    std::size_t get_offset(int Key) {
        std::size_t offset = 0;
        for(int i = 0; i < current_measurements.size(); i++) {
            if(current_measurements(i) == Key) {
                break;
            }

            offset += Detail::GetFieldSize<Fields...>(current_measurements(i));
        }

        return offset;
    }
};

}

#endif