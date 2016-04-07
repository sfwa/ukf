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

    /* Get size of measurement vector. */
    static constexpr std::size_t size() {
        return Detail::GetCompositeVectorDimension<typename Fields::type...>();
    }

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

private:

};

/*
This class provides a dynamic measurement vector, to be used when not all
measurements are available every time step.
*/
template <typename... Fields>
class DynamicMeasurementVector : public MeasurementVector<MeasurementVectorDynamicBaseType, Fields...> {
public:
    using Base = MeasurementVector<MeasurementVectorDynamicBaseType, Fields...>;

    /* Get maximum size of dynamic measurement vector. */
    static constexpr std::size_t max_size() {
        return Detail::GetCompositeVectorDimension<typename Fields::type...>();
    }

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