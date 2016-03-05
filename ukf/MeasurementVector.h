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
#include <cstddef>
#include <utility>
#include <Eigen/Core>
#include "Config.h"
#include "StateVector.h"

namespace UKF {

/* Alias for the Eigen type from which MeasurementVector inherits. */
template <typename... Fields>
using MeasurementVectorBaseType = Eigen::Matrix<
	real_t,
	Eigen::Dynamic,
	1,
	0,
	detail::GetCompositeVectorDimension<Fields...>(),
	1>;

/*
Templated measurement vector class. A particular UKF implementation should
specialise this class in a similar way to the StateVector class, but with a
list of measurements which are intended to be provided to the filter.
*/
template <typename... Fields>
class MeasurementVector : public MeasurementVectorBaseType<typename Fields::type...> {
private:
    using Base = MeasurementVectorBaseType<typename Fields::type...>;
    using field_types = std::tuple<typename Fields::type...>;

public:
    /* Inherit Eigen::Matrix constructors and assignment operators. */
    using Base::Base;
    using Base::operator=;

    /* Get maximum size of measurement vector. */
    static constexpr std::size_t max_size() {
        return detail::GetCompositeVectorDimension<typename Fields::type...>();
    }
};

}

#endif