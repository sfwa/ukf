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

#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Core>
#include "UKF/Config.h"

namespace UKF {

template <int Rows, int Columns>
using Matrix = Eigen::Matrix<real_t, Rows, Columns>;

template <int Rows, int Columns>
using MatrixDynamic = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, 0, Rows, Columns>;

template <int Rows, int Columns>
using Array = Eigen::Array<real_t, Rows, Columns>;

template <int Length>
using Vector = Eigen::Matrix<real_t, Length, 1>;

template <int Length>
using VectorDynamic = Eigen::Matrix<real_t, Eigen::Dynamic, 1, 0, Length, 1>;

using Quaternion = Eigen::Quaternion<real_t>;

class FieldVector : public Vector<3> {
public:
	/* Inherit Eigen::Matrix constructors and assignment operators. */
    using Base = Vector<3>;
    using Base::Base;
    using Base::operator=;
};

}

#endif