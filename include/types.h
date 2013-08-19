/*
Copyright (C) 2013 Daniel Dyer

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

#ifndef TYPES_H_
#define TYPES_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "config.h"

#ifdef SINGLE_PRECISION
typedef float real_t;
#endif

#ifdef DOUBLE_PRECISION
typedef double real_t;
#endif

typedef Eigen::Matrix<real_t, 3, 1> Vector3r;
typedef Eigen::Quaternion<real_t> Quaternionr;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1> VectorXr;

#define G_ACCEL (9.80665)
#define RHO (1.225)

#endif
