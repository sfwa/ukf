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

#ifndef CORE_H
#define CORE_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Config.h"
#include "StateVector.h"

namespace UKF {

/*
UKF core class. This class contains the intermediate values used in the filter
iteration, and the methods used to carry out a filter iteration itself.

The algorithm used here is based largely on the work of Edgar Kraft, described
in the paper "A Quaternion-based Unscented Kalman Filter for Orientation
Tracking", retrieved from:
http://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf
Comments in the code will occasionally refer to equations or sections in the
paper.

The attitude-related code makes use of the MRP method described in the paper
"Unscented Filtering for Spacecraft Attitude Estimation" by John L. Crassidis
and F. Landis Markley.
*/
template <typename StateVectorType, typename MeasurementVectorType>
class Core {
public:
    /* Top-level function used to carry out a filter step. */
    void iterate(const MeasurementVectorType &m) {
        /* Add process noise covariance to the state covariance and scale. */

        /* Sigma point distribution. */
        StateVectorType::SigmaPointDistribution sigma_points;
    }
};

private:
    /* State vector. */
    StateVectorType state;

    /* Covariance matrix. */
    StateVectorType::CovarianceMatrix covariance;

}

#endif