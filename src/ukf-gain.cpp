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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

#include "types.h"
#include "ukf.h"
#include "state.h"
#include "dynamics.h"
#include "debug.h"

/*
This implementation is split into three files, ukf.cpp, ukf-estimates.cpp and
ukf-gain.cpp in order to work around memory limitations in the TI CCS
compiler.

See ukf.cpp for details.
*/

/* The Kalman gain is calculated using sections 3.5.3 and 3.6. */
void UnscentedKalmanFilter::calculate_kalman_gain() {
    /*
    First calculate the cross-correlation matrix between the a priori estimate
    and the predicted measurement as shown in equation 70.
    */
    cross_correlation =
        UKF_SIGMA_WC0 * (w_prime.col(0) * z_prime.col(0).transpose());
    for(uint32_t i = 1; i < UKF_NUM_SIGMA; i++) {
        cross_correlation +=
            UKF_SIGMA_WCI * (w_prime.col(i) * z_prime.col(i).transpose());
    }

    /* Next, calculate the Kalman gain as described in equation 72. */
    kalman_gain = cross_correlation * innovation_covariance.inverse();
}
