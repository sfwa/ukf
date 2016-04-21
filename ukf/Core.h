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
#include "Types.h"
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
    /* Aliases needed during filter iteration. */
    template <S, M>
    using CrossCorrelation = Eigen::Matrix<
        real_t,
        S::covariance_size(),
        M::RowsAtCompileTime,
        0,
        S::covariance_size(),
        M::MaxRowsAtCompileTime>;

    /* Top-level function used to carry out a filter step. */
    void step(real_t rime_step, const ControlVectorType &u, const MeasurementVectorType &z) {
        a_priori_step(time_step, u);
        innovation_step(z);
        a_posteriori_step();
    }

    /*
    The a priori step calculates the sigma point distribution using the
    previous state covariance and previous state estimate, and then
    propagates all sigma points through the process model.
    From the transformed sigma point distribution, the a priori mean and
    covariance are calculated.
    */
    void a_priori_step(real_t time_step, const ControlVectorType &u) {
        /*
        Add process noise covariance to the state covariance and calculate
        the sigma point distribution.
        */
        sigma_points = state.calculate_sigma_point_distribution(
            covariance + /* FIXME: Insert process noise covariance here. */);

        /* Propagate the sigma points through the process model. */
        /* FIXME: Insert process model here. */

        /* Calculate the a priori estimate mean, deltas and covariance. */
        a_priori_mean = StateVectorType::calculate_sigma_point_mean(sigma_points);
        w_prime = a_priori_mean.calculate_sigma_point_deltas(sigma_points);
        a_priori_covariance = StateVectorType::calculate_sigma_point_covariance(w_prime);
    }

    /*
    In the innovation step, the a priori sigma point distribution is further
    propagated using the measurement model. For the measurement sigma point
    distribution, the estimated measurement mean and covariance are
    calculated. A measurement vector is then supplied, and the innovation
    and innovation covariance are calculated.
    */
    void innovation_step(const MeasurementVectorType &z) {
        /* Propagate the sigma points through the measurement model. */
        MeasurementVectorType::SigmaPointDistribution<StateVectorType::num_sigma> measurement_sigma_points =
            z.calculate_sigma_point_distribution<StateVectorType>(sigma_points);

        /* Calculate the measurement prediction mean, deltas and covariance. */
        MeasurementVectorType z_pred = z.calculate_sigma_point_mean<StateVectorType>(measurement_sigma_points);
        z_prime = z_pred.calculate_sigma_point_deltas<StateVectorType>(measurement_sigma_points);
        innovation_covariance = z_pred.calculate_sigma_point_covariance<StateVectorType>(z_prime);

        /*
        Calculate innovation and innovation covariance. Innovation is simply
        the difference between the measurement and the predicted measurement,
        and innovation covariance is the sum of the predicted measurement
        covariance and the measurement covariance.
        See equations 44 and 45 from the Kraft paper for details.
        */
        innovation = z - z_pred;
        innovation_covariance += /* FIXME: Insert measurement covariance here. */;
    }

    /*
    In the a posteriori step, the cross-correlation matrix and innovation
    covariance inverse are calculated and used to calculate the Kalman gain.
    This step is split from the innovation step to allow the user to carry
    out innovation consistency checks and identify failed sensors, if
    desired.
    */
    void a_posteriori_step() {
        /*
        Calculate the cross-correlation matrix described in equations 70 and
        71 from from the Kraft paper.
        */
        CrossCorrelation<StateVectorType, MeasurementVectorType> cross_correlation =
            Parameters::Sigma_WC0<S> * (w_prime.col(0) * z_prime.col(0).transpose());
        for(int i = 1; i < S::num_sigma(); i++) {
            cross_correlation += Parameters::Sigma_WCI<S> * (w_prime.col(i) * z_prime.col(i).transpose());
        }

        /*
        Calculate the Kalman gain as described in equation 72 from the Kraft
        paper.
        */
        CrossCorrelation<StateVectorType, MeasurementVectorType> kalman_gain =
            cross_correlation * innovation_covariance.inverse();

        /* Do the final update step. */
    }
};

private:
    StateVectorType state;
    StateVectorType::CovarianceMatrix covariance;
    StateVectorType::SigmaPointDistribution sigma_points;

    StateVectorType a_priori_mean;
    StateVectorType::CovarianceMatrix a_priori_covariance;
    StateVectorType::SigmaPointDeltas w_prime;

    MeasurementVectorType innovation;
    MeasurementVectorType::CovarianceMatrix innovation_covariance;
    MeasurementVectorType::SigmaPointDeltas<StateVectorType::num_sigma> z_prime;
}

#endif