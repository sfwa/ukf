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
The algorithm used here is based largely on the work of Edgar Kraft, described
in the paper "A Quaternion-based Unscented Kalman Filter for Orientation
Tracking", retrieved from:
http://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf
Comments in the code will occasionally refer to equations or sections in the
paper.
The attitude-related code makes use of the MRP method described in the paper
"Unscented Filtering for Spacecraft Attitude Estimation" by John L. Crassidis
and F. Landis Markley.

This implementation is split into two files, ukf.cpp and ukf-estimates.cpp,
in order to work around memory limitations in the TI CCS compiler.
*/

UnscentedKalmanFilter::UnscentedKalmanFilter(SensorModel &sensor_model) :
sensor(sensor_model) {
#ifdef INTEGRATOR_RK4
    integrator = IntegratorRK4();
#else
#ifdef INTEGRATOR_HEUN
    integrator = IntegratorHeun();
#else
#ifdef INTEGRATOR_EULER
    integrator = IntegratorEuler();
#endif
#endif
#endif

    state = StateVector::Zero();
    state.attitude() << 0, 0, 0, 1;

    process_noise_covariance <<
        1e-6, 1e-6, 1e-6,
        1e-6, 1e-6, 1e-6,
        1e-6, 1e-6, 1e-6,
        1e-6, 1e-6, 1e-6,
        1e-6, 1e-6, 1e-6,
        1e-6, 1e-6, 1e-6,
        1e-6, 1e-6, 1e-6,
        1e-6, 1e-6, 1e-6;

    state_covariance = StateVectorCovariance::Zero();
    state_covariance.diagonal() <<
        M_PI * M_PI * 0.0625, M_PI * M_PI * 0.0625, 1000,
        50, 50, 50,
        10, 10, 10,
        M_PI * 0.25, M_PI * 0.25, M_PI * 0.25,
        2, 2, 2,
        5, 5, 5,
        20, 20, 20,
        0, 0, 0;
    dynamics = NULL;
}

/* Follows section 3.1, however we're using the scaled unscented transform. */
void UnscentedKalmanFilter::create_sigma_points() {
    Eigen::Matrix<real_t, UKF_DIM, UKF_DIM> S;

    //AssertNormalized(Quaternionr(state.attitude()));

    /* Add the process noise before calculating the square root. */
    state_covariance.diagonal() += process_noise_covariance;
    state_covariance *= UKF_DIM_PLUS_LAMBDA;

    //AssertPositiveDefinite(state_covariance);

    /* Compute the 'square root' of the covariance matrix. */
    S = state_covariance.llt().matrixL();

    /* Create centre sigma point. */
    sigma_points.col(0) = state;

    /* For each column in S, create the two complementary sigma points. */
    for(uint8_t i = 0; i < UKF_DIM; i++) {
        /*
        Construct error quaternions using the MRP method, equation 34 from the
        Markley paper.
        */
        Vector3r d_p = S.col(i).segment<3>(9);
        real_t x_2 = d_p.squaredNorm();
        real_t err_w = (-UKF_MRP_A * x_2 +
            UKF_MRP_F * std::sqrt(UKF_MRP_F_2 + (1.0 - UKF_MRP_A_2) * x_2)) /
            (UKF_MRP_F_2 + x_2);
        Vector3r err_xyz = ((1.0 / UKF_MRP_F) * (UKF_MRP_A + err_w)) * d_p;
        Quaternionr noise;
        noise.vec() = err_xyz;
        noise.w() = err_w;

        Quaternionr temp;

        /* Create positive sigma point. */
        sigma_points.col(i+1).segment<9>(0) =
            state.segment<9>(0) + S.col(i).segment<9>(0);
        temp = noise * Quaternionr(state.attitude());
        //AssertNormalized(temp);
        sigma_points.col(i+1).segment<4>(9) << temp.vec(), temp.w();
        sigma_points.col(i+1).segment<12>(13) =
            state.segment<12>(13) + S.col(i).segment<12>(12);

        /* Create negative sigma point. */
        sigma_points.col(i+1 + UKF_DIM).segment<9>(0) =
            state.segment<9>(0) - S.col(i).segment<9>(0);
        temp = noise.conjugate() * Quaternionr(state.attitude());
        //AssertNormalized(temp);
        sigma_points.col(i+1 + UKF_DIM).segment<4>(9) << temp.vec(), temp.w();
        sigma_points.col(i+1 + UKF_DIM).segment<12>(13) =
            state.segment<12>(13) - S.col(i).segment<12>(12);
    }
}

/* Calculating the innovation follows equations 44 and 45. */
void UnscentedKalmanFilter::calculate_innovation() {
    /*
    The innovation is just the predicted measurement subtracted from the
    actual measurement. The actual measurement comes from the sensor model.
    */
    innovation = sensor.collate() - measurement_estimate_mean;

    /*
    Innovation covariance is just the measurement estimate covariance added to
    the measurement noise covariance.
    */
    innovation_covariance = measurement_estimate_covariance;
    innovation_covariance.diagonal() += sensor.get_covariance();
}

/* This just does all the steps in order to complete a filter iteration. */
void UnscentedKalmanFilter::iterate(real_t dt, ControlVector c) {
    create_sigma_points();
    apriori_estimate(dt, c);

    /* Check to make sure there's a need to run the measurement step. */
    if(sensor.size() > 0) {
        measurement_estimate();
        calculate_innovation();
        calculate_kalman_gain();
        aposteriori_estimate();
    } else {
        state = sigma_points.col(0);
        state_covariance = apriori_covariance;
    }

    if(state.attitude()(3) < 0) {
        state.attitude() *= -1.0;
    }
    state.attitude().normalize();
}
