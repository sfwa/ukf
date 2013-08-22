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

#ifndef UKF_H
#define UKF_H

#include "types.h"
#include "state.h"
#include "sensors.h"
#include "integrator.h"
#include "dynamics.h"

#define UKF_NUM_SIGMA (2*UKF_STATE_DIM + 1)

/*
Definitions for parameters of the Scaled Unscented Transform.
See "The Unscented Kalman Filter for Nonlinear Estimation" by Eric A. Wan and
Rudolph van der Merwe for details.

Note that alpha^2 here should be small for correct operation of the filter.
Most literature seems to quote about 1e-3 for alpha (1e-6 for alpha^2), but
the parameters suggested in "Gaussian Processes for State Space Models and
Change Point Detection" by Ryan Tuner (2011) provide a more stable filter.
*/
#define UKF_ALPHA_2 (1.0)
#define UKF_BETA (0.0)
#define UKF_KAPPA (3.0)
#define UKF_LAMBDA (UKF_ALPHA_2*(UKF_STATE_DIM + UKF_KAPPA) - UKF_STATE_DIM)
#define UKF_DIM_PLUS_LAMBDA (UKF_ALPHA_2*(UKF_STATE_DIM + UKF_KAPPA))

/*
Definitions for parameters used to calculated MRP vectors.
See the Markley paper for further details.
*/
#define UKF_MRP_A (1.0)
#define UKF_MRP_A_2 (UKF_MRP_A*UKF_MRP_A)
#define UKF_MRP_F (2.0*(UKF_MRP_A + 1))
#define UKF_MRP_F_2 (UKF_MRP_F*UKF_MRP_F)

/*
Definitions for sigma point weights. The naming convention follows that used
in in the paper given above.
*/
#define UKF_SIGMA_WM0 (UKF_LAMBDA/(UKF_DIM_PLUS_LAMBDA))
#define UKF_SIGMA_WC0 (UKF_SIGMA_WM0 + (1.0 - UKF_ALPHA_2 + UKF_BETA))
#define UKF_SIGMA_WMI (1.0/(2.0*(UKF_DIM_PLUS_LAMBDA)))
#define UKF_SIGMA_WCI (UKF_SIGMA_WMI)

/*
Unscented Kalman Filter object.
*/
class UnscentedKalmanFilter {
    /* State information and parameters. */
    State state;
    StateVectorCovariance state_covariance;
    ProcessCovariance process_noise_covariance;

    /*
    Reference to sensor model and pointer to dynamics model.
    A dynamics model isn't required – if not provided, the filter will assume
    constant angular and linear velocity for the process model.
    */
    SensorModel &sensor;
    DynamicsModel *dynamics;

#ifdef UKF_USE_EIGEN
    /* Intermediate variables used during filter iteration. */
    Eigen::Matrix<
        real_t,
        StateVector::RowsAtCompileTime,
        UKF_NUM_SIGMA> sigma_points;
    State apriori_mean;
    StateVectorCovariance apriori_covariance;

    Eigen::Matrix<
        real_t,
        UKF_STATE_DIM,
        UKF_NUM_SIGMA> w_prime;

    Eigen::Matrix<
        real_t,
        Eigen::Dynamic,
        UKF_NUM_SIGMA,
        0,
        UKF_MEASUREMENT_MAX_DIM> z_prime;

    MeasurementVector measurement_estimate_mean;
    Eigen::Matrix<
        real_t,
        Eigen::Dynamic,
        Eigen::Dynamic,
        0,
        UKF_MEASUREMENT_MAX_DIM,
        UKF_MEASUREMENT_MAX_DIM> measurement_estimate_covariance;

    MeasurementVector innovation;
    Eigen::Matrix<
        real_t,
        Eigen::Dynamic,
        Eigen::Dynamic,
        0,
        UKF_MEASUREMENT_MAX_DIM,
        UKF_MEASUREMENT_MAX_DIM> innovation_covariance;

    Eigen::Matrix<
        real_t,
        UKF_STATE_DIM,
        Eigen::Dynamic,
        0,
        UKF_STATE_DIM,
        UKF_MEASUREMENT_MAX_DIM> cross_correlation;

    Eigen::Matrix<
        real_t,
        UKF_STATE_DIM,
        Eigen::Dynamic,
        0,
        UKF_STATE_DIM,
        UKF_MEASUREMENT_MAX_DIM> kalman_gain;

    /* Integrator object, depends on selection in `config.h`. */
#ifdef UKF_INTEGRATOR_RK4
    IntegratorRK4 integrator;
#else
#ifdef UKF_INTEGRATOR_HEUN
    IntegratorHeun integrator;
#else
#ifdef UKF_INTEGRATOR_EULER
    IntegratorEuler integrator;
#endif
#endif
#endif

    void create_sigma_points();
    void apriori_estimate(real_t dt, ControlVector c);
    void measurement_estimate();
    void calculate_innovation();
    void calculate_kalman_gain();
    void aposteriori_estimate();
#endif

public:
    UnscentedKalmanFilter(SensorModel &sensor_model);
    const State& get_state() const { return state; }
    void set_state(const State &in) { state = in; }
    const StateVectorCovariance& get_state_covariance() const {
        return state_covariance;
    }
    SensorModel* get_sensor_model() { return &sensor; }
    void set_process_noise(ProcessCovariance in) {
        process_noise_covariance = in;
    }
    void set_dynamics_model(DynamicsModel *in) { dynamics = in; }
    void iterate(real_t dt, ControlVector c);
};

#endif
