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

#include <cmath>

#include "types.h"
#include "dynamics.h"

DynamicsModel::~DynamicsModel() {}

/*
Runs the dynamics model and calculates the expected linear and angular
accelerations for this timestep. Eventually this will need to take control
surface inputs as a parameter.

For now, it calculates the expected acceleration as the centripetal
acceleration expected to be felt for the current velocity and angular
velocity, and doesn't attempt to calculate angular acceleration at all.
*/
AccelerationVector CentripetalModel::evaluate(
const State &in, const ControlVector &control) const {
    #pragma unused(control)

    AccelerationVector output;

    /* First convert velocity to body frame. */
    Eigen::Matrix<real_t, 3, 1> velocity_body;
    velocity_body = Quaternionr(in.attitude()) * in.velocity();

    /* Calculate centripetal acceleration. */
    output.segment<3>(0) = in.angular_velocity().cross(velocity_body);

    /* Clear angular acceleration. */
    output.segment<3>(3) << 0, 0, 0;

    return output;
}

/*
Runs a basic fixed-wing flight dynamics model, calculating:
- Lift and drag as a quartic function of alpha
- Side force (see header for details)
- Thrust as a function of prop area, prop exit velocity coefficient, forward
  velocity, and throttle
- Pitch moment as a linear function of alpha, a quadratic function of pitch
  rate, and a linear function of control surface position
- Roll moment as a linear function of roll rate and control surface position
- Yaw moment as a linear function of alpha, yaw rate, and control surface
  position
*/
AccelerationVector FixedWingFlightDynamicsModel::evaluate(
const State &in, const ControlVector &control) const {
    /* Cache state data for convenience */
    Quaternionr attitude = Quaternionr(in.attitude());
    real_t yaw_rate = in.angular_velocity()[2],
           pitch_rate = in.angular_velocity()[1],
           roll_rate = in.angular_velocity()[0];

    /* External axes */
    Vector3r airflow;
    real_t v, v_inv, horizontal_v2, vertical_v2;

    airflow = attitude * (in.wind_velocity() - in.velocity());
    v = airflow.norm();
    horizontal_v2 = airflow.y() * airflow.y() + airflow.x() * airflow.x();
    vertical_v2 = airflow.z() * airflow.z() + airflow.x() * airflow.x();

    if (horizontal_v2 < UKF_DYNAMICS_MIN_V) {
        /*
        Airflow too slow for any of this to work; just return acceleration
        due to gravity
        */
        AccelerationVector output;
        output.segment<3>(0) = attitude * Vector3r(0, 0, G_ACCEL);

        /* No angular acceleration */
        output.segment<3>(3) << 0, 0, 0;
        return output;
    }

    v_inv = (real_t)1.0 / v;

    /* Determine alpha and beta: alpha = atan(wz/wx), beta = atan(wy/|wxz|) */
    real_t alpha, beta, qbar, sin_alpha, cos_alpha, sin_beta, cos_beta,
           lift, drag, side_force, roll_moment, pitch_moment, yaw_moment,
           a2;
    qbar = RHO * horizontal_v2 * (real_t)0.5;
    alpha = std::atan2(-airflow.z(), -airflow.x());
    beta = std::asin(airflow.y() * v_inv);

    sin_alpha = std::sin(alpha);
    cos_alpha = std::cos(alpha);
    sin_beta = std::sin(beta);
    cos_beta = std::cos(beta);

    a2 = alpha * alpha;

    lift = -5 * a2 * alpha + a2 + 2.5 * alpha + 0.12;
    if (alpha < 0.0) {
        lift = std::min(lift, 0.8 * sin_alpha * cos_alpha);
    } else {
        lift = std::max(lift, 0.8 * sin_alpha * cos_alpha);
    }

    drag = 0.05 + 0.7 * sin_alpha * sin_alpha + 0.2 * sin_beta * sin_beta;
    side_force = 0.3 * sin_beta;

    pitch_moment = 0.001 - 0.1 * sin_alpha * cos_alpha - 0.003 * pitch_rate -
                   0.01 * control[1] - 0.01 * control[2];
    roll_moment = -0.015 * roll_rate + 0.025 * control[1] - 0.025 * control[2];
    yaw_moment = -0.02 * sin_beta - 0.05 * yaw_rate -
                 0.01 * std::abs(control[1]) + 0.01 * std::abs(control[2]);

    /*
    Determine motor thrust and torque:
    https://www.grc.nasa.gov/WWW/Wright/airplane/propth.html

    Propeller thrust =
    F = 0.5 * rho * A * (Ve^2 - V0^2)
    where A = propeller disc area, Ve = exit velocity, V0 = air velocity

    In this formulation,
    F = 0.5 * rho * A * ((k * rpm)^2 - V0^2)

    Presumably the relationship between thrust and torque on airframe is
    constant? If so, they're related by prop_ct.
    */
    real_t thrust = 0.0;
    if (motor_idx < UKF_CONTROL_DIM) {
        real_t ve = prop_cve * control[motor_idx], v0 = airflow.x();
        thrust = (real_t)0.5 * RHO * prop_area * (ve * ve - v0 * v0);
        if (thrust < 0.0) {
            /* Folding prop, so assume no drag */
            thrust = 0.0;
        }
    }

    /*
    Sum and apply forces and moments
    */
    Vector3r sum_force, sum_torque;
    sum_force << thrust + qbar * (lift * sin_alpha - drag * cos_alpha - side_force * sin_beta),
                 qbar * side_force * cos_beta,
                 -qbar * (lift * cos_alpha + drag * sin_alpha);
    sum_torque = qbar * Vector3r(roll_moment, pitch_moment, yaw_moment);

    /* Calculate linear acceleration (F / m) */
    AccelerationVector output;
    output.segment<3>(0) = sum_force * mass_inv +
                           attitude * Vector3r(0, 0, G_ACCEL);

    /* Calculate angular acceleration (tau / inertia tensor) */
    output.segment<3>(3) = inertia_tensor_inv * sum_torque;

    return output;
}
