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
    Vector3r airflow, lift_axis, drag_axis, side_axis;
    real_t v, v_inv;

    airflow = attitude * (in.wind_velocity() - in.velocity());
    v = airflow.norm();

    if (v < UKF_DYNAMICS_MIN_V || v > UKF_DYNAMICS_MAX_V ||
            std::abs(yaw_rate) > UKF_DYNAMICS_MAX_RATE ||
            std::abs(roll_rate) > UKF_DYNAMICS_MAX_RATE ||
            std::abs(pitch_rate) > UKF_DYNAMICS_MAX_RATE) {
        /* Airflow too slow or too fast for any of this to work */
        return AccelerationVector();
    }

    v_inv = 1.0 / v;

    /*
    Lift is always perpendicular to airflow, drag is always parallel, and
    side is always towards the starboard wing of the aircraft.
    */
    drag_axis = airflow * v_inv;
    lift_axis = drag_axis.cross(Vector3r(0, 1, 0));
    side_axis = drag_axis.cross(lift_axis);

    /* Determine alpha and beta: alpha = atan(wz/wx), beta = atan(wy/|wxz|) */
    real_t alpha, beta, qbar, alpha2, beta2;
    qbar = RHO * v * v * 0.5;
    alpha = std::atan2(-airflow.z(), -airflow.x());
    beta = std::asin(airflow.y() * v_inv);

    if (std::abs(alpha) > M_PI * 0.63 || std::abs(beta) > M_PI * 0.63) {
        /* Alpha or beta completely out of range */
        return AccelerationVector();
    }

    alpha2 = alpha * alpha;
    beta2 = beta * beta;

    /*
    Set up a^4, a^3, a^2, a, 1 so that we can use dot product for the
    polynomial evaluation.
    */
    Eigen::Matrix<real_t, 5, 1> alpha_coeffs;
    alpha_coeffs << alpha2 * alpha2, alpha2 * alpha, alpha2, alpha, 1.0;

    /* Evaluate quartics in alpha to determine lift and drag */
    real_t lift = alpha_coeffs.dot(c_lift_alpha),
           drag = alpha_coeffs.dot(c_drag_alpha);

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
    if (motor_idx < control.rows()) {
        real_t ve = prop_cve * control[motor_idx], v0 = airflow.x();
        thrust = 0.5 * RHO * prop_area * (ve * ve - v0 * v0);
        if (thrust < 0.0) {
            /* Folding prop, so assume no drag */
            thrust = 0.0;
        }
    }

    /*
    Set up the side force coefficient vector, and calculate side force
    */
    Eigen::Matrix<real_t, 8 + UKF_CONTROL_DIM, 1> side_coeffs;
    side_coeffs.segment<8>(0) << alpha2, alpha, beta2, beta,
                                 alpha2 * beta, alpha * beta,
                                 yaw_rate, roll_rate;
    side_coeffs.segment<UKF_CONTROL_DIM>(8) << control;

    real_t side_force = side_coeffs.dot(c_side_force);

    /*
    Calculate pitch moment
    */
    Eigen::Matrix<real_t, 2 + UKF_CONTROL_DIM, 1> pitch_coeffs;
    pitch_coeffs.segment<2>(0) <<
        alpha, pitch_rate * pitch_rate * (pitch_rate < 0.0 ? -1.0 : 1.0);
    pitch_coeffs.segment<UKF_CONTROL_DIM>(2) << control;

    real_t pitch_moment = pitch_coeffs.dot(c_pitch_moment);

    /*
    Roll moment
    */
    Eigen::Matrix<real_t, 1 + UKF_CONTROL_DIM, 1> roll_coeffs;
    roll_coeffs[0] = roll_rate;
    roll_coeffs.segment<UKF_CONTROL_DIM>(1) << control;

    real_t roll_moment = roll_coeffs.dot(c_roll_moment);

    /*
    Yaw moment
    */
    Eigen::Matrix<real_t, 2 + UKF_CONTROL_DIM, 1> yaw_coeffs;
    yaw_coeffs.segment<2>(0) << beta, yaw_rate;
    yaw_coeffs.segment<UKF_CONTROL_DIM>(2) << control;

    real_t yaw_moment = yaw_coeffs.dot(c_yaw_moment);

    /*
    Sum and apply forces and moments
    */
    Vector3r sum_force, sum_torque;

    sum_force = qbar * (lift * lift_axis + drag * drag_axis +
                side_force * side_axis) + thrust * motor_thrust;
    sum_torque = qbar * Vector3r(roll_moment, pitch_moment, yaw_moment);

    /* Calculate linear acceleration (F / m) */
    AccelerationVector output;
    output.segment<3>(0) = sum_force * mass_inv +
                           attitude * Vector3r(0, 0, G_ACCEL);

    /* Calculate angular acceleration (tau / inertia tensor) */
    output.segment<3>(3) = inertia_tensor_inv * sum_torque;

    return output;
}
