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

    if (v < UKF_DYNAMICS_MIN_V) {
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
    real_t alpha, beta, qbar, alpha2, beta2;
    qbar = RHO * v * v * (real_t)0.5;
    alpha = std::atan2(-airflow.z(), -airflow.x());
    beta = std::asin(airflow.y() * v_inv);

    alpha2 = alpha * alpha;
    beta2 = beta * beta;

    /*
    Set up a^4, a^3, a^2, a, 1 so that we can use dot product for the
    polynomial evaluation.
    */
    Vector5r alpha_coeffs;
    alpha_coeffs << alpha2 * alpha2, alpha2 * alpha, alpha2, alpha, 1.0;

    /* Evaluate quartics in alpha to determine lift and drag */
    real_t lift = alpha_coeffs.dot(c_lift_alpha),
           drag = alpha_coeffs.dot(c_drag_alpha);

    /*
    Assume lift is somewhat well-behaved for alpha in the range [-0.25, 0.25].
    If outside that range, clamp it to 0 so that the polynomial doesn't have
    to model the full possible range.
    */
    if (std::abs(alpha) > 1.0) {
        lift = 0.0;
    } else if (alpha > 0.25) {
        lift = std::max(lift, 0.0);
    } else if (alpha < -0.25) {
        lift = std::min(lift, 0.0);
    }

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
    Set up the side force coefficient vector, and calculate side force
    */
    Vector4r side_coeffs;
    side_coeffs << beta2, beta, yaw_rate, roll_rate;

    real_t side_force = side_coeffs.dot(c_side_force) +
        control.dot(c_side_force_control);

    /*
    Calculate pitch moment
    */
    Vector2r pitch_coeffs;
    pitch_coeffs << alpha,
        pitch_rate * (pitch_rate < 0.0 ? -pitch_rate : pitch_rate);

    real_t pitch_moment = pitch_coeffs.dot(c_pitch_moment) +
        control.dot(c_pitch_moment_control);

    /*
    Roll moment
    */
    Vector1r roll_coeffs;
    roll_coeffs << roll_rate;
    real_t roll_moment = roll_coeffs.dot(c_roll_moment) +
        control.dot(c_roll_moment_control);

    /*
    Yaw moment
    */
    Vector2r yaw_coeffs;
    yaw_coeffs << beta, yaw_rate;
    real_t yaw_moment = yaw_coeffs.dot(c_yaw_moment) +
        control.dot(c_yaw_moment_control);

    /*
    Sum and apply forces and moments
    */
    Vector3r sum_force, sum_torque;

    sum_force = qbar * (
            (lift * std::cos(alpha) + drag * std::sin(alpha)) * Vector3r(0, 0, -1) +
            (lift * std::sin(alpha) - drag * std::cos(alpha) - side_force * std::sin(beta)) * Vector3r(1, 0, 0) +
            side_force * std::cos(beta) * Vector3r(0, 1, 0)
        ) +
        thrust * Vector3r(1, 0, 0);
    sum_torque = qbar * Vector3r(roll_moment, pitch_moment, yaw_moment);

    /* Calculate linear acceleration (F / m) */
    AccelerationVector output;
    output.segment<3>(0) = sum_force * mass_inv +
                           attitude * Vector3r(0, 0, G_ACCEL);

    /* Calculate angular acceleration (tau / inertia tensor) */
    output.segment<3>(3) = inertia_tensor_inv * sum_torque;

    return output;
}

/*
Runs a dynamics model with hard-coded coefficients for the X8.
*/
AccelerationVector X8DynamicsModel::evaluate(
const State &in, const ControlVector &control) const {
    /* Cache state data for convenience */
    Quaternionr attitude = Quaternionr(in.attitude());
    real_t yaw_rate = in.angular_velocity()[2],
           pitch_rate = in.angular_velocity()[1],
           roll_rate = in.angular_velocity()[0];

    /* External axes */
    Vector3r airflow;
    real_t v, v_inv, horizontal_v2, vertical_v2, vertical_v, vertical_v_inv;

    airflow = attitude * (in.wind_velocity() - in.velocity());
    v = airflow.norm();
    horizontal_v2 = airflow.y() * airflow.y() + airflow.x() * airflow.x();
    vertical_v2 = airflow.z() * airflow.z() + airflow.x() * airflow.x();

    v_inv = (real_t)1.0 / std::max(v, 1.0);
    vertical_v = std::sqrt(vertical_v2);
    vertical_v_inv = (real_t)1.0 / std::max(vertical_v, 1.0);

    /* Determine alpha and beta: alpha = atan(wz/wx), beta = atan(wy/|wxz|) */
    real_t alpha, sin_alpha, cos_alpha, sin_beta, cos_beta, a2, sin_cos_alpha;
    alpha = std::atan2(-airflow.z(), -airflow.x());

    sin_alpha = -airflow.z() * vertical_v_inv;
    cos_alpha = -airflow.x() * vertical_v_inv;
    sin_cos_alpha = sin_alpha * cos_alpha;
    sin_beta = airflow.y() * v_inv;
    cos_beta = vertical_v * v_inv;

    a2 = alpha * alpha;

    real_t lift, drag, side_force, roll_moment, pitch_moment, yaw_moment;
    lift = -5 * a2 * alpha + a2 + 2.5 * alpha + 0.12;
    if (alpha < 0.0) {
        lift = std::min(lift, 0.8 * sin_cos_alpha);
    } else {
        lift = std::max(lift, 0.8 * sin_cos_alpha);
    }

    drag = 0.05 + 0.7 * sin_alpha * sin_alpha;
    side_force = 0.3 * sin_beta * cos_beta;

    pitch_moment = 0.001 - 0.1 * sin_cos_alpha - 0.003 * pitch_rate -
                   0.01 * control[1] - 0.01 * control[2];
    roll_moment = -0.03 * sin_beta - 0.015 * roll_rate +
                  0.025 * control[1] - 0.025 * control[2];
    yaw_moment = -0.02 * sin_beta - 0.05 * yaw_rate -
                 0.01 * std::abs(control[1]) + 0.01 * std::abs(control[2]);

    /*
    Determine motor thrust and torque.
    */
    real_t thrust, ve = 0.0025 * control[0], v0 = airflow.x();
    thrust = (real_t)0.5 * RHO * 0.025 * (ve * ve - v0 * v0);
    if (thrust < 0.0) {
        /* Folding prop, so assume no drag */
        thrust = 0.0;
    }

    /*
    Sum and apply forces and moments
    */
    Vector3r sum_force, sum_torque;
    real_t qbar = RHO * horizontal_v2 * (real_t)0.5;
    sum_force << thrust + qbar * (lift * sin_alpha - drag * cos_alpha - side_force * sin_beta),
                 qbar * side_force * cos_beta,
                 -qbar * (lift * cos_alpha + drag * sin_alpha);
    sum_torque = qbar * Vector3r(roll_moment, pitch_moment, yaw_moment);

    /* Calculate linear acceleration (F / m) */
    AccelerationVector output;
    output.segment<3>(0) = sum_force * mass_inv +
                           attitude * Vector3r(0, 0, G_ACCEL);

    /* Calculate angular acceleration (tau / inertia tensor) */
    /*output.segment<3>(3) = inertia_tensor_inv * sum_torque;*/
    output.segment<3>(3) <<
        qbar * (3.364222 * roll_moment + 0.27744448 * yaw_moment),
        qbar * 5.8823528 * pitch_moment,
        qbar * (0.27744448 * roll_moment + 2.4920163 * yaw_moment);

    return output;
}
