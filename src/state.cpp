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
#include "state.h"

/*
Runs the kinematics model on the state vector and returns a vector with the
derivative of each components (except for the accelerations, which must be
calculated directly using a dynamics model).
Contents are as follows:
    - Rate of change in position:
        (3-vector, latitude (rad/s), longitude (rad/s), altitude (m/s))
    - Rate of change in linear velocity (3-vector, m/s^2, NED frame)
    - Rate of change in linear acceleration (3-vector, m/s^3, body frame)
    - Rate of change in attitude (quaternion (x, y, z, w), 1/s, body frame)
    - Rate of change in angular velocity (3-vector, rad/s^2, body frame)
    - Rate of change in angular acceleration (3-vector, rad/s^3, body frame)
    - Rate of change in wind velocity (3-vector, m/s, NED frame)
    - Rate of change in gyro bias (3-vector, rad/s^2, body frame)
*/
const StateVectorDerivative State::model() {
    StateVectorDerivative output;

    /* Calculate the normal and meridional radii of curvature. */
    Vector3r pos = position();
    real_t tempA = WGS84_A*std::cos(pos[0]), tempB = WGS84_B*std::sin(pos[0]),
           temp = tempA * tempA + tempB * tempB,
           temp_sqrt = std::sqrt(temp);
    real_t M = WGS84_AB2 / (temp_sqrt * temp);
    real_t N = WGS84_A2 / temp_sqrt;

    /*
    Calculate change in position. Using the small angle approximation, this
    becomes very simple – no trig required for latitude derivative, and one
    cosine function for longitude derivative.
    */
    Vector3r vel = velocity();
    output.segment<3>(0) <<
        vel[0] / (M + pos[2]),
        (vel[1] / (N + pos[2])) * std::cos(pos[0]),
        -vel[2];

    /* Calculate change in velocity. */
    Quaternionr attitude_q = Quaternionr(attitude());
    output.segment<3>(3) = attitude_q.conjugate() * acceleration();

    /* Change in linear acceleration is zero. */
    output.segment<3>(6) << 0, 0, 0;

    /* Calculate change in attitude. */
    Eigen::Matrix<real_t, 4, 1> omega_q;
    omega_q << angular_velocity(), 0;

    attitude_q = Quaternionr(omega_q).conjugate() * attitude_q;
    output.segment<4>(9) << attitude_q.vec(), attitude_q.w();
    output.segment<4>(9) *= 0.5;

    /* Calculate change in angular velocity (just angular acceleration). */
    output.segment<3>(13) = angular_acceleration();

    /* Change in angular acceleration is zero. */
    output.segment<3>(16) << 0, 0, 0;

    /* Change in wind velocity is zero. */
    output.segment<3>(19) << 0, 0, 0;

    /* Change in gyro bias is zero. */
    output.segment<3>(22) << 0, 0, 0;

    return output;
}
