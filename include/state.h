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

#ifndef STATE_H
#define STATE_H

#include "types.h"

/*
Definition for filter state vector.
Contents are as follows:
    - Position (3-vector, latitude (rad), longitude (rad), altitude (m))
    - Linear Velocity (3-vector, m/s, NED frame)
    - Linear Acceleration (3-vector, m/s^2, body frame)
    - Attitude (quaternion (x, y, z, w), describes rotation from local NED
      frame to body frame.)
    - Angular Velocity (3-vector, rad/s, body frame)
    - Angular Acceleration (3-vector, rad/s^2, body frame)
    - Wind Velocity (3-vector, m/s, NED frame)
    - Gyro bias (3-vector, rad/s, body frame)
*/
class State: public StateVector {
public:
    State(void) : StateVector() {}

    template<typename OtherDerived>
    State(const Eigen::MatrixBase<OtherDerived>& other)
        : StateVector(other) { }

    template<typename OtherDerived>
    State & operator= (const Eigen::MatrixBase<OtherDerived>& other)
    {
        StateVector::operator=(other);
        return *this;
    }

    const StateVectorDerivative model();

    /* Read-only accessors */
    const Vector3r position() const {
        return segment<3>(0);
    }

    const Vector3r velocity() const {
        return segment<3>(3);
    }

    const Vector3r acceleration() const {
        return segment<3>(6);
    }

    const Vector4r attitude() const {
        return segment<4>(9);
    }

    const Vector3r angular_velocity() const {
        return segment<3>(13);
    }

    const Vector3r angular_acceleration() const {
        return segment<3>(16);
    }

    const Vector3r wind_velocity() const {
        return segment<3>(19);
    }

    const Vector3r gyro_bias() const {
        return segment<3>(22);
    }

    /* Mutable accessors */
    Eigen::VectorBlock<StateVector, 3> position() {
        return segment<3>(0);
    }

    Eigen::VectorBlock<StateVector, 3> velocity() {
        return segment<3>(3);
    }

    Eigen::VectorBlock<StateVector, 3> acceleration() {
        return segment<3>(6);
    }

    Eigen::VectorBlock<StateVector, 4> attitude() {
        return segment<4>(9);
    }

    Eigen::VectorBlock<StateVector, 3> angular_velocity() {
        return segment<3>(13);
    }

    Eigen::VectorBlock<StateVector, 3> angular_acceleration() {
        return segment<3>(16);
    }

    Eigen::VectorBlock<StateVector, 3> wind_velocity() {
        return segment<3>(19);
    }

    Eigen::VectorBlock<StateVector, 3> gyro_bias() {
        return segment<3>(22);
    }
};

#endif
