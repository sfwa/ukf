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

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <stdint.h>
#include <cmath>

#include "types.h"
#include "state.h"

/* Disable dynamics model if velocity is less than 1m/s */
#define UKF_DYNAMICS_MIN_V 1.0
/* Airframe minimums */
#define UKF_AIRFRAME_MIN_MASS 0.1
#define UKF_AIRFRAME_MIN_MOMENT 1e-6

/*
Dynamics model base class. The public interface is via the evaluate() method,
which returns a 6-dimensional column vector containing linear acceleration
and angular acceleration. In the future, this will need to accept a control
input vector as a parameter.
*/
class DynamicsModel {
public:
    virtual ~DynamicsModel();
    virtual AccelerationVector evaluate(
    const State &in, const ControlVector &control) const = 0;
};

/*
Just a basic dynamics model which predicts only acceleration using a
centripetal acceleration model.
*/
class CentripetalModel: public DynamicsModel {
public:
    AccelerationVector evaluate(
    const State &in, const ControlVector &control) const;
};

/*
Dynamics model tuned for the X8
*/
class X8DynamicsModel: public DynamicsModel {
    /* Use reciprocal of mass for performance */
    real_t mass_inv;

    /* Store inverse inertia tensor for performance */
    Matrix3x3r inertia_tensor_inv;

public:
    X8DynamicsModel(void) {
        mass_inv = (real_t)1.0 / 3.8;

        Matrix3x3r inertia_tensor;
        inertia_tensor <<
            3.0e-1, 0, -0.334e-1,
            0, 1.7e-1, 0,
            -0.334e-1, 0, 4.05e-1;
        inertia_tensor_inv = inertia_tensor.inverse();
    }

    AccelerationVector evaluate(
    const State &in, const ControlVector &control) const;
};

/*
A custom dynamics model wrapper, which calls a pointer to a function accepting
state vector, control vector and output vector parameters (as C arrays).
*/
class CustomDynamicsModel: public DynamicsModel {
    ModelFunction dynamics_model;

public:
    CustomDynamicsModel(void) {
        dynamics_model = NULL;
    }

    void set_function(ModelFunction func) {
        assert(func);

        dynamics_model = func;
    }

    AccelerationVector evaluate(
    const State &in, const ControlVector &control) const;
};

#endif
