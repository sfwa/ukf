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

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "UKF/Config.h"

namespace UKF {

/* Fourth-order integrator. */
class IntegratorRK4 {
public:
    template <typename S, typename... U>
    static S integrate(real_t delta, const S& state, const U&... input) {
        S a = state.derivative(input...);
        S b = S(state + real_t(0.5) * delta * a).derivative(input...);
        S c = S(state + real_t(0.5) * delta * b).derivative(input...);
        S d = S(state + delta * c).derivative(input...);
        return state + (delta / real_t(6.0)) * (a + (b * real_t(2.0)) + (c * real_t(2.0)) + d);
    }
};

/* Second-order integrator. */
class IntegratorHeun {
public:
    template <typename S, typename... U>
    static S integrate(real_t delta, const S& state, const U&... input) {
        S initial = state.derivative(input...);
        S predictor = state + delta * initial;
        return state + (delta * real_t(0.5)) * (initial + predictor.derivative(input...));
    }
};

/* First-order integrator. */
class IntegratorEuler {
    public:
    template <typename S, typename... U>
    static S integrate(real_t delta, const S& state, const U&... input) {
        return state + delta * state.derivative(input...);
    }
};

}

#endif
