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

#include "Config.h"

/*
Integrator base class. The public interface is via the integrate() method,
which takes a template parameter that at a minimum must support addition,
subtraction and scalar multiplication. It also must have a public method
"derivative" which takes no arguments and returns a type which can be added
to the template parameter and also supports scalar multiplication.
*/
template <typename Derived>
class Integrator {
public:
    template <typename S>
    static S integrate(real_t delta, const S &state) {
        return Derived::integrate(delta, state);
    }
};

class IntegratorRK4: Integrator<IntegratorRK4> {
public:
    template <typename S>
    static S integrate(real_t delta, const S &state) {
        S a = state.model();
        S b = static_cast<S>(state + 0.5f * delta * a).model();
        S c = static_cast<S>(state + 0.5f * delta * b).model();
        S d = static_cast<S>(state + delta * c).model();
        return state + (delta / 6.0f) * (a + (b * 2.0f) + (c * 2.0f) + d);
    }
};

class IntegratorHeun: Integrator<IntegratorHeun> {
public:
    template <typename S>
    static S integrate(real_t delta, const S &state) {
        S initial = state.model();
        S predictor = state + delta * initial;
        return state + (delta * 0.5f) * (initial + predictor.model());
    }
};

class IntegratorEuler: Integrator<IntegratorEuler> {
public:
    template <typename S>
    static S integrate(real_t delta, const S &state) {
        return state + delta * state.model();
    }
};

#endif
