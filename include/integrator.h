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

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

/*
Integrator base class. The public interface is via the integrate() method,
which takes a template parameter that at a minimum must support addition,
subtraction and scalar multiplication. It also must have a public method
"model" which takes no arguments and returns a type which can be added to
the template parameter and also supports scalar multiplication.
*/
template<typename Derived>
class Integrator {
public:
    template<typename StateModel>
    StateModel integrate(StateModel in, real_t delta) const {
        static_cast<Derived *>(this)->integrate(in, delta);
    }
};

class IntegratorRK4: Integrator<IntegratorRK4> {
public:
    template<typename StateModel>
    const StateModel integrate(StateModel in, real_t delta) const {
        StateModel a = in.model();
        StateModel b = static_cast<StateModel>(in + (real_t)0.5 * delta * a).model();
        StateModel c = static_cast<StateModel>(in + (real_t)0.5 * delta * b).model();
        StateModel d = static_cast<StateModel>(in + delta * c).model();
        return in + (delta / (real_t)6.0) * (a + (b * (real_t)2.0) + (c * (real_t)2.0) + d);
    }
};

class IntegratorHeun: Integrator<IntegratorHeun> {
public:
    template<typename StateModel>
    const StateModel integrate(StateModel in, real_t delta) const {
        StateModel initial = in.model();
        StateModel predictor = in + delta * initial;
        return in + (delta * (real_t)0.5) * (initial + predictor.model());
    }
};

class IntegratorEuler: Integrator<IntegratorEuler> {
public:
    template<typename StateModel>
    const StateModel integrate(StateModel in, real_t delta) const {
        return in + delta * in.model();
    }
};

#endif
