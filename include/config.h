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

#ifndef CONFIG_H_
#define CONFIG_H_

/* Whether to use Eigen3 or architecture-specific implementations */
#define UKF_USE_EIGEN

/* Uncomment one of the following to select the float-point precision */
/* #define UKF_SINGLE_PRECISION */
#define UKF_DOUBLE_PRECISION

/*
Choose the integration method here:
    - UKF_INTEGRATOR_RK4: 4th-order integration method. Requires most CPU time.
    - UKF_INTEGRATOR_HEUN: 2nd-order Heun's method. Requires middle CPU time.
    - UKF_INTEGRATOR_EULER: 1st-orer Euler's method. Requires least CPU time.
*/

#define UKF_INTEGRATOR_RK4
/* #define UKF_INTEGRATOR_HEUN */
/* #define UKF_INTEGRATOR_EULER */

/*
UKF vector dimensioning: note that the UKF state vector is one element larger
than UKF_STATE_DIM, as the attitude is stored as a quaternion.
*/

#define UKF_CONTROL_DIM 4
#define UKF_STATE_DIM 24
#define UKF_MEASUREMENT_MAX_DIM 32

#endif
