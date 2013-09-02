# ukf

Unscented Kalman filter library.


## Structure

`ukf/include` contains the header files, divided into the following:
* `dynamics.h`: Dynamics models. A simple centripetal model and a fixed-wing
aerodynamic model are provided.
* `integrator.h`: Integration routines. A 4th-order Runge-Kutta,
2nd-order Heun method, or 1st-order Euler method can be selected.
* `sensors.h`: Defines the `SensorModel` interface and a model of the I/O
board's sensors.
* `state.h`: Declares the `State` type, a 22-dimensional column vector
used to represent the state of the Kalman filter. Individual components of
the `State` type can be accessed using the accessor methods provided.
* `types.h`: Defines the `real_t` type (either a single- or double-precision
floating-point, depending on the precision selected in `config.h`), and a few
required constants.
* `ukf.h`: Defines the interface to the Unscented Kalman Filter, along with
sigma point scaling constants, MRP parameters and sigma point weights.

`ukf/src` contains source files, hopefully divided into logical components:
* `dynamics.cpp`: Functions relating a state vector to the linear and angular
acceleration.
* `sensors.cpp`: The I/O board sensor model, including measurement prediction
and mean finding functions.
* `state.cpp`: Kinematic state transition function.
* `ukf.cpp`: The Kalman filter proper.
* `ukf-estimates.cpp`: Kalman filter estimation functions, broken out into a
separate file to work around issues in the Texas Instruments CCS compiler.

`ukf/test` contains unit tests, built using the googletest framework.

`ukf/c` is the source to `libcukf`, a C interface to the C++ UKF static library.

`ukf/python` contains a ctypes-based Python wrapper for `libcukf`.

`ukf/ccs-c66x` contains a project for Texas Instruments Code Composer Studio 5,
targeting the Keystone DSP platform (C66x cores). Includes a modified version
of Eigen.


## Configuration

`src/config.h` is currently used for configuration. The parameters which can
be configured here are as follows:
* The precision (single or double) of the floating-point values used by the
library;
* The integration method used (RK4, Heun or Euler).


TODO: Use cmake to automatically generate `config.h`


## Building

Requires `cmake` version 2.8.7 or higher.

Create a build directory outside the source tree, then use cmake to generate
the makefile.

`mkdir ukf_build`

`cd ukf_build`

`cmake /path/to/ukf`

Now, build the library using the `make` command. An appropriate version of
Eigen will be downloaded automatically.

To build the dynamic library, run `make cukf`. A dynamic library appropriate
for the host platform should be built.


## Testing

The `googletest` library is used for unit testing.

After creating the build directory, `make check` will automatically download
and build googletest, build the unit tests and then run them.

To build the unit tests without running them, use `make unittest`. The unit
tests can then manually be run (with more detailed reporting) by running
`test/unittest` in the build directory.


## Python module installation

Requires `cmake` version 2.8.7 or higher.

Run `python setup.py install` to build the C shared library and install the
Python interface (the `ukf` module) in your `site-packages` directory.

Alternatively, just run `pip install https://github.com/sfwa/ukf/archive/master.zip#egg=ukf-1.0.0`
to download and install.


## Compiling with CCS5

