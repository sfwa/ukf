# ukf

Unscented Kalman filter library.


## Structure

`ukf/include` contains the header files, divided into the following:
* `state.h` declares the `State` type, a 22-dimensional column vector
used to represent the state of the Kalman filter. Individual components of
the `State` type can be accessed using the accessor methods provided;
* `types.h` defines the `real_t` type (either a single- or double-precision
floating-point, depending on the precision selected in `config.h`).

`ukf/src` contains source files, hopefully divided into logical components:
* `integrator.h` contains the integration routines. A 4th-order Runge-Kutta,
2nd-order Heun method, or 1st-order Euler method can be selected.
* `state.cpp` contains the kinematic state transition function, expressed as
a series of ordinary differential equations.
* `dynamics.cpp` defines a number of dynamics model. Currently, the only
model implemented simply estimates centripetal acceleration.
* `ukf.cpp` contains the various filter steps.

`ukf/test` contains unit tests, built using the googletest framework.

`ukf/c` provides a dynamic library with a C interface to the UKF library.

`ukf/python` contains ctypes-based Python code for using the C interface.


## Configuration

`src/config.h` is currently used for configuration. The parameters which can
be configured here are as follows:
* The precision (single or double) of the floating-point values used by the
library;
* The integration method used (RK4, Heun or Euler).


TODO: Use cmake to automatically generate `config.h`


## Building

Required packages:
* `cmake` version 2.8.7 or higher;
* `eigen` version 3.0 or higher, installed to `/usr/local/include`.

Create a build directory outside the source tree, then use cmake to generate
the makefile.

`mkdir ukf_build`

`cd ukf_build`

`cmake /path/to/ukf`

Now, build the library using the `make` command.

To build the dynamic library, run `make cukf`. A dynamic library appropriate
for the host platform should be built.


## Testing

The `googletest` library is used for unit testing.

After creating the build directory, `make check` will automatically download
and build googletest, build the unit tests and then run them.

To build the unit tests without running them, use `make unittest`. The unit
tests can then manually be run (with more detailed reporting) by running
`test/unittest` in the build directory.
