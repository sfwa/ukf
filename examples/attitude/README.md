# Attitude estimation

This small square-root filter estimates orientation and body angular velocity
from gyroscope readings and two nonparallel reference vectors. Its input data
are generated deterministically, so it runs without sensor hardware or data
files:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target ukf_attitude_example --parallel 2
./build/examples/attitude/ukf_attitude_example
```

The executable prints the estimated quaternion, angular velocity, and rotation
error. It exits unsuccessfully if estimates become nonfinite, the quaternion
loses unit norm, or tracking accuracy exceeds conservative bounds. Debug builds
enable `EIGEN_NO_MALLOC`.

## Field mapping and dynamics

`Orientation` and `BodyRate` are state keys. `Gyroscope`, `GravityDirection`, and
`MagneticDirection` are different observation keys, with different numeric IDs.
The `expected_measurement<AttitudeState, Key>` specializations explicitly define
each relationship; matching key names or numbers are not required.

The quaternion rotates navigation-frame vectors into body coordinates. Angular
velocity is a physical three-vector expressed in the body frame, in radians per
second. With that convention, the user-supplied process model is
`q_dot = -0.5 * [0, omega_body] * q`, with constant angular velocity. RK4
integrates this derivative. Angular acceleration, biases, and other dynamics
would also be supplied by the application; the library does not infer them by
differentiating observations.

The two reference vectors are normalized, synthetic gravity and magnetic
directions. Their Cartesian components use ordinary `Vector<3>` measurements,
as in the larger AHRS example. The quaternion is part of the **state**;
native quaternion-valued **measurements** are not currently supported.

## Covariance conventions

The state has seven stored coefficients but six covariance dimensions: three
orientation-error coordinates and three angular-velocity coordinates. A
square-root filter stores factors `S` with `P = S * S.transpose()`. Its diagonal
measurement-noise entries therefore contain standard deviations, not variances.

This example explicitly assumes diagonal random-walk process increments with
`Q(dt) = diag(q_attitude, q_rate) * dt` and supplies its square root. General
models, especially coupled position/velocity or attitude/rate uncertainty, need
their own discrete-time covariance. The synthetic disturbances and tuning here
are illustrative, not a sensor calibration or a production AHRS model.
