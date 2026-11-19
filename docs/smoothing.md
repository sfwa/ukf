# Fixed-lag and offline smoothing

The ordinary filters estimate the current state given the samples received so
far. An unscented Rauch–Tung–Striebel (RTS) smoother uses later observations to
revise an earlier estimate. A lag of two updates means that the estimate at
time k incorporates measurements through k+2. The live filter operates at its
original rate. Smoothing is experimental in 1.0 (see
[RELEASING.md](../RELEASING.md#experimental-interfaces)).

The [random-walk example](../examples/smoothing/smoothing.cpp) starts with
mean 0 and variance 1. With process variance 1, sensor variance 1 and successive
readings of 1 and 2, the two-step smoother revises the initial mean to 0.5
and its variance to 0.625. Build and run it from the repository root.

```sh
cmake --build build --target ukf_smoothing_example
./build/examples/smoothing/ukf_smoothing_example
```

Include `UKF/Smoother.h`. Call `reset({state, covariance})` with the initial
filtered estimate. Repeat the following steps for each subsequent update.

1. Set process noise for the interval. Call
   `a_priori_step_with_cross_covariance(cross, dt, inputs...)` on `Core` or
   `SquareRootCore` and save its predicted state and full covariance.
2. Run the usual innovation and correction stages. When there is no accepted
   observation, use the prediction as the filtered estimate.
3. Call `smoother.push(prediction, cross, filtered_estimate)` after deciding
   which estimate to keep. For a square-root filter, convert its lower factor
   to full covariance with `L * L.transpose()` at both snapshots.
4. Call `smoother.estimate(steps_back, output)`. Zero returns the latest
   filtered estimate. Valid ages run from zero to `smoother.size()-1`.
   A full window contains `Lag+1` estimates. Query its oldest retained state
   before pushing a new update if every delayed output is required.

The cross covariance pairs each sigma point at the previous filtered state
with the same point after propagation. Independent additive process noise
enters the predicted covariance only. The exposed measurement sigma points
have already been redrawn, so they cannot recover those pairings.

For offline use, retain the filtered estimates, predictions and crosses.
Initialise the last smoothed estimate from the last filtered one. Working
backwards, call `unscented_rts_step(filtered_k, predicted_next, cross,
smoothed_next, smoothed_k)` for every transition. Both smoothing APIs use the
state's left-perturbation MRP convention for quaternion differences and
updates. They transport covariance between tangent coordinates with local
Jacobians. This approximation assumes uncertainty is local about the supplied
unit quaternions.

The ring holds at most `Lag` records with two states and three covariance-sized
matrices each, plus the latest estimate. Prediction tracing, the ring and
backward steps use fixed-size storage. A backward query runs one RTS step per
requested lag. Singular positive-semidefinite covariances can represent exact
observations. The smoothing gain is solved in unit-variance coordinates,
so changing a state's physical units preserves its numerical rank. Directions
of the correlation matrix indistinguishable from zero within roundoff receive
zero gain. Check each Boolean result. A failed reset or push preserves the
existing window and a failed query or offline step preserves its output.
