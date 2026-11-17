# Models, noise and filter updates

The [README quick start](../README.md#quick-start) is a complete small program.
The [tracking](../examples/tracking/) and [attitude](../examples/attitude/)
examples are built and run by CI. Code fragments below use application-defined
`State`, `Measurement` and filter aliases unless they explicitly define them.

## Fields and measurement mapping

A field key identifies a component within its own state or measurement type.
Keys must be unique within that type. State and measurement keys do not need to
have matching names or integer values. The specialisation
`Measurement::expected_measurement<State, SensorKey>(state, inputs...)` computes
what that sensor should report for a state; this function connects the two.
For example, a `PositionSensor` field can read a state's `Position` field, while
a gyroscope field reads its `AngularVelocity` field.

Use `real_t` for a scalar and `UKF::Vector<N>` for an ordinary vector. Scalar
states, scalar observations and multiple observations of a scalar state are
supported by all three filter variants. State vectors also support
`UKF::Quaternion`, and so do measurement vectors, for an orientation reported
directly by a sensor such as a camera, motion-capture system or another
attitude estimator. A quaternion measurement has four stored coefficients but a
three-dimensional uncertainty. Its innovation is the rotation vector of
`z * conjugate(z_predicted)`, the shortest rotation from the predicted to the
measured orientation, in the same normalised MRP convention as quaternion
states. It does not depend on the sign of either quaternion. Its noise is three
rotation variances (standard deviations for square-root filters), or a 3 × 3
block of [dense noise](dense-measurement-covariance.md), in those rotation
coordinates. A measurement type's `covariance_size()` counts three for each
quaternion field, and `Measurement::InnovationVector` (also the cores'
`InnovationVectorType`) holds the innovation in those coordinates: each
quaternion field becomes a `UKF::Vector<3>` there. For measurements without
quaternion fields `InnovationVector` is the measurement type itself, as before.

`UKF::FieldVector` is a directional measurement such as gravity or a magnetic
field. Its innovation uses a rotation vector, and the noise Jacobian maps raw
Cartesian sensor uncertainty into those rotation coordinates. Rotation and
Jacobian calculations use the same normalized MRP convention. The rotation is
taken between unit vectors, so it does not depend on the field's magnitude or
units: a magnetometer in teslas behaves as one in gauss, in single precision
too. The noise Jacobian accounts for the measured magnitude. At exactly
opposite directions there is no unique rotation axis; the implementation chooses
a deterministic perpendicular axis. That convention does not create additional
physical information or make the antiparallel derivative unique.

A fixed measurement requires every declared field on every update. A dynamic
measurement contains only fields explicitly set for that update. Construct a
fresh dynamic vector each iteration so a missing sensor does not retain an old
value. Fields can be inserted in any order; retrieve them by key rather than by
a hard-coded component index. Noise storage, however, always follows the
**declared** field order, including fields absent from a dynamic observation.

## Process models and inputs

`UKF::Core<State, Measurement>` defaults to a direct model. Specialise
`State::process_model<>(dt, inputs...)` to return the next state, as in the
tracking example. This is a discrete transition, not a derivative.

For continuous-time dynamics, specialise `State::derivative<>(inputs...)` and
select `UKF::IntegratorEuler`, `UKF::IntegratorHeun` or `UKF::IntegratorRK4`:

```cpp
using IntegratedFilter = UKF::Core<State, Measurement, UKF::IntegratorRK4>;
using RootFilter = UKF::SquareRootCore<State, Measurement, UKF::IntegratorRK4>;
```

If the integrator is omitted (or is `UKF::NotUsedIntegrator`), both `Core`
and `SquareRootCore` call the direct `process_model<>(dt, inputs...)`
specialisation; with an integrator they integrate `derivative<>(inputs...)`.
`SquareRootParameterEstimationCore<State, Measurement>` has no
state-transition model or time-step argument: prediction adds the configured
parameter process noise.

Pass optional process/measurement inputs after the measurement:
`filter.step(dt, measurement, input1, input2)`. Both models receive the same
inputs as const references. Provide matching specialisations for the exact
input types, including measurement fields that ignore an input. The parameter
estimator instead uses `filter.step(measurement, inputs...)`.

For attitude, angular velocity and angular acceleration are physical
three-vectors, not quaternion coefficient derivatives or Euler-angle rates.
The application supplies quaternion kinematics consistent with its frame and
multiplication convention. The attitude example declares its convention and
implements that derivative. The filter does not infer a dynamics model by
differentiating a sequence of quaternion observations.

Sigma-point scaling is configured through specialisations in `UKF::Parameters`,
for example `AlphaSquared<State>`, as described in the
[state-vector header](../include/UKF/StateVector.h) and the
[algorithm references](../README.md#references). Set these before using the
corresponding filter type. The weights and spread derive from `AlphaSquared`,
`Kappa` and `Beta` (through `ScaledDimension`, which is n + lambda), so adjust
those rather than `Lambda`. A small `AlphaSquared` makes the central weight
negative, which a strongly nonlinear model can turn into an indefinite
predicted covariance (the step then fails, as described below). In single
precision keep `AlphaSquared` well above 1e-6: the outer weights,
1 / (2 alpha² (n + kappa)), multiply the rounding of every sigma point, so with
1e-6 a mean is good to only a few tenths of a percent of its value. Parameters such as `MRP_A<Measurement>` and
`MRP_F<Measurement>` determine the directional-measurement rotation coordinates;
the same convention is used for their innovations and noise Jacobians.

## Initialisation and covariance

Initialise every state field, the state uncertainty, process noise and
measurement noise before the first update. Normalise quaternion state values.
The uncertainty of a quaternion state is three-dimensional, so use
`State::CovarianceMatrix` and `State::covariance_size()` rather than the number
of stored state coefficients when allocating uncertainty matrices.

| Quantity | `Core` | Square-root filters |
| --- | --- | --- |
| State uncertainty | `covariance = P` | `root_covariance = S`, with `P = S * S.transpose()` |
| Process noise | `process_noise_covariance = Q` | `process_noise_root_covariance = LQ`, with `Q = LQ * LQ.transpose()` |
| Default measurement noise | `measurement_covariance`: variances | `measurement_root_covariance`: standard deviations |
| Dense measurement noise | Full covariance `R` | Full factor `LR`, with `R = LR * LR.transpose()` |

All noise settings belong to the **filter instance**, so filters with the same
template types can have different noise and can change it for each observation.
The default measurement storage contains only diagonal values. See
[dense measurement covariance](dense-measurement-covariance.md) for correlations
within or between sensor fields. Square-root factors are not elementwise square
roots of a covariance matrix; a Cholesky factor is a suitable choice for a
positive-definite covariance.

Both state-estimation filters use discrete additive process noise:

```text
x_next = process_model(x, dt, inputs) + noise,  covariance(noise) = Q(dt)
```

`Q` is added after propagating the prior, and the measurement update includes
that added uncertainty. The parameter estimator similarly combines full state
and process covariances, preserving correlations. Set the noise for the actual
time interval; the library does not convert continuous-time noise densities or
rescale a previously assigned `Q` when `dt` changes.

For a one-dimensional position/velocity state driven by continuous white
acceleration with spectral density `q`, one possible model is:

```text
Q(dt) = q * [ dt³/3   dt²/2 ]
            [ dt²/2   dt    ]
```

This differs from a model with one constant random acceleration per step.
Choose the discretisation that matches the physical assumptions of the model.

Covariances must be finite, symmetric and positive semidefinite. The standard
filter can generate sigma points from singular or zero state covariance; it
uses a checked factorisation, tolerates roundoff-sized negative eigenvalues and
stores the repaired covariance. Square-root correction supports semidefinite
posteriors, including exact observations. Supply valid measurement noise as
well; the filters do not validate every user input.

## Error handling

The library does not use C++ exceptions and can be built with them disabled.
Each of `step()`, `a_priori_step()`, `innovation_step()` and
`a_posteriori_step()` returns `true` on success. On failure it returns `false`
and leaves the public state, (root) covariance, innovation and innovation
(root) covariance exactly as they were. A step fails when a covariance (or a
square-root downdate) is indefinite beyond roundoff, or when a result is not
finite. A step that depends on a failed (or missing) earlier step also returns
`false` until the next successful `a_priori_step()`.

`step()` is all-or-nothing: if its correction fails after a successful
prediction, the prediction is rolled back too. To keep the prediction when a
measurement cannot be applied, call the split stages described below and check
each result. Always check the result before using an update:

```cpp
if (!filter.step(dt, measurement)) {
    // The filter is unchanged: skip this measurement or reinitialise.
}
```

## Confidence and outliers

A detector confidence score is not automatically an inverse variance. Estimate
`R(score)` using representative sensor errors, including units and correlations.
The tracking example uses a declared synthetic noise model and assigns exactly
that covariance, so the relationship is known for the demonstration.

A heuristic such as `R = baseline_R / (score * score)` needs calibration. For
the same heuristic, a square-root factor scales as `LR = baseline_LR / score`,
not by `1 / score²`. Reject invalid, nonfinite or zero scores before computing
noise or entering innovation calculation; guard application-specific upper and
lower bounds too. Start each calculation from the baseline rather than
repeatedly scaling the last observation's covariance.

Confidence weighting and [innovation gating](innovation-gating.md) serve
different purposes. Covariance determines how much an accepted observation
influences the estimate. A gate checks its residual against the predicted
uncertainty and can reject the whole update. A gate runs after innovation
calculation, so it cannot repair invalid noise values supplied beforehand.

## Prediction, correction and missing observations

A full state-filter update is `filter.step(dt, measurement, inputs...)`. The
three stages are also public:

```cpp
if (filter.a_priori_step(dt)) {           // state and uncertainty are now predictions.
    if (filter.innovation_step(measurement) &&
            filter.a_posteriori_step()) {     // state and uncertainty are now corrected.
        // Corrected.
    } // Otherwise the prediction is kept.
}
```

Forward the same model inputs to the first two calls when using inputs.
After a successful `a_priori_step`, `filter.get_sigma_points()` returns a const
reference to its sigma points, one per column in the state vector's layout:
the points drawn from the predicted state and uncertainty, process noise
included, which the innovation step passes through the measurement model. The
a posteriori step leaves them unchanged, and they are unspecified after a step
that fails. To apply a gate, use `step_with_gate(...)` or the gated posterior overload described
in the gating guide. Do not apply both the gated and ungated posterior to the
same innovation.

When no observation is available, call only `a_priori_step`. An empty dynamic
measurement also yields prediction-only behaviour with `step` on any of the
three filters (it returns `true`). With `step_with_gate`, an empty measurement
returns `false` without invoking the callback, retains the prediction and
clears its innovation diagnostics.

Prediction mutates the filter. To preview a future state without advancing the
live estimator, copy it and configure process noise for the forecast interval:

```cpp
auto forecast = filter;
forecast.process_noise_covariance = forecast_Q; // Q for horizon, supplied by the application.
if (forecast.a_priori_step(horizon)) {
    // Read forecast.state; filter is unchanged.
}
```

For a square-root filter assign `process_noise_root_covariance` instead. Long
forecasts can be broken into steps, using the correct process noise and inputs
for each step. The parameter estimator's `a_priori_step()` has no `dt`; it adds
its configured noise once per call.

## Read the estimate, not the observation

After a successful update, read `filter.state.get_field<StateKey>()`.
`measurement.get_field<SensorKey>()` still contains the raw observation. The
tracking example prints raw and estimated position errors separately and uses
consistent ground-truth motion, process dynamics and measurement units.

Tune process and measurement noise against measured error and model mismatch.
A large process covariance or a small measurement covariance makes an estimate
follow measurements more closely. A filter need not visibly smooth every
trajectory, and comparing raw output alone cannot show whether filtering worked.
The library performs causal filtering; it does not currently provide an RTS or
fixed-lag smoother that revises past estimates using future measurements.
