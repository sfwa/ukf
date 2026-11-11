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
`UKF::Quaternion`. Native quaternion measurement fields are not currently
supported: an orientation has four stored coefficients but needs a
three-dimensional rotation residual. Do not treat four quaternion coefficients
as independent linear measurements.

`UKF::FieldVector` is a directional measurement such as gravity or a magnetic
field. Its innovation uses a rotation vector, and the noise Jacobian maps raw
Cartesian sensor uncertainty into those rotation coordinates. Rotation and
Jacobian calculations use the same normalized MRP convention. At exactly
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

The square-root state filter requires an explicit integrator. Its ordinary
integration path calls `process_model<IntegratorType>`; do not assume that
omitting the integrator or selecting `NotUsedIntegrator` enables the direct
model dispatch used by `Core`. `SquareRootParameterEstimationCore<State,
Measurement>` has no state-transition model or time-step argument: prediction
adds the configured parameter process noise.

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
corresponding filter type. Parameters such as `MRP_A<Measurement>` and
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
uses a checked factorisation and tolerates roundoff-sized negative eigenvalues.
Genuinely indefinite or nonfinite state covariance raises `std::domain_error`
when factorisation is required. Supply valid measurement noise as well; these
checks do not validate every user input.

Square-root correction supports semidefinite posteriors, including exact
observations. If a downdate would produce an indefinite posterior beyond
roundoff, it throws `std::domain_error` before applying that correction to the
predicted state or root covariance. An entire `step()` is not a transaction:
prediction may already have completed when an exception occurs.

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
filter.a_priori_step(dt);       // state and uncertainty are now predictions.
filter.innovation_step(measurement);
filter.a_posteriori_step();     // state and uncertainty are now corrected.
```

Forward the same model inputs to the first two calls when using inputs. To
apply a gate, use `step_with_gate(...)` or the gated posterior overload described
in the gating guide. Do not apply both the gated and ungated posterior to the
same innovation.

When no observation is available, call only `a_priori_step`. An empty dynamic
measurement also yields prediction-only behaviour with `Core::step`, or with
`step_with_gate` on any of the three filters. The gated call returns `false`
without invoking the callback and clears its innovation diagnostics. For
ordinary square-root `step` calls, use prediction-only explicitly when the
measurement is empty.

Prediction mutates the filter. To preview a future state without advancing the
live estimator, copy it and configure process noise for the forecast interval:

```cpp
auto forecast = filter;
forecast.process_noise_covariance = forecast_Q; // Q for horizon, supplied by the application.
forecast.a_priori_step(horizon);
// Read forecast.state; filter is unchanged.
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
