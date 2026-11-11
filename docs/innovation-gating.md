# Innovation gating

All three filters accept a templated callable that decides whether to apply the
measurement correction. Gating is opt-in through `step_with_gate` or a posterior
overload. The [tracking example](../examples/tracking/) runs a complete gate
with dynamic, correlated measurements and calibrated synthetic noise.

The call forms below use an already initialised filter, application-defined
measurement and model inputs. For state estimation, use:

```cpp
bool applied = filter.step_with_gate(dt, measurement, gate, model_inputs...);
```

For `SquareRootParameterEstimationCore`, omit `dt`:

```cpp
bool applied = filter.step_with_gate(measurement, gate, model_inputs...);
```

After prediction and innovation calculation, the filter invokes
`gate(const Filter&, const Measurement&, const Inputs&...)` once. The filter
reference exposes the predicted state, innovation and `innovation_covariance`
(`Core`) or `innovation_root_covariance` (both square-root filters). The
measurement and model inputs are the same objects supplied to the step.

Returning `true` applies the usual correction. Returning `false` rejects the
whole measurement vector: the predicted state and covariance remain intact,
while innovation diagnostics remain available for logging. The step returns
whether a correction was applied. An empty dynamic measurement returns `false`
without invoking the gate. If the gate throws, no correction has been applied;
the predicted state and innovation diagnostics remain available. A later
correction can also throw if its covariance is invalid: `true` is returned only
after that correction succeeds.

The callable is invoked synchronously and is neither copied nor stored. A lambda
can capture a threshold, counters, or other application state. No type erasure or
heap allocation is added by the gating interface. The const filter argument
prevents accidental updates through that argument; avoid mutating the filter
through a separately captured reference inside the gate.

## Normalized innovation squared

For a positive-definite innovation covariance, a Mahalanobis-distance gate for
`Core` can be written as follows:

```cpp
real_t threshold = 9;
auto gate = [threshold](const auto& current, const auto&, const auto&...) {
    // Materialize the solve before the temporary decomposition is destroyed.
    const auto weighted = current.innovation_covariance.llt()
        .solve(current.innovation).eval();
    const real_t nis = current.innovation.dot(weighted);
    return std::isfinite(nis) && nis <= threshold;
};
bool applied = filter.step_with_gate(dt, measurement, gate, model_inputs...);
```

Include `<cmath>` for `std::isfinite`. This example assumes the covariance is
positive definite; applications accepting wider covariance inputs should check
factorisation success and define their singular-covariance policy. Choose the threshold for the number of
active innovation components and the application's desired rejection rate;
`9` is only an example, not a dimension-independent default. For singular
innovation covariances, define a deliberate policy for residuals in the null
space instead of using this positive-definite example unchanged.

Square-root filters store a lower factor `L` with innovation covariance
`L * L.transpose()`. Whiten the residual by solving `L * w = innovation`:

```cpp
auto root_gate = [threshold](const auto& current, const auto&, const auto&...) {
    const auto whitened = current.innovation_root_covariance
        .template triangularView<Eigen::Lower>().solve(current.innovation).eval();
    const real_t nis = whitened.squaredNorm();
    return std::isfinite(nis) && nis <= threshold;
};
bool applied = root_filter.step_with_gate(dt, measurement, root_gate, model_inputs...);
// For a parameter estimator: parameter_filter.step_with_gate(measurement, root_gate, model_inputs...);
```

This computes the same statistic without forming a covariance matrix or an
inverse. It assumes a nonsingular innovation root.

## Split-step use

To inspect prediction or innovation independently, use the posterior overload:

```cpp
filter.a_priori_step(dt, model_inputs...);
filter.innovation_step(measurement, model_inputs...);
bool applied = filter.a_posteriori_step(gate, measurement, model_inputs...);
```

Only the const filter argument is implicit. Any other context is passed
explicitly, allowing a separate gating policy from the process/measurement
model inputs:

```cpp
bool applied = filter.a_posteriori_step(
    [](const auto& current, const real_t& limit) {
        return current.innovation.squaredNorm() <= limit;
    }, residual_squared_limit);
```

This overload must follow `innovation_step`, just like the original posterior
method. It returns `false` for an empty innovation without calling the gate.
Gating selects whether to correct with the entire measurement vector. An
application needing per-sensor rejection must select its observation set or
implement a separate update policy; returning a mask is not part of this API.
Do not call an ungated posterior after a rejected gated posterior unless you
intend to override that rejection.

## Confidence is a separate input

Set valid per-observation measurement noise before calling the filter. A gate
runs after innovation calculation, so invalid confidence values or division by
zero must be handled earlier. Weighting an observation by its calibrated
covariance and rejecting an inconsistent residual are complementary choices.
See [confidence and outliers](usage-guide.md#confidence-and-outliers).
