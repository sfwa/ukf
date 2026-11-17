# Innovation gating

All three filters accept a templated callable that decides whether to apply the
measurement correction. Existing `step(...)` and `a_posteriori_step()` calls keep
their behavior; gating is opt-in.

For state estimation, use:

```cpp
bool applied = filter.step_with_gate(dt, measurement, gate, model_inputs...);
```

For `SquareRootParameterEstimationCore`, omit `dt`:

```cpp
bool applied = filter.step_with_gate(measurement, gate, model_inputs...);
```

After successful prediction and innovation calculation, the filter invokes
`gate(const Filter&, const Measurement&, const Inputs&...)` once. The filter
reference exposes the predicted state, innovation and `innovation_covariance`
(`Core`) or `innovation_root_covariance` (both square-root filters). The
measurement and model inputs are the same objects supplied to the step.

Returning `true` applies the usual correction. Returning `false` rejects the
whole measurement vector: the predicted state and covariance remain intact,
while innovation diagnostics remain available for logging. The step returns
whether a correction was applied. An empty dynamic measurement returns `false`
without invoking the gate, retaining the prediction and clearing the innovation
diagnostics.

Like ordinary `step`, a gated full step restores the pre-step state, covariance,
innovation and innovation covariance on numerical failure. A failed prediction
or innovation never invokes the gate, including when diagnostics remain from an
earlier successful update. If an accepting gate is followed by a failed
correction, the step returns `false` and rolls back the prediction as well.
Thus `false` can mean rejection, an empty observation, or numerical failure;
use the checked split stages below and record the gate's decision when the
application needs to distinguish those outcomes.

The library reports numerical failure without exceptions and works with
exceptions disabled. If an application enables exceptions and its gate throws,
the exception propagates; no correction has been applied, and the prediction
and innovation diagnostics remain available.

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

Include `<cmath>` for `std::isfinite`. Choose the threshold for the number of
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
if (filter.a_priori_step(dt, model_inputs...) &&
        filter.innovation_step(measurement, model_inputs...)) {
    bool applied = filter.a_posteriori_step(gate, measurement, model_inputs...);
    // Handle acceptance, rejection, or a failed correction here.
}
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
method. It returns `false` for a missing, failed, or empty innovation without
calling the gate. If the gate accepts, the return value reports the actual
correction result. A rejected or failed split posterior leaves the predicted
state, covariance and innovation diagnostics unchanged; unlike a failed full
step, it does not roll back preceding successful split stages.
Gating selects whether to correct with the entire measurement vector; an
application needing separate decisions per sensor can perform separate
measurement updates with its dynamic measurement vector.
