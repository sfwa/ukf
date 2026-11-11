# Tracking, confidence, gating, and forecasting

Build from the repository root and run:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target ukf_tracking_example --parallel 2
./build/examples/tracking/ukf_tracking_example
```

The executable tracks one-dimensional position and velocity. Its state fields
and measurement fields use different keys: `expected_measurement` explicitly
maps each sensor to the state. The process model is a direct, discrete-time
constant-velocity function, so it uses `NotUsedIntegrator` and needs no
`derivative` specialization.

Synthetic measurements have a known confidence-dependent noise scale and
correlated position/velocity errors. The example updates a dense covariance
matrix before each observation. It includes scalar position-only observations,
reordered position/velocity observations, missing samples, invalid confidence,
and injected outliers. Dense covariance entries retain declaration order even
when a dynamic measurement is populated in another order.

A normalized-innovation-squared gate rejects whole measurement vectors that
are inconsistent with the prediction. Its thresholds depend on the number of
active observations and illustrate a Gaussian-noise policy; the deterministic
sinusoids do not establish a statistical false-rejection rate. Missing or
invalid observations leave the measurement
empty; the filter predicts without calling the gate. Rejection retains that
prediction, while the innovation remains available for diagnosis.

Confidence calibration is specific to the sensor. Here, the sensor generator
uses the same standard deviations as the covariance model. In a real system,
estimate `R(score)` from residuals rather than treating a detector probability
as inverse variance. If a chosen calibration uses `R / score²`, a square-root
filter needs its factor scaled by `1 / score`, not by `1 / score²`. Check
invalid or out-of-range scores before innovation calculation.

The program prints raw and estimated position RMSE after initialization,
including injected outliers, followed by the actual `filter.state` and a
two-second forecast. Forecasting uses `a_priori_step` on a copy so the live
filter remains unchanged. Process noise is discretized for the forecast
interval as well as for each regular update. This is forward filtering and
prediction; it does not perform retrospective smoothing using future samples.

The deterministic signals make the example repeatable. A nonzero exit status
reports failed consistency checks, including ineffective filtering, missed
outliers, changes to rejected predictions, or an incorrect forecast. This is a
demonstration of a calibrated toy model, not a universal performance guarantee.
