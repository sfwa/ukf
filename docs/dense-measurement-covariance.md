# Dense measurement covariance

Diagonal measurement noise is the default: variances for `Core`, standard
deviations for the square-root filters. Select the measurement type's
`DenseCovarianceMatrix` as the final core template argument to include
correlations:

```cpp
using DenseCore = UKF::Core<MyStateVector, MyMeasurementVector,
    UKF::IntegratorRK4, MyMeasurementVector::DenseCovarianceMatrix>;
using DenseSquareRootCore = UKF::SquareRootCore<MyStateVector, MyMeasurementVector,
    UKF::IntegratorRK4, MyMeasurementVector::DenseCovarianceMatrix>;
using DenseParameterCore = UKF::SquareRootParameterEstimationCore<MyStateVector,
    MyMeasurementVector, MyMeasurementVector::DenseCovarianceMatrix>;
```

These aliases assume application-defined state, measurement and model
specialisations. To use a direct process model with `DenseCore`, select
`UKF::NotUsedIntegrator` in place of `UKF::IntegratorRK4`. The
[tracking example](../examples/tracking/) is a complete implementation.

## Covariance versus root factor

For a two-component measurement with correlated noise:

```cpp
DenseCore filter;
filter.measurement_covariance << real_t(4), real_t(0.6),
                                 real_t(0.6), real_t(1);

DenseSquareRootCore root_filter;
root_filter.measurement_root_covariance =
    filter.measurement_covariance.llt().matrixL();
```

Initialise state and process noise too, as described in the
[usage guide](usage-guide.md#initialisation-and-covariance). `Core` accepts the
full symmetric positive-semidefinite covariance `R`. The square-root filters
accept a factor `L` satisfying `R = L * L.transpose()`; passing `R` itself would
use the wrong noise. A lower Cholesky factor is conventional for
positive-definite `R`, but a measurement-noise factor need not be triangular.
If `R` is singular, provide a valid factor rather than assuming ordinary LLT
succeeded. Do not take the elementwise square root of a correlated covariance.

The same interface supports a single scalar observation: its dense covariance
and root are `1 × 1` matrices, containing the variance and standard deviation
respectively. The default diagonal vector interface remains available for every
filter variant.

## Dynamic fields and ordering

For both fixed and dynamic measurements, dense matrices contain **all declared
fields in declaration order**, with vector components contiguous. Keep the
matrix size and ordering even when only some dynamic fields are populated or
when they are inserted in another order. For example:

```cpp
enum SensorKey { PositionSensor = 100, VelocitySensor = 101 };
using Sensors = UKF::DynamicMeasurementVector<
    UKF::Field<PositionSensor, real_t>,
    UKF::Field<VelocitySensor, real_t>>;
Sensors::DenseCovarianceMatrix R;
R << real_t(4), real_t(0.6),
     real_t(0.6), real_t(1); // Declared order: position, velocity.

Sensors sample;
sample.set_field<VelocitySensor>(real_t(2)); // Only velocity is present.
// Assign the entire R to the matching dense filter, not R(1,1) or a resized R.
```

For this sample, the filter uses velocity variance `R(1,1)`. If position is
subsequently inserted, the active measurement order is velocity then position;
the filter also reorders covariance and cross-correlation terms. This differs
from the permanent declaration order of noise storage.

A square-root factor requires particular care: selection removes inactive
**rows**, while retaining every factor **column**. Otherwise shared noise from
an omitted field could incorrectly disappear. Pass the full factor to the
filter and let it perform that selection. Storage is fixed or bounded, and the
diagonal default continues to use its original efficient path.

## Directional fields and changing confidence

Noise for a `UKF::FieldVector` is supplied in raw Cartesian sensor coordinates.
The library applies the rotation-vector Jacobian to both sides of `R`, including
correlations with other fields; for a factor it applies that Jacobian on the
left. Supply sensor noise, not an already transformed innovation covariance.
See the usage guide for [directional and quaternion conventions](usage-guide.md#fields-and-measurement-mapping).

Noise is per filter instance and may change for each observation. Estimate the
relationship between confidence and covariance from sensor errors; the
[confidence guide](usage-guide.md#confidence-and-outliers) explains covariance
versus root scaling. After assigning valid noise, an
[innovation gate](innovation-gating.md) can separately decide whether the
observation is consistent with the prediction.
