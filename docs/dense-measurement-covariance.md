# Dense measurement covariance

Diagonal measurement noise remains the default, storing only variances (or
standard deviations for a square-root filter). To include correlations, select
`MyMeasurementVector::DenseCovarianceMatrix` as the final core template argument:

```C++
using DenseCore = UKF::Core<MyStateVector, MyMeasurementVector,
    UKF::IntegratorRK4, MyMeasurementVector::DenseCovarianceMatrix>;
using DenseSquareRootCore = UKF::SquareRootCore<MyStateVector, MyMeasurementVector,
    UKF::IntegratorRK4, MyMeasurementVector::DenseCovarianceMatrix>;
using DenseParameterCore = UKF::SquareRootParameterEstimationCore<MyStateVector,
    MyMeasurementVector, MyMeasurementVector::DenseCovarianceMatrix>;

DenseCore dense_filter;
dense_filter.measurement_covariance.setIdentity();
dense_filter.measurement_covariance(0, 1) = real_t(0.25);
dense_filter.measurement_covariance(1, 0) = real_t(0.25);

DenseSquareRootCore square_root_filter;
// For positive-definite noise covariance R, obtain a factor L with R = L * L.transpose().
square_root_filter.measurement_root_covariance =
    dense_filter.measurement_covariance.llt().matrixL();
```

Initialise state and process noise as in the [README examples](../README.md#initialisation). Standard filters
accept the full symmetric positive-semidefinite noise covariance `R`;
square-root filters accept a full factor `L`, **not** `R`. A lower Cholesky
factor is conventional, but the factor need not be triangular.

For both fixed and dynamic measurements, dense matrices contain every declared
field in declaration order, with vector components contiguous. Keep this matrix
size and ordering when only some dynamic fields are populated, even if they
are inserted in another order. The filter selects and reorders the covariance
for the fields present in that measurement. For a dense root it selects rows
and retains all factor columns, so absent sensors do not remove shared noise
contributions. Storage is fixed or bounded; the default diagonal interface
continues to use its original efficient path.

Noise for a `FieldVector` is specified in raw Cartesian sensor coordinates.
The library applies its rotation-vector Jacobian to both sides of the
covariance, including correlations with other fields (or to the left of a
square-root factor), before using it as innovation noise.

