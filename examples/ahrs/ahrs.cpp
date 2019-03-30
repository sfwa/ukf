#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
#include "ahrs.h"

/*
This is an implementation of an Unscented Kalman filter for a 9-axis AHRS,
using accelerometer, gyroscope and magnetometer to estimate attitude and
angular velocity.
*/

/* Value of g in m/s^2. */
#define G_ACCEL (9.80665)

/* Default magnetic field norm in Gauss. */
#define MAG_NORM (0.45)

enum AHRS_Keys {
    /* AHRS filter fields. */
    Attitude,
    AngularVelocity,

    /* Parameter estimation filter fields. */
    AccelerometerBias,
    GyroscopeBias,
    MagnetometerBias,
    MagnetometerScaleFactor,
    MagneticFieldNorm,
    MagneticFieldInclination,

    /* AHRS measurement vector fields. */
    Accelerometer,
    Gyroscope,
    Magnetometer
};

/*
The AHRS state vector contains the following:
- Attitude as a quaternion (NED frame to body frame)
- Angular velocity (body frame, rad/s)
*/
using AHRS_StateVector = UKF::StateVector<
    UKF::Field<Attitude, UKF::Quaternion>,
    UKF::Field<AngularVelocity, UKF::Vector<3>>
>;

namespace UKF {
namespace Parameters {
template <> constexpr real_t AlphaSquared<AHRS_StateVector> = 1e-2;
template <> constexpr real_t Kappa<AHRS_StateVector> = 3.0;
}
}

/*
In addition to the AHRS filter, an online parameter estimation filter is also
implemented in order to calculate biases in each of the sensors.
The magnetometer scale factor is represented as a direction cosine matrix
with no normalisation constraint.
*/

using AHRS_SensorErrorVector = UKF::StateVector<
    UKF::Field<AccelerometerBias, UKF::Vector<3>>,
    UKF::Field<GyroscopeBias, UKF::Vector<3>>,
    UKF::Field<MagnetometerBias, UKF::Vector<3>>,
    UKF::Field<MagnetometerScaleFactor, UKF::Vector<3>>,
    UKF::Field<MagneticFieldNorm, real_t>,
    UKF::Field<MagneticFieldInclination, real_t>
>;

namespace UKF {
namespace Parameters {
template <> constexpr real_t AlphaSquared<AHRS_SensorErrorVector> = 1.0;
template <> constexpr real_t Kappa<AHRS_SensorErrorVector> = 3.0;
}

/* AHRS process model. */
template <> template <>
AHRS_StateVector AHRS_StateVector::derivative<>() const {
    AHRS_StateVector output;

    /* Calculate change in attitude. */
    UKF::Quaternion omega_q;
    omega_q.vec() = get_field<AngularVelocity>() * 0.5;
    omega_q.w() = 0;
    output.set_field<Attitude>(omega_q.conjugate() * get_field<Attitude>());

    /* Assume constant angular velocity. */
    output.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));

    return output;
}
}

using AHRS_MeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<Accelerometer, UKF::Vector<3>>,
    UKF::Field<Gyroscope, UKF::Vector<3>>,
    UKF::Field<Magnetometer, UKF::Vector<3>>
>;

namespace UKF {
/*
This is the measurement model that's actually used in the filter, because
it's the one which takes the parameter estimation filter state as an input.
*/
template <> template <>
UKF::Vector<3> AHRS_MeasurementVector::expected_measurement
<AHRS_StateVector, Accelerometer, AHRS_SensorErrorVector>(
        const AHRS_StateVector& state, const AHRS_SensorErrorVector& input) {
    return input.get_field<AccelerometerBias>() + state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -G_ACCEL);
}

template <> template <>
UKF::Vector<3> AHRS_MeasurementVector::expected_measurement
<AHRS_StateVector, Gyroscope, AHRS_SensorErrorVector>(
        const AHRS_StateVector& state, const AHRS_SensorErrorVector& input) {
    return input.get_field<GyroscopeBias>() + state.get_field<AngularVelocity>();
}

template <> template <>
UKF::Vector<3> AHRS_MeasurementVector::expected_measurement
<AHRS_StateVector, Magnetometer, AHRS_SensorErrorVector>(
        const AHRS_StateVector& state, const AHRS_SensorErrorVector& input) {
    return input.get_field<MagnetometerBias>().array() + input.get_field<MagnetometerScaleFactor>().array() *
        (state.get_field<Attitude>() * UKF::Vector<3>(
            input.get_field<MagneticFieldNorm>() *
                std::cos(std::atan(input.get_field<MagneticFieldInclination>())),
            0.0,
            -input.get_field<MagneticFieldNorm>() *
                std::sin(std::atan(input.get_field<MagneticFieldInclination>())))).array();
}

}

using AHRS_Filter = UKF::SquareRootCore<
    AHRS_StateVector,
    AHRS_MeasurementVector,
    UKF::IntegratorHeun
>;

namespace UKF {
/*
AHRS parameter estimation filter process model. Since the evolution of sensor
errors is by definition unpredictable, this does nothing.
*/
template <> template <>
AHRS_SensorErrorVector AHRS_SensorErrorVector::derivative<>() const {
    return AHRS_SensorErrorVector::Zero();
}

/*
AHRS parameter estimation filter measurement model. These take in the current
state estimate and sensor scale factor and bias estimates, and use them to
calculate predicted measurements.
These functions are just the same as the state measurement model, but with
their arguments flipped.
*/
template <> template <>
UKF::Vector<3> AHRS_MeasurementVector::expected_measurement
<AHRS_SensorErrorVector, Accelerometer, AHRS_StateVector>(
        const AHRS_SensorErrorVector& state, const AHRS_StateVector& input) {
    return state.get_field<AccelerometerBias>() + input.get_field<Attitude>() * UKF::Vector<3>(0, 0, -G_ACCEL);
}

template <> template <>
UKF::Vector<3> AHRS_MeasurementVector::expected_measurement
<AHRS_SensorErrorVector, Gyroscope, AHRS_StateVector>(
        const AHRS_SensorErrorVector& state, const AHRS_StateVector& input) {
    return state.get_field<GyroscopeBias>() + input.get_field<AngularVelocity>();
}

template <> template <>
UKF::Vector<3> AHRS_MeasurementVector::expected_measurement
<AHRS_SensorErrorVector, Magnetometer, AHRS_StateVector>(
        const AHRS_SensorErrorVector& state, const AHRS_StateVector& input) {
    return state.get_field<MagnetometerBias>().array() + state.get_field<MagnetometerScaleFactor>().array() *
        (input.get_field<Attitude>() * UKF::Vector<3>(
            state.get_field<MagneticFieldNorm>() *
                std::cos(std::atan(state.get_field<MagneticFieldInclination>())),
            0.0,
            -state.get_field<MagneticFieldNorm>() *
                std::sin(std::atan(state.get_field<MagneticFieldInclination>())))).array();
}

}

/* Just use the Euler integrator since there's no process model. */
using AHRS_ParameterEstimationFilter = UKF::SquareRootParameterEstimationCore<
    AHRS_SensorErrorVector,
    AHRS_MeasurementVector
>;

static AHRS_Filter ahrs;
static AHRS_ParameterEstimationFilter ahrs_errors;
static AHRS_MeasurementVector meas;
static UKF::Vector<3> acceleration;

/*
The following functions provide a ctypes-compatible interface for ease of
testing.
*/

void ukf_init() {
    /* Initialise state vector and covariance. */
    ahrs.state.set_field<Attitude>(UKF::Quaternion(1, 0, 0, 0));
    ahrs.state.set_field<AngularVelocity>(UKF::Vector<3>(0, 0, 0));
    acceleration << UKF::Vector<3>(0, 0, 0);
    ahrs.root_covariance = AHRS_StateVector::CovarianceMatrix::Zero();
    ahrs.root_covariance.diagonal() <<
        1e0, 1e0, 3.2e0,
        1e-3 * UKF::Vector<3>::Ones();

    /* Set measurement noise covariance. */
    ahrs.measurement_root_covariance <<
        0.5 * UKF::Vector<3>::Ones(),
        0.004 * UKF::Vector<3>::Ones(),
        0.1 * UKF::Vector<3>::Ones();

    /* Set process noise covariance. */
    ahrs.process_noise_root_covariance = AHRS_StateVector::CovarianceMatrix::Zero();
    ahrs.process_noise_root_covariance.diagonal() <<
        5e-5 * UKF::Vector<3>::Ones(),
        5e-3 * UKF::Vector<3>::Ones();

    /* Initialise scale factor and bias errors. */
    ahrs_errors.state.set_field<AccelerometerBias>(UKF::Vector<3>(0, 0, 0));
    ahrs_errors.state.set_field<GyroscopeBias>(UKF::Vector<3>(0, 0, 0));
    ahrs_errors.state.set_field<MagnetometerBias>(UKF::Vector<3>(0, 0, 0));
    ahrs_errors.state.set_field<MagnetometerScaleFactor>(UKF::Vector<3>(1, 1, 1));
    ahrs_errors.state.set_field<MagneticFieldNorm>(MAG_NORM);
    ahrs_errors.state.set_field<MagneticFieldInclination>(0.0);

    /* Initialise scale factor and bias error covariance. */
    ahrs_errors.root_covariance = AHRS_SensorErrorVector::CovarianceMatrix::Zero();
    ahrs_errors.root_covariance.diagonal() <<
        0.8 * UKF::Vector<3>::Ones(),
        0.02 * UKF::Vector<3>::Ones(),
        5.0e-2 * UKF::Vector<3>::Ones(), 1.0e-1 * UKF::Vector<3>::Ones(),
        0.4, 0.7;

    /* Set measurement noise covariance. */
    ahrs_errors.measurement_root_covariance << ahrs.measurement_root_covariance;

    /*
    Set bias error process noise – this is derived from bias instability.

    Bias instability is actually characterised as a 1/f flicker noise rather
    than the white noise (which is what we're specifying using the process
    noise covariance), so these are tuned by hand to values which allow the
    filter to track biases over time, but not change too quickly.
    */
    ahrs_errors.process_noise_root_covariance = AHRS_SensorErrorVector::CovarianceMatrix::Zero();
    ahrs_errors.process_noise_root_covariance.diagonal() <<
        1e-5f * UKF::Vector<3>::Ones(),
        1e-7f * UKF::Vector<3>::Ones(),
        1e-5f * UKF::Vector<3>::Ones(), 1e-6f * UKF::Vector<3>::Ones(),
        1e-7f, 1e-7f;
}

void ukf_set_attitude(real_t w, real_t x, real_t y, real_t z) {
    ahrs.state.set_field<Attitude>(UKF::Quaternion(w, x, y, z));
}

void ukf_set_angular_velocity(real_t x, real_t y, real_t z) {
    ahrs.state.set_field<AngularVelocity>(UKF::Vector<3>(x, y, z));
}

void ukf_get_state(struct ukf_state_t *in) {
    in->attitude[0] = ahrs.state.get_field<Attitude>().x();
    in->attitude[1] = ahrs.state.get_field<Attitude>().y();
    in->attitude[2] = ahrs.state.get_field<Attitude>().z();
    in->attitude[3] = ahrs.state.get_field<Attitude>().w();
    in->angular_velocity[0] = ahrs.state.get_field<AngularVelocity>()[0];
    in->angular_velocity[1] = ahrs.state.get_field<AngularVelocity>()[1];
    in->angular_velocity[2] = ahrs.state.get_field<AngularVelocity>()[2];
    in->acceleration[0] = acceleration[0];
    in->acceleration[1] = acceleration[1];
    in->acceleration[2] = acceleration[2];
}

void ukf_set_state(struct ukf_state_t *in) {
    ahrs.state.set_field<Attitude>(
        UKF::Quaternion(in->attitude[3], in->attitude[0], in->attitude[1], in->attitude[2]));
    ahrs.state.set_field<AngularVelocity>(
        UKF::Vector<3>(in->angular_velocity[0], in->angular_velocity[1], in->angular_velocity[2]));
}

void ukf_get_state_covariance(
        real_t state_covariance[AHRS_StateVector::covariance_size()*AHRS_StateVector::covariance_size()]) {
    Eigen::Map<typename AHRS_StateVector::CovarianceMatrix> covariance_map(state_covariance);
    covariance_map = ahrs.root_covariance * ahrs.root_covariance.transpose();
}

void ukf_get_state_covariance_diagonal(
        real_t state_covariance_diagonal[AHRS_StateVector::covariance_size()]) {
    Eigen::Map<UKF::Vector<AHRS_StateVector::covariance_size()>> covariance_map(state_covariance_diagonal);
    covariance_map = (ahrs.root_covariance * ahrs.root_covariance.transpose()).diagonal();
}

void ukf_get_state_error(struct ukf_state_error_t *in) {
    AHRS_StateVector::StateVectorDelta state_error;
    state_error = (ahrs.root_covariance * ahrs.root_covariance.transpose()).cwiseAbs().rowwise().sum().cwiseSqrt();

    in->attitude[0] = state_error[0];
    in->attitude[1] = state_error[1];
    in->attitude[2] = state_error[2];
    in->angular_velocity[0] = state_error[3];
    in->angular_velocity[1] = state_error[4];
    in->angular_velocity[2] = state_error[5];
}

/*
This assumes accelerometer, gyroscope and magnetometer measurements are
present and in that order.
*/
void ukf_get_innovation(struct ukf_innovation_t *in) {
    in->accel[0] = ahrs.innovation[0];
    in->accel[1] = ahrs.innovation[1];
    in->accel[2] = ahrs.innovation[2];
    in->gyro[0] = ahrs.innovation[3];
    in->gyro[1] = ahrs.innovation[4];
    in->gyro[2] = ahrs.innovation[5];
    in->mag[0] = ahrs.innovation[6];
    in->mag[1] = ahrs.innovation[7];
    in->mag[2] = ahrs.innovation[8];
}

void ukf_sensor_clear() {
    meas = AHRS_MeasurementVector();
}

void ukf_sensor_set_accelerometer(real_t x, real_t y, real_t z) {
    meas.set_field<Accelerometer>(UKF::Vector<3>(x, y, z));
}

void ukf_sensor_set_gyroscope(real_t x, real_t y, real_t z) {
    meas.set_field<Gyroscope>(UKF::Vector<3>(x, y, z));
}

void ukf_sensor_set_magnetometer(real_t x, real_t y, real_t z) {
    meas.set_field<Magnetometer>(UKF::Vector<3>(x, y, z));
}

void ukf_set_params(struct ukf_sensor_params_t *in) {
    ahrs.measurement_root_covariance <<
        std::sqrt(in->accel_covariance[0]), std::sqrt(in->accel_covariance[1]), std::sqrt(in->accel_covariance[2]),
        std::sqrt(in->gyro_covariance[0]), std::sqrt(in->gyro_covariance[1]), std::sqrt(in->gyro_covariance[2]),
        std::sqrt(in->mag_covariance[0]), std::sqrt(in->mag_covariance[1]), std::sqrt(in->mag_covariance[2]);
    ahrs_errors.measurement_root_covariance << ahrs.measurement_root_covariance;
}

void ukf_iterate(float dt) {
    /*
    Split the parameter estimation filter into a priori and a posteriori
    steps, to reduce the CPU load each iteration.
    */
    static int step = 0;

    switch(step++) {
        case 0:
            /*
            The time delta is not used by the parameter estimation filter, so
            there's no need to adjust it.
            */
            ahrs_errors.a_priori_step();
            ahrs_errors.innovation_step(meas, ahrs.state);
            break;
        case 1:
            ahrs_errors.a_posteriori_step();

            /* Clip parameters to physically reasonable values. */
            ahrs_errors.state.set_field<MagneticFieldNorm>(
                std::max(0.2f, std::min(0.7f, ahrs_errors.state.get_field<MagneticFieldNorm>())));

            UKF::Vector<3> temp;
            temp = ahrs_errors.state.get_field<AccelerometerBias>();
            temp[0] = std::max(real_t(-G_ACCEL/4.0), std::min(real_t(G_ACCEL/4.0), temp[0]));
            temp[1] = std::max(real_t(-G_ACCEL/4.0), std::min(real_t(G_ACCEL/4.0), temp[1]));
            temp[2] = std::max(real_t(-G_ACCEL/4.0), std::min(real_t(G_ACCEL/4.0), temp[2]));
            ahrs_errors.state.set_field<AccelerometerBias>(temp);

            temp = ahrs_errors.state.get_field<MagnetometerScaleFactor>();
            temp[0] = std::max(real_t(0.5f), std::min(2.0f, temp[0]));
            temp[1] = std::max(real_t(0.5f), std::min(2.0f, temp[1]));
            temp[2] = std::max(real_t(0.5f), std::min(2.0f, temp[2]));
            ahrs_errors.state.set_field<MagnetometerScaleFactor>(temp);

            step = 0;
            break;
    }

    /*
    Do a normal iteration for the AHRS filter, with the current state of
    the parameter estimation filter as the measurement input.
    */
    ahrs.a_priori_step(dt);
    ahrs.innovation_step(meas, ahrs_errors.state);
    ahrs.a_posteriori_step();

    acceleration << meas.get_field<Accelerometer>() - ahrs_errors.state.get_field<AccelerometerBias>() -
        (ahrs.state.get_field<Attitude>() * UKF::Vector<3>(0, 0, -G_ACCEL));
}

void ukf_set_process_noise(real_t process_noise_covariance[AHRS_StateVector::covariance_size()]) {
    Eigen::Map<typename AHRS_StateVector::StateVectorDelta> covariance_map(process_noise_covariance);
    ahrs.process_noise_root_covariance = AHRS_StateVector::CovarianceMatrix::Zero();
    ahrs.process_noise_root_covariance.diagonal() << covariance_map;
    ahrs.process_noise_root_covariance = ahrs.process_noise_root_covariance.llt().matrixU();
}

void ukf_get_parameters(struct ukf_sensor_errors_t *in) {
    in->accel_bias[0] = ahrs_errors.state.get_field<AccelerometerBias>()[0];
    in->accel_bias[1] = ahrs_errors.state.get_field<AccelerometerBias>()[1];
    in->accel_bias[2] = ahrs_errors.state.get_field<AccelerometerBias>()[2];
    in->gyro_bias[0] = ahrs_errors.state.get_field<GyroscopeBias>()[0];
    in->gyro_bias[1] = ahrs_errors.state.get_field<GyroscopeBias>()[1];
    in->gyro_bias[2] = ahrs_errors.state.get_field<GyroscopeBias>()[2];
    in->mag_bias[0] = ahrs_errors.state.get_field<MagnetometerBias>()[0];
    in->mag_bias[1] = ahrs_errors.state.get_field<MagnetometerBias>()[1];
    in->mag_bias[2] = ahrs_errors.state.get_field<MagnetometerBias>()[2];
    in->mag_scale[0] = ahrs_errors.state.get_field<MagnetometerScaleFactor>()[0];
    in->mag_scale[1] = ahrs_errors.state.get_field<MagnetometerScaleFactor>()[1];
    in->mag_scale[2] = ahrs_errors.state.get_field<MagnetometerScaleFactor>()[2];
    in->mag_field_norm = ahrs_errors.state.get_field<MagneticFieldNorm>();
    in->mag_field_inclination = std::atan(ahrs_errors.state.get_field<MagneticFieldInclination>());
}

void ukf_get_parameters_error(struct ukf_sensor_errors_t *in) {
    AHRS_SensorErrorVector::StateVectorDelta parameters_error;
    parameters_error =
        (ahrs_errors.root_covariance * ahrs_errors.root_covariance.transpose()).cwiseAbs().rowwise().sum().cwiseSqrt();

    in->accel_bias[0] = parameters_error[0];
    in->accel_bias[1] = parameters_error[1];
    in->accel_bias[2] = parameters_error[2];
    in->gyro_bias[0] = parameters_error[3];
    in->gyro_bias[1] = parameters_error[4];
    in->gyro_bias[2] = parameters_error[5];
    in->mag_bias[0] = parameters_error[6];
    in->mag_bias[1] = parameters_error[7];
    in->mag_bias[2] = parameters_error[8];
    in->mag_scale[0] = parameters_error[9];
    in->mag_scale[1] = parameters_error[10];
    in->mag_scale[2] = parameters_error[11];
    in->mag_field_norm = parameters_error[12];
    in->mag_field_inclination = parameters_error[13];
}

uint32_t ukf_config_get_state_dim() {
    return AHRS_StateVector::covariance_size();
}

uint32_t ukf_config_get_measurement_dim() {
    return AHRS_MeasurementVector::max_size();
}

enum ukf_precision_t ukf_config_get_precision() {
    if(sizeof(real_t) == 8) {
        return UKF_PRECISION_DOUBLE;
    } else {
        return UKF_PRECISION_FLOAT;
    }
}
