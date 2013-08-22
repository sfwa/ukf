/*
Copyright (C) 2013 Daniel Dyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef INTERFACE_H
#define INTERFACE_H

#ifdef UKF_SINGLE_PRECISION
typedef float real_t;
#endif

#ifdef UKF_DOUBLE_PRECISION
typedef double real_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Dynamics model types. */
enum ukf_model_t {
    UKF_MODEL_NONE = 0,
    UKF_MODEL_CENTRIPETAL = 1,
    UKF_MODEL_FIXED_WING = 2
};

/*
Parameters for ioboard sensor model.

Members are as follows:
    - accel_orientation: Quaternion from body frame to accelerometer frame.
    - accel_offset: Accelerometer offset from CoG in metres.
    - gyro_orientation: Quaternion from body frame to gyroscope frame.
    - mag_orientation: Quaternion from body frame to magnetometer frame.
    - mag_field: Expected magnetic field in NED frame, µT.

    - accel_covariance: covariance of accelerometer readings (XYZ)
    - gyro_covariance: covariance of gyro readings (XYZ)
    - mag_covariance: covariance of magnetometer readings (XYZ)
    - gps_position_covariance: covariance of GPS position readings (LLH)
    - gps_velocity_covariance: covariance of GPS velocity readings (NED)
    - pitot_covariance: covariance of pitot sensor readings
    - barometer_amsl_covariance: covariance of barometric pressure sensor
        readings
*/
struct ukf_ioboard_params_t {
    real_t accel_orientation[4];
    real_t accel_offset[3];
    real_t gyro_orientation[4];
    real_t mag_orientation[4];
    real_t mag_field[3];

    real_t accel_covariance[3];
    real_t gyro_covariance[3];
    real_t mag_covariance[3];
    real_t gps_position_covariance[3];
    real_t gps_velocity_covariance[3];
    real_t pitot_covariance;
    real_t barometer_amsl_covariance;
};

struct ukf_state_t {
    real_t position[3];
    real_t velocity[3];
    real_t acceleration[3];
    real_t attitude[4]; /* w, x, y, z */
    real_t angular_velocity[3];
    real_t angular_acceleration[3];
    real_t wind_velocity[3];
    real_t gyro_bias[3];
};

/* Functions for setting different parts of the state vector. */
void ukf_set_position(real_t lat, real_t lon, real_t alt);
void ukf_set_velocity(real_t x, real_t y, real_t z);
void ukf_set_acceleration(real_t x, real_t y, real_t z);
void ukf_set_attitude(real_t w, real_t x, real_t y, real_t z);
void ukf_set_angular_velocity(real_t x, real_t y, real_t z);
void ukf_set_angular_acceleration(real_t x, real_t y, real_t z);
void ukf_set_wind_velocity(real_t x, real_t y, real_t z);
void ukf_set_gyro_bias(real_t x, real_t y, real_t z);

/* Functions for getting the state vector and covariance. */
void ukf_set_state(struct ukf_state_t *in);
void ukf_get_state(struct ukf_state_t *in);
void ukf_get_state_covariance(
    real_t state_covariance[UKF_STATE_DIM * UKF_STATE_DIM]);

/*
Functions for setting sensor data. Before each frame, call the sensor_clear()
function to clear old sensor data.
*/
void ukf_sensor_clear();
void ukf_sensor_set_accelerometer(real_t x, real_t y, real_t z);
void ukf_sensor_set_gyroscope(real_t x, real_t y, real_t z);
void ukf_sensor_set_magnetometer(real_t x, real_t y, real_t z);
void ukf_sensor_set_gps_position(real_t lat, real_t lon, real_t alt);
void ukf_sensor_set_gps_velocity(real_t x, real_t y, real_t z);
void ukf_sensor_set_pitot_tas(real_t tas);
void ukf_sensor_set_barometer_amsl(real_t amsl);

/*
UKF-related functions.
*/
void ukf_set_params(struct ukf_ioboard_params_t *in);
void ukf_set_process_noise(real_t process_noise_covariance[UKF_STATE_DIM]);
void ukf_choose_dynamics(enum ukf_model_t t);
void ukf_iterate(float dt, real_t control_vector[UKF_CONTROL_DIM]);

/*
Functions to set airframe properties and coefficients for the fixed-wing
dynamics model.
*/
void ukf_fixedwingdynamics_set_mass(real_t mass);
void ukf_fixedwingdynamics_set_inertia_tensor(real_t inertia_tensor[9]);
void ukf_fixedwingdynamics_set_prop_coeffs(real_t in_prop_area,
    real_t in_prop_cve);
void ukf_fixedwingdynamics_set_drag_coeffs(real_t coeffs[5]);
void ukf_fixedwingdynamics_set_lift_coeffs(real_t coeffs[5]);
void ukf_fixedwingdynamics_set_side_coeffs(real_t coeffs[8],
    real_t control[UKF_CONTROL_DIM]);
void ukf_fixedwingdynamics_set_pitch_moment_coeffs(real_t coeffs[2],
    real_t control[UKF_CONTROL_DIM]);
void ukf_fixedwingdynamics_set_roll_moment_coeffs(real_t coeffs[1],
    real_t control[UKF_CONTROL_DIM]);
void ukf_fixedwingdynamics_set_yaw_moment_coeffs(real_t coeffs[2],
    real_t control[UKF_CONTROL_DIM]);

/*
Functions to access the compiled configuration
*/
enum ukf_precision_t {
    UKF_PRECISION_FLOAT = 0,
    UKF_PRECISION_DOUBLE = 1
};

uint32_t ukf_config_get_state_dim();
uint32_t ukf_config_get_measurement_dim();
uint32_t ukf_config_get_control_dim();
enum ukf_precision_t ukf_config_get_precision();

#ifdef __cplusplus
}
#endif

#endif
