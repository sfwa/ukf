/*
Copyright (C) 2016 Daniel Dyer

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

#ifdef __cplusplus
extern "C" {
#endif

#define UKF_STATE_DIM 24
#define UKF_MEASUREMENT_DIM 17
#define UKF_CONTROL_DIM 4

/* Dynamics model types. */
enum ukf_model_t {
    UKF_MODEL_NONE = 0,
    UKF_MODEL_CENTRIPETAL = 1,
    UKF_MODEL_CUSTOM = 2,
    UKF_MODEL_X8 = 3
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
    real_t accel_orientation[4]; /* x, y, z, W */
    real_t accel_offset[3]; /* forwards, starboard, down from CoG (m) */
    real_t gyro_orientation[4]; /* x, y, z, W */
    real_t mag_orientation[4]; /* x, y, z, W */
    real_t mag_field[3]; /* North, East, Down (µT) */

    real_t accel_covariance[3];
    real_t gyro_covariance[3];
    real_t mag_covariance[3];
    real_t gps_position_covariance[3];
    real_t gps_velocity_covariance[3];
    real_t pitot_covariance;
    real_t barometer_amsl_covariance;
};

struct ukf_state_t {
    real_t position[3]; /* lat (rad), lon (rad), alt (m above ellipsoid) */
    real_t velocity[3]; /* North (m/s), East (m/s), Down (m/s) */
    real_t acceleration[3]; /* forwards (m/s), starboard (m/s), down (m/s) */
    real_t attitude[4]; /* x, y, z, W */
    real_t angular_velocity[3]; /* rolling (rad/s), pitching (rad/s),
                                   yawing (rad/s) */
    real_t angular_acceleration[3];
    real_t wind_velocity[3]; /* North (m/s), East (m/s), Down (m/s) */
    real_t gyro_bias[3]; /* X (rad/s), Y (rad/s), Z (rad/s) */
};

/*
Dynamics model function, for custom model support. These functions have to
be compatible between the C++ and C versions of the UKF, so they take pointers
to C arrays representing the state vector, the control vector, and the output
vector.
*/
typedef void (*ukf_model_function_t)(const real_t *, const real_t *,
                                     real_t *);


void ukf_init(void);

/* Functions for setting different parts of the state vector. */
void ukf_set_position(real_t lat, real_t lon, real_t alt);
void ukf_set_velocity(real_t x, real_t y, real_t z);
void ukf_set_acceleration(real_t x, real_t y, real_t z);
/*
Note: W, x, y, z in the parameters for ukf_set_attitude differs to the stored
attitude representation in struct ukf_state_t, which is x, y, z, W
*/
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
void ukf_get_state_covariance_diagonal(
    real_t state_covariance_diagonal[UKF_STATE_DIM]);
void ukf_get_state_error(real_t state_error[UKF_STATE_DIM]);

/*
Functions for setting sensor data. Before each frame, call the sensor_clear()
function to clear old sensor data.
*/
void ukf_sensor_clear(void);
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
void ukf_set_custom_dynamics_model(ukf_model_function_t func);
/* dt is the time delta in seconds */
void ukf_iterate(float dt, real_t control_vector[UKF_CONTROL_DIM]);

/*
Functions to access the compiled configuration
*/
enum ukf_precision_t {
    UKF_PRECISION_FLOAT = 0,
    UKF_PRECISION_DOUBLE = 1
};

uint32_t ukf_config_get_state_dim(void);
uint32_t ukf_config_get_measurement_dim(void);
uint32_t ukf_config_get_control_dim(void);
enum ukf_precision_t ukf_config_get_precision(void);

#ifdef __cplusplus
}
#endif

#endif
