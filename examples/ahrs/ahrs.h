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

#ifndef AHRS_H
#define AHRS_H

#ifdef __cplusplus
extern "C" {
#endif

#define UKF_STATE_DIM 6
#define UKF_MEASUREMENT_DIM 9

/*
Parameters for sensor model.

Members are as follows:
    - accel_covariance: covariance of accelerometer readings (XYZ)
    - gyro_covariance: covariance of gyro readings (XYZ)
    - mag_covariance: covariance of magnetometer readings (XYZ)
*/
struct ukf_sensor_params_t {
    real_t accel_covariance[3];
    real_t gyro_covariance[3];
    real_t mag_covariance[3];
};

struct ukf_state_t {
    real_t attitude[4]; /* x, y, z, W */
    real_t angular_velocity[3]; /* rolling (rad/s), pitching (rad/s),
                                   yawing (rad/s) */
    real_t acceleration[3]; /* x, y, z (m/s^2) */
};

struct ukf_state_error_t {
    real_t attitude[3]; /* roll, pitch, yaw */
    real_t angular_velocity[3]; /* rolling (rad/s), pitching (rad/s),
                                   yawing (rad/s) */
};

struct ukf_sensor_errors_t {
    real_t accel_bias[3];
    real_t gyro_bias[3];
    real_t mag_bias[3];
    real_t mag_scale[3];
    real_t mag_field_norm;
    real_t mag_field_inclination;
};

struct ukf_innovation_t {
    real_t accel[3];
    real_t gyro[3];
    real_t mag[3];
};

void ukf_init(void);

/*
Note: W, x, y, z in the parameters for ukf_set_attitude differs to the stored
attitude representation in struct ukf_state_t, which is x, y, z, W
*/
void ukf_set_attitude(real_t w, real_t x, real_t y, real_t z);
void ukf_set_angular_velocity(real_t x, real_t y, real_t z);

/* Functions for getting the state vector and covariance. */
void ukf_set_state(struct ukf_state_t *in);
void ukf_get_state(struct ukf_state_t *in);
void ukf_get_state_covariance(
    real_t state_covariance[UKF_STATE_DIM * UKF_STATE_DIM]);
void ukf_get_state_covariance_diagonal(
    real_t state_covariance_diagonal[UKF_STATE_DIM]);
void ukf_get_state_error(struct ukf_state_error_t *in);

void ukf_get_innovation(struct ukf_innovation_t *in);

/*
Functions for setting sensor data. Before each frame, call the sensor_clear()
function to clear old sensor data.
*/
void ukf_sensor_clear(void);
void ukf_sensor_set_accelerometer(real_t x, real_t y, real_t z);
void ukf_sensor_set_gyroscope(real_t x, real_t y, real_t z);
void ukf_sensor_set_magnetometer(real_t x, real_t y, real_t z);

/*
UKF-related functions.
*/
void ukf_set_params(struct ukf_sensor_params_t *in);
void ukf_set_process_noise(real_t process_noise_covariance[UKF_STATE_DIM]);
/* dt is the time delta in seconds */
void ukf_iterate(float dt);

/*
Functions to get the current bias and scale factor estimates from the
parameter estimation filter.
*/
void ukf_get_parameters(struct ukf_sensor_errors_t *in);
void ukf_get_parameters_error(struct ukf_sensor_errors_t *in);

/*
Functions to access the compiled configuration
*/
enum ukf_precision_t {
    UKF_PRECISION_FLOAT = 0,
    UKF_PRECISION_DOUBLE = 1
};

uint32_t ukf_config_get_state_dim(void);
uint32_t ukf_config_get_measurement_dim(void);
enum ukf_precision_t ukf_config_get_precision(void);

#ifdef __cplusplus
}
#endif

#endif
