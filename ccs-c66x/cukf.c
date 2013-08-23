#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

#include "config.h"
#include "cukf.h"
#include "cukfmath.h"

/* Disable dynamics model if velocity is less than 1m/s or greater than 100m/s */
#define UKF_DYNAMICS_MIN_V 1.0
#define UKF_DYNAMICS_MAX_V 100.0
/* Disable dynamics model if angular velocity exceeds ~120deg/s */
#define UKF_DYNAMICS_MAX_RATE (M_PI*0.63)
/* Airframe minimums */
#define UKF_AIRFRAME_MIN_MASS 0.1
#define UKF_AIRFRAME_MIN_MOMENT 1e-6

#define G_ACCEL ((real_t)9.80665)
#define RHO ((real_t)1.225)

/* WGS84 reference ellipsoid constants */
#define WGS84_A (6378137.0)
#define WGS84_B (6356752.314245)
#define WGS84_E2 (0.0066943799901975848)
#define WGS84_A2 (WGS84_A*WGS84_A)
#define WGS84_B2 (WGS84_B*WGS84_B)
#define WGS84_AB2 (WGS84_A2*WGS84_B2)

struct _ukf_ioboard_model_t {
    /* See include/sensors.h */
    struct ukf_ioboard_params_t configuration;

    real_t accelerometer[3];
    real_t gyroscope[3];
    real_t magnetometer[3];
    real_t gps_position[3];
    real_t gps_velocity[3];
    real_t pitot_tas;
    real_t barometer_amsl;

    /* True if sensor has data this timestep */
    struct {
        bool accelerometer : 1;
        bool gyroscope : 1;
        bool magnetometer : 1;
        bool gps_position : 1;
        bool gps_velocity : 1;
        bool pitot_tas : 1;
        bool barometer_amsl : 1;
    } flags;
};

struct _fixedwingdynamics_params_t {
    /* See include/dynamics.h */
    real_t prop_area, prop_cve;

    real_t mass_inv;
    real_t inertia_tensor_inv[9];

    real_t c_side_force[8];
    real_t c_side_force_control[UKF_CONTROL_DIM];

    real_t c_pitch_moment[2];
    real_t c_pitch_moment_control[UKF_CONTROL_DIM];

    real_t c_yaw_moment[2];
    real_t c_yaw_moment_control[UKF_CONTROL_DIM];

    real_t c_roll_moment[1];
    real_t c_roll_moment_control[UKF_CONTROL_DIM];

    real_t c_drag_alpha[5], c_lift_alpha[5];

    int8_t motor_idx;
};

/* Sensor and process configuration */
static struct _ukf_ioboard_model_t sensor_model;
static real_t process_noise[UKF_STATE_DIM];

/* UKF state */
static struct ukf_state_t state;
static real_t state_covariance[UKF_STATE_DIM * UKF_STATE_DIM];

/* UKF inputs */
static real_t control[UKF_CONTROL_DIM];

/* Dynamics model configuration */
static enum ukf_model_t dynamics_model = UKF_MODEL_NONE;
static struct _fixedwingdynamics_params_t fixedwing_params;

/* Private interface */
void _ukf_state_model(struct ukf_state_t *in) {
    assert(in);
    /* See src/state.cpp */

    /* Change in position */
    real_t cos_lat = cos(in->position[0]),
           tempA = WGS84_A * cos_lat,
           tempB = WGS84_B * sin(in->position[0]),
           temp = (tempA * tempA + tempB * tempB),
           temp_sqrt_inv = sqrt_inv(temp);
    real_t M = WGS84_AB2 * temp_sqrt_inv / temp;
    real_t N = WGS84_A2 * temp_sqrt_inv;

    in->position[0] = in->velocity[X] / (M + in->position[2]);
    in->position[1] = (in->velocity[Y] / (N + in->position[2])) * cos_lat;
    in->position[2] = -in->velocity[Z];

    /* Change in velocity */
    real_t a[4];
    a[X] = -in->attitude[X];
    a[Y] = -in->attitude[Y];
    a[Z] = -in->attitude[Z];
    a[W] = in->attitude[W];
    _mul_quat_vec3(in->velocity, a, in->acceleration);

    /* No change in acceleration */
    memset(in->acceleration, 0, sizeof(in->acceleration));

    /*
    Change in attitude (XYZW): delta_att = 0.5 * omega_v.conj() * att
    (note a is conjugate of attitude)
    */
    #define Q (in->angular_velocity)
    in->attitude[X] = 0.5 * (-Q[X]*a[W] + Q[Y]*a[Z] - Q[Z]*a[Y]);
    in->attitude[Y] = 0.5 * (-Q[Y]*a[W] + Q[Z]*a[X] - Q[X]*a[Z]);
    in->attitude[Z] = 0.5 * (-Q[Z]*a[W] + Q[X]*a[Y] - Q[Y]*a[X]);
    in->attitude[W] = 0.5 * (-Q[X]*a[X] - Q[Y]*a[Y] - Q[Z]*a[Z]);
    #undef Q

    /* Change in angular velocity */
    memcpy(in->angular_velocity, in->angular_acceleration,
        sizeof(in->angular_velocity));

    /* No change in angular acceleration, wind velocity or gyro bias */
    memset(in->angular_acceleration, 0, sizeof(in->angular_acceleration));
    memset(in->wind_velocity, 0, sizeof(in->wind_velocity));
    memset(in->gyro_bias, 0, sizeof(in->gyro_bias));
}

void _ukf_state_integrate_rk4(struct ukf_state_t *in, real_t delta) {
    assert(in);
    assert(delta >= 0.0);
    if (delta < 1e-6) {
        return;
    }

    /* See include/integrator.h */
    struct ukf_state_t a, b, c, d;

    /* a = in.model() */
    memcpy(&a, in, sizeof(a));
    _ukf_state_model(&a);

    /* b = (in + 0.5 * delta * a).model() */
    memcpy(&b, &a, sizeof(b));
    _mul_vec_scalar_add_vec((real_t*)&b, delta * 0.5, (real_t*)in,
        UKF_STATE_DIM + 1);
    _ukf_state_model(&b);

    /* c = (in + 0.5 * delta * b).model() */
    memcpy(&c, &b, sizeof(c));
    _mul_vec_scalar_add_vec((real_t*)&c, delta * 0.5, (real_t*)in,
        UKF_STATE_DIM + 1);
    _ukf_state_model(&c);

    /* d = (in + delta * c).model */
    memcpy(&d, &c, sizeof(d));
    _mul_vec_scalar_add_vec((real_t*)&d, delta, (real_t*)in,
        UKF_STATE_DIM + 1);
    _ukf_state_model(&d);

    /* in = in + (delta / 6.0) * (a + (b * 2.0) + (c * 2.0) + d) */
    _mul_vec_scalar_add_vec((real_t*)&c, 2.0, (real_t*)&d, UKF_STATE_DIM + 1);
    _mul_vec_scalar_add_vec((real_t*)&b, 2.0, (real_t*)&a, UKF_STATE_DIM + 1);
    _add_vec_vec((real_t*)&b, (real_t*)&c, UKF_STATE_DIM + 1);
    /*
    UKF_STATE_DIM + 1 because this time round we want to copy in's angular
    acceleration, wind velocity and gyro bias to b
    */
    _mul_vec_scalar_add_vec((real_t*)&b, delta / 6.0, (real_t*)in,
        UKF_STATE_DIM + 1);

    memcpy(in, &b, sizeof(b));
}

void _ukf_state_fixed_wing_dynamics(struct ukf_state_t *in) {
    assert(in);
    /* FIXME */
    assert(UKF_CONTROL_DIM == 4);

    /* See src/dynamics.cpp */
    memset(in->acceleration, 0, sizeof(in->acceleration));
    memset(in->angular_acceleration, 0, sizeof(in->angular_acceleration));

    real_t yaw_rate = in->angular_velocity[Z],
           pitch_rate = in->angular_velocity[Y],
           roll_rate = in->angular_velocity[X];

    if (fabs(yaw_rate) > UKF_DYNAMICS_MAX_RATE ||
            fabs(pitch_rate) > UKF_DYNAMICS_MAX_RATE ||
            fabs(roll_rate) > UKF_DYNAMICS_MAX_RATE) {
        return;
    }

    /* External axes */
    real_t ned_airflow[3], airflow[3];
    real_t v2, v_inv;

    ned_airflow[X] = in->wind_velocity[X] - in->velocity[X];
    ned_airflow[Y] = in->wind_velocity[Y] - in->velocity[Y];
    ned_airflow[Z] = in->wind_velocity[Z] - in->velocity[Z];
    _mul_quat_vec3(airflow, in->attitude, ned_airflow);

    v2 = airflow[X]*airflow[X] + airflow[Y]*airflow[Y] + airflow[Z]*airflow[Z];
    if (v2 < UKF_DYNAMICS_MIN_V * UKF_DYNAMICS_MIN_V ||
            v2 > UKF_DYNAMICS_MAX_V * UKF_DYNAMICS_MAX_V) {
        return;
    }

    v_inv = sqrt_inv(v2);

    /* Determine alpha and beta: alpha = atan(wz/wx), beta = atan(wy/|wxz|) */
    real_t alpha, beta, qbar, alpha2, beta2;
    qbar = RHO * v2 * 0.5;
    alpha = atan2(-airflow[Z], -airflow[X]);
    beta = asin(airflow[Y] * v_inv);

    if (fabs(alpha) > M_PI * 0.63 || fabs(beta) > M_PI * 0.63) {
        return;
    }

    alpha2 = alpha * alpha;
    beta2 = beta * beta;

    /* Evaluate quartics in alpha to determine lift and drag */
    real_t lift = fixedwing_params.c_lift_alpha[0]*alpha2*alpha2 +
                  fixedwing_params.c_lift_alpha[1]*alpha2*alpha +
                  fixedwing_params.c_lift_alpha[2]*alpha2 +
                  fixedwing_params.c_lift_alpha[3]*alpha +
                  fixedwing_params.c_lift_alpha[4];
    real_t drag = fixedwing_params.c_drag_alpha[0]*alpha2*alpha2 +
                  fixedwing_params.c_drag_alpha[1]*alpha2*alpha +
                  fixedwing_params.c_drag_alpha[2]*alpha2 +
                  fixedwing_params.c_drag_alpha[3]*alpha +
                  fixedwing_params.c_drag_alpha[4];

    real_t thrust = 0.0;
    if (fixedwing_params.motor_idx < UKF_CONTROL_DIM) {
        real_t ve = fixedwing_params.prop_cve *
                        control[fixedwing_params.motor_idx],
               v0 = airflow[X];
        thrust = 0.5 * RHO * fixedwing_params.prop_area * (ve * ve - v0 * v0);
        if (thrust < 0.0) {
            thrust = 0.0;
        }
    }

    real_t side_force;
    side_force = fixedwing_params.c_side_force[0]*alpha2 +
                 fixedwing_params.c_side_force[1]*alpha +
                 fixedwing_params.c_side_force[2]*beta2 +
                 fixedwing_params.c_side_force[3]*beta +
                 fixedwing_params.c_side_force[4]*alpha2*beta +
                 fixedwing_params.c_side_force[5]*alpha*beta +
                 fixedwing_params.c_side_force[6]*yaw_rate +
                 fixedwing_params.c_side_force[7]*roll_rate +
                 fixedwing_params.c_side_force_control[0]*control[0] +
                 fixedwing_params.c_side_force_control[1]*control[1] +
                 fixedwing_params.c_side_force_control[2]*control[2] +
                 fixedwing_params.c_side_force_control[3]*control[3];

    real_t pitch_moment;
    pitch_moment = fixedwing_params.c_pitch_moment[0]*alpha +
                   fixedwing_params.c_pitch_moment[1]*pitch_rate*pitch_rate*
                        (pitch_rate < 0.0 ? -1.0 : 1.0) +
                   fixedwing_params.c_pitch_moment_control[0]*control[0] +
                   fixedwing_params.c_pitch_moment_control[1]*control[1] +
                   fixedwing_params.c_pitch_moment_control[2]*control[2] +
                   fixedwing_params.c_pitch_moment_control[3]*control[3];

    real_t roll_moment;
    roll_moment = fixedwing_params.c_roll_moment[0]*roll_rate +
                  fixedwing_params.c_roll_moment_control[0]*control[0] +
                  fixedwing_params.c_roll_moment_control[1]*control[1] +
                  fixedwing_params.c_roll_moment_control[2]*control[2] +
                  fixedwing_params.c_roll_moment_control[3]*control[3];

    real_t yaw_moment;
    yaw_moment = fixedwing_params.c_yaw_moment[0]*beta +
                 fixedwing_params.c_yaw_moment[1]*yaw_rate +
                 fixedwing_params.c_yaw_moment_control[0]*control[0] +
                 fixedwing_params.c_yaw_moment_control[1]*control[1] +
                 fixedwing_params.c_yaw_moment_control[2]*control[2] +
                 fixedwing_params.c_yaw_moment_control[3]*control[3];

    /*
    Rotate G_ACCEL by current attitude and add to body-frame force components
    */
    real_t tx = 2.0 * in->attitude[Y] * G_ACCEL,
           ty = 2.0 * -in->attitude[X] * G_ACCEL;

    in->acceleration[X] = (qbar * drag - thrust) * fixedwing_params.mass_inv +
                            (in->attitude[W]*tx - in->attitude[Z]*ty);
    in->acceleration[Y] = (qbar * side_force) * fixedwing_params.mass_inv +
                            (in->attitude[W]*ty + in->attitude[Z]*tx);
    in->acceleration[Z] = -(qbar * lift) * fixedwing_params.mass_inv +
                            (G_ACCEL + in->attitude[X]*ty - in->attitude[Y]*tx);

    /* Calculate angular acceleration (tau / inertia tensor) */
    #define M(x) (fixedwing_params.inertia_tensor_inv[x])
    in->angular_acceleration[X] = qbar *
        (M(0)*roll_moment + M(1)*pitch_moment + M(2)*yaw_moment);
    in->angular_acceleration[Y] = qbar *
        (M(3)*roll_moment + M(4)*pitch_moment + M(5)*yaw_moment);
    in->angular_acceleration[Z] = qbar *
        (M(6)*roll_moment + M(7)*pitch_moment + M(8)*yaw_moment);
    #undef M
}

/* Public interface */
void ukf_set_position(real_t lat, real_t lon, real_t alt) {
    state.position[0] = lat;
    state.position[1] = lon;
    state.position[2] = alt;
}

void ukf_set_velocity(real_t x, real_t y, real_t z) {
    state.velocity[X] = x;
    state.velocity[Y] = y;
    state.velocity[Z] = z;
}

void ukf_set_acceleration(real_t x, real_t y, real_t z) {
    state.acceleration[X] = x;
    state.acceleration[Y] = y;
    state.acceleration[Z] = z;
}

void ukf_set_attitude(real_t w, real_t x, real_t y, real_t z) {
    state.attitude[X] = x;
    state.attitude[Y] = y;
    state.attitude[Z] = z;
    state.attitude[W] = w;
}

void ukf_set_angular_velocity(real_t x, real_t y, real_t z) {
    state.angular_velocity[X] = x;
    state.angular_velocity[Y] = y;
    state.angular_velocity[Z] = z;
}

void ukf_set_angular_acceleration(real_t x, real_t y, real_t z) {
    state.angular_acceleration[X] = x;
    state.angular_acceleration[Y] = y;
    state.angular_acceleration[Z] = z;
}

void ukf_set_wind_velocity(real_t x, real_t y, real_t z) {
    state.wind_velocity[X] = x;
    state.wind_velocity[Y] = y;
    state.wind_velocity[Z] = z;
}

void ukf_set_gyro_bias(real_t x, real_t y, real_t z) {
    state.gyro_bias[X] = x;
    state.gyro_bias[Y] = y;
    state.gyro_bias[Z] = z;
}

void ukf_get_state(struct ukf_state_t *in) {
    assert(in);
    memcpy(in, &state, sizeof(state));
}

void ukf_set_state(struct ukf_state_t *in) {
    assert(in);
    memcpy(&state, in, sizeof(state));
}

void ukf_get_state_covariance(real_t in[UKF_STATE_DIM * UKF_STATE_DIM]) {
    assert(in);
    memcpy(in, state_covariance, sizeof(state_covariance));
}

void ukf_sensor_clear() {
    memset(&sensor_model.flags, 0, sizeof(sensor_model.flags));
}

void ukf_sensor_set_accelerometer(real_t x, real_t y, real_t z) {
    sensor_model.accelerometer[X] = x;
    sensor_model.accelerometer[Y] = y;
    sensor_model.accelerometer[Z] = z;
    sensor_model.flags.accelerometer = true;
}

void ukf_sensor_set_gyroscope(real_t x, real_t y, real_t z) {
    sensor_model.gyroscope[X] = x;
    sensor_model.gyroscope[Y] = y;
    sensor_model.gyroscope[Z] = z;
    sensor_model.flags.gyroscope = true;
}

void ukf_sensor_set_magnetometer(real_t x, real_t y, real_t z) {
    sensor_model.magnetometer[X] = x;
    sensor_model.magnetometer[Y] = y;
    sensor_model.magnetometer[Z] = z;
    sensor_model.flags.magnetometer = true;
}

void ukf_sensor_set_gps_position(real_t lat, real_t lon, real_t alt) {
    sensor_model.gps_position[X] = lat;
    sensor_model.gps_position[Y] = lon;
    sensor_model.gps_position[Z] = alt;
    sensor_model.flags.gps_position = true;
}

void ukf_sensor_set_gps_velocity(real_t x, real_t y, real_t z) {
    sensor_model.gps_velocity[X] = x;
    sensor_model.gps_velocity[Y] = y;
    sensor_model.gps_velocity[Z] = z;
    sensor_model.flags.gps_velocity = true;
}

void ukf_sensor_set_pitot_tas(real_t tas) {
    sensor_model.pitot_tas = tas;
    sensor_model.flags.pitot_tas = true;
}

void ukf_sensor_set_barometer_amsl(real_t amsl) {
    sensor_model.barometer_amsl = amsl;
    sensor_model.flags.barometer_amsl = true;
}

void ukf_set_params(struct ukf_ioboard_params_t *in) {
    assert(in);
    memcpy(&sensor_model.configuration, in,
        sizeof(sensor_model.configuration));
}

void ukf_choose_dynamics(enum ukf_model_t t) {
    assert(t == UKF_MODEL_NONE || t == UKF_MODEL_CENTRIPETAL ||
        t == UKF_MODEL_FIXED_WING);
    dynamics_model = t;
}

void ukf_iterate(float dt, real_t in_control[UKF_CONTROL_DIM]) {
    assert(in_control);
    memcpy(&control, in_control, sizeof(control));

    /* TODO */
    /*
    create_sigma_points();
    apriori_estimate(dt, c);

    if(sensor.size() > 0) {
        measurement_estimate();
        calculate_innovation();
        calculate_kalman_gain();
        aposteriori_estimate();
    } else {
        state = sigma_points.col(0);
        state_covariance = apriori_covariance;
    }

    if(state.attitude()(3) < 0) {
        state.attitude() *= -1.0;
    }
    state.attitude().normalize();
    */
}

void ukf_set_process_noise(real_t in[UKF_STATE_DIM]) {
    assert(in);
    memcpy(&process_noise, in, sizeof(process_noise));
}

void ukf_fixedwingdynamics_set_mass(real_t mass) {
    assert(mass > 0.1);
    fixedwing_params.mass_inv = 1.0 / mass;
}

void ukf_fixedwingdynamics_set_inertia_tensor(real_t in[9]) {
    assert(in);
    _inv_mat3x3(fixedwing_params.inertia_tensor_inv, in);
}

void ukf_fixedwingdynamics_set_prop_coeffs(real_t in_area, real_t in_cve){
    fixedwing_params.prop_area = in_area;
    fixedwing_params.prop_cve = in_cve;
}

void ukf_fixedwingdynamics_set_drag_coeffs(real_t coeffs[5]) {
    assert(coeffs);
    memcpy(fixedwing_params.c_drag_alpha, coeffs,
        sizeof(fixedwing_params.c_drag_alpha));
}

void ukf_fixedwingdynamics_set_lift_coeffs(real_t coeffs[5]) {
    assert(coeffs);
    memcpy(fixedwing_params.c_lift_alpha, coeffs,
        sizeof(fixedwing_params.c_lift_alpha));
}

void ukf_fixedwingdynamics_set_side_coeffs(real_t coeffs[8],
real_t control[UKF_CONTROL_DIM]) {
    assert(coeffs);
    assert(control);
    memcpy(fixedwing_params.c_side_force, coeffs,
        sizeof(fixedwing_params.c_side_force));
    memcpy(fixedwing_params.c_side_force_control, control,
        sizeof(fixedwing_params.c_side_force_control));
}

void ukf_fixedwingdynamics_set_pitch_moment_coeffs(real_t coeffs[2],
real_t control[UKF_CONTROL_DIM]) {
    assert(coeffs);
    assert(control);
    memcpy(fixedwing_params.c_pitch_moment, coeffs,
        sizeof(fixedwing_params.c_pitch_moment));
    memcpy(fixedwing_params.c_pitch_moment_control, control,
        sizeof(fixedwing_params.c_pitch_moment_control));
}

void ukf_fixedwingdynamics_set_roll_moment_coeffs(real_t coeffs[1],
real_t control[UKF_CONTROL_DIM]) {
    assert(coeffs);
    assert(control);
    memcpy(fixedwing_params.c_roll_moment, coeffs,
        sizeof(fixedwing_params.c_roll_moment));
    memcpy(fixedwing_params.c_roll_moment_control, control,
        sizeof(fixedwing_params.c_roll_moment_control));
}

void ukf_fixedwingdynamics_set_yaw_moment_coeffs(real_t coeffs[2],
real_t control[UKF_CONTROL_DIM]) {
    assert(coeffs);
    assert(control);
    memcpy(fixedwing_params.c_yaw_moment, coeffs,
        sizeof(fixedwing_params.c_yaw_moment));
    memcpy(fixedwing_params.c_yaw_moment_control, control,
        sizeof(fixedwing_params.c_yaw_moment_control));
}

uint32_t ukf_config_get_state_dim() {
    return UKF_STATE_DIM;
}

uint32_t ukf_config_get_measurement_dim() {
    return UKF_MEASUREMENT_DIM;
}

uint32_t ukf_config_get_control_dim() {
    return UKF_CONTROL_DIM;
}

enum ukf_precision_t ukf_config_get_precision() {
#ifdef UKF_SINGLE_PRECISION
    return UKF_PRECISION_FLOAT;
#else
    return UKF_PRECISION_DOUBLE;
#endif
}
