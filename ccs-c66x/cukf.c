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

/* See include/ukf.h */
#define UKF_NUM_SIGMA (2*UKF_STATE_DIM + 1)

#define UKF_ALPHA_2 (1.0)
#define UKF_BETA (0.0)
#define UKF_KAPPA (3.0)
#define UKF_LAMBDA (UKF_ALPHA_2*(UKF_STATE_DIM + UKF_KAPPA) - UKF_STATE_DIM)
#define UKF_DIM_PLUS_LAMBDA (UKF_ALPHA_2*(UKF_STATE_DIM + UKF_KAPPA))

#define UKF_MRP_A (1.0)
#define UKF_MRP_A_2 (UKF_MRP_A*UKF_MRP_A)
#define UKF_MRP_F (2.0*(UKF_MRP_A + 1))
#define UKF_MRP_F_2 (UKF_MRP_F*UKF_MRP_F)

#define UKF_SIGMA_WM0 (UKF_LAMBDA/(UKF_DIM_PLUS_LAMBDA))
#define UKF_SIGMA_WC0 (UKF_SIGMA_WM0 + (1.0 - UKF_ALPHA_2 + UKF_BETA))
#define UKF_SIGMA_WMI (1.0/(2.0*(UKF_DIM_PLUS_LAMBDA)))
#define UKF_SIGMA_WCI (UKF_SIGMA_WMI)

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
static real_t sensor_model_mag_field_norm;

/* UKF state */
static struct ukf_state_t state;
static real_t state_covariance[UKF_STATE_DIM * UKF_STATE_DIM]; /* 4608B */

/*
UKF temporaries -- global to avoid blowing up function stack.
w_prime is 9408 bytes, measurement_estimate_sigma is 7840 bytes,
measurement_estimate_covariance is 3200 bytes, cross_correlation is 3840 bytes
*/
real_t w_prime[UKF_STATE_DIM * UKF_NUM_SIGMA];
real_t measurement_estimate_sigma[UKF_MEASUREMENT_DIM * UKF_NUM_SIGMA];
real_t measurement_estimate_covariance[UKF_MEASUREMENT_DIM *
                                       UKF_MEASUREMENT_DIM];
real_t cross_correlation[UKF_STATE_DIM * UKF_MEASUREMENT_DIM];

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
           temp = tempA * tempA + tempB * tempB,
           temp_sqrt_inv = sqrt_inv(temp);
    real_t M = WGS84_AB2 * (temp_sqrt_inv / temp);
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
    Change in attitude (XYZW): delta_att = 0.5 * (omega_v.conj() * att)
    */
    a[X] = -a[X];
    a[Y] = -a[Y];
    a[Z] = -a[Z];
    real_t omega_q_conj[4] = {
        -in->angular_velocity[X],
        -in->angular_velocity[Y],
        -in->angular_velocity[Z],
        0
    };
    _mul_quat_quat(in->attitude, omega_q_conj, a);
    in->attitude[X] *= 0.5;
    in->attitude[Y] *= 0.5;
    in->attitude[Z] *= 0.5;
    in->attitude[W] *= 0.5;

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
    if (delta < 1e-6) {
        return;
    }

    /* See include/integrator.h */
    struct ukf_state_t a, b, c, d;

    /* a = in.model() */
    memcpy(&a, in, sizeof(a));
    _ukf_state_model(&a);

    /* b = (in + 0.5 * delta * a).model() */
    _mul_state_scalar_add_state(&b, &a, delta * 0.5, in);
    _ukf_state_model(&b);

    /* c = (in + 0.5 * delta * b).model() */
    _mul_state_scalar_add_state(&c, &b, delta * 0.5, in);
    _ukf_state_model(&c);

    /* d = (in + delta * c).model */
    _mul_state_scalar_add_state(&d, &c, delta, in);
    _ukf_state_model(&d);

    /* in = in + (delta / 6.0) * (a + (b * 2.0) + (c * 2.0) + d) */
    delta = delta / 6.0;
    in->position[0] += delta * (a.position[0] + b.position[0] * 2.0 +
        c.position[0] * 2.0 + d.position[0]);
    in->position[1] += delta * (a.position[1] + b.position[1] * 2.0 +
        c.position[1] * 2.0 + d.position[1]);
    in->position[2] += delta * (a.position[2] + b.position[2] * 2.0 +
        c.position[2] * 2.0 + d.position[2]);
    in->velocity[0] += delta * (a.velocity[0] + b.velocity[0] * 2.0 +
        c.velocity[0] * 2.0 + d.velocity[0]);
    in->velocity[1] += delta * (a.velocity[1] + b.velocity[1] * 2.0 +
        c.velocity[1] * 2.0 + d.velocity[1]);
    in->velocity[2] += delta * (a.velocity[2] + b.velocity[2] * 2.0 +
        c.velocity[2] * 2.0 + d.velocity[2]);
    in->attitude[0] += delta * (a.attitude[0] + b.attitude[0] * 2.0 +
        c.attitude[0] * 2.0 + d.attitude[0]);
    in->attitude[1] += delta * (a.attitude[1] + b.attitude[1] * 2.0 +
        c.attitude[1] * 2.0 + d.attitude[1]);
    in->attitude[2] += delta * (a.attitude[2] + b.attitude[2] * 2.0 +
        c.attitude[2] * 2.0 + d.attitude[2]);
    in->attitude[3] += delta * (a.attitude[3] + b.attitude[3] * 2.0 +
        c.attitude[3] * 2.0 + d.attitude[3]);
    in->angular_velocity[0] += delta * (a.angular_velocity[0] +
        b.angular_velocity[0] * 2.0 + c.angular_velocity[0] * 2.0 +
        d.angular_velocity[0]);
    in->angular_velocity[1] += delta * (a.angular_velocity[1] +
        b.angular_velocity[1] * 2.0 + c.angular_velocity[1] * 2.0 +
        d.angular_velocity[1]);
    in->angular_velocity[2] += delta * (a.angular_velocity[2] +
        b.angular_velocity[2] * 2.0 + c.angular_velocity[2] * 2.0 +
        d.angular_velocity[2]);
}

void _ukf_state_centripetal_dynamics(struct ukf_state_t *in,
real_t control[UKF_CONTROL_DIM]) {
    assert(in && control);

    real_t velocity_body[3];
    _mul_quat_vec3(velocity_body, in->attitude, in->velocity);
    _cross_vec3(in->acceleration, in->angular_velocity, velocity_body);
    memset(in->angular_acceleration, 0, sizeof(in->angular_acceleration));
}

void _ukf_state_fixed_wing_dynamics(struct ukf_state_t *in,
real_t control[UKF_CONTROL_DIM]) {
    assert(in && control);
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
    real_t g_accel[3], g[3] = { 0, 0, G_ACCEL };
    _mul_quat_vec3(g_accel, in->attitude, g);

    in->acceleration[X] = (thrust - qbar * drag) * fixedwing_params.mass_inv +
                            g_accel[X];
    in->acceleration[Y] = (qbar * side_force) * fixedwing_params.mass_inv +
                            g_accel[Y];
    in->acceleration[Z] = -(qbar * lift) * fixedwing_params.mass_inv +
                            g_accel[Z];

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

    #define B sensor_model.configuration.mag_field
    sensor_model_mag_field_norm = sqrt(B[X]*B[X] + B[Y]*B[Y] + B[Z]*B[Z]);
    #undef B
}

void ukf_choose_dynamics(enum ukf_model_t t) {
    assert(t == UKF_MODEL_NONE || t == UKF_MODEL_CENTRIPETAL ||
        t == UKF_MODEL_FIXED_WING);
    dynamics_model = t;
}

void _ukf_sensor_predict(real_t measurement_estimate[UKF_MEASUREMENT_DIM],
struct ukf_state_t *sigma) {
    assert(measurement_estimate && sigma);
    assert(fabs(sigma->attitude[X]*sigma->attitude[X] +
                sigma->attitude[Y]*sigma->attitude[Y] +
                sigma->attitude[Z]*sigma->attitude[Z] +
                sigma->attitude[W]*sigma->attitude[W] - 1.0) < 1e-6);

    /* see src/sensors.cpp line 123 */

    size_t i = 0;

    if (sensor_model.flags.accelerometer) {
        real_t temp1[3], temp2[3];

        _cross_vec3(temp1, sigma->angular_velocity,
                    sensor_model.configuration.accel_offset);
        _cross_vec3(temp2, sigma->angular_acceleration,
                    sensor_model.configuration.accel_offset);
        temp1[X] += temp2[X] + sigma->acceleration[X];
        temp1[Y] += temp2[Y] + sigma->acceleration[Y];
        temp1[Z] += temp2[Z] + sigma->acceleration[Z];

        _mul_quat_vec3(&measurement_estimate[i],
                       sensor_model.configuration.accel_orientation,
                       temp1);
        i += 3;

        real_t g[3] = { 0, 0, -G_ACCEL };
        _mul_quat_vec3(temp1, sigma->attitude, g);
        _mul_quat_vec3(&measurement_estimate[i],
                       sensor_model.configuration.accel_orientation,
                       temp1);
        i += 3;
    }

    if (sensor_model.flags.gyroscope) {
        real_t net_av[3] = {
            sigma->angular_velocity[X] + sigma->gyro_bias[X],
            sigma->angular_velocity[Y] + sigma->gyro_bias[Y],
            sigma->angular_velocity[Z] + sigma->gyro_bias[Z]
        };
        _mul_quat_vec3(&measurement_estimate[i],
                       sensor_model.configuration.gyro_orientation,
                       net_av);
        i += 3;
    }

    if (sensor_model.flags.magnetometer) {
        real_t temp[3];
        _mul_quat_vec3(temp, sigma->attitude,
                       sensor_model.configuration.mag_field);
        _mul_quat_vec3(&measurement_estimate[i],
                       sensor_model.configuration.mag_orientation,
                       temp);
        i += 3;
    }

    if (sensor_model.flags.gps_position) {
        measurement_estimate[i++] = sigma->position[0];
        measurement_estimate[i++] = sigma->position[1];
        measurement_estimate[i++] = sigma->position[2];
    }

    if (sensor_model.flags.gps_velocity) {
        measurement_estimate[i++] = sigma->velocity[0];
        measurement_estimate[i++] = sigma->velocity[1];
        measurement_estimate[i++] = sigma->velocity[2];
    }

    if (sensor_model.flags.pitot_tas) {
        /* predicted = (attitude * (in.velocity() - in.wind_velocity()))[X] */
        real_t airflow[3] = {
            sigma->velocity[X] - sigma->wind_velocity[X],
            sigma->velocity[Y] - sigma->wind_velocity[Y],
            sigma->velocity[Z] - sigma->wind_velocity[Z]
        }, temp[3];

        _mul_quat_vec3(temp, sigma->attitude, airflow);

        measurement_estimate[i++] = temp[X];
    }

    if (sensor_model.flags.barometer_amsl) {
        measurement_estimate[i++] = sigma->position[2];
    }
}

size_t _ukf_sensor_collate(real_t measurement_estimate[UKF_MEASUREMENT_DIM]) {
    assert(measurement_estimate);

    size_t i = 0;

    if (sensor_model.flags.accelerometer) {
        measurement_estimate[i++] = sensor_model.accelerometer[X];
        measurement_estimate[i++] = sensor_model.accelerometer[Y];
        measurement_estimate[i++] = sensor_model.accelerometer[Z];
    }

    if (sensor_model.flags.gyroscope) {
        measurement_estimate[i++] = sensor_model.gyroscope[X];
        measurement_estimate[i++] = sensor_model.gyroscope[Y];
        measurement_estimate[i++] = sensor_model.gyroscope[Z];
    }

    if (sensor_model.flags.magnetometer) {
        measurement_estimate[i++] = sensor_model.magnetometer[X];
        measurement_estimate[i++] = sensor_model.magnetometer[Y];
        measurement_estimate[i++] = sensor_model.magnetometer[Z];
    }

    if (sensor_model.flags.gps_position) {
        measurement_estimate[i++] = sensor_model.gps_position[X];
        measurement_estimate[i++] = sensor_model.gps_position[Y];
        measurement_estimate[i++] = sensor_model.gps_position[Z];
    }

    if (sensor_model.flags.gps_velocity) {
        measurement_estimate[i++] = sensor_model.gps_velocity[X];
        measurement_estimate[i++] = sensor_model.gps_velocity[Y];
        measurement_estimate[i++] = sensor_model.gps_velocity[Z];
    }

    if (sensor_model.flags.pitot_tas) {
        measurement_estimate[i++] = sensor_model.pitot_tas;
    }

    if (sensor_model.flags.barometer_amsl) {
        measurement_estimate[i++] = sensor_model.barometer_amsl;
    }

    return i;
}

void _ukf_sensor_get_covariance(real_t covariance[UKF_MEASUREMENT_DIM]) {
    assert(covariance);

    size_t i = 0;

    #define C sensor_model.configuration
    if (sensor_model.flags.accelerometer) {
        covariance[i++] = C.accel_covariance[X];
        covariance[i++] = C.accel_covariance[Y];
        covariance[i++] = C.accel_covariance[Z];
    }

    if (sensor_model.flags.gyroscope) {
        covariance[i++] = C.gyro_covariance[X];
        covariance[i++] = C.gyro_covariance[Y];
        covariance[i++] = C.gyro_covariance[Z];
    }

    if (sensor_model.flags.magnetometer) {
        covariance[i++] = C.mag_covariance[X];
        covariance[i++] = C.mag_covariance[Y];
        covariance[i++] = C.mag_covariance[Z];
    }

    if (sensor_model.flags.gps_position) {
        covariance[i++] = C.gps_position_covariance[X];
        covariance[i++] = C.gps_position_covariance[Y];
        covariance[i++] = C.gps_position_covariance[Z];
    }

    if (sensor_model.flags.gps_velocity) {
        covariance[i++] = C.gps_velocity_covariance[X];
        covariance[i++] = C.gps_velocity_covariance[Y];
        covariance[i++] = C.gps_velocity_covariance[Z];
    }

    if (sensor_model.flags.pitot_tas) {
        covariance[i++] = C.pitot_covariance;
    }

    if (sensor_model.flags.barometer_amsl) {
        covariance[i++] = C.barometer_amsl_covariance;
    }
    #undef C
}

void _ukf_sensor_calculate_deltas(
real_t measurement_estimate[UKF_MEASUREMENT_DIM],
real_t measurement_estimate_mean[UKF_MEASUREMENT_DIM - 3],
uint32_t sigma_idx) {
    assert(measurement_estimate && measurement_estimate_mean &&
           sigma_idx < UKF_NUM_SIGMA);

    #define E measurement_estimate
    #define M measurement_estimate_mean
    #define sigma (&measurement_estimate_sigma[sigma_idx*UKF_MEASUREMENT_DIM])
    if (sensor_model.flags.accelerometer) {
        E[0] = sigma[0] + sigma[3] - M[0];
        E[1] = sigma[1] + sigma[4] - M[1];
        E[2] = sigma[2] + sigma[5] - M[2];
        E[3] = sigma[6] - M[3];
        E[4] = sigma[7] - M[4];
        E[5] = sigma[8] - M[5];
        E[6] = sigma[9] - M[6];
        E[7] = sigma[10] - M[7];
        E[8] = sigma[11] - M[8];
        E[9] = sigma[12] - M[9];
        E[10] = sigma[13] - M[10];
        E[11] = sigma[14] - M[11];
        E[12] = sigma[15] - M[12];
        E[13] = sigma[16] - M[13];
        E[14] = sigma[17] - M[14];
        E[15] = sigma[18] - M[15];
        E[16] = sigma[19] - M[16];
    } else {
        E[0] = sigma[0] - M[0];
        E[1] = sigma[1] - M[1];
        E[2] = sigma[2] - M[2];
        E[3] = sigma[3] - M[3];
        E[4] = sigma[4] - M[4];
        E[5] = sigma[5] - M[5];
        E[6] = sigma[6] - M[6];
        E[7] = sigma[7] - M[7];
        E[8] = sigma[8] - M[8];
        E[9] = sigma[9] - M[9];
        E[10] = sigma[10] - M[10];
        E[11] = sigma[11] - M[11];
        E[12] = sigma[12] - M[12];
        E[13] = sigma[13] - M[13];
        E[14] = sigma[14] - M[14];
        E[15] = sigma[15] - M[15];
        E[16] = sigma[16] - M[16];
        E[17] = sigma[17] - M[17];
        E[18] = sigma[18] - M[18];
        E[19] = sigma[19] - M[19];
    }
    #undef E
    #undef M
    #undef sigma
}

void _ukf_process_sigma(struct ukf_state_t *sigma, uint32_t idx, real_t dt,
real_t control[4], struct ukf_state_t *apriori_mean,
real_t w_prime[UKF_STATE_DIM * UKF_NUM_SIGMA],
real_t measurement_estimate_mean[UKF_MEASUREMENT_DIM],
real_t measurement_estimate[UKF_MEASUREMENT_DIM],
struct ukf_state_t *central_sigma) {
    assert(sigma && control && apriori_mean && w_prime &&
        measurement_estimate_mean && measurement_estimate && control);

    /*
    This function handles a number of sequential steps from ukf.cpp in an
    iterative manner, taking each sigma point, applying the kinematics and
    process models, saving the difference between the sigma point and the
    central point to w_prime, and predicting and averaging the measurement
    estimates
    */

    /* apriori_estimate -- src/ukf.cpp line 149 */
    _ukf_state_integrate_rk4(sigma, dt);
    _normalize_quat(sigma->attitude, sigma->attitude, false);

    /* src/ukf.cpp line 165 */
    if (dynamics_model == UKF_MODEL_FIXED_WING) {
        _ukf_state_fixed_wing_dynamics(sigma, control);
    } else if (dynamics_model == UKF_MODEL_CENTRIPETAL) {
        _ukf_state_centripetal_dynamics(sigma, control);
    } else {
        sigma->angular_acceleration[X] = 0;
        sigma->angular_acceleration[Y] = 0;
        sigma->angular_acceleration[Z] = 0;
    }

    /* save offset from central point to w_prime, src/ukf.cpp line 192 */
    _add_state_accum(apriori_mean, sigma);

    if (idx == 0) {
        /* Processing the centre point, so w_prime is 0 */
        memset(w_prime, 0, sizeof(real_t) * UKF_STATE_DIM);
    } else {
        #define Sd(x) (sigma->x - central_sigma->x)
        w_prime[idx*UKF_STATE_DIM + 0] = Sd(position[0]);
        w_prime[idx*UKF_STATE_DIM + 1] = Sd(position[1]);
        w_prime[idx*UKF_STATE_DIM + 2] = Sd(position[2]);
        w_prime[idx*UKF_STATE_DIM + 3] = Sd(velocity[X]);
        w_prime[idx*UKF_STATE_DIM + 4] = Sd(velocity[Y]);
        w_prime[idx*UKF_STATE_DIM + 5] = Sd(velocity[Z]);
        w_prime[idx*UKF_STATE_DIM + 6] = Sd(acceleration[X]);
        w_prime[idx*UKF_STATE_DIM + 7] = Sd(acceleration[Y]);
        w_prime[idx*UKF_STATE_DIM + 8] = Sd(acceleration[Z]);

        w_prime[idx*UKF_STATE_DIM + 12] = Sd(angular_velocity[X]);
        w_prime[idx*UKF_STATE_DIM + 13] = Sd(angular_velocity[Y]);
        w_prime[idx*UKF_STATE_DIM + 14] = Sd(angular_velocity[Z]);
        w_prime[idx*UKF_STATE_DIM + 15] = Sd(angular_acceleration[X]);
        w_prime[idx*UKF_STATE_DIM + 16] = Sd(angular_acceleration[Y]);
        w_prime[idx*UKF_STATE_DIM + 17] = Sd(angular_acceleration[Z]);
        w_prime[idx*UKF_STATE_DIM + 18] = Sd(wind_velocity[X]);
        w_prime[idx*UKF_STATE_DIM + 19] = Sd(wind_velocity[Y]);
        w_prime[idx*UKF_STATE_DIM + 20] = Sd(wind_velocity[Z]);
        w_prime[idx*UKF_STATE_DIM + 21] = Sd(gyro_bias[X]);
        w_prime[idx*UKF_STATE_DIM + 22] = Sd(gyro_bias[Y]);
        w_prime[idx*UKF_STATE_DIM + 23] = Sd(gyro_bias[Z]);
        #undef Sd
    }

    /* src/ukf.cpp line 213 */
    if (idx == 0) {
        w_prime[idx*UKF_STATE_DIM + 9] = 0;
        w_prime[idx*UKF_STATE_DIM + 10] = 0;
        w_prime[idx*UKF_STATE_DIM + 11] = 0;
    } else {
        real_t err_q[4], cs_conj[4] = {
            -central_sigma->attitude[X],
            -central_sigma->attitude[Y],
            -central_sigma->attitude[Z],
            central_sigma->attitude[W]
        };
        _mul_quat_quat(err_q, sigma->attitude, cs_conj);

        real_t qw_recip = 1.0 / (UKF_MRP_A + err_q[W]);
        w_prime[idx*UKF_STATE_DIM + 9] = UKF_MRP_F * err_q[X] * qw_recip;
        w_prime[idx*UKF_STATE_DIM + 10] = UKF_MRP_F * err_q[Y] * qw_recip;
        w_prime[idx*UKF_STATE_DIM + 11] = UKF_MRP_F * err_q[Z] * qw_recip;
    }

    /* measurement_estimate, src/ukf.cpp line 250 */
    _ukf_sensor_predict(measurement_estimate, sigma);

    /*
    mean calculation -- straight weighted arithmetic mean for most sensors,
    but magnetometer and acceleration due to gravity get normalised based on
    expected field strength

    src/ukf.cpp line 258; src/sensors.cpp line 265
    */
    measurement_estimate_mean[0] += measurement_estimate[0];
    measurement_estimate_mean[1] += measurement_estimate[1];
    measurement_estimate_mean[2] += measurement_estimate[2];
    measurement_estimate_mean[3] += measurement_estimate[3];
    measurement_estimate_mean[4] += measurement_estimate[4];
    measurement_estimate_mean[5] += measurement_estimate[5];
    measurement_estimate_mean[6] += measurement_estimate[6];
    measurement_estimate_mean[7] += measurement_estimate[7];
    measurement_estimate_mean[8] += measurement_estimate[8];
    measurement_estimate_mean[9] += measurement_estimate[9];
    measurement_estimate_mean[10] += measurement_estimate[10];
    measurement_estimate_mean[11] += measurement_estimate[11];
    measurement_estimate_mean[12] += measurement_estimate[12];
    measurement_estimate_mean[13] += measurement_estimate[13];
    measurement_estimate_mean[14] += measurement_estimate[14];
    measurement_estimate_mean[15] += measurement_estimate[15];
    measurement_estimate_mean[16] += measurement_estimate[16];
    measurement_estimate_mean[17] += measurement_estimate[17];
    measurement_estimate_mean[18] += measurement_estimate[18];
    measurement_estimate_mean[19] += measurement_estimate[19];
}

void ukf_iterate(float dt, real_t control[UKF_CONTROL_DIM]) {
    assert(control);
    assert(UKF_STATE_DIM == 24);

/*
    uint64_t t = rdtsc();
*/

    /* See src/ukf.cpp, line 85 */
    state_covariance[0] += process_noise[0];
    state_covariance[25] += process_noise[1];
    state_covariance[50] += process_noise[2];
    state_covariance[75] += process_noise[3];
    state_covariance[100] += process_noise[4];
    state_covariance[125] += process_noise[5];
    state_covariance[150] += process_noise[6];
    state_covariance[175] += process_noise[7];
    state_covariance[200] += process_noise[8];
    state_covariance[225] += process_noise[9];
    state_covariance[250] += process_noise[10];
    state_covariance[275] += process_noise[11];
    state_covariance[300] += process_noise[12];
    state_covariance[325] += process_noise[13];
    state_covariance[350] += process_noise[14];
    state_covariance[375] += process_noise[15];
    state_covariance[400] += process_noise[16];
    state_covariance[425] += process_noise[17];
    state_covariance[450] += process_noise[18];
    state_covariance[475] += process_noise[19];
    state_covariance[500] += process_noise[20];
    state_covariance[525] += process_noise[21];
    state_covariance[550] += process_noise[22];
    state_covariance[575] += process_noise[23];

    /* In-place LLT on state covariance matrix times UKF_DIM_PLUS_LAMBDA */
    _cholesky_mat_mul(state_covariance, state_covariance, UKF_DIM_PLUS_LAMBDA,
        24);
    /* Clear out the upper triangle */
    uint32_t i, j;
    for (i = 0; i < UKF_STATE_DIM; i++) {
        for (j = i + 1; j < UKF_STATE_DIM; j++) {
            state_covariance[j*UKF_STATE_DIM + i] = 0.0;
        }
    }

/*
    printf("S:\n");
    _print_matrix(state_covariance, UKF_STATE_DIM, UKF_STATE_DIM);
*/

    /*
    Keep central sigma point means separate from the rest, since the high
    weight will otherwise reduce precision
    */
    struct ukf_state_t sigma, central_sigma, apriori_mean, apriori_central;
    real_t measurement_estimate_mean[UKF_MEASUREMENT_DIM],
           measurement_estimate_central[UKF_MEASUREMENT_DIM];

    memset(&apriori_mean, 0, sizeof(apriori_mean));
    memset(&apriori_central, 0, sizeof(apriori_central));
    memset(measurement_estimate_mean, 0, sizeof(measurement_estimate_mean));
    memset(measurement_estimate_central, 0,
           sizeof(measurement_estimate_central));

    /* Central point */
    memcpy(&central_sigma, &state, sizeof(central_sigma));
    _ukf_process_sigma(&central_sigma, 0, dt, control, &apriori_central,
        w_prime, measurement_estimate_central,
        &measurement_estimate_sigma[0*UKF_MEASUREMENT_DIM], &central_sigma);

    /* Other points */
    uint32_t col;
    for (i = 0, col = 0; i < 24; i++, col = i * 24) {
        /* src/ukf.cpp line 103  */
        real_t d_p[3] = {
            state_covariance[col + 9],
            state_covariance[col + 10],
            state_covariance[col + 11]
        };
        real_t x_2 = (d_p[X]*d_p[X] + d_p[Y]*d_p[Y] + d_p[Z]*d_p[Z]);
        real_t err_w = (-UKF_MRP_A * x_2 +
            UKF_MRP_F * sqrt(UKF_MRP_F_2 + (1.0 - UKF_MRP_A_2) * x_2)) /
            (UKF_MRP_F_2 + x_2);
        real_t err_scale = (1.0 / UKF_MRP_F) * (UKF_MRP_A + err_w);
        real_t noise_q[4] = {
            err_scale * d_p[X], err_scale * d_p[Y], err_scale * d_p[Z], err_w
        };

        /* Create positive sigma point. */
        sigma.position[0] = state.position[0] + state_covariance[col + 0];
        sigma.position[1] = state.position[1] + state_covariance[col + 1];
        sigma.position[2] = state.position[2] + state_covariance[col + 2];
        sigma.velocity[0] = state.velocity[0] + state_covariance[col + 3];
        sigma.velocity[1] = state.velocity[1] + state_covariance[col + 4];
        sigma.velocity[2] = state.velocity[2] + state_covariance[col + 5];
        sigma.acceleration[0] = state.acceleration[0] +
            state_covariance[col + 6];
        sigma.acceleration[1] = state.acceleration[1] +
            state_covariance[col + 7];
        sigma.acceleration[2] = state.acceleration[2] +
            state_covariance[col + 8];

        _mul_quat_quat(sigma.attitude, noise_q, state.attitude);

        sigma.angular_velocity[0] = state.angular_velocity[0] +
            state_covariance[col + 12];
        sigma.angular_velocity[1] = state.angular_velocity[1] +
            state_covariance[col + 13];
        sigma.angular_velocity[2] = state.angular_velocity[2] +
            state_covariance[col + 14];
        sigma.angular_acceleration[0] = state.angular_acceleration[0] +
            state_covariance[col + 15];
        sigma.angular_acceleration[1] = state.angular_acceleration[1] +
            state_covariance[col + 16];
        sigma.angular_acceleration[2] = state.angular_acceleration[2] +
            state_covariance[col + 17];
        sigma.wind_velocity[0] = state.wind_velocity[0] +
            state_covariance[col + 18];
        sigma.wind_velocity[1] = state.wind_velocity[1] +
            state_covariance[col + 19];
        sigma.wind_velocity[2] = state.wind_velocity[2] +
            state_covariance[col + 20];
        sigma.gyro_bias[0] = state.gyro_bias[0] + state_covariance[col + 21];
        sigma.gyro_bias[1] = state.gyro_bias[1] + state_covariance[col + 22];
        sigma.gyro_bias[2] = state.gyro_bias[2] + state_covariance[col + 23];

        _ukf_process_sigma(&sigma, i + 1, dt, control, &apriori_mean,
            w_prime, measurement_estimate_mean,
            &measurement_estimate_sigma[(i + 1)*UKF_MEASUREMENT_DIM],
            &central_sigma);

        /* Create negative sigma point. */
        sigma.position[0] = state.position[0] - state_covariance[col + 0];
        sigma.position[1] = state.position[1] - state_covariance[col + 1];
        sigma.position[2] = state.position[2] - state_covariance[col + 2];
        sigma.velocity[0] = state.velocity[0] - state_covariance[col + 3];
        sigma.velocity[1] = state.velocity[1] - state_covariance[col + 4];
        sigma.velocity[2] = state.velocity[2] - state_covariance[col + 5];
        sigma.acceleration[0] = state.acceleration[0] -
            state_covariance[col + 6];
        sigma.acceleration[1] = state.acceleration[1] -
            state_covariance[col + 7];
        sigma.acceleration[2] = state.acceleration[2] -
            state_covariance[col + 8];

        noise_q[X] *= -1.0;
        noise_q[Y] *= -1.0;
        noise_q[Z] *= -1.0;
        _mul_quat_quat(sigma.attitude, noise_q, state.attitude);

        sigma.angular_velocity[0] = state.angular_velocity[0] -
            state_covariance[col + 12];
        sigma.angular_velocity[1] = state.angular_velocity[1] -
            state_covariance[col + 13];
        sigma.angular_velocity[2] = state.angular_velocity[2] -
            state_covariance[col + 14];
        sigma.angular_acceleration[0] = state.angular_acceleration[0] -
            state_covariance[col + 15];
        sigma.angular_acceleration[1] = state.angular_acceleration[1] -
            state_covariance[col + 16];
        sigma.angular_acceleration[2] = state.angular_acceleration[2] -
            state_covariance[col + 17];
        sigma.wind_velocity[0] = state.wind_velocity[0] -
            state_covariance[col + 18];
        sigma.wind_velocity[1] = state.wind_velocity[1] -
            state_covariance[col + 19];
        sigma.wind_velocity[2] = state.wind_velocity[2] -
            state_covariance[col + 20];
        sigma.gyro_bias[0] = state.gyro_bias[0] - state_covariance[col + 21];
        sigma.gyro_bias[1] = state.gyro_bias[1] - state_covariance[col + 22];
        sigma.gyro_bias[2] = state.gyro_bias[2] - state_covariance[col + 23];

        _ukf_process_sigma(&sigma, i + 1 + UKF_STATE_DIM, dt, control,
            &apriori_mean, w_prime, measurement_estimate_mean,
            &measurement_estimate_sigma[
                (i + 1 + UKF_STATE_DIM)*UKF_MEASUREMENT_DIM],
            &central_sigma);
    }

    /* Add central sigma points to the means */
    for (i = 0; i < UKF_MEASUREMENT_DIM; i++) {
        measurement_estimate_mean[i] =
            measurement_estimate_mean[i] * UKF_SIGMA_WMI
            + measurement_estimate_central[i] * UKF_SIGMA_WM0;
    }

    apriori_central.position[0] *= UKF_SIGMA_WM0;
    apriori_central.position[1] *= UKF_SIGMA_WM0;
    apriori_central.position[2] *= UKF_SIGMA_WM0;
    apriori_central.velocity[0] *= UKF_SIGMA_WM0;
    apriori_central.velocity[1] *= UKF_SIGMA_WM0;
    apriori_central.velocity[2] *= UKF_SIGMA_WM0;
    apriori_central.acceleration[0] *= UKF_SIGMA_WM0;
    apriori_central.acceleration[1] *= UKF_SIGMA_WM0;
    apriori_central.acceleration[2] *= UKF_SIGMA_WM0;
    apriori_central.attitude[0] *= UKF_SIGMA_WM0;
    apriori_central.attitude[1] *= UKF_SIGMA_WM0;
    apriori_central.attitude[2] *= UKF_SIGMA_WM0;
    apriori_central.attitude[3] *= UKF_SIGMA_WM0;
    apriori_central.angular_velocity[0] *= UKF_SIGMA_WM0;
    apriori_central.angular_velocity[1] *= UKF_SIGMA_WM0;
    apriori_central.angular_velocity[2] *= UKF_SIGMA_WM0;
    apriori_central.angular_acceleration[0] *= UKF_SIGMA_WM0;
    apriori_central.angular_acceleration[1] *= UKF_SIGMA_WM0;
    apriori_central.angular_acceleration[2] *= UKF_SIGMA_WM0;
    apriori_central.wind_velocity[0] *= UKF_SIGMA_WM0;
    apriori_central.wind_velocity[1] *= UKF_SIGMA_WM0;
    apriori_central.wind_velocity[2] *= UKF_SIGMA_WM0;
    apriori_central.gyro_bias[0] *= UKF_SIGMA_WM0;
    apriori_central.gyro_bias[1] *= UKF_SIGMA_WM0;
    apriori_central.gyro_bias[2] *= UKF_SIGMA_WM0;

    _mul_state_scalar_add_state(&apriori_mean, &apriori_mean, UKF_SIGMA_WMI,
                                &apriori_central);

/*
    printf("Apriori mean:\n");
    _print_matrix(apriori_mean.position, UKF_STATE_DIM, 1);
*/

    /*
    The initial mean estimate includes acceleration due to gravity separately,
    and both it and field strength need to be normalised to the expected
    field strengths.

    src/sensors.cpp line 265
    */
    real_t z_prime_col[UKF_MEASUREMENT_DIM];
    size_t measurement_dim = 0;
    j = 0;
    if (sensor_model.flags.accelerometer) {
        #define G(x) measurement_estimate_mean[3 + x]
        real_t g_norm = G_ACCEL /
            sqrt(G(X)*G(X) + G(Y)*G(Y) + G(Z)*G(Z));
        measurement_estimate_mean[X] += G(X) * g_norm;
        measurement_estimate_mean[Y] += G(Y) * g_norm;
        measurement_estimate_mean[Z] += G(Z) * g_norm;
        #undef G

        measurement_dim += 3;
        j += 6;
    }
    if (sensor_model.flags.gyroscope) {
        measurement_dim += 3;
        j += 3;
    }
    if (sensor_model.flags.magnetometer) {
        #define B(x) measurement_estimate_mean[j + x]
        real_t mag_norm = sensor_model_mag_field_norm /
            sqrt(B(X)*B(X) + B(Y)*B(Y) + B(Z)*B(Z));
        B(X) *= mag_norm;
        B(Y) *= mag_norm;
        B(Z) *= mag_norm;
        #undef B

        measurement_dim += 3;
        j += 3;
    }
    if (sensor_model.flags.gps_position) {
        measurement_dim += 3;
        j += 3;
    }
    if (sensor_model.flags.gps_velocity) {
        measurement_dim += 3;
        j += 3;
    }
    if (sensor_model.flags.pitot_tas) {
        measurement_dim += 1;
        j += 1;
    }
    if (sensor_model.flags.barometer_amsl) {
        measurement_dim += 1;
        j += 1;
    }

    if (j != measurement_dim) {
        /*
        Remove separate accelerometer reading due to gravity values; these
        will always be measurement_estimate_mean[3:6]
        */
        measurement_estimate_mean[3] = measurement_estimate_mean[6];
        measurement_estimate_mean[4] = measurement_estimate_mean[7];
        measurement_estimate_mean[5] = measurement_estimate_mean[8];
        measurement_estimate_mean[6] = measurement_estimate_mean[9];
        measurement_estimate_mean[7] = measurement_estimate_mean[10];
        measurement_estimate_mean[8] = measurement_estimate_mean[11];
        measurement_estimate_mean[9] = measurement_estimate_mean[12];
        measurement_estimate_mean[10] = measurement_estimate_mean[13];
        measurement_estimate_mean[11] = measurement_estimate_mean[14];
        measurement_estimate_mean[12] = measurement_estimate_mean[15];
        measurement_estimate_mean[13] = measurement_estimate_mean[16];
        measurement_estimate_mean[14] = measurement_estimate_mean[17];
        measurement_estimate_mean[15] = measurement_estimate_mean[18];
        measurement_estimate_mean[16] = measurement_estimate_mean[19];
    }

/*
    printf("Measurement estimate mean:\n");
    _print_matrix(measurement_estimate_mean, measurement_dim, 1);
*/

    /*
    Calculate z_prime columns for each sigma point in sequence;
    src/ukf.cpp line 265

    These are only used to calculate measurement_estimate_covariance and
    cross_correlation (src/ukf.cpp line 298).

    In each case, sum together all sigma points before adding the central
    point, since the weighting of the central point will cause loss of
    precision otherwise.
    */
    _ukf_sensor_calculate_deltas(z_prime_col, measurement_estimate_mean, 1);

    _mul_mat(cross_correlation, &w_prime[1 * UKF_STATE_DIM], z_prime_col, 1,
             measurement_dim, UKF_STATE_DIM, 1, UKF_SIGMA_WCI);
    _mul_mat(measurement_estimate_covariance, z_prime_col, z_prime_col, 1,
             measurement_dim, measurement_dim, 1, UKF_SIGMA_WCI);
    _mul_mat(state_covariance, &w_prime[1 * UKF_STATE_DIM],
             &w_prime[1 * UKF_STATE_DIM], 1, UKF_STATE_DIM,
             UKF_STATE_DIM, 1, UKF_SIGMA_WCI);

    for (i = 2; i < UKF_NUM_SIGMA; i++) {
        _ukf_sensor_calculate_deltas(z_prime_col, measurement_estimate_mean,
                                     i);

        _mul_mat_accum(cross_correlation, &w_prime[i * UKF_STATE_DIM],
                       z_prime_col, 1, measurement_dim, UKF_STATE_DIM, 1,
                       UKF_SIGMA_WCI);
        _mul_mat_accum(measurement_estimate_covariance, z_prime_col,
                       z_prime_col, 1, measurement_dim, measurement_dim, 1,
                       UKF_SIGMA_WCI);
        _mul_mat_accum(state_covariance, &w_prime[i * UKF_STATE_DIM],
                       &w_prime[i * UKF_STATE_DIM], 1, UKF_STATE_DIM,
                       UKF_STATE_DIM, 1, UKF_SIGMA_WCI);
    }

    _ukf_sensor_calculate_deltas(z_prime_col, measurement_estimate_mean, 0);
    _mul_mat_accum(cross_correlation, w_prime, z_prime_col, 1,
                   measurement_dim, UKF_STATE_DIM, 1, UKF_SIGMA_WC0);
    _mul_mat_accum(measurement_estimate_covariance, z_prime_col, z_prime_col,
                   1, measurement_dim, measurement_dim, 1, UKF_SIGMA_WC0);
    _mul_mat_accum(state_covariance, w_prime, w_prime, 1, UKF_STATE_DIM,
                   UKF_STATE_DIM, 1, UKF_SIGMA_WC0);

    /* Easy case if no measurements */
    if (measurement_dim == 0) {
        memcpy(&state, &central_sigma, sizeof(state));
        _normalize_quat(state.attitude, state.attitude, true);
        return;
    }

    /*
    Done with w_prime and measurement_sigma_points, so we can use them for
    other variables as necessary
    */

    assert(measurement_dim <= UKF_MEASUREMENT_DIM - 3);
    /*
    MEASUREMENT_DIM - 3 because we always remove the extra 3 accelerometer
    readings for the gravity component above

    src/ukf.cpp line 282
    */
    real_t innovation[UKF_MEASUREMENT_DIM - 3],
           sensor_covariance[UKF_MEASUREMENT_DIM - 3];
    _ukf_sensor_collate(innovation);
    _ukf_sensor_get_covariance(sensor_covariance);

    for (i = 0; i < measurement_dim; i++) {
        innovation[i] -= measurement_estimate_mean[i];
        measurement_estimate_covariance[i*measurement_dim + i] +=
            sensor_covariance[i];
    }

/*
    printf("Apriori covariance:\n");
    _print_matrix(state_covariance, UKF_STATE_DIM, UKF_STATE_DIM);

    printf("Measurement estimate covariance:\n");
    _print_matrix(measurement_estimate_covariance, measurement_dim,
                  measurement_dim);

    printf("Cross correlation:\n");
    _print_matrix(cross_correlation, measurement_dim, UKF_STATE_DIM);
*/

    /*
    Use w_prime (9408 bytes) and measurement_estimate_sigma (7840 bytes)
    for temporary variable storage.
    update_temp is 17 bytes, measurement_estimate_covariance_i is 2312 bytes,
    kalman_gain is 3264 bytes (5593 bytes used from w_prime);
    kalman_gain_t is 3264 bytes, state_temp1 is 3264 bytes (6528 bytes used
    from measurement_estimate_sigma).

    src/ukf.cpp line 306, 313
    */
    real_t *update_temp = w_prime,
           *measurement_estimate_covariance_i = &w_prime[UKF_MEASUREMENT_DIM - 3],
           *kalman_gain = &measurement_estimate_covariance_i[
                (UKF_MEASUREMENT_DIM - 3)*(UKF_MEASUREMENT_DIM - 3)],
           *kalman_gain_t = measurement_estimate_sigma,
           *state_temp1 = &measurement_estimate_sigma[
                UKF_STATE_DIM*(UKF_MEASUREMENT_DIM - 3)];

    /*
    use kalman_gain as temporary working space required by the invert function
    */
    _inv_mat(measurement_estimate_covariance_i,
        measurement_estimate_covariance, measurement_dim,
        kalman_gain);

/*
    printf("Measurement estimate covariance inverse:\n");
    _print_matrix(measurement_estimate_covariance_i, measurement_dim,
                  measurement_dim);

    printf("Innovation:\n");
    _print_matrix(innovation, measurement_dim, 1);
*/

    _mul_mat(kalman_gain, cross_correlation,
        measurement_estimate_covariance_i, measurement_dim,
        measurement_dim, UKF_STATE_DIM, measurement_dim, 1.0);

/*
    printf("Kalman gain:\n");
    _print_matrix(kalman_gain, measurement_dim, UKF_STATE_DIM);
    printf("Measurement estimate covariance:\n");
    _print_matrix(measurement_estimate_covariance, measurement_dim,
                  measurement_dim);
*/

    _mul_mat(update_temp, kalman_gain, innovation, measurement_dim,
        1, UKF_STATE_DIM, measurement_dim, 1.0);

/*
    printf("Update temp:\n");
    _print_matrix(update_temp, UKF_STATE_DIM, 1);

    printf("Apriori mean:\n");
    _print_matrix(apriori_mean.position, UKF_STATE_DIM, 1);
*/

    /* Update the state */
    state.position[0] = apriori_mean.position[0] + update_temp[0];
    state.position[1] = apriori_mean.position[1] + update_temp[1];
    state.position[2] = apriori_mean.position[2] + update_temp[2];
    state.velocity[X] = apriori_mean.velocity[X] + update_temp[3];
    state.velocity[Y] = apriori_mean.velocity[Y] + update_temp[4];
    state.velocity[Z] = apriori_mean.velocity[Z] + update_temp[5];
    state.acceleration[X] = apriori_mean.acceleration[X] + update_temp[6];
    state.acceleration[Y] = apriori_mean.acceleration[Y] + update_temp[7];
    state.acceleration[Z] = apriori_mean.acceleration[Z] + update_temp[8];

    real_t d_p[3] = { update_temp[9], update_temp[10], update_temp[11] };
    real_t x_2 = d_p[X]*d_p[X] + d_p[Y]*d_p[Y] + d_p[Z]*d_p[Z];
    real_t d_q_w = (-UKF_MRP_A * x_2 + UKF_MRP_F *
                sqrt(UKF_MRP_F_2 + ((real_t)1.0 - UKF_MRP_A_2) * x_2)) /
                (UKF_MRP_F_2 + x_2);
    real_t d_q_scale = (1.0 / UKF_MRP_F) * (UKF_MRP_A + d_q_w);
    real_t d_q[4] = {
        d_q_scale * d_p[X], d_q_scale * d_p[Y], d_q_scale * d_p[Z], d_q_w
    };

    _mul_quat_quat(state.attitude, d_q, central_sigma.attitude);

/*
    printf("Attitude:\n");
    _print_matrix(state.attitude, 4, 1);
*/

    _normalize_quat(state.attitude, state.attitude, true);

    state.angular_velocity[X] = apriori_mean.angular_velocity[X] +
        update_temp[12];
    state.angular_velocity[Y] = apriori_mean.angular_velocity[Y] +
        update_temp[13];
    state.angular_velocity[Z] = apriori_mean.angular_velocity[Z] +
        update_temp[14];
    state.angular_acceleration[X] = apriori_mean.angular_acceleration[X] +
        update_temp[15];
    state.angular_acceleration[Y] = apriori_mean.angular_acceleration[Y] +
        update_temp[16];
    state.angular_acceleration[Z] = apriori_mean.angular_acceleration[Z] +
        update_temp[17];
    state.wind_velocity[X] = apriori_mean.wind_velocity[X] + update_temp[18];
    state.wind_velocity[Y] = apriori_mean.wind_velocity[Y] + update_temp[19];
    state.wind_velocity[Z] = apriori_mean.wind_velocity[Z] + update_temp[20];
    state.gyro_bias[X] = apriori_mean.gyro_bias[X] + update_temp[21];
    state.gyro_bias[Y] = apriori_mean.gyro_bias[Y] + update_temp[22];
    state.gyro_bias[Z] = apriori_mean.gyro_bias[Z] + update_temp[23];

    /* Update the state covariance; src/ukf.cpp line 352 */
    _mul_mat(state_temp1, kalman_gain, measurement_estimate_covariance,
        measurement_dim, measurement_dim, UKF_STATE_DIM, measurement_dim,
        1.0);
/*
    printf("Update matrix:\n");
    _print_matrix(state_temp1, measurement_dim, UKF_STATE_DIM);
*/

    _transpose_mat(kalman_gain_t, kalman_gain, measurement_dim,
        UKF_STATE_DIM);

    /*
    All we need now is kalman_gain_t and state_temp1, both of which are stored
    in measurement_estimate_sigma's space. Thus, we can re-use the w_prime
    allocation
    */
    real_t *state_temp2 = w_prime;
    _mul_mat(state_temp2, state_temp1, kalman_gain_t, measurement_dim,
             UKF_STATE_DIM, UKF_STATE_DIM, measurement_dim, -1.0);
    _add_mat_accum(state_covariance, state_temp2,
                   UKF_STATE_DIM * UKF_STATE_DIM);

/*
    printf("State covariance:\n");
    _print_matrix(state_covariance, UKF_STATE_DIM, UKF_STATE_DIM);

    printf("Delta cycles: %llu\n", rdtsc() - t);

    printf("State:\n");
    _print_matrix(state.position, UKF_STATE_DIM, 1);
*/
}

void ukf_init(void) {
    real_t state_covariance_diag[UKF_STATE_DIM] = {
        (real_t)M_PI * (real_t)M_PI * (real_t)0.0625,
            (real_t)M_PI * (real_t)M_PI * (real_t)0.0625, 1000,
        50, 50, 50,
        10, 10, 10,
        (real_t)M_PI * (real_t)0.25, (real_t)M_PI * (real_t)0.25,
            (real_t)M_PI * (real_t)0.25,
        2, 2, 2,
        5, 5, 5,
        20, 20, 20,
        0, 0, 0
    };

    memset(state_covariance, 0, sizeof(state_covariance));

    uint32_t i;
    for (i = 0; i < UKF_STATE_DIM; i++) {
        process_noise[i] = (real_t)1e-6;
        state_covariance[i + i*UKF_STATE_DIM] = state_covariance_diag[i];
    }

    struct ukf_state_t base_state = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };
    ukf_set_state(&base_state);

    struct ukf_ioboard_params_t base_config = {
        /* sensor offsets/orientations */
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0, 1},
        {1, 0, 0},
        /* sensor covariance */
        {1, 1, 1},
        {1, 1, 1},
        {1, 1, 1},
        {1, 1, 1},
        {1, 1, 1},
        1,
        1
    };
    ukf_set_params(&base_config);
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
