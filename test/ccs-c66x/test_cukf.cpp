#include <gtest/gtest.h>
#include <math.h>
#include <stdint.h>
#include "cukf.h"
#include "cukfmath.h"

#include <xmmintrin.h>

/* Prototypes for internal functions */
extern "C" {
void _ukf_state_integrate_rk4(struct ukf_state_t *in, real_t delta);
void _ukf_state_model(struct ukf_state_t *in);
}

TEST(C66xMathTest, QuaternionVectorMultiplication) {
    real_t q[4], v[3], result[3];

    /* Quaternion identity */
    q[X] = 0;
    q[Y] = 0;
    q[Z] = 0;
    q[W] = 1;

    v[X] = 2;
    v[Y] = 3;
    v[Z] = 4;

    _mul_quat_vec3(result, q, v);

    EXPECT_FLOAT_EQ(2, result[X]);
    EXPECT_FLOAT_EQ(3, result[Y]);
    EXPECT_FLOAT_EQ(4, result[Z]);

    /* Conjugate of identity */
    q[W] = -1;

    _mul_quat_vec3(result, q, v);

    EXPECT_FLOAT_EQ(2, result[X]);
    EXPECT_FLOAT_EQ(3, result[Y]);
    EXPECT_FLOAT_EQ(4, result[Z]);

    /* Rotation about X */
    /* TODO */
}

TEST(C66xMathTest, VectorAxPlusB) {
    struct ukf_state_t a, b, c;

    a.position[X] = 1;
    a.position[Y] = 2;
    a.position[Z] = 3;
    a.velocity[X] = 4;
    a.velocity[Y] = 5;
    a.velocity[Z] = 6;
    a.acceleration[X] = 7;
    a.acceleration[Y] = 8;
    a.acceleration[Z] = 9;
    a.attitude[X] = 10;
    a.attitude[Y] = 11;
    a.attitude[Z] = 12;
    a.attitude[W] = 13;
    a.angular_velocity[X] = 14;
    a.angular_velocity[Y] = 15;
    a.angular_velocity[Z] = 16;
    a.angular_acceleration[X] = 17;
    a.angular_acceleration[Y] = 18;
    a.angular_acceleration[Z] = 19;
    a.wind_velocity[X] = 20;
    a.wind_velocity[Y] = 21;
    a.wind_velocity[Z] = 22;
    a.gyro_bias[X] = 23;
    a.gyro_bias[Y] = 24;
    a.gyro_bias[Z] = 25;

    b.position[X] = 26;
    b.position[Y] = 27;
    b.position[Z] = 28;
    b.velocity[X] = 29;
    b.velocity[Y] = 30;
    b.velocity[Z] = 31;
    b.acceleration[X] = 32;
    b.acceleration[Y] = 33;
    b.acceleration[Z] = 34;
    b.attitude[X] = 35;
    b.attitude[Y] = 36;
    b.attitude[Z] = 37;
    b.attitude[W] = 38;
    b.angular_velocity[X] = 39;
    b.angular_velocity[Y] = 40;
    b.angular_velocity[Z] = 41;
    b.angular_acceleration[X] = 42;
    b.angular_acceleration[Y] = 43;
    b.angular_acceleration[Z] = 44;
    b.wind_velocity[X] = 45;
    b.wind_velocity[Y] = 46;
    b.wind_velocity[Z] = 47;
    b.gyro_bias[X] = 48;
    b.gyro_bias[Y] = 49;
    b.gyro_bias[Z] = 50;

    memcpy(&c, &a, sizeof(c));

    _mul_vec_scalar_add_vec((real_t*)&c, -1.0, (real_t*)&a, UKF_STATE_DIM + 1);
    EXPECT_FLOAT_EQ(0, c.position[X]);
    EXPECT_FLOAT_EQ(0, c.position[Y]);
    EXPECT_FLOAT_EQ(0, c.position[Z]);
    EXPECT_FLOAT_EQ(0, c.velocity[X]);
    EXPECT_FLOAT_EQ(0, c.velocity[Y]);
    EXPECT_FLOAT_EQ(0, c.velocity[Z]);
    EXPECT_FLOAT_EQ(0, c.acceleration[X]);
    EXPECT_FLOAT_EQ(0, c.acceleration[Y]);
    EXPECT_FLOAT_EQ(0, c.acceleration[Z]);
    EXPECT_FLOAT_EQ(0, c.attitude[X]);
    EXPECT_FLOAT_EQ(0, c.attitude[Y]);
    EXPECT_FLOAT_EQ(0, c.attitude[Z]);
    EXPECT_FLOAT_EQ(0, c.attitude[W]);
    EXPECT_FLOAT_EQ(0, c.angular_velocity[X]);
    EXPECT_FLOAT_EQ(0, c.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(0, c.angular_velocity[Z]);
    EXPECT_FLOAT_EQ(0, c.angular_acceleration[X]);
    EXPECT_FLOAT_EQ(0, c.angular_acceleration[Y]);
    EXPECT_FLOAT_EQ(0, c.angular_acceleration[Z]);
    EXPECT_FLOAT_EQ(0, c.wind_velocity[X]);
    EXPECT_FLOAT_EQ(0, c.wind_velocity[Y]);
    EXPECT_FLOAT_EQ(0, c.wind_velocity[Z]);
    EXPECT_FLOAT_EQ(0, c.gyro_bias[X]);
    EXPECT_FLOAT_EQ(0, c.gyro_bias[Y]);
    EXPECT_FLOAT_EQ(0, c.gyro_bias[Z]);

    _mul_vec_scalar_add_vec((real_t*)&a, 0.5, (real_t*)&b, UKF_STATE_DIM + 1);
    EXPECT_FLOAT_EQ(26.5, a.position[X]);
    EXPECT_FLOAT_EQ(28.0, a.position[Y]);
    EXPECT_FLOAT_EQ(29.5, a.position[Z]);
    EXPECT_FLOAT_EQ(31.0, a.velocity[X]);
    EXPECT_FLOAT_EQ(32.5, a.velocity[Y]);
    EXPECT_FLOAT_EQ(34.0, a.velocity[Z]);
    EXPECT_FLOAT_EQ(35.5, a.acceleration[X]);
    EXPECT_FLOAT_EQ(37.0, a.acceleration[Y]);
    EXPECT_FLOAT_EQ(38.5, a.acceleration[Z]);
    EXPECT_FLOAT_EQ(40.0, a.attitude[X]);
    EXPECT_FLOAT_EQ(41.5, a.attitude[Y]);
    EXPECT_FLOAT_EQ(43.0, a.attitude[Z]);
    EXPECT_FLOAT_EQ(44.5, a.attitude[W]);
    EXPECT_FLOAT_EQ(46.0, a.angular_velocity[X]);
    EXPECT_FLOAT_EQ(47.5, a.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(49.0, a.angular_velocity[Z]);
    EXPECT_FLOAT_EQ(50.5, a.angular_acceleration[X]);
    EXPECT_FLOAT_EQ(52.0, a.angular_acceleration[Y]);
    EXPECT_FLOAT_EQ(53.5, a.angular_acceleration[Z]);
    EXPECT_FLOAT_EQ(55.0, a.wind_velocity[X]);
    EXPECT_FLOAT_EQ(56.5, a.wind_velocity[Y]);
    EXPECT_FLOAT_EQ(58.0, a.wind_velocity[Z]);
    EXPECT_FLOAT_EQ(59.5, a.gyro_bias[X]);
    EXPECT_FLOAT_EQ(61.0, a.gyro_bias[Y]);
    EXPECT_FLOAT_EQ(62.5, a.gyro_bias[Z]);
}

TEST(C66xMathTest, VectorAdd) {
    real_t v1[5] = { 1, 2, 3, 4, 5 }, v2[5] = { 6, 7, 8, 9, 10 },
           v3[5] = { -7, -9, -11, -13, -15 };

    _add_vec_vec(v1, v2, 5);
    EXPECT_FLOAT_EQ(7, v1[0]);
    EXPECT_FLOAT_EQ(9, v1[1]);
    EXPECT_FLOAT_EQ(11, v1[2]);
    EXPECT_FLOAT_EQ(13, v1[3]);
    EXPECT_FLOAT_EQ(15, v1[4]);

    _add_vec_vec(v3, v1, 5);
    EXPECT_FLOAT_EQ(0, v3[0]);
    EXPECT_FLOAT_EQ(0, v3[1]);
    EXPECT_FLOAT_EQ(0, v3[2]);
    EXPECT_FLOAT_EQ(0, v3[3]);
    EXPECT_FLOAT_EQ(0, v3[4]);
}

TEST(C66xMathTest, Matrix3x3Inverse) {
    real_t m1[9] = {
                1.1, 0.2, 0.3,
                0.4, 1.5, 0.6,
                0.7, 0.8, 1.9
           },
           result[9];

    _inv_mat3x3(result, m1);
    EXPECT_FLOAT_EQ(1.021551724138, result[0]);
    EXPECT_FLOAT_EQ(-0.060344827586, result[1]);
    EXPECT_FLOAT_EQ(-0.142241379310, result[2]);
    EXPECT_FLOAT_EQ(-0.146551724138, result[3]);
    EXPECT_FLOAT_EQ(0.810344827586, result[4]);
    EXPECT_FLOAT_EQ(-0.232758620690, result[5]);
    EXPECT_FLOAT_EQ(-0.314655172414, result[6]);
    EXPECT_FLOAT_EQ(-0.318965517241, result[7]);
    EXPECT_FLOAT_EQ(0.676724137931, result[8]);
}

/* State model tests */

TEST(C66xStateTest, KinematicsModelBasic) {
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {0, 1, 2},
        {3, 4, 5},
        {0, 0, 0, 1},
        {1, 0, 0},
        {6, 7, 8},
        {9, 10, 11},
        {12, 13, 14}
    };

    _ukf_state_model(&test_state);

    EXPECT_FLOAT_EQ(0, test_state.position[X]);
    EXPECT_FLOAT_EQ(1.5678559e-07, test_state.position[Y]);
    EXPECT_FLOAT_EQ(-2, test_state.position[Z]);

    EXPECT_FLOAT_EQ(3, test_state.velocity[X]);
    EXPECT_FLOAT_EQ(4, test_state.velocity[Y]);
    EXPECT_FLOAT_EQ(5, test_state.velocity[Z]);

    EXPECT_FLOAT_EQ(6, test_state.angular_velocity[X]);
    EXPECT_FLOAT_EQ(7, test_state.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(8, test_state.angular_velocity[Z]);

    EXPECT_FLOAT_EQ(-0.5, test_state.attitude[X]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[Y]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[Z]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[W]);
}

TEST(C66xStateTest, KinematicsModelRotated) {
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {0, 1, 2},
        {3, 4, 5},
        {0.70710678118655, 0, 0, 0.70710678118655},
        {1, 0, 0},
        {6, 7, 8},
        {9, 10, 11},
        {12, 13, 14}
    };

    _ukf_state_model(&test_state);

    EXPECT_FLOAT_EQ(0, test_state.position[X]);
    EXPECT_FLOAT_EQ(1.5678559e-07, test_state.position[Y]);
    EXPECT_FLOAT_EQ(-2, test_state.position[Z]);

    EXPECT_FLOAT_EQ(3, test_state.velocity[X]);
    EXPECT_FLOAT_EQ(5, test_state.velocity[Y]);
    EXPECT_FLOAT_EQ(-4, test_state.velocity[Z]);

    EXPECT_FLOAT_EQ(6, test_state.angular_velocity[X]);
    EXPECT_FLOAT_EQ(7, test_state.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(8, test_state.angular_velocity[Z]);

    EXPECT_FLOAT_EQ(-0.35355339059328, test_state.attitude[X]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[Y]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[Z]);
    EXPECT_FLOAT_EQ(0.35355339059328, test_state.attitude[W]);
}

/* Integrator tests */

TEST(C66xIntegratorTest, ConstantAcceleration) {
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {3, 4, 5},
        {6, 7, 8},
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    for (float i = 0; i < 10; i += 0.01) {
        _ukf_state_integrate_rk4(&test_state, 0.01);
    }

    EXPECT_FLOAT_EQ(5.20898e-05, test_state.position[X]);
    EXPECT_FLOAT_EQ(6.1148545e-05, test_state.position[Y]);
    EXPECT_FLOAT_EQ(-450, test_state.position[Z]);
    EXPECT_FLOAT_EQ(63, test_state.velocity[X]);
    EXPECT_FLOAT_EQ(74, test_state.velocity[Y]);
    EXPECT_FLOAT_EQ(85, test_state.velocity[Z]);
}

TEST(C66xIntegratorTest, ConstantAngularAcceleration) {
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 1},
        {0, 0, 0},
        {0, 0, 0}
    };

    for(float i = 0; i < 10; i += 0.01) {
        _ukf_state_integrate_rk4(&test_state, 0.01);
    }

    EXPECT_FLOAT_EQ(0, test_state.attitude[X]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[Y]);
    EXPECT_FLOAT_EQ(0.13235217, test_state.attitude[Z]);
    EXPECT_FLOAT_EQ(0.99120271, test_state.attitude[W]);
    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(10, test_state.angular_velocity[Z]);
}

TEST(C66xIntegratorTest, CircularMotion) {
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {10, 0, 0},
        {0, 6.28318530717959, 0},
        {0, 0, 0, 1},
        {0, 0, 0.628318530717959},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    for(float i = 0; i < 10; i += 0.01) {
        _ukf_state_integrate_rk4(&test_state, 0.01);
    }

    EXPECT_NEAR(0, test_state.position[X], 1e-6);
    EXPECT_NEAR(0, test_state.position[Y], 1e-6);
    EXPECT_NEAR(0, test_state.position[Z], 1e-6);
    EXPECT_NEAR(10, test_state.velocity[X], 1e-6);
    EXPECT_NEAR(0, test_state.velocity[Y], 1e-6);
    EXPECT_NEAR(0, test_state.velocity[Z], 1e-6);
    EXPECT_FLOAT_EQ(0, test_state.acceleration[X]);
    EXPECT_FLOAT_EQ(6.28318530717959, test_state.acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.acceleration[Z]);

    EXPECT_NEAR(0, test_state.attitude[X], 1e-6);
    EXPECT_NEAR(0, test_state.attitude[Y], 1e-6);
    EXPECT_NEAR(0, test_state.attitude[Z], 1e-6);
    EXPECT_NEAR(1, fabs(test_state.attitude[W]), 1e-6);
    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(0.628318530717959, test_state.angular_velocity[Z]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Z]);
}
