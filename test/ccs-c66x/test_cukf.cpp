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

TEST(C66xMathTest, StateAxPlusB) {
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

    _mul_state_scalar_add_state(&c, &c, -1.0, &a);
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

    _mul_state_scalar_add_state(&a, &a, 0.5, &b);
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

TEST(C66xMathTest, MatrixCholeskyLLT) {
    real_t m1[16] = {
        18, 22,  54,  42,
        22, 70,  86,  62,
        54, 86, 174, 134,
        42, 62, 134, 106
    };
    /* lower triangle, column major */
    real_t expected[16] = {
        4.24264, 5.18545, 12.72792, 9.8994951,
        0, 6.5659051, 3.0460384, 1.6245539,
        0, 0, 1.6497422, 1.8497111,
        0, 0, 0, 1.3926213
    };

    _cholesky_mat_mul(m1, m1, 1.0, 4);
    EXPECT_FLOAT_EQ(expected[0], m1[0]);
    EXPECT_FLOAT_EQ(expected[1], m1[1]);
    EXPECT_FLOAT_EQ(expected[2], m1[2]);
    EXPECT_FLOAT_EQ(expected[3], m1[3]);
    EXPECT_FLOAT_EQ(expected[5], m1[5]);
    EXPECT_FLOAT_EQ(expected[6], m1[6]);
    EXPECT_FLOAT_EQ(expected[7], m1[7]);
    EXPECT_FLOAT_EQ(expected[10], m1[10]);
    EXPECT_FLOAT_EQ(expected[11], m1[11]);
    EXPECT_FLOAT_EQ(expected[15], m1[15]);
}

TEST(C66xMathTest, MatrixMultiply) {
    real_t a[16] = {
        1, 1, 1, 1,
        2, 4, 8, 16,
        3, 9, 27, 81,
        4, 16, 64, 256
    };
    real_t b[16] = {
        4.0, -3.0, 4.0/3.0, -1.0/4.0,
        -13.0/3.0, 19.0/4.0, -7.0/3.0, 11.0/24.0,
        3.0/2.0, -2.0, 7.0/6.0, -1.0/4.0,
        -1.0/6.0, 1.0/4.0, -1.0/6.0, 1.0/24.0
    };
    real_t c[16];

    _mul_mat(c, a, b, 4, 4, 4, 4, 1.0);
    EXPECT_NEAR(1.0, c[0], 1e-12);
    EXPECT_NEAR(0.0, c[1], 1e-12);
    EXPECT_NEAR(0.0, c[2], 1e-12);
    EXPECT_NEAR(0.0, c[3], 1e-12);
    EXPECT_NEAR(0.0, c[4], 1e-12);
    EXPECT_NEAR(1.0, c[5], 1e-12);
    EXPECT_NEAR(0.0, c[6], 1e-12);
    EXPECT_NEAR(0.0, c[7], 1e-12);
    EXPECT_NEAR(0.0, c[8], 1e-12);
    EXPECT_NEAR(0.0, c[9], 1e-12);
    EXPECT_NEAR(1.0, c[10], 1e-12);
    EXPECT_NEAR(0.0, c[11], 1e-12);
    EXPECT_NEAR(0.0, c[12], 1e-12);
    EXPECT_NEAR(0.0, c[13], 1e-12);
    EXPECT_NEAR(0.0, c[14], 1e-12);
    EXPECT_NEAR(1.0, c[15], 1e-12);

    real_t d[4] = { 1, 2, 3, 4}, e[4] = { 5, 6, 7, 8 };

    _mul_mat(c, d, e, 1, 4, 4, 1, 1.0);
    EXPECT_NEAR(5.0, c[0], 1e-12);
    EXPECT_NEAR(10.0, c[1], 1e-12);
    EXPECT_NEAR(15.0, c[2], 1e-12);
    EXPECT_NEAR(20.0, c[3], 1e-12);
    EXPECT_NEAR(6.0, c[4], 1e-12);
    EXPECT_NEAR(12.0, c[5], 1e-12);
    EXPECT_NEAR(18.0, c[6], 1e-12);
    EXPECT_NEAR(24.0, c[7], 1e-12);
    EXPECT_NEAR(7.0, c[8], 1e-12);
    EXPECT_NEAR(14.0, c[9], 1e-12);
    EXPECT_NEAR(21.0, c[10], 1e-12);
    EXPECT_NEAR(28.0, c[11], 1e-12);
    EXPECT_NEAR(8.0, c[12], 1e-12);
    EXPECT_NEAR(16.0, c[13], 1e-12);
    EXPECT_NEAR(24.0, c[14], 1e-12);
    EXPECT_NEAR(32.0, c[15], 1e-12);
}

TEST(C66xMathTest, MatrixMultiplyAccum) {
    real_t c[16], d[4] = { 1, 2, 3, 4}, e[4] = { 5, 6, 7, 8 };

    memset(c, 0, sizeof(c));
    c[0] = c[5] = c[10] = c[15] = 1.0;

    _mul_mat_accum(c, d, e, 1, 4, 4, 1, 1.0);
    EXPECT_NEAR(6.0, c[0], 1e-12);
    EXPECT_NEAR(10.0, c[1], 1e-12);
    EXPECT_NEAR(15.0, c[2], 1e-12);
    EXPECT_NEAR(20.0, c[3], 1e-12);
    EXPECT_NEAR(6.0, c[4], 1e-12);
    EXPECT_NEAR(13.0, c[5], 1e-12);
    EXPECT_NEAR(18.0, c[6], 1e-12);
    EXPECT_NEAR(24.0, c[7], 1e-12);
    EXPECT_NEAR(7.0, c[8], 1e-12);
    EXPECT_NEAR(14.0, c[9], 1e-12);
    EXPECT_NEAR(22.0, c[10], 1e-12);
    EXPECT_NEAR(28.0, c[11], 1e-12);
    EXPECT_NEAR(8.0, c[12], 1e-12);
    EXPECT_NEAR(16.0, c[13], 1e-12);
    EXPECT_NEAR(24.0, c[14], 1e-12);
    EXPECT_NEAR(33.0, c[15], 1e-12);
}

TEST(C66xMathTest, VectorCross) {
    real_t v1[3] = { 1, 2, 3 }, v2[3] = { 6, 5, 4}, result[3];

    _cross_vec3(result, v1, v2);
    EXPECT_FLOAT_EQ(-7, result[X]);
    EXPECT_FLOAT_EQ(14, result[Y]);
    EXPECT_FLOAT_EQ(-7, result[Z]);

    real_t v3[3] = { 0.1, 1.0, 0.5 }, v4[3] = { -5.0, 0.1, 3.0 };
    _cross_vec3(result, v3, v4);
    EXPECT_FLOAT_EQ(2.95, result[X]);
    EXPECT_FLOAT_EQ(-2.8, result[Y]);
    EXPECT_FLOAT_EQ(5.01, result[Z]);
}

TEST(C66xMathTest, QuaternionQuaternionMultiply) {
    real_t q1[4] = { 1, -2, 5, 2 }, q2[4] = { 2, 3, 4, 1 }, result[4];

    _mul_quat_quat(result, q1, q2);
    EXPECT_FLOAT_EQ(-18, result[X]);
    EXPECT_FLOAT_EQ(10, result[Y]);
    EXPECT_FLOAT_EQ(20, result[Z]);
    EXPECT_FLOAT_EQ(-14, result[W]);
}

TEST(C66xMathTest, QuaternionNormalize) {
    real_t q1[4] = { -18, 10, 20, -14 }, result[4];

    _normalize_quat(result, q1, true);
    EXPECT_FLOAT_EQ(0.56360185, result[X]);
    EXPECT_FLOAT_EQ(-0.31311214, result[Y]);
    EXPECT_FLOAT_EQ(-0.62622428, result[Z]);
    EXPECT_FLOAT_EQ(0.438357, result[W]);

    _normalize_quat(q1, q1, false);
    EXPECT_FLOAT_EQ(-0.56360185, q1[X]);
    EXPECT_FLOAT_EQ(0.31311214, q1[Y]);
    EXPECT_FLOAT_EQ(0.62622428, q1[Z]);
    EXPECT_FLOAT_EQ(-0.438357, q1[W]);
}

TEST(C66xMathTest, MatrixTranspose) {
    real_t M1[25] = {
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15,
        16, 17, 18, 19, 20,
        21, 22, 23, 24, 25
    }, result[25];

    _transpose_mat(result, M1, 5, 5);
    EXPECT_FLOAT_EQ(1, result[0]);
    EXPECT_FLOAT_EQ(6, result[1]);
    EXPECT_FLOAT_EQ(11, result[2]);
    EXPECT_FLOAT_EQ(16, result[3]);
    EXPECT_FLOAT_EQ(21, result[4]);
    EXPECT_FLOAT_EQ(2, result[5]);
    EXPECT_FLOAT_EQ(7, result[6]);
    EXPECT_FLOAT_EQ(12, result[7]);
    EXPECT_FLOAT_EQ(17, result[8]);
    EXPECT_FLOAT_EQ(22, result[9]);
    EXPECT_FLOAT_EQ(3, result[10]);
    EXPECT_FLOAT_EQ(8, result[11]);
    EXPECT_FLOAT_EQ(13, result[12]);
    EXPECT_FLOAT_EQ(18, result[13]);
    EXPECT_FLOAT_EQ(23, result[14]);
    EXPECT_FLOAT_EQ(4, result[15]);
    EXPECT_FLOAT_EQ(9, result[16]);
    EXPECT_FLOAT_EQ(14, result[17]);
    EXPECT_FLOAT_EQ(19, result[18]);
    EXPECT_FLOAT_EQ(24, result[19]);
    EXPECT_FLOAT_EQ(5, result[20]);
    EXPECT_FLOAT_EQ(10, result[21]);
    EXPECT_FLOAT_EQ(15, result[22]);
    EXPECT_FLOAT_EQ(20, result[23]);
    EXPECT_FLOAT_EQ(25, result[24]);

    real_t M2[8] = {
        1, 2, 3, 4,
        5, 6, 7, 8
    };

    _transpose_mat(result, M2, 2, 4);
    EXPECT_FLOAT_EQ(1, result[0]);
    EXPECT_FLOAT_EQ(5, result[1]);
    EXPECT_FLOAT_EQ(2, result[2]);
    EXPECT_FLOAT_EQ(6, result[3]);
    EXPECT_FLOAT_EQ(3, result[4]);
    EXPECT_FLOAT_EQ(7, result[5]);
    EXPECT_FLOAT_EQ(4, result[6]);
    EXPECT_FLOAT_EQ(8, result[7]);

    _transpose_mat(result, M2, 4, 2);
    EXPECT_FLOAT_EQ(1, result[0]);
    EXPECT_FLOAT_EQ(3, result[1]);
    EXPECT_FLOAT_EQ(5, result[2]);
    EXPECT_FLOAT_EQ(7, result[3]);
    EXPECT_FLOAT_EQ(2, result[4]);
    EXPECT_FLOAT_EQ(4, result[5]);
    EXPECT_FLOAT_EQ(6, result[6]);
    EXPECT_FLOAT_EQ(8, result[7]);
}

TEST(C66xMathTest, MatrixInvertSPD) {
    real_t M1[16] = {
        18, 22,  54,  42,
        22, 70,  86,  62,
        54, 86, 174, 134,
        42, 62, 134, 106
    };
    real_t expected[16] = {
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0
    };
    real_t temp[16], result[16];

    _inv_mat(result, M1, 4, temp);
    EXPECT_FLOAT_EQ(161/64.0, result[0]);
    EXPECT_FLOAT_EQ(31/64.0, result[1]);
    EXPECT_FLOAT_EQ(-83/64.0, result[2]);
    EXPECT_FLOAT_EQ(23/64.0, result[3]);
    EXPECT_FLOAT_EQ(31/64.0, result[4]);
    EXPECT_FLOAT_EQ(9/64.0, result[5]);
    EXPECT_FLOAT_EQ(-21/64.0, result[6]);
    EXPECT_FLOAT_EQ(9/64.0, result[7]);
    EXPECT_FLOAT_EQ(-83/64.0, result[8]);
    EXPECT_FLOAT_EQ(-21/64.0, result[9]);
    EXPECT_FLOAT_EQ(65/64.0, result[10]);
    EXPECT_FLOAT_EQ(-37/64.0, result[11]);
    EXPECT_FLOAT_EQ(23/64.0, result[12]);
    EXPECT_FLOAT_EQ(9/64.0, result[13]);
    EXPECT_FLOAT_EQ(-37/64.0, result[14]);
    EXPECT_FLOAT_EQ(33/64.0, result[15]);
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

TEST(C66xUKFTest, NoSensorConstantVelocity) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {10, 0, 0},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 1},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };
    real_t control[4] = { 0, 0, 0, 0 };

    ukf_init();
    ukf_set_state(&test_state);

    for(real_t i = 0; i < 10; i += 0.001) {
        ukf_iterate(0.001, control);
    }

    ukf_get_state(&test_state);

    EXPECT_FLOAT_EQ(1.5785803e-05, test_state.position[0]);
    EXPECT_FLOAT_EQ(0, test_state.position[1]);
    EXPECT_FLOAT_EQ(0, test_state.position[2]);

    EXPECT_FLOAT_EQ(10, test_state.velocity[X]);
    EXPECT_FLOAT_EQ(0, test_state.velocity[Y]);
    EXPECT_FLOAT_EQ(0, test_state.velocity[Z]);

    EXPECT_FLOAT_EQ(0, test_state.acceleration[X]);
    EXPECT_FLOAT_EQ(0, test_state.acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.acceleration[Z]);

    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(1, test_state.angular_velocity[Z]);

    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Z]);

    EXPECT_FLOAT_EQ(0, test_state.attitude[X]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[Y]);
    EXPECT_FLOAT_EQ(0.95878226, test_state.attitude[Z]);
    EXPECT_FLOAT_EQ(0.28414184, test_state.attitude[W]);
}

TEST(C66xUKFTest, NoSensorCircularMotion) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
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
    real_t control[4] = { 0, 0, 0, 0};

    ukf_init();
    ukf_set_state(&test_state);

    for(real_t i = 0; i < 10; i += 0.001) {
        ukf_iterate(0.001, control);
    }

    ukf_get_state(&test_state);

    EXPECT_NEAR(0, test_state.position[0], 1e-8);
    EXPECT_NEAR(0, test_state.position[1], 1e-12);
    EXPECT_FLOAT_EQ(0, test_state.position[2]);

    EXPECT_FLOAT_EQ(10, test_state.velocity[X]);
    EXPECT_NEAR(0, test_state.velocity[Y], 1e-2);
    EXPECT_FLOAT_EQ(0, test_state.velocity[Z]);

    EXPECT_FLOAT_EQ(0, test_state.acceleration[X]);
    EXPECT_FLOAT_EQ(6.28318530717959, test_state.acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.acceleration[Z]);

    EXPECT_FLOAT_EQ(0, test_state.attitude[X]);
    EXPECT_FLOAT_EQ(0, test_state.attitude[Y]);
    EXPECT_NEAR(0, test_state.attitude[Z], 1e-3);
    EXPECT_NEAR(1, test_state.attitude[W], 1e-6);

    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_velocity[Y]);
    EXPECT_FLOAT_EQ(0.628318530717959, test_state.angular_velocity[Z]);

    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Z]);
}

TEST(C66xUKFTest, AngularSensorsConstantAngularVelocity) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 1, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };
    struct ukf_ioboard_params_t test_config = {
        /* sensor offsets/orientations */
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0, 1},
        {1, 0, 0},
        /* sensor covariance */
        {0.001, 0.001, 0.001},
        {0.01, 0.01, 0.01},
        {0.01, 0.01, 0.01},
        {5, 5, 5},
        {5, 5, 5},
        1,
        1
    };
    real_t control[4] = { 0, 0, 0, 0 };

    ukf_init();
    ukf_set_state(&test_state);
    ukf_set_params(&test_config);

    real_t ref[4] = { 0, 0, 0, 1 };

    for (real_t i = 0; i < 10.0; i += 0.01) {
        real_t q[4] = { 0, -1, 0, 0 }, temp_q[4];
        _mul_quat_quat(temp_q, q, ref);
        ref[X] += 0.01 * 0.5 * temp_q[X];
        ref[Y] += 0.01 * 0.5 * temp_q[Y];
        ref[Z] += 0.01 * 0.5 * temp_q[Z];
        ref[W] += 0.01 * 0.5 * temp_q[W];
        _normalize_quat(ref, ref, true);

        real_t accel[3] = { 0, 0, -9.80665 }, accel_res[3];
        real_t mag[3] = { 1, 0, 0 }, mag_res[3];
        real_t gyro[3] = { 0, 1, 0 };

        _mul_quat_vec3(accel_res, ref, accel);
        _mul_quat_vec3(mag_res, ref, mag);

        ukf_sensor_clear();
        ukf_sensor_set_gyroscope(gyro[X], gyro[Y], gyro[Z]);
        ukf_sensor_set_accelerometer(
            accel_res[X], accel_res[Y], accel_res[Z]);
        ukf_sensor_set_magnetometer(mag_res[X], mag_res[Y], mag_res[Z]);

        ukf_iterate(0.01, control);
    }

    ukf_get_state(&test_state);

    EXPECT_NEAR(0, test_state.position[0], 1e-5);
    EXPECT_NEAR(0, test_state.position[1], 1e-18);
    EXPECT_NEAR(0, test_state.position[2], 3);

    EXPECT_NEAR(0, test_state.velocity[X], 2);
    EXPECT_NEAR(0, test_state.velocity[Y], 1e-12);
    EXPECT_NEAR(0, test_state.velocity[Z], 1);

    EXPECT_NEAR(0, test_state.acceleration[X], 1e-5);
    EXPECT_NEAR(0, test_state.acceleration[Y], 1e-16);
    EXPECT_NEAR(0, test_state.acceleration[Z], 1e-5);

    EXPECT_NEAR(0, test_state.attitude[X], 1e-6);
    EXPECT_NEAR(0.958924, test_state.attitude[Y], 2e-3);
    EXPECT_NEAR(0, test_state.attitude[Z], 1e-6);
    EXPECT_NEAR(0.283662, test_state.attitude[W], 5e-3);

    EXPECT_NEAR(0, test_state.angular_velocity[X], 1e-16);
    EXPECT_NEAR(1, test_state.angular_velocity[Y], 1e-4);
    EXPECT_NEAR(0, test_state.angular_velocity[Z], 1e-18);

    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Z]);

    EXPECT_FLOAT_EQ(0, test_state.wind_velocity[X]);
    EXPECT_FLOAT_EQ(0, test_state.wind_velocity[Y]);
    EXPECT_FLOAT_EQ(0, test_state.wind_velocity[Z]);

    EXPECT_NEAR(0, test_state.gyro_bias[X], 1e-18);
    EXPECT_NEAR(0, test_state.gyro_bias[Y], 1e-4);
    EXPECT_NEAR(0, test_state.gyro_bias[Z], 1e-18);
}

TEST(C66xUKFTest, MagnetometerAccelerometerAtRest) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };
    struct ukf_ioboard_params_t test_config = {
        /* sensor offsets/orientations */
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0, 1},
        {21.2578, 4.4132, -55.9578},
        /* sensor covariance */
        {10, 10, 10},
        {0.01, 0.01, 0.01},
        {25, 25, 25},
        {5, 5, 5},
        {5, 5, 5},
        1,
        1
    };
    real_t control[4] = { 0, 0, 0, 0 };

    ukf_init();
    ukf_set_state(&test_state);
    ukf_set_params(&test_config);
    ukf_choose_dynamics(UKF_MODEL_CENTRIPETAL);

    for(real_t i = 0; i < 10.0; i += 0.001) {
        ukf_sensor_clear();
        ukf_sensor_set_accelerometer(0, 0, -9.80665);
        ukf_sensor_set_magnetometer(-4.4132, 21.2578, -55.9578);
        ukf_sensor_set_gps_velocity(0, 0, 0);

        ukf_iterate(0.001, control);
    }

    ukf_get_state(&test_state);

    EXPECT_NEAR(0, test_state.position[0], 1e-7);
    EXPECT_NEAR(0, test_state.position[1], 1e-7);
    EXPECT_NEAR(0, test_state.position[2], 0.5);

    EXPECT_NEAR(0, test_state.velocity[X], 1e-2);
    EXPECT_NEAR(0, test_state.velocity[Y], 1e-3);
    EXPECT_NEAR(0, test_state.velocity[Z], 1e-2);

    EXPECT_NEAR(0, test_state.acceleration[X], 1e-5);
    EXPECT_NEAR(0, test_state.acceleration[Y], 1e-5);
    EXPECT_NEAR(0, test_state.acceleration[Z], 1e-5);

    EXPECT_NEAR(0, test_state.attitude[X], 1e-5);
    EXPECT_NEAR(0, test_state.attitude[Y], 1e-5);
    EXPECT_NEAR(0.707107, test_state.attitude[Z], 1e-5);
    EXPECT_NEAR(0.707107, test_state.attitude[W], 1e-5);

    EXPECT_NEAR(0, test_state.angular_velocity[X], 1e-5);
    EXPECT_NEAR(0, test_state.angular_velocity[Y], 1e-5);
    EXPECT_NEAR(0, test_state.angular_velocity[Z], 1e-4);

    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[X]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Y]);
    EXPECT_FLOAT_EQ(0, test_state.angular_acceleration[Z]);

    EXPECT_FLOAT_EQ(0, test_state.wind_velocity[X]);
    EXPECT_FLOAT_EQ(0, test_state.wind_velocity[Y]);
    EXPECT_FLOAT_EQ(0, test_state.wind_velocity[Z]);

    EXPECT_NEAR(0, test_state.gyro_bias[X], 1e-20);
    EXPECT_NEAR(0, test_state.gyro_bias[Y], 1e-20);
    EXPECT_NEAR(0, test_state.gyro_bias[Z], 1e-20);
}

TEST(C66xUKFTest, AllSensorsAtRest) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    struct ukf_state_t test_state = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };
    struct ukf_ioboard_params_t test_config = {
        /* sensor offsets/orientations */
        {0, 0, 0, 1},
        {0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0, 1},
        {21.2578, 4.4132, -55.9578},
        /* sensor covariance */
        {5, 5, 5},
        {0.1, 0.1, 0.1},
        {1, 1, 1},
        {5, 5, 5},
        {5, 5, 5},
        1,
        1
    };
    real_t control[4] = { 0, 0, 0, 0 };

    ukf_init();
    ukf_set_state(&test_state);
    ukf_set_params(&test_config);
    ukf_choose_dynamics(UKF_MODEL_CENTRIPETAL);

    for(real_t i = 0; i < 10.0; i += 0.001) {
        ukf_sensor_clear();
        ukf_sensor_set_gyroscope(0, 0, 0);
        ukf_sensor_set_accelerometer(0, 0, -9.80665);
        ukf_sensor_set_magnetometer(-4.4132, 21.2578, -55.9578);
        ukf_sensor_set_gps_velocity(0, 0, 0);
        ukf_sensor_set_gps_position(145.0 / 180.0 * M_PI, 37.0 / 180.0 * M_PI,
                                    50);
        ukf_sensor_set_barometer_amsl(53);
        ukf_sensor_set_pitot_tas(3);

        ukf_iterate(0.001, control);
    }

    ukf_get_state(&test_state);

    EXPECT_NEAR(0, test_state.attitude[X], 1e-3);
    EXPECT_NEAR(0, test_state.attitude[Y], 1e-3);
    EXPECT_NEAR(0.707107, test_state.attitude[Z], 1e-3);
    EXPECT_NEAR(0.707107, test_state.attitude[W], 1e-3);
}
