#ifndef CUKFMATH_H_
#define CUKFMATH_H_

#define X 0
#define Y 1
#define Z 2
#define W 3

#ifdef UKF_USE_DSP_INTRINSICS
#undef cos
#define cos Cosdp
#undef sin
#define sin Sindp
#undef atan2
#define atan2 Atan2dp
#undef sqrt
#define sqrt Sqrtdp
#define sqrt_inv(x) Rsqrtdp(x)
#else
#define sqrt_inv(x) (1.0 / sqrt((x)))
#endif

/* DEBUG */
#ifdef UKF_DEBUG
#include <stdio.h>

void _print_matrix(real_t mat[], size_t rows, size_t cols) {
    for (size_t i = 0; i < cols; i++) {
        for (size_t j = 0; j < rows; j++) {
            printf("%12.6g ", mat[j*cols + i]);
        }
        printf("\n");
    }
}

typedef uint64_t cycles_t;
static inline cycles_t rdtsc() {
    cycles_t result;
    __asm__ __volatile__ ("rdtsc" : "=A" (result));
    return result;
}
#endif
/* END DEBUG */

#ifndef M_PI
#define M_PI ((real_t)3.14159265358979323846)
#define M_PI_2 (M_PI / 2.0)
#define M_PI_4 (M_PI / 4.0)
#endif

static inline void _cross_vec3(real_t res[3], real_t v1[3], real_t v2[3]) {
    assert(res && v1 && v2 && res != v1 && res != v2);
    res[X] = v1[Y]*v2[Z] - v1[Z]*v2[Y];
    res[Y] = v1[Z]*v2[X] - v1[X]*v2[Z];
    res[Z] = v1[X]*v2[Y] - v1[Y]*v2[X];
}

static inline void _mul_quat_vec3(real_t res[3], real_t q[4], real_t v[3]) {
    /*
    Multiply a quaternion by a vector (i.e. transform a vectory by a
    quaternion)

    v' = q * v * conjugate(q), or:
    t = 2 * cross(q.xyz, v)
    v' = v + q.w * t + cross(q.xyz, t)

    http://molecularmusings.wordpress.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    */

    assert(res && q && v && res != v && res != q);

    real_t t[3] = {
        2.0 * (q[Y]*v[Z] - q[Z]*v[Y]),
        2.0 * (q[Z]*v[X] - q[X]*v[Z]),
        2.0 * (q[X]*v[Y] - q[Y]*v[X])
    };

    res[X] = v[X] + q[W]*t[X] + (q[Y]*t[Z] - q[Z]*t[Y]);
    res[Y] = v[Y] + q[W]*t[Y] + (q[Z]*t[X] - q[X]*t[Z]);
    res[Z] = v[Z] + q[W]*t[Z] + (q[X]*t[Y] - q[Y]*t[X]);
}

static inline void _mul_quat_quat(real_t res[4], real_t q1[4], real_t q2[4]) {
    assert(res && q1 && q2 && res != q1 && res != q2);

    res[W] = q1[W]*q2[W] - q1[X]*q2[X] - q1[Y]*q2[Y] - q1[Z]*q2[Z];
    res[X] = q1[W]*q2[X] + q1[X]*q2[W] + q1[Y]*q2[Z] - q1[Z]*q2[Y];
    res[Y] = q1[W]*q2[Y] - q1[X]*q2[Z] + q1[Y]*q2[W] + q1[Z]*q2[X];
    res[Z] = q1[W]*q2[Z] + q1[X]*q2[Y] - q1[Y]*q2[X] + q1[Z]*q2[W];
}

static inline void _normalize_quat(real_t res[4], real_t q[4],
bool force_pos) {
    assert(res && q);

    real_t q_norm = 1.0 / sqrt(q[X]*q[X] + q[Y]*q[Y] + q[Z]*q[Z] + q[W]*q[W]);
    if (force_pos && q[W] < 0.0) {
        q_norm = -q_norm;
    }

    res[X] = q[X] * q_norm;
    res[Y] = q[Y] * q_norm;
    res[Z] = q[Z] * q_norm;
    res[W] = q[W] * q_norm;
}

static inline void _mul_state_scalar_add_state(struct ukf_state_t *res,
struct ukf_state_t *s1, real_t a, struct ukf_state_t *s2) {
    assert(res && s1 && s2 && s1 != s2);

    res->position[0] = s1->position[0] * a + s2->position[0];
    res->position[1] = s1->position[1] * a + s2->position[1];
    res->position[2] = s1->position[2] * a + s2->position[2];

    res->velocity[0] = s1->velocity[0] * a + s2->velocity[0];
    res->velocity[1] = s1->velocity[1] * a + s2->velocity[1];
    res->velocity[2] = s1->velocity[2] * a + s2->velocity[2];

    res->acceleration[0] = s1->acceleration[0] * a + s2->acceleration[0];
    res->acceleration[1] = s1->acceleration[1] * a + s2->acceleration[1];
    res->acceleration[2] = s1->acceleration[2] * a + s2->acceleration[2];

    res->attitude[0] = s1->attitude[0] * a + s2->attitude[0];
    res->attitude[1] = s1->attitude[1] * a + s2->attitude[1];
    res->attitude[2] = s1->attitude[2] * a + s2->attitude[2];
    res->attitude[3] = s1->attitude[3] * a + s2->attitude[3];

    res->angular_velocity[0] = s1->angular_velocity[0] * a +
        s2->angular_velocity[0];
    res->angular_velocity[1] = s1->angular_velocity[1] * a +
        s2->angular_velocity[1];
    res->angular_velocity[2] = s1->angular_velocity[2] * a +
        s2->angular_velocity[2];

    res->angular_acceleration[0] = s1->angular_acceleration[0] * a +
        s2->angular_acceleration[0];
    res->angular_acceleration[1] = s1->angular_acceleration[1] * a +
        s2->angular_acceleration[1];
    res->angular_acceleration[2] = s1->angular_acceleration[2] * a +
        s2->angular_acceleration[2];

    res->wind_velocity[0] = s1->wind_velocity[0] * a + s2->wind_velocity[0];
    res->wind_velocity[1] = s1->wind_velocity[1] * a + s2->wind_velocity[1];
    res->wind_velocity[2] = s1->wind_velocity[2] * a + s2->wind_velocity[2];

    res->gyro_bias[0] = s1->gyro_bias[0] * a + s2->gyro_bias[0];
    res->gyro_bias[1] = s1->gyro_bias[1] * a + s2->gyro_bias[1];
    res->gyro_bias[2] = s1->gyro_bias[2] * a + s2->gyro_bias[2];
}

static inline void _add_state_accum(struct ukf_state_t *res,
struct ukf_state_t *s1) {
    assert(res && s1);

    res->position[0] += s1->position[0];
    res->position[1] += s1->position[1];
    res->position[2] += s1->position[2];

    res->velocity[0] += s1->velocity[0];
    res->velocity[1] += s1->velocity[1];
    res->velocity[2] += s1->velocity[2];

    res->acceleration[0] += s1->acceleration[0];
    res->acceleration[1] += s1->acceleration[1];
    res->acceleration[2] += s1->acceleration[2];

    res->attitude[0] += s1->attitude[0];
    res->attitude[1] += s1->attitude[1];
    res->attitude[2] += s1->attitude[2];
    res->attitude[3] += s1->attitude[3];

    res->angular_velocity[0] += s1->angular_velocity[0];
    res->angular_velocity[1] += s1->angular_velocity[1];
    res->angular_velocity[2] += s1->angular_velocity[2];

    res->angular_acceleration[0] += s1->angular_acceleration[0];
    res->angular_acceleration[1] += s1->angular_acceleration[1];
    res->angular_acceleration[2] += s1->angular_acceleration[2];

    res->wind_velocity[0] += s1->wind_velocity[0];
    res->wind_velocity[1] += s1->wind_velocity[1];
    res->wind_velocity[2] += s1->wind_velocity[2];

    res->gyro_bias[0] += s1->gyro_bias[0];
    res->gyro_bias[1] += s1->gyro_bias[1];
    res->gyro_bias[2] += s1->gyro_bias[2];
}

static inline void _inv_mat3x3(real_t out[9], real_t M[9]) {
    /*
    M = 0 1 2
        3 4 5
        6 7 8
    */
    assert(M && out);
    real_t det = M[0] * (M[8]*M[4] - M[5]*M[7]) -
                 M[1] * (M[8]*M[3] - M[5]*M[6]) +
                 M[2] * (M[7]*M[3] - M[4]*M[6]);

    assert(fabs(det) > 1e-6);
    det = 1.0 / det;

    out[0] =  (M[8]*M[4] - M[5]*M[7]) * det;
    out[3] = -(M[8]*M[3] - M[5]*M[6]) * det;
    out[6] =  (M[7]*M[3] - M[4]*M[6]) * det;
    out[1] = -(M[8]*M[1] - M[2]*M[7]) * det;
    out[4] =  (M[8]*M[0] - M[2]*M[6]) * det;
    out[7] = -(M[7]*M[0] - M[1]*M[6]) * det;
    out[2] =  (M[5]*M[1] - M[2]*M[4]) * det;
    out[5] = -(M[5]*M[0] - M[2]*M[3]) * det;
    out[8] =  (M[4]*M[0] - M[1]*M[3]) * det;
}

static inline void _mul_mat(real_t C[], real_t B[], real_t A[], size_t AR,
size_t AC, size_t BR, size_t BC, real_t mul) {
    /*
    Calculate C = AB, where A has AN rows and AM columns, and B has BN rows
    and BM columns. C must be AN x BM.
    */
    assert(A && B && C && AR == BC);
    uint32_t i, j, k;
    for (i = 0; i < AC; i++) {
        for (j = 0; j < BR; j++) {
            real_t t = A[i * AR] * B[j];
            for (k = 1; k < AR; k++) {
                t += A[i * AR + k] * B[k * BR + j];
            }

            C[i * BR + j] = t * mul;
        }
    }
}

static inline void _mul_mat_accum(real_t C[], real_t B[], real_t A[],
size_t AR, size_t AC, size_t BR, size_t BC, real_t mul) {
    /* Calculate C = C + AB */
    assert(A && B && C && AR == BC);
    uint32_t i, j, k, ci;
    for (i = 0; i < AC; i++) {
        for (j = 0; j < BR; j++) {
            ci = i * BR + j;
            for (k = 0; k < AR; k++) {
                C[ci] += A[i * AR + k] * B[k * BR + j] * mul;
            }
        }
    }
}

static inline void _add_mat_accum(real_t B[], real_t A[], size_t n) {
    /* Calculate B = B + A */
    assert(A && B);
    uint32_t i;
    for (i = 0; i < n; i++) {
        B[i] += A[i];
    }
}

static inline void _transpose_mat(real_t B[], real_t A[], size_t rows,
size_t cols) {
    /* Transpose the matrix A */
    assert(A && B && A != B && rows && cols);
    uint32_t i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            B[j*rows + i] = A[i*cols + j];
        }
    }
}

static inline void _cholesky_mat_mul(real_t L[], real_t A[], real_t mul,
size_t n) {
    assert(L && A && n);

    /*
    24x24:
    7464 mult
    552 div
    24 sqrt
    ~23,000 cycles, compared with ~6000 for Eigen
    */

    uint32_t i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            real_t s = 0;
            for (k = 0; k < j; k++) {
                s += L[i + k*n] * L[j + k*n];
            }

            L[i + j*n] = (i == j) ? sqrt(A[i + i*n]*mul - s) :
                (1.0 / L[j + j*n] * (A[i + j*n]*mul - s));
        }
    }
}

static inline void _inv_mat(real_t out[], real_t M[], size_t n,
real_t temp[]) {
    /*
    Matrix inverse via Cholesky decomposition. Requires M is symmetric
    positive-definite.

    Follows algorithm suggested in
    http://adorio-research.org/wordpress/?p=4560
    */
    assert(M && out && n && M != out && temp && M != temp && out != temp);

    _cholesky_mat_mul(temp, M, 1.0, n);

    /* TODO: check if this is required */
    memset(out, 0, sizeof(real_t) * n * n);

    int32_t i, j, k;
    for (j = n - 1; j >= 0; j--) {
        real_t tjj_recip = 1.0 / temp[j*n + j], s = 0.0;

        for (k = j + 1; k < n; k++) {
            s += temp[j*n + k] * out[j*n + k];
        }

        out[j*n + j] = tjj_recip*tjj_recip - s * tjj_recip;

        for (i = j - 1; i >= 0; i--) {
            s = 0.0;
            for (k = i + 1; k < n; k++) {
                s += temp[i*n + k] * out[k*n + j];
            }
            out[j*n + i] = out[i*n + j] = -s / temp[i*n + i];
        }
    }
}

#endif
