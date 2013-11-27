#ifndef CUKFMATH_H_
#define CUKFMATH_H_

#define X 0
#define Y 1
#define Z 2
#define W 3

#ifndef absval
#define absval(x) ((x) < 0 ? -x : x)
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifdef UKF_USE_DSP_INTRINSICS

static inline double sqrt_inv(double a) {
    double x, half_a;

    if (a < 0) {
        x = _lltod(0x7FEFFFFFFFFFFFFF);
    } else {
        half_a = a * 0.5;
        x = _rsqrdp(a);
        x = x * (1.5 - half_a*x*x);
        x = x * (1.5 - half_a*x*x);
        x = x * (1.5 - half_a*x*x);
    }

    return x;
}

static inline double fsqrt(double a) {
    double  y, X0, X1, X2, X4;
    int     upper;

    upper = _clr(_hi(a), 31, 31);
    y = _itod(upper, _lo(a));

    X0 = _rsqrdp(y);
    X1 = X0 * (1.5 - (y*X0*X0*0.5));
    X2 = X1 * (1.5 - (y*X1*X1*0.5));
    X4 = y * X2 * (1.5 - (y*X2*X2*0.5));

    if (a <= 0.0) {
        X4 = 0.0;
    }
    if (a > 1.7976931348623157E+308) {
        X4 = 1.7976931348623157E+308;
    }

    return X4;
}

static inline double divide(double a, double b) {
    double  x;
    if (a == 0.0) {
        x = 0.0;
    } else {
        x = _rcpdp(b);
        x = x * (2.0 - b*x);
        x = x * (2.0 - b*x);
        x = x * (2.0 - b*x);
        x = x * a;
    }

    return x;
}

static inline double recip(double b) {
    double  x;
    x = _rcpdp(b);
    x = x * (2.0 - b*x);
    x = x * (2.0 - b*x);
    x = x * (2.0 - b*x);

    return x;
}

#else

#define sqrt_inv(x) (1.0 / sqrt((x)))
#define divide(a, b) ((a) / (b))
#define recip(a) (1.0 / (a))
#define fsqrt(a) sqrt((a))

#endif

static inline double fatan2(double y, double x) {
#define PI_FLOAT     3.14159265
#define PIBY2_FLOAT  1.5707963
// |error| < 0.005
    if (x == 0.0) {
        if (y > 0.0) return PIBY2_FLOAT;
        if (y == 0.0) return 0.0;
        return -PIBY2_FLOAT;
    }
    double atan;
    double z = divide(y, x);
    if (absval(z) < 1.0) {
        atan = divide(z, (1.0 + 0.28 * z * z));
        if (x < 0.0) {
            if (y < 0.0) return atan - PI_FLOAT;
            return atan + PI_FLOAT;
        }
    } else {
        atan = PIBY2_FLOAT - z / (z * z + 0.28);
        if (y < 0.0f) return atan - PI_FLOAT;
    }
    return atan;
#undef PI_FLOAT
#undef PIBY2_FLOAT
}

/* DEBUG */
#ifdef UKF_DEBUG
#include <stdio.h>

void _print_matrix(const char *label, real_t mat[], size_t rows,
size_t cols) {
    printf("%s", label);
    for (size_t i = 0; i < cols; i++) {
        for (size_t j = 0; j < rows; j++) {
            printf("%12.6g ", mat[j*cols + i]);
        }
        printf("\n");
    }
}

#ifndef __TI_COMPILER_VERSION__
static inline uint64_t rdtsc() {
    uint64_t result;
    __asm__ __volatile__ (
        "rdtsc" \
        : "=A" (result) \
        : /* no input parameters*/ \
        : "%eax", "%edx", "memory");
    return result;
}
#else
#define rdtsc() 0
#endif

#else

#define rdtsc() 0
#define _print_matrix(a, b, c, d)

#undef assert
#define assert(x)

#endif
/* END DEBUG */

/* Non-TI compatibility */
#ifndef __TI_COMPILER_VERSION__
#define _nassert(x)
#endif

#ifndef M_PI
#define M_PI ((real_t)3.14159265358979323846)
#define M_PI_2 (M_PI * 0.5)
#define M_PI_4 (M_PI * 0.25)
#endif

static inline void _cross_vec3(real_t *restrict res,
const real_t *restrict v1, const real_t *restrict v2) {
    assert(res && v1 && v2 && res != v1 && res != v2 && v1 != v2);
    _nassert((size_t)res % 8 == 0);
    _nassert((size_t)v1 % 8 == 0);
    _nassert((size_t)v2 % 8 == 0);

    real_t r0, r1, r2, r3, r4, r5;

    r0 = v1[Y]*v2[Z];
    r1 = v1[Z]*v2[Y];
    r2 = v1[Z]*v2[X];
    r3 = v1[X]*v2[Z];
    r4 = v1[X]*v2[Y];
    r5 = v1[Y]*v2[X];

    res[X] = r0 - r1;
    res[Y] = r2 - r3;
    res[Z] = r4 - r5;
}

static inline void _mul_quat_vec3(real_t *restrict res,
const real_t *restrict q, const real_t *restrict v) {
    /*
    Multiply a quaternion by a vector (i.e. transform a vectory by a
    quaternion)

    v' = q * v * conjugate(q), or:
    t = 2 * cross(q.xyz, v)
    v' = v + q.w * t + cross(q.xyz, t)

    http://molecularmusings.wordpress.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    */

    assert(res && q && v && res != v && res != q);
    _nassert((size_t)res % 8 == 0);
    _nassert((size_t)q % 8 == 0);
    _nassert((size_t)v % 8 == 0);

    register real_t rx, ry, rz, tx, ty, tz;

    tx = q[Y]*v[Z];
    ty = q[Z]*v[X];
    tx -= q[Z]*v[Y];
    ty -= q[X]*v[Z];
    tz = q[X]*v[Y];
    ty *= 2.0;
    tz -= q[Y]*v[X];
    tx *= 2.0;
    tz *= 2.0;

    rx = v[X];
    rx += q[W]*tx;
    rx += q[Y]*tz;
    rx -= q[Z]*ty;
    res[X] = rx;

    ry = v[Y];
    ry += q[W]*ty;
    ry += q[Z]*tx;
    ry -= q[X]*tz;
    res[Y] = ry;

    rz = v[Z];
    rz += q[W]*tz;
    rz -= q[Y]*tx;
    rz += q[X]*ty;
    res[Z] = rz;
}

static inline void _mul_quat_quat(real_t res[4],
const real_t q1[4], const real_t q2[4]) {
    assert(res && q1 && q2 && res != q1 && res != q2);
    _nassert((size_t)res % 8 == 0);
    _nassert((size_t)q1 % 8 == 0);
    _nassert((size_t)q2 % 8 == 0);

    res[W] = q1[W]*q2[W] - q1[X]*q2[X] - q1[Y]*q2[Y] - q1[Z]*q2[Z];
    res[X] = q1[W]*q2[X] + q1[X]*q2[W] + q1[Y]*q2[Z] - q1[Z]*q2[Y];
    res[Y] = q1[W]*q2[Y] - q1[X]*q2[Z] + q1[Y]*q2[W] + q1[Z]*q2[X];
    res[Z] = q1[W]*q2[Z] + q1[X]*q2[Y] - q1[Y]*q2[X] + q1[Z]*q2[W];
}

static inline void _normalize_quat(real_t res[4], const real_t q[4],
bool force_pos) {
    assert(res && q);
    _nassert((size_t)res % 8 == 0);
    _nassert((size_t)q % 8 == 0);

    real_t q_norm = sqrt_inv(q[X]*q[X] + q[Y]*q[Y] + q[Z]*q[Z] + q[W]*q[W]);
    if (force_pos && q[W] < 0.0) {
        q_norm = -q_norm;
    }

    res[X] = q[X] * q_norm;
    res[Y] = q[Y] * q_norm;
    res[Z] = q[Z] * q_norm;
    res[W] = q[W] * q_norm;
}

static inline void _mul_state_scalar_add_state(
struct ukf_state_t *restrict res,
const struct ukf_state_t *restrict s1, const real_t a,
const struct ukf_state_t *restrict s2) {
    assert(res && s1 && s2 && s1 != s2 && res != s1 && res != s2);
    _nassert((size_t)res % 8 == 0);
    _nassert((size_t)s1 % 8 == 0);
    _nassert((size_t)s2 % 8 == 0);

    uint32_t i;
    real_t *const restrict s1ptr = (real_t*)s1,
           *const restrict s2ptr = (real_t*)s2,
           *const restrict rptr = (real_t*)res;

    #pragma MUST_ITERATE(24)
    #pragma UNROLL(2)
    for (i = 0; i < UKF_STATE_DIM; i++) {
        rptr[i] = s2ptr[i] + s1ptr[i] * a;
    }
    rptr[i] = s1ptr[i] * a + s2ptr[i];
}

static inline void _mul_state_inplace(struct ukf_state_t *restrict res,
real_t a) {
    assert(res);
    _nassert((size_t)res % 8 == 0);

    uint32_t i;
    real_t *const restrict rptr = (real_t*)res;

    #pragma MUST_ITERATE(24)
    #pragma UNROLL(2)
    for (i = 0; i < UKF_STATE_DIM; i++) {
        rptr[i] *= a;
    }
    rptr[i] *= a;
}

static inline void _add_state_inplace(struct ukf_state_t *restrict res,
struct ukf_state_t *restrict s1) {
    assert(res && s1 && res != s1);
    _nassert((size_t)res % 8 == 0);
    _nassert((size_t)s1 % 8 == 0);

    uint32_t i;
    real_t *const restrict s1ptr = (real_t*)s1,
           *const restrict rptr = (real_t*)res;

    #pragma MUST_ITERATE(24)
    #pragma UNROLL(2)
    for (i = 0; i < UKF_STATE_DIM; i++) {
        rptr[i] += s1ptr[i];
    }
    rptr[i] += s1ptr[i];
}

static inline void _mul_wprime(real_t *restrict out,
const real_t *restrict in,
const real_t weight) {
    /*
    Multiply a column of wprime by its transpose, and add it to the state
    */
    assert(out && in && out != in);
    _nassert((size_t)out % 8 == 0);
    _nassert((size_t)in % 8 == 0);

    real_t *restrict outcol1, *restrict outcol2, *restrict outcol3,
           *restrict outcol4;
    register real_t colv1, colv2, colv3, colv4;
    uint32_t i, j;

    for (i = 0; i < UKF_STATE_DIM; i += 4) {
        colv1 = in[i] * weight;
        colv2 = in[i+1] * weight;
        colv3 = in[i+2] * weight;
        colv4 = in[i+3] * weight;
        outcol1 = &out[i * UKF_STATE_DIM];
        outcol2 = &out[(i+1) * UKF_STATE_DIM];
        outcol3 = &out[(i+2) * UKF_STATE_DIM];
        outcol4 = &out[(i+3) * UKF_STATE_DIM];

        #pragma MUST_ITERATE(UKF_STATE_DIM)
        for (j = 0; j < UKF_STATE_DIM; j++) {
            outcol1[j] = colv1 * in[j];
            outcol2[j] = colv2 * in[j];
            outcol3[j] = colv3 * in[j];
            outcol4[j] = colv4 * in[j];
        }
    }
}

static inline void _mul_wprime_accum(real_t *restrict out,
const real_t *restrict in,
const real_t weight) {
    /*
    Multiply a column of wprime by its transpose, and add it to the state
    */
    assert(out && in && out != in);
    _nassert((size_t)out % 8 == 0);
    _nassert((size_t)in % 8 == 0);

    real_t *restrict outcol1, *restrict outcol2, *restrict outcol3,
           *restrict outcol4;
    register real_t colv1, colv2, colv3, colv4;
    uint32_t i, j;

    for (i = 0; i < UKF_STATE_DIM; i += 4) {
        colv1 = in[i] * weight;
        colv2 = in[i+1] * weight;
        colv3 = in[i+2] * weight;
        colv4 = in[i+3] * weight;
        outcol1 = &out[i * UKF_STATE_DIM];
        outcol2 = &out[(i+1) * UKF_STATE_DIM];
        outcol3 = &out[(i+2) * UKF_STATE_DIM];
        outcol4 = &out[(i+3) * UKF_STATE_DIM];

        #pragma MUST_ITERATE(UKF_STATE_DIM)
        for (j = 0; j < UKF_STATE_DIM; j++) {
            outcol1[j] += colv1 * in[j];
            outcol2[j] += colv2 * in[j];
            outcol3[j] += colv3 * in[j];
            outcol4[j] += colv4 * in[j];
        }
    }
}

static inline void _mul_vec_outer(real_t *restrict out,
const real_t *restrict in1, const real_t *restrict in2, const size_t r1,
const size_t c2, const real_t weight) {
    /*
    Multiply a column of wprime by its transpose, and add it to the state
    */
    assert(out && in1 && in2 && out != in1 && out != in2);
    _nassert((size_t)out % 8 == 0);
    _nassert((size_t)in1 % 8 == 0);
    _nassert((size_t)in2 % 8 == 0);

    real_t *restrict outcol;
    register real_t colv;
    uint32_t i, j;

    for (i = 0; i < r1; i++) {
        colv = in1[i] * weight;
        outcol = &out[i * c2];

        #pragma MUST_ITERATE(1,24)
        #pragma UNROLL(2)
        for (j = 0; j < c2; j++) {
            outcol[j] = colv * in2[j];
        }
    }
}

static inline void _mul_vec_outer_accum(real_t *restrict out,
const real_t *restrict in1, const real_t *restrict in2, const size_t r1,
const size_t c2, const real_t weight) {
    /*
    Multiply a column of wprime by its transpose, and add it to the state
    */
    assert(out && in1 && in2 && out != in1 && out != in2);
    _nassert((size_t)out % 8 == 0);
    _nassert((size_t)in1 % 8 == 0);
    _nassert((size_t)in2 % 8 == 0);

    real_t *restrict outcol;
    register real_t colv;
    uint32_t i, j;

    for (i = 0; i < r1; i++) {
        colv = in1[i] * weight;
        outcol = &out[i * c2];

        #pragma MUST_ITERATE(1,24)
        #pragma UNROLL(2)
        for (j = 0; j < c2; j++) {
            outcol[j] += colv * in2[j];
        }
    }
}

static void _mul_mat(real_t *restrict C, const real_t B[],
const real_t A[], const size_t AR, const size_t AC, const size_t BR,
const size_t BC, const real_t mul) {
#pragma unused(BC)
    /*
    Calculate C = AB, where A has AN rows and AM columns, and B has BN rows
    and BM columns. C must be AN x BM.
    */
    assert(A && B && C && AR == BC);
    _nassert((size_t)C % 8 == 0);
    _nassert((size_t)B % 8 == 0);
    _nassert((size_t)A % 8 == 0);

    size_t i, j, k;
    #pragma MUST_ITERATE(1,24)
    for (i = 0; i < AC; i++) {
        #pragma MUST_ITERATE(1,24)
        for (j = 0; j < BR; j++) {
            real_t t = 0.0;
            #pragma MUST_ITERATE(1,24)
            for (k = 0; k < AR; k++) {
                t += A[i * AR + k] * B[k * BR + j];
            }

            C[i * BR + j] = t * mul;
        }
    }
}

static void _add_mat_accum(real_t B[], const real_t A[], const size_t n) {
    /* Calculate B = B + A */
    assert(A && B);
    _nassert((size_t)B % 8 == 0);
    _nassert((size_t)A % 8 == 0);

    uint32_t i;
    #pragma MUST_ITERATE(1,24*24)
    for (i = 0; i < n; i++) {
        B[i] += A[i];
    }
}

static void _transpose_mat(real_t *restrict B, const real_t *restrict A,
const size_t rows, const size_t cols) {
    /* Transpose the matrix A */
    assert(A && B && A != B && rows && cols);
    _nassert((size_t)B % 8 == 0);
    _nassert((size_t)A % 8 == 0);

    uint32_t i, j, ic, jr;
    #pragma MUST_ITERATE(1,24)
    for (i = 0, ic = 0; i < rows; i++, ic += cols) {
        #pragma MUST_ITERATE(1,24)
        for (j = 0, jr = 0; j < cols; j++, jr += rows) {
            B[jr + i] = A[ic + j];
        }
    }
}

static void _cholesky_mat_mul(real_t L[], const real_t A[], const real_t mul,
const size_t n) {
    assert(L && A && n);
    _nassert((size_t)L % 8 == 0);
    _nassert((size_t)A % 8 == 0);

    /*
    24x24:
    7464 mult
    552 div
    24 sqrt
    ~23,000 cycles, compared with ~6000 for Eigen
    ~ 110,000 cycles on C66x
    */

    size_t i, j, kn, in, jn;
    for (i = 0, in = 0; i < n; i++, in += n) {
        L[i + 0] = (i == 0) ? fsqrt(A[i + in]*mul) : recip(L[0]) * (A[i]*mul);

        for (j = 1, jn = n; j <= i; j++, jn += n) {
            real_t s = 0;
            #pragma MUST_ITERATE(1,24)
            for (kn = 0; kn < j*n; kn += n) {
                s += L[i + kn] * L[j + kn];
            }

            L[i + jn] = (i == j) ? fsqrt(A[i + in]*mul - s) :
                recip(L[j + jn]) * (A[i + jn]*mul - s);
        }
    }
}

static void _inv_mat(real_t *restrict out, const real_t *restrict M,
const size_t n, real_t *restrict temp) {
    /*
    Matrix inverse via Cholesky decomposition. Requires M is symmetric
    positive-definite, and temp is a scratch space the same size as M.

    Follows algorithm suggested in
    http://adorio-research.org/wordpress/?p=4560
    */
    assert(M && out && n && M != out && temp && M != temp && out != temp);
    _nassert((size_t)out % 8 == 0);
    _nassert((size_t)M % 8 == 0);
    _nassert((size_t)temp % 8 == 0);

    _cholesky_mat_mul(temp, M, 1.0, n);

    int32_t i, j;
    size_t k, jn, in;
    real_t tjj_recip, s;
    for (j = (int32_t)n - 1; j >= 0; j--) {
        jn = (size_t)j*n;
        tjj_recip = recip(temp[jn + (size_t)j]);
        s = 0.0;

        for (k = (size_t)j + 1; k < n; k++) {
            s += temp[jn + k] * out[jn + k];
        }

        out[jn + (size_t)j] = tjj_recip*tjj_recip - s * tjj_recip;

        for (i = j - 1; i >= 0; i--) {
            in = (size_t)i*n;
            s = 0.0;
            for (k = (size_t)i + 1; k < n; k++) {
                s += temp[in + k] * out[k*n+ (size_t)j];
            }
            out[jn + (size_t)i] = out[in + (size_t)j] =
                divide(-s, temp[in + (size_t)i]);
        }
    }
}

#endif
