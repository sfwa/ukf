/*
Small maths helpers for the figures: vectors, quaternions, 4x4 matrices,
bell curves and repeatable random numbers.

Vectors are plain arrays. Quaternions are [w, x, y, z] and turn vectors from
a body frame into the world frame, as in the estimator's output. Matrices
are 16 numbers, column by column, as WebGL takes them.
*/

export const DEG = Math.PI / 180;
export const clamp = (value, low, high) => Math.min(high, Math.max(low, value));
export const lerp = (a, b, t) => a + (b - a) * t;
export const smooth = t => t * t * (3 - 2 * t);

/* Vectors. */
export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const norm = a => Math.hypot(a[0], a[1], a[2]);
export const normalize = a => {
    const length = norm(a);
    return length > 0 ? scale(a, 1 / length) : [0, 0, 0];
};
export const mix = (a, b, t) => [lerp(a[0], b[0], t), lerp(a[1], b[1], t), lerp(a[2], b[2], t)];

/* Quaternions [w, x, y, z]. */
export const qIdentity = () => [1, 0, 0, 0];
export const qConj = q => [q[0], -q[1], -q[2], -q[3]];
export const qMul = (a, b) => [
    a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
    a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
    a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
    a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
];
export const qNormalize = q => {
    const length = Math.hypot(q[0], q[1], q[2], q[3]);
    return [q[0] / length, q[1] / length, q[2] / length, q[3] / length];
};
export const qFromAxisAngle = (axis, angle) => {
    const unit = normalize(axis);
    const s = Math.sin(angle / 2);
    return [Math.cos(angle / 2), unit[0] * s, unit[1] * s, unit[2] * s];
};
/* A rotation vector: its direction is the axis, its length the angle. */
export const qFromRotationVector = r => {
    const angle = norm(r);
    return angle < 1e-12 ? [1, r[0] / 2, r[1] / 2, r[2] / 2] : qFromAxisAngle(r, angle);
};
export const qRotate = (q, v) => {
    /* v + 2 w (u x v) + 2 u x (u x v), with u the vector part. */
    const u = [q[1], q[2], q[3]];
    const t = scale(cross(u, v), 2);
    return add(add(v, scale(t, q[0])), cross(u, t));
};
export const qSlerp = (a, b, t) => {
    let cosine = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    const sign = cosine < 0 ? -1 : 1;
    cosine *= sign;
    if (cosine > 0.9995) {
        return qNormalize([0, 1, 2, 3].map(i => lerp(a[i], sign * b[i], t)));
    }
    const angle = Math.acos(cosine);
    const wa = Math.sin((1 - t) * angle) / Math.sin(angle);
    const wb = sign * Math.sin(t * angle) / Math.sin(angle);
    return [0, 1, 2, 3].map(i => wa * a[i] + wb * b[i]);
};
/* The shortest rotation which takes direction a to direction b. */
export const qBetween = (a, b) => {
    const from = normalize(a), to = normalize(b);
    const cosine = dot(from, to);
    if (cosine < -0.999999) {
        const axis = Math.abs(from[0]) < 0.9 ? cross(from, [1, 0, 0]) : cross(from, [0, 1, 0]);
        return qFromAxisAngle(axis, Math.PI);
    }
    const axis = cross(from, to);
    return qNormalize([1 + cosine, axis[0], axis[1], axis[2]]);
};
/* Read quaternion number "index" from a flat array of w x y z. */
export const qAt = (array, index) => [array[4 * index], array[4 * index + 1], array[4 * index + 2], array[4 * index + 3]];
export const vAt = (array, index) => [array[3 * index], array[3 * index + 1], array[3 * index + 2]];

/* 4x4 matrices, column by column. */
export const mIdentity = () => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
export const mMul = (a, b) => {
    const out = new Array(16);
    for (let column = 0; column < 4; column++) {
        for (let row = 0; row < 4; row++) {
            let sum = 0;
            for (let k = 0; k < 4; k++) {
                sum += a[k * 4 + row] * b[column * 4 + k];
            }
            out[column * 4 + row] = sum;
        }
    }
    return out;
};
export const mTranslate = v => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, v[0], v[1], v[2], 1];
export const mScale = s => {
    const v = typeof s === 'number' ? [s, s, s] : s;
    return [v[0], 0, 0, 0, 0, v[1], 0, 0, 0, 0, v[2], 0, 0, 0, 0, 1];
};
export const mFromQuat = q => {
    const [w, x, y, z] = q;
    return [
        1 - 2 * (y * y + z * z), 2 * (x * y + w * z), 2 * (x * z - w * y), 0,
        2 * (x * y - w * z), 1 - 2 * (x * x + z * z), 2 * (y * z + w * x), 0,
        2 * (x * z + w * y), 2 * (y * z - w * x), 1 - 2 * (x * x + y * y), 0,
        0, 0, 0, 1,
    ];
};
/* A matrix whose columns are the given axes, at the given origin. */
export const mFromAxes = (x, y, z, origin = [0, 0, 0]) =>
    [x[0], x[1], x[2], 0, y[0], y[1], y[2], 0, z[0], z[1], z[2], 0, origin[0], origin[1], origin[2], 1];
export const mPose = (position, q, size = 1) => mMul(mTranslate(position), mMul(mFromQuat(q), mScale(size)));
export const mPerspective = (fovY, aspect, near, far) => {
    const f = 1 / Math.tan(fovY / 2);
    return [f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) / (near - far), -1, 0, 0, 2 * far * near / (near - far), 0];
};
export const mLookAt = (eye, target, up) => {
    const z = normalize(sub(eye, target));
    const x = normalize(cross(up, z));
    const y = cross(z, x);
    return [x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0, -dot(x, eye), -dot(y, eye), -dot(z, eye), 1];
};
export const mTransformPoint = (m, v) => {
    const x = m[0] * v[0] + m[4] * v[1] + m[8] * v[2] + m[12];
    const y = m[1] * v[0] + m[5] * v[1] + m[9] * v[2] + m[13];
    const z = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14];
    const w = m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15];
    return [x, y, z, w];
};

/* Bell curves. */
export const gaussian = (x, mean, sd) => Math.exp(-0.5 * ((x - mean) / sd) ** 2) / (sd * Math.sqrt(2 * Math.PI));

/*
The ellipse which holds a two-dimensional bell curve out to "sigmas"
standard deviations: the half-lengths of its axes and the angle of the long
one. The covariance is [[a, b], [b, c]].
*/
export const covarianceEllipse = (a, b, c, sigmas = 1) => {
    const mean = (a + c) / 2;
    const spread = Math.sqrt(((a - c) / 2) ** 2 + b * b);
    const major = Math.sqrt(Math.max(mean + spread, 0)) * sigmas;
    const minor = Math.sqrt(Math.max(mean - spread, 0)) * sigmas;
    const angle = Math.abs(b) < 1e-300 && a >= c ? 0 : Math.atan2(mean + spread - a, b);
    return {major, minor, angle};
};

/* Lower-triangular square root of a small symmetric matrix (array of rows). */
export const cholesky = matrix => {
    const n = matrix.length;
    const root = matrix.map(row => row.map(() => 0));
    for (let i = 0; i < n; i++) {
        for (let j = 0; j <= i; j++) {
            let sum = matrix[i][j];
            for (let k = 0; k < j; k++) {
                sum -= root[i][k] * root[j][k];
            }
            root[i][j] = i === j ? Math.sqrt(Math.max(sum, 0)) : (root[j][j] > 0 ? sum / root[j][j] : 0);
        }
    }
    return root;
};

/* Repeatable random numbers, so that a figure looks the same on every visit. */
export class Random {
    constructor(seed = 1) {
        this.state = seed >>> 0;
        this.spare = null;
    }

    /* Uniform in [0, 1) (the "mulberry32" generator). */
    uniform() {
        this.state = (this.state + 0x6d2b79f5) >>> 0;
        let t = this.state;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    }

    /* Bell-shaped, mean 0, standard deviation 1. */
    normal() {
        if (this.spare !== null) {
            const value = this.spare;
            this.spare = null;
            return value;
        }
        let u = 0;
        while (u === 0) {
            u = this.uniform();
        }
        const radius = Math.sqrt(-2 * Math.log(u));
        const angle = 2 * Math.PI * this.uniform();
        this.spare = radius * Math.sin(angle);
        return radius * Math.cos(angle);
    }
}

/* Round a tick spacing to 1, 2 or 5 times a power of ten. */
export const niceStep = (range, targetTicks) => {
    const rough = range / Math.max(targetTicks, 1);
    const power = 10 ** Math.floor(Math.log10(rough));
    const fraction = rough / power;
    return (fraction < 1.5 ? 1 : fraction < 3.5 ? 2 : fraction < 7.5 ? 5 : 10) * power;
};
