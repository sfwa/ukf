/*
Attitude errors, computed as the benchmark computes them
(examples/attitude_benchmark/metrics.py, compute_errors).

The error between an estimate and the reference is split into two parts:
the tilt error (inclination: how far "up" is wrong) and the heading error
(a turn about the vertical). A six-axis IMU cannot know its heading, so the
estimate is first turned about the vertical by one fixed angle, measured in
the first second of the reference. Scoring starts five seconds into the
recording.

Quaternions are w, x, y, z; they turn IMU-frame vectors into a world frame
whose +Z axis points up.
*/

const DEGREES = 180 / Math.PI;

function wrap(angle) {
    const turn = 2 * Math.PI;
    return ((angle + Math.PI) % turn + turn) % turn - Math.PI;
}

/*
Returns per-sample errors in degrees (NaN where there is no estimate or no
reference), the fixed heading offset in radians, and the RMS errors over the
scored samples.
*/
export function attitudeErrors(time, estimate, reference, valid, count = time.length,
        burnIn = 5, alignmentSeconds = 1) {
    const samples = time.length;
    const inclination = new Float32Array(samples).fill(NaN);
    const heading = new Float32Array(samples).fill(NaN);
    const total = new Float32Array(samples).fill(NaN);

    /* delta = estimate * conjugate(reference), for the samples which have both. */
    const delta = new Float64Array(4 * samples);
    const rawHeading = new Float64Array(samples).fill(NaN);
    let first = -1;
    for (let i = 0; i < count; i++) {
        if (!valid[i]) {
            continue;
        }
        if (first < 0) {
            first = i;
        }
        const aw = estimate[4 * i], ax = estimate[4 * i + 1], ay = estimate[4 * i + 2], az = estimate[4 * i + 3];
        const bw = reference[4 * i], bx = -reference[4 * i + 1], by = -reference[4 * i + 2], bz = -reference[4 * i + 3];
        let w = aw * bw - ax * bx - ay * by - az * bz;
        let x = aw * bx + ax * bw + ay * bz - az * by;
        let y = aw * by - ax * bz + ay * bw + az * bx;
        let z = aw * bz + ax * by - ay * bx + az * bw;
        const norm = Math.hypot(w, x, y, z);
        w /= norm; x /= norm; y /= norm; z /= norm;
        delta[4 * i] = w; delta[4 * i + 1] = x; delta[4 * i + 2] = y; delta[4 * i + 3] = z;
        const twist = Math.hypot(w, z);
        inclination[i] = 2 * Math.acos(Math.min(twist, 1)) * DEGREES;
        if (twist >= 1e-10) {
            rawHeading[i] = wrap(2 * Math.atan2(z, w));
        }
    }
    if (first < 0) {
        return {inclination, heading, total, offset: 0, rms: null, scored: 0, scoreStart: samples};
    }

    /* One fixed turn about the vertical, from the first second of reference. */
    const alignmentEnd = time[first] + alignmentSeconds;
    let sumCos = 0, sumSin = 0;
    for (let i = first; i < count && time[i] < alignmentEnd; i++) {
        if (valid[i] && Number.isFinite(rawHeading[i])) {
            sumCos += Math.cos(rawHeading[i]);
            sumSin += Math.sin(rawHeading[i]);
        }
    }
    const offset = Math.atan2(sumSin, sumCos);
    const cw = Math.cos(offset / 2), cz = -Math.sin(offset / 2);

    const scoreFrom = Math.max(time[0] + burnIn, alignmentEnd);
    let scoreStart = samples;
    const sums = {inclination: 0, heading: 0, total: 0};
    let scored = 0, scoredHeading = 0;
    for (let i = first; i < count; i++) {
        if (!valid[i]) {
            continue;
        }
        /* correction * delta, where correction = (cw, 0, 0, cz). Only w and z are needed. */
        const w = delta[4 * i], z = delta[4 * i + 3];
        const aw = cw * w - cz * z;
        const az = cw * z + cz * w;
        const twist = Math.hypot(aw, az);
        total[i] = 2 * Math.acos(Math.min(Math.abs(aw), 1)) * DEGREES;
        if (twist >= 1e-10) {
            heading[i] = wrap(2 * Math.atan2(az, aw)) * DEGREES;
        }
        if (time[i] >= scoreFrom) {
            if (scoreStart === samples) {
                scoreStart = i;
            }
            scored++;
            sums.inclination += inclination[i] ** 2;
            sums.total += total[i] ** 2;
            if (Number.isFinite(heading[i])) {
                sums.heading += heading[i] ** 2;
                scoredHeading++;
            }
        }
    }

    const rms = scored > 0 ? {
        inclination: Math.sqrt(sums.inclination / scored),
        heading: scoredHeading > 0 ? Math.sqrt(sums.heading / scoredHeading) : NaN,
        total: Math.sqrt(sums.total / scored),
    } : null;
    return {inclination, heading, total, offset, rms, scored, scoreStart};
}
