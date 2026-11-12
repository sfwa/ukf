/*
Run a packed recording through the estimator module in Node, with the same
JavaScript the page uses, and write the results for test_explainer.py to
compare with the native benchmark program.

Usage:
    node run_wasm.mjs MODULE.wasm RECORDING.pack OUTPUT_DIRECTORY

Writes, as little-endian doubles:
    imu.f64          time, gyroscope x y z, accelerometer x y z per sample
    reference.f64    reference attitude w x y z per sample
    mav-srukf.f64, mav-srukf-gravity.f64
                     estimated attitude w x y z per sample
and summary.json with the page's error figures and the step checks.
*/

import fs from 'node:fs';
import path from 'node:path';
import {decodePack} from '../web/js/data/pack.js';
import {EstimatorModule} from '../web/js/data/wasm.js';
import {attitudeErrors} from '../web/js/data/metrics.js';

const [modulePath, packPath, outputDirectory] = process.argv.slice(2);
if (!modulePath || !packPath || !outputDirectory) {
    console.error('Usage: node run_wasm.mjs MODULE.wasm RECORDING.pack OUTPUT_DIRECTORY');
    process.exit(2);
}

const write = (name, values) => fs.writeFileSync(path.join(outputDirectory, name),
    Buffer.from(values.buffer, values.byteOffset, values.byteLength));

const packBytes = fs.readFileSync(packPath);
const recording = await decodePack(packBytes.buffer.slice(packBytes.byteOffset, packBytes.byteOffset + packBytes.byteLength));
const estimator = await EstimatorModule.create(fs.readFileSync(modulePath));
const samples = recording.samples;

const imu = new Float64Array(7 * samples);
for (let i = 0; i < samples; i++) {
    imu[7 * i] = recording.time[i];
    for (let axis = 0; axis < 3; axis++) {
        imu[7 * i + 1 + axis] = recording.gyroscope[3 * i + axis];
        imu[7 * i + 4 + axis] = recording.accelerometer[3 * i + axis];
    }
}
write('imu.f64', imu);
write('reference.f64', recording.reference.attitude);

estimator.load(recording.time, recording.gyroscope, recording.accelerometer);
const summary = {
    sequence: recording.name,
    samples,
    constants: estimator.constants,
    defaultTuning: estimator.defaultTuning,
    filters: {},
};
const attitudeField = estimator.traceFields.attitude;
for (const [filter, tuning] of [['mav-srukf', {}], ['mav-srukf-gravity', {use_rotor_model: 0}]]) {
    const started = performance.now();
    const result = estimator.run(tuning);
    const elapsed = performance.now() - started;
    const attitude = new Float64Array(4 * samples);
    for (let i = 0; i < samples; i++) {
        for (let k = 0; k < 4; k++) {
            attitude[4 * i + k] = result.trace[i * estimator.traceStride + attitudeField.offset + k];
        }
    }
    write(`${filter}.f64`, attitude);

    /* Take steps apart all through the recording, and keep the worst check. */
    let worstCheck = 0, steps = 0;
    for (let i = 1; i < samples; i += 211) {
        const step = estimator.explain(i);
        if (step === null) {
            throw new Error(`Cannot take step ${i} apart`);
        }
        worstCheck = Math.max(worstCheck, step[estimator.stepFields.check.offset]);
        steps++;
    }

    const errors = attitudeErrors(recording.time, attitude, recording.reference.attitude,
        recording.reference.valid, result.valid);
    summary.filters[filter] = {
        valid: result.valid,
        elapsed_ms: elapsed,
        explained_steps: steps,
        worst_step_difference: worstCheck,
        heading_offset_rad: errors.offset,
        scored_samples: errors.scored,
        rms_deg: errors.rms,
    };
}
fs.writeFileSync(path.join(outputDirectory, 'summary.json'), JSON.stringify(summary, null, 2) + '\n');
console.log(JSON.stringify(Object.fromEntries(Object.entries(summary.filters).map(
    ([name, value]) => [name, {valid: value.valid, rms_deg: value.rms_deg, elapsed_ms: Math.round(value.elapsed_ms)}]))));
