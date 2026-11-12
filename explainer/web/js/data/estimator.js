/*
The page's side of the estimator: start it, give it a recording, run it,
and read the results.

The estimator is the benchmark's MAV_Estimator, compiled to WebAssembly
(explainer/wasm). It runs in a worker when the browser allows, and on the
main thread otherwise.
*/

import {EstimatorModule} from './wasm.js';
import {attitudeErrors} from './metrics.js';

let compiled = null;

function compileModule(url) {
    compiled ??= (async () => {
        try {
            return await WebAssembly.compileStreaming(fetch(url));
        } catch {
            /* Some servers do not label .wasm files; compile from the bytes instead. */
            const response = await fetch(url);
            if (!response.ok) {
                throw new Error(`Cannot load the estimator (${response.status})`);
            }
            return WebAssembly.compile(await response.arrayBuffer());
        }
    })();
    return compiled;
}

/* The result of one run of the estimator over a recording. */
export class Run {
    constructor(description, recording, result, tuning) {
        this.recording = recording;
        this.tuning = tuning;
        this.trace = result.trace;
        this.valid = result.valid;                 /* samples which have an estimate */
        this.failed = result.valid < result.samples;
        this.elapsed = result.elapsed ?? 0;
        this.stride = description.traceStride;
        this.fields = description.traceFields;
        this.errorsCache = null;
    }

    /* Component k of the named value at sample i. */
    value(name, i, k = 0) {
        return this.trace[i * this.stride + this.fields[name].offset + k];
    }

    /* All components of the named value at sample i. */
    vector(name, i) {
        const start = i * this.stride + this.fields[name].offset;
        return Array.from(this.trace.subarray(start, start + this.fields[name].size));
    }

    /* The estimated attitude at sample i: w, x, y, z; IMU frame to world. */
    attitude(i) {
        return this.vector('attitude', i);
    }

    /* One component of a value for every sample, as its own array. */
    column(name, k = 0) {
        const out = new Float32Array(this.recording.samples);
        const offset = this.fields[name].offset + k;
        for (let i = 0; i < out.length; i++) {
            out[i] = this.trace[i * this.stride + offset];
        }
        return out;
    }

    /* Attitude errors against the reference, as the benchmark scores them. */
    get errors() {
        if (this.errorsCache === null) {
            const samples = this.recording.samples;
            const estimate = new Float64Array(4 * samples);
            const offset = this.fields.attitude.offset;
            for (let i = 0; i < samples; i++) {
                for (let k = 0; k < 4; k++) {
                    estimate[4 * i + k] = this.trace[i * this.stride + offset + k];
                }
            }
            this.errorsCache = attitudeErrors(this.recording.time, estimate, this.recording.reference.attitude,
                this.recording.reference.valid, this.valid);
        }
        return this.errorsCache;
    }
}

/* One filter step taken apart: named blocks of numbers (see explainer_api.cpp). */
export class Step {
    constructor(description, data) {
        this.data = data;
        this.fields = description.stepFields;
    }

    /* A single number. */
    number(name) {
        return this.data[this.fields[name].offset];
    }

    /* Column j of a block (for a block of sigma points, point j). */
    column(name, j = 0) {
        const field = this.fields[name];
        const start = field.offset + j * field.rows;
        return this.data.subarray(start, start + field.rows);
    }

    /* Entry (row, col) of a matrix block. */
    at(name, row, col) {
        const field = this.fields[name];
        return this.data[field.offset + col * field.rows + row];
    }
}

export class Estimator {
    /* Start an estimator. "url" is the address of ukf_explainer.wasm. */
    static async create(url) {
        const module = await compileModule(url);
        const estimator = new Estimator();
        try {
            estimator.worker = new Worker(new URL('./worker.js', import.meta.url), {type: 'module'});
            estimator.worker.onmessage = ({data}) => estimator.#settle(data);
            estimator.description = await estimator.#ask({type: 'init', module}, [], 15000);
        } catch (error) {
            console.warn('The estimator worker did not start; running on the main thread.', error);
            estimator.worker?.terminate();
            estimator.worker = null;
            estimator.local = await EstimatorModule.create(module);
            estimator.description = estimator.local.describe();
        }
        return estimator;
    }

    constructor() {
        this.worker = null;
        this.local = null;
        this.pending = new Map();
        this.nextId = 1;
        this.recording = null;
        this.queued = null;
        this.busy = false;
    }

    get constants() {
        return this.description.constants;
    }

    get defaultTuning() {
        return this.description.defaultTuning;
    }

    #ask(message, transfer = [], timeout = 0) {
        return new Promise((resolve, reject) => {
            const id = this.nextId++;
            this.pending.set(id, {resolve, reject});
            this.worker.onerror = event => {
                this.pending.delete(id);
                reject(new Error(event.message || 'The estimator worker failed'));
            };
            if (timeout > 0) {
                setTimeout(() => {
                    if (this.pending.delete(id)) {
                        reject(new Error('The estimator worker did not answer'));
                    }
                }, timeout);
            }
            this.worker.postMessage({...message, id}, transfer);
        });
    }

    #settle({id, result, error}) {
        const waiting = this.pending.get(id);
        if (!waiting) {
            return;
        }
        this.pending.delete(id);
        if (error) {
            waiting.reject(new Error(error));
        } else {
            waiting.resolve(result);
        }
    }

    /*
    Give the estimator a recording. "sensors" may replace the gyroscope and
    accelerometer arrays, to try what happens when a sensor goes wrong.
    */
    async load(recording, sensors = {}) {
        this.recording = recording;
        const gyroscope = sensors.gyroscope ?? recording.gyroscope;
        const accelerometer = sensors.accelerometer ?? recording.accelerometer;
        if (this.local) {
            this.local.load(recording.time, gyroscope, accelerometer);
        } else {
            await this.#ask({type: 'load', time: recording.time, gyroscope, accelerometer});
        }
    }

    /* Run the estimator. "tuning" holds only the values which differ from the benchmark's. */
    async run(tuning = {}) {
        const recording = this.recording;
        let result;
        if (this.local) {
            const started = performance.now();
            result = this.local.run(tuning);
            result.elapsed = performance.now() - started;
        } else {
            result = await this.#ask({type: 'run', tuning});
        }
        return new Run(this.description, recording, result, tuning);
    }

    /*
    Like run(), for a control which changes quickly: if a run is already
    under way, only the newest request is kept. Resolves to null for a
    request which was overtaken.
    */
    async runLatest(prepare) {
        if (this.busy) {
            this.queued?.resolve(null);
            return new Promise(resolve => {
                this.queued = {prepare, resolve};
            });
        }
        this.busy = true;
        try {
            return await prepare();
        } finally {
            this.busy = false;
            if (this.queued) {
                const next = this.queued;
                this.queued = null;
                this.runLatest(next.prepare).then(next.resolve);
            }
        }
    }

    /*
    Look ahead from sample "index" for "steps" samples without corrections.
    Returns {snapshots, points, size, data}: data holds, for each snapshot
    and each sigma point, a state of "size" numbers (attitude w x y z,
    velocity, gyro bias). Null if the index is out of range.
    */
    async predict(index, steps) {
        const data = this.local ? this.local.predict(index, steps) : await this.#ask({type: 'predict', index, steps});
        if (!data) {
            return null;
        }
        const points = this.constants.sigma_count, size = this.constants.state_size;
        return {snapshots: data.length / (points * size), points, size, data};
    }

    /* Take the step which ends at sample "index" of the last run apart. */
    async explain(index) {
        const data = this.local ? this.local.explain(index) : await this.#ask({type: 'explain', index});
        return data ? new Step(this.description, data) : null;
    }
}
