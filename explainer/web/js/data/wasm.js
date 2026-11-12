/*
The estimator module: a thin wrapper over the C interface declared in
explainer/wasm/explainer_api.h. It runs in a worker (see worker.js) and in
Node (see explainer/test), and has no dependency on the page.

The module describes its own result layouts, so nothing here knows an offset.
*/

const SAMPLE_STRIDE = 7;
const ABI_VERSION = 1;

export class EstimatorModule {
    /* "source" is the bytes of ukf_explainer.wasm, or a compiled WebAssembly.Module. */
    static async create(source) {
        const module = source instanceof WebAssembly.Module ? source : await WebAssembly.compile(source);
        const instance = await WebAssembly.instantiate(module, {
            env: {emscripten_notify_memory_growth() {}},
        });
        return new EstimatorModule(instance);
    }

    constructor(instance) {
        this.api = instance.exports;
        /* Runs the C++ static initialisers. */
        this.api._initialize();
        if (this.api.ex_abi_version() !== ABI_VERSION) {
            throw new Error('The estimator module and the page are different versions');
        }

        this.constants = {};
        for (let i = 0; i < this.api.ex_constant_count(); i++) {
            this.constants[this.#string(this.api.ex_constant_name(i))] = this.api.ex_constant(i);
        }

        this.tuningNames = [];
        for (let i = 0; i < this.api.ex_tuning_count(); i++) {
            this.tuningNames.push(this.#string(this.api.ex_tuning_name(i)));
        }
        this.api.ex_tuning_reset();
        this.defaultTuning = this.#readTuning();

        this.traceStride = this.api.ex_trace_stride();
        this.traceFields = {};
        for (let i = 0; i < this.api.ex_trace_field_count(); i++) {
            this.traceFields[this.#string(this.api.ex_trace_field_name(i))] = {
                offset: this.api.ex_trace_field_offset(i),
                size: this.api.ex_trace_field_size(i),
            };
        }

        this.stepSize = this.api.ex_step_size();
        this.stepFields = {};
        for (let i = 0; i < this.api.ex_step_field_count(); i++) {
            this.stepFields[this.#string(this.api.ex_step_field_name(i))] = {
                offset: this.api.ex_step_field_offset(i),
                rows: this.api.ex_step_field_rows(i),
                cols: this.api.ex_step_field_cols(i),
            };
        }
        this.samples = 0;
    }

    /* Everything the page needs to know about the module, as plain data. */
    describe() {
        return {
            constants: this.constants,
            tuningNames: this.tuningNames,
            defaultTuning: this.defaultTuning,
            traceStride: this.traceStride,
            traceFields: this.traceFields,
            stepSize: this.stepSize,
            stepFields: this.stepFields,
        };
    }

    #string(pointer) {
        const bytes = new Uint8Array(this.api.memory.buffer);
        let end = pointer;
        while (bytes[end] !== 0) {
            end++;
        }
        return new TextDecoder().decode(bytes.subarray(pointer, end));
    }

    #readTuning() {
        const values = new Float64Array(this.api.memory.buffer, this.api.ex_tuning(), this.tuningNames.length);
        return Object.fromEntries(this.tuningNames.map((name, i) => [name, values[i]]));
    }

    /*
    Give the module a recording: time (n), gyroscope (3n) and accelerometer
    (3n). The values are copied.
    */
    load(time, gyroscope, accelerometer) {
        const count = time.length;
        const pointer = this.api.ex_samples(count);
        if (pointer === 0) {
            throw new Error('A recording needs at least two samples');
        }
        const values = new Float64Array(this.api.memory.buffer, pointer, count * SAMPLE_STRIDE);
        for (let i = 0; i < count; i++) {
            values[SAMPLE_STRIDE * i] = time[i];
            for (let axis = 0; axis < 3; axis++) {
                values[SAMPLE_STRIDE * i + 1 + axis] = gyroscope[3 * i + axis];
                values[SAMPLE_STRIDE * i + 4 + axis] = accelerometer[3 * i + axis];
            }
        }
        this.samples = count;
    }

    /*
    Run the estimator over the loaded recording. "tuning" holds the values
    to change; every other value is the benchmark's. Returns a copy of the
    trace, and how many samples have an estimate.
    */
    run(tuning = {}) {
        this.api.ex_tuning_reset();
        const values = new Float64Array(this.api.memory.buffer, this.api.ex_tuning(), this.tuningNames.length);
        for (const [name, value] of Object.entries(tuning)) {
            const index = this.tuningNames.indexOf(name);
            if (index < 0) {
                throw new Error(`The estimator has no tuning value called ${name}`);
            }
            values[index] = value;
        }
        const valid = this.api.ex_run();
        if (valid < 0) {
            throw new Error('There is no recording to run');
        }
        const trace = new Float64Array(this.api.memory.buffer, this.api.ex_trace(), this.samples * this.traceStride);
        return {valid, samples: this.samples, trace: trace.slice()};
    }

    /*
    Move the sigma points of the estimate at sample "index" through the
    process model for "steps" samples, with no corrections. Returns the
    snapshots one after another: (steps + 1) x sigma points x state size.
    */
    predict(index, steps) {
        const snapshots = this.api.ex_predict(index, steps);
        if (snapshots === 0) {
            return null;
        }
        const size = snapshots * this.constants.sigma_count * this.constants.state_size;
        return new Float64Array(this.api.memory.buffer, this.api.ex_prediction(), size).slice();
    }

    /* Take the step which ends at sample "index" of the last run apart. */
    explain(index) {
        const code = this.api.ex_explain(index);
        if (code !== 0) {
            return null;
        }
        return new Float64Array(this.api.memory.buffer, this.api.ex_step(), this.stepSize).slice();
    }
}
