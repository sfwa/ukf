/*
The playground: run the estimator over a whole recording with your own
tuning, and compare the result with the benchmark's tuning.

Every change runs the real estimator again in a worker. The sensor faults
are made by changing the recorded gyroscope samples before they are given
to the estimator; the estimator itself is never changed.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addToggle, addChoice, addButton, addRow, addReadout} from '../core/controls.js';
import {Axes, text, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {clamp, DEG} from '../core/math.js';
import {stack, seconds} from './story-kit.js';
import {loadIndex, loadRecording, loadRun, createEstimator} from '../data/story.js';

/* The tuning values on show: the estimator's own names, with plain labels and units for display. */
const TUNING = [
    {name: 'gyroscope_noise', label: 'Gyroscope noise', min: 0.0003, max: 0.03, unit: '°/√s', scale: 1 / DEG, digits: 3},
    {name: 'rotor_model_noise', label: 'Rotor model noise', min: 0.02, max: 2, unit: 'm/s/√s', scale: 1, digits: 3},
    {name: 'rotor_drag_noise', label: 'Drag reading noise', min: 0.0067, max: 0.67, unit: 'm/s²·√s', scale: 1, digits: 3},
    {name: 'gyro_bias_walk', label: 'Bias drift', min: 0.0000175, max: 0.00175, unit: '°/s/√s', scale: 1 / DEG, digits: 4},
    {name: 'heading_sigma_limit', label: 'Heading uncertainty limit', min: 0.03, max: 3, unit: '°', scale: 1 / DEG, digits: 1},
];

const VIEWS = [
    {value: 'tilt', label: 'Tilt error'},
    {value: 'heading', label: 'Heading error'},
    {value: 'velocity', label: 'Velocity'},
    {value: 'bias', label: 'Gyroscope bias'},
];

const FAULT_AXIS = 2;           /* the gyroscope axis which the faults act on */
const KICK_AT = 0.35;           /* where in the recording the knock happens (fraction) */
const KICK_SECONDS = 0.1;
const JUMP_AT = 0.6;            /* where the bias jump happens (fraction) */
const VELOCITY_AXIS = 1;        /* the IMU axis plotted in the velocity view (sideways) */

class PlaygroundFigure extends Figure {
    setup() {
        this.recording = null;
        this.run = null;
        this.baseline = null;
        this.failure = null;
        this.loaded = null;
        this.working = false;
        this.params.t = 20;

        const top = addRow(this);
        this.recordingSlot = top;
        addChoice(this, 'view', {label: 'Plot', options: VIEWS, value: 'tilt', parent: addRow(this)});
        addToggle(this, 'use_rotor_model', {label: 'Use the rotor model', value: true});
        for (const item of TUNING) {
            addSlider(this, item.name, {
                label: item.label, min: item.min, max: item.max, value: Math.sqrt(item.min * item.max), log: true,
                format: v => `${(v * item.scale).toFixed(item.digits)} ${item.unit}`,
            });
        }
        addSlider(this, 'kick', {label: 'Knock the estimate by', min: 0, max: 40, value: 0, step: 1, format: v => `${v.toFixed(0)}°`});
        addSlider(this, 'jump', {label: 'Gyroscope bias jump', min: 0, max: 3, value: 0, step: 0.05, format: v => `${v.toFixed(2)} °/s`});
        const buttons = addRow(this);
        addButton(this, {label: "Back to the benchmark's tuning", parent: buttons, onClick: () => this.reset()});
        addButton(this, {label: 'Save this run as CSV', parent: buttons, onClick: () => this.save()});
        this.table = document.createElement('table');
        this.table.className = 'numbers';
        this.controls.append(this.table);
        this.status = addReadout(this);
        this.readout = addReadout(this);

        this.canvas.classList.add('scrub');
        const move = at => {
            if (this.axes && this.recording) {
                this.set({t: clamp(this.axes.fromX(at.x), 0, this.recording.duration)}, true);
            }
        };
        this.onPointer({down: move, drag: move});

        this.estimator = createEstimator();
        Promise.all([loadIndex(), this.estimator]).then(([list, estimator]) => {
            this.defaults = estimator.defaultTuning;
            addChoice(this, 'recording', {
                label: 'Recording',
                options: list.map(item => ({value: item.name, label: item.name.replace(/_/g, ' ')})),
                value: list[0].name,
                parent: this.recordingSlot,
            });
            this.reset();
        }).catch(error => this.fail(error));
    }

    fail(error) {
        console.error(error);
        this.failure = error.message;
        this.invalidate();
    }

    reset() {
        const values = {use_rotor_model: true, kick: 0, jump: 0};
        for (const item of TUNING) {
            values[item.name] = this.defaults[item.name];
        }
        this.set(values, true);
    }

    /* The tuning values which differ from the benchmark's. */
    get tuning() {
        const tuning = {};
        if (!this.params.use_rotor_model) {
            tuning.use_rotor_model = 0;
        }
        for (const item of TUNING) {
            if (Math.abs(this.params[item.name] / this.defaults[item.name] - 1) > 1e-9) {
                tuning[item.name] = this.params[item.name];
            }
        }
        return tuning;
    }

    changed() {
        /* While a value glides to its new setting, wait: one run at the end is enough. */
        if (this.defaults && this.params.recording && this.transitions.size === 0) {
            this.refresh();
        }
    }

    /* The recorded gyroscope samples with the chosen faults added. */
    faulty(recording) {
        const {kick, jump} = this.params;
        if (kick === 0 && jump === 0) {
            return {};
        }
        const gyroscope = recording.gyroscope.slice();
        const kickStart = recording.duration * KICK_AT, jumpStart = recording.duration * JUMP_AT;
        for (let i = 0; i < recording.samples; i++) {
            const t = recording.time[i];
            if (kick > 0 && t >= kickStart && t < kickStart + KICK_SECONDS) {
                gyroscope[3 * i + FAULT_AXIS] += kick * DEG / KICK_SECONDS;
            }
            if (jump > 0 && t >= jumpStart) {
                gyroscope[3 * i + FAULT_AXIS] += jump * DEG;
            }
        }
        return {gyroscope};
    }

    async refresh() {
        try {
            const name = this.params.recording;
            const estimator = await this.estimator;
            if (!this.recording || this.recording.name !== name) {
                this.status.textContent = 'Loading the recording…';
                this.baseline = null;
                const recording = await loadRecording(name);
                if (this.params.recording !== name) {
                    return;
                }
                this.recording = recording;
                this.set({t: clamp(this.params.t, 0, recording.duration)}, true);
                loadRun(name).then(result => {
                    if (this.params.recording === name) {
                        this.baseline = result.run;
                        this.fillTable();
                        this.invalidate();
                    }
                });
            }
            const recording = this.recording;
            const tuning = this.tuning;
            const faults = `${name}|${this.params.kick}|${this.params.jump}`;
            this.working = true;
            this.invalidate();
            const run = await estimator.runLatest(async () => {
                if (this.loaded !== faults) {
                    await estimator.load(recording, this.faulty(recording));
                    this.loaded = faults;
                }
                return estimator.run(tuning);
            });
            if (run && run.recording.name === this.params.recording) {
                this.run = run;
                this.working = false;
                this.fillTable();
                this.invalidate();
            }
        } catch (error) {
            this.fail(error);
        }
    }

    fillTable() {
        const run = this.run;
        if (!run) {
            return;
        }
        const number = value => (Number.isFinite(value) ? `${value.toFixed(3)}°` : 'none');
        const mine = run.errors.rms, base = this.baseline?.errors.rms;
        const rows = [
            ['', 'This run', 'Benchmark tuning'],
            ['Tilt error (RMS)', mine ? number(mine.inclination) : 'none', base ? number(base.inclination) : '…'],
            ['Heading error (RMS)', mine ? number(mine.heading) : 'none', base ? number(base.heading) : '…'],
        ];
        this.table.replaceChildren(...rows.map((cells, r) => {
            const row = document.createElement('tr');
            cells.forEach((value, k) => {
                const cell = document.createElement(r === 0 || k === 0 ? 'th' : 'td');
                cell.textContent = value;
                row.append(cell);
            });
            return row;
        }));
        const recording = run.recording;
        this.status.textContent = run.failed
            ? `The filter stopped at ${recording.time[Math.max(run.valid - 1, 0)].toFixed(1)} s: a step could not be completed with this tuning.`
            : `${recording.duration.toFixed(0)} s of flight, ${recording.samples} filter steps, run in ${(run.elapsed / 1000).toFixed(2)} s.`;
    }

    /* Write the run to a CSV file. */
    save() {
        const run = this.run;
        if (!run) {
            return;
        }
        const names = ['attitude', 'velocity', 'gyro_bias', 'sd', 'rotor_vibration', 'resting', 'rotor_model_active'];
        const header = ['time', 'qw', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz', 'bias_x', 'bias_y', 'bias_z',
            'sd_tilt_x', 'sd_tilt_y', 'sd_heading', 'sd_vx', 'sd_vy', 'sd_vz', 'sd_bias_x', 'sd_bias_y', 'sd_bias_z',
            'rotor_vibration', 'rest', 'rotor_model_active', 'tilt_error_deg', 'heading_error_deg'];
        const lines = [header.join(',')];
        const errors = run.errors;
        for (let i = 0; i < run.valid; i++) {
            const row = [run.recording.time[i]];
            for (const name of names) {
                row.push(...run.vector(name, i));
            }
            row.push(errors.inclination[i], errors.heading[i]);
            lines.push(row.join(','));
        }
        const link = document.createElement('a');
        link.href = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], {type: 'text/csv'}));
        link.download = `${run.recording.name}-estimate.csv`;
        link.click();
        setTimeout(() => URL.revokeObjectURL(link.href), 1000);
    }

    draw(ctx, width, height) {
        const c = this.colors;
        if (this.failure || !this.run) {
            text(ctx, this.failure ?? 'Loading the recording…', width / 2, height / 2,
                {align: 'center', baseline: 'middle', color: c.muted, size: 14});
            return;
        }
        const run = this.run, recording = run.recording, baseline = this.baseline;
        const time = recording.time, reference = recording.reference;
        const view = this.params.view;
        const [main, strip] = stack(width, height, [1, 0.07], {top: 56, gap: 26, bottom: 26});
        const end = run.valid;
        const sd = (i, k) => run.value('sd', i, k);
        const jump = i => (this.params.jump > 0 && time[i] >= recording.duration * JUMP_AT ? this.params.jump : 0);

        let axes, entries;
        const opacity = this.working ? 0.45 : 1;
        if (view === 'tilt') {
            /* Tall enough to show a knock, if one was asked for. */
            axes = new Axes(0, recording.duration, 0, Math.max(6, Math.ceil(this.params.kick * 1.15)));
            axes.place(main.left, main.top, main.width, main.height).grid(ctx, {yLabel: 'Tilt error (degrees)', xLabels: false});
            ctx.save();
            axes.clip(ctx);
            axes.band(ctx, time, () => 0, i => Math.hypot(sd(i, 0), sd(i, 1)) / DEG, {end, color: c.blue, opacity: 0.2 * opacity});
            if (baseline) {
                axes.line(ctx, time, baseline.errors.inclination, {color: alpha(c.muted, 0.9), width: 1.5});
            }
            axes.line(ctx, time, run.errors.inclination, {end, color: alpha(c.blue, opacity)});
            ctx.restore();
            entries = [{label: 'this run', color: c.blue}, {label: 'what the filter expects its error to be', color: c.blue, band: true},
                {label: 'benchmark tuning', color: c.muted}];
        } else if (view === 'heading') {
            axes = new Axes(0, recording.duration, -30, 30);
            axes.place(main.left, main.top, main.width, main.height).grid(ctx, {yLabel: 'Heading error (degrees)', xLabels: false});
            ctx.save();
            axes.clip(ctx);
            axes.band(ctx, time, i => -sd(i, 2) / DEG, i => sd(i, 2) / DEG, {end, color: c.blue, opacity: 0.2 * opacity});
            if (baseline) {
                axes.line(ctx, time, baseline.errors.heading, {color: alpha(c.muted, 0.9), width: 1.5});
            }
            axes.line(ctx, time, run.errors.heading, {end, color: alpha(c.blue, opacity)});
            ctx.restore();
            entries = [{label: 'this run', color: c.blue}, {label: "the filter's own uncertainty", color: c.blue, band: true},
                {label: 'benchmark tuning', color: c.muted}];
        } else if (view === 'velocity') {
            axes = new Axes(0, recording.duration, -3, 3);
            axes.place(main.left, main.top, main.width, main.height).grid(ctx, {yLabel: 'Sideways velocity (m/s)', xLabels: false});
            if (!this.bodyVelocity || this.bodyVelocity.name !== recording.name) {
                this.bodyVelocity = {name: recording.name, values: sidewaysVelocity(recording)};
            }
            const velocity = i => run.value('velocity', i, VELOCITY_AXIS);
            ctx.save();
            axes.clip(ctx);
            axes.band(ctx, time, i => velocity(i) - 2 * sd(i, 3 + VELOCITY_AXIS), i => velocity(i) + 2 * sd(i, 3 + VELOCITY_AXIS),
                {end, color: c.blue, opacity: 0.2 * opacity});
            axes.line(ctx, time, this.bodyVelocity.values, {color: c.ink});
            axes.line(ctx, time, run.trace, {end, color: alpha(c.blue, opacity), stride: run.stride, offset: run.fields.velocity.offset + VELOCITY_AXIS});
            ctx.restore();
            entries = [{label: 'this run', color: c.blue}, {label: "the filter's own uncertainty", color: c.blue, band: true},
                {label: 'reference', color: c.ink}];
        } else {
            axes = new Axes(0, recording.duration, -1, 8);
            axes.place(main.left, main.top, main.width, main.height).grid(ctx, {yLabel: 'Gyroscope bias, one axis (degrees per second)', xLabels: false});
            const bias = i => run.value('gyro_bias', i, FAULT_AXIS) / DEG;
            ctx.save();
            axes.clip(ctx);
            axes.band(ctx, time, i => bias(i) - 2 * sd(i, 6 + FAULT_AXIS) / DEG, i => bias(i) + 2 * sd(i, 6 + FAULT_AXIS) / DEG,
                {end, color: c.blue, opacity: 0.2 * opacity});
            /* The true bias is the reference's value, plus the jump if one was added. */
            ctx.strokeStyle = c.ink;
            ctx.lineWidth = 2;
            ctx.beginPath();
            let pen = false;
            for (let i = 0; i < recording.samples; i += 20) {
                if (!reference.valid[i]) {
                    pen = false;
                    continue;
                }
                const px = axes.x(time[i]), py = axes.y(reference.gyroBias[3 * i + FAULT_AXIS] / DEG + jump(i));
                if (pen) {
                    ctx.lineTo(px, py);
                } else {
                    ctx.moveTo(px, py);
                    pen = true;
                }
            }
            ctx.stroke();
            axes.line(ctx, time, run.trace, {end, color: alpha(c.blue, opacity), stride: run.stride,
                offset: run.fields.gyro_bias.offset + FAULT_AXIS, scale: 1 / DEG});
            ctx.restore();
            entries = [{label: 'this run', color: c.blue}, {label: "the filter's own uncertainty", color: c.blue, band: true},
                {label: 'reference', color: c.ink}];
        }
        this.axes = axes;
        legend(ctx, entries, axes.left, 11);

        /* Marks for the faults. */
        const marks = [];
        if (this.params.kick > 0) {
            marks.push({t: recording.duration * KICK_AT, label: 'knock'});
        }
        if (this.params.jump > 0) {
            marks.push({t: recording.duration * JUMP_AT, label: 'bias jump'});
        }
        for (const mark of marks) {
            const px = Math.round(axes.x(mark.t)) + 0.5;
            ctx.strokeStyle = c.orange;
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.moveTo(px, axes.top);
            ctx.lineTo(px, axes.top + axes.height);
            ctx.stroke();
            text(ctx, mark.label, px + 5, axes.top + 4, {baseline: 'top', color: c.ink2, size: 11.5, halo: true});
        }

        /* Below: when the rotor model was in use. */
        const stripAxes = new Axes(0, recording.duration, 0, 1).place(strip.left, strip.top, strip.width, strip.height);
        ctx.fillStyle = c.surface2;
        ctx.fillRect(strip.left, strip.top, strip.width, strip.height);
        stripAxes.spans(ctx, time, i => run.value('rotor_model_active', i) > 0, {end, color: c.ink2, opacity: 0.4});
        stripAxes.grid(ctx, {yTicks: 0, xGrid: false, xFormat: seconds});
        text(ctx, 'rotor model in use', strip.left, strip.top - 5, {baseline: 'bottom', color: c.muted, size: 11});

        /* The cursor. */
        const t = clamp(this.params.t, 0, recording.duration);
        const i = clamp(Math.round(t / recording.duration * (recording.samples - 1)), 0, Math.max(end - 1, 0));
        axes.cursor(ctx, t);
        const parts = [`${t.toFixed(1)} s`];
        if (view === 'tilt' && Number.isFinite(run.errors.inclination[i])) {
            parts.push(`tilt error ${run.errors.inclination[i].toFixed(2)}°`, `the filter expects about ${(Math.hypot(sd(i, 0), sd(i, 1)) / DEG).toFixed(2)}°`);
        } else if (view === 'heading' && Number.isFinite(run.errors.heading[i])) {
            parts.push(`heading error ${run.errors.heading[i].toFixed(1)}°`, `uncertainty ${(sd(i, 2) / DEG).toFixed(1)}°`);
        } else if (view === 'velocity') {
            parts.push(`filter ${run.value('velocity', i, VELOCITY_AXIS).toFixed(2)} m/s`, `reference ${this.bodyVelocity.values[i].toFixed(2)} m/s`);
        } else if (view === 'bias') {
            parts.push(`filter ${(run.value('gyro_bias', i, FAULT_AXIS) / DEG).toFixed(2)} °/s`);
        }
        this.readout.textContent = parts.join('  ·  ');
    }
}

/* The reference velocity along one IMU axis, for every sample (NaN without a reference). */
function sidewaysVelocity(recording) {
    const reference = recording.reference;
    const out = new Float32Array(recording.samples).fill(NaN);
    for (let i = 0; i < recording.samples; i++) {
        if (reference.valid[i]) {
            /* Turn the world-frame velocity into the IMU frame: v_body = conj(q) v q. */
            const w = reference.attitude[4 * i], x = -reference.attitude[4 * i + 1], y = -reference.attitude[4 * i + 2], z = -reference.attitude[4 * i + 3];
            const vx = reference.velocity[3 * i], vy = reference.velocity[3 * i + 1], vz = reference.velocity[3 * i + 2];
            const tx = 2 * (y * vz - z * vy), ty = 2 * (z * vx - x * vz), tz = 2 * (x * vy - y * vx);
            const body = [vx + w * tx + (y * tz - z * ty), vy + w * ty + (z * tx - x * tz), vz + w * tz + (x * ty - y * tx)];
            out[i] = body[VELOCITY_AXIS];
        }
    }
    return out;
}

registerFigure('playground', PlaygroundFigure);
