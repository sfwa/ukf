/*
Which measurements the filter uses, and when. This recording starts with
the aircraft carried by hand, then standing on the ground while its rotors
are tried, then flying.

Top: the vibration level, which tells the filter whether the rotors are
turning. Middle: the measurements in use. Bottom: the speed the filter
believes, against the reference.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addReadout} from '../core/controls.js';
import {Axes, text, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {clamp} from '../core/math.js';
import {stack, seconds} from './story-kit.js';
import {loadRun} from '../data/story.js';

const RECORDING = 'MH_01_easy';
const WINDOW = 60;
const ROWS = [
    {key: 'rest', label: 'gyroscope reading at rest'},
    {key: 'zero', label: 'velocity is about zero'},
    {key: 'drag', label: 'rotor drag'},
];

class ModesFigure extends Figure {
    setup() {
        this.params.t = 32;
        this.data = null;
        this.failure = null;
        this.readout = addReadout(this);
        this.canvas.classList.add('scrub');
        const move = at => {
            if (this.axes) {
                this.set({t: clamp(this.axes.fromX(at.x), 0, WINDOW)}, true);
            }
        };
        this.onPointer({down: move, drag: move});
        loadRun(RECORDING).then(result => this.ready(result)).catch(error => {
            this.failure = error.message;
            this.invalidate();
        });
    }

    ready({recording, run, tuning}) {
        const count = Math.min(recording.samples, Math.round(WINDOW / recording.duration * (recording.samples - 1)) + 1);
        const data = {
            recording, run, tuning, count,
            vibration: new Float32Array(count),
            speed: new Float32Array(count),
            referenceSpeed: new Float32Array(count).fill(NaN),
        };
        for (let i = 0; i < count; i++) {
            data.vibration[i] = run.value('rotor_vibration', i);
            data.speed[i] = Math.hypot(run.value('velocity', i, 0), run.value('velocity', i, 1), run.value('velocity', i, 2));
            if (recording.reference.valid[i]) {
                const v = recording.reference.velocity;
                data.referenceSpeed[i] = Math.hypot(v[3 * i], v[3 * i + 1], v[3 * i + 2]);
            }
        }
        this.data = data;
        this.invalidate();
    }

    /* Which measurements the filter used at sample i. */
    used(i) {
        const {run} = this.data;
        const powered = run.value('rotor_model_active', i) > 0;
        return {rest: run.value('resting', i) > 0, zero: !powered, drag: powered};
    }

    draw(ctx, width, height) {
        const c = this.colors;
        if (!this.data) {
            text(ctx, this.failure ?? 'Loading the recording…', width / 2, height / 2,
                {align: 'center', baseline: 'middle', color: c.muted, size: 14});
            return;
        }
        const {recording, tuning, count, vibration, speed, referenceSpeed} = this.data;
        const time = recording.time;
        const [top, middle, bottom] = stack(width, height, [1.15, 0.62, 1], {top: 26, gap: 36, bottom: 24});
        const vibrationAxes = new Axes(0, WINDOW, 0, 2.5).place(top.left, top.top, top.width, top.height);
        const speedAxes = new Axes(0, WINDOW, 0, 3).place(bottom.left, bottom.top, bottom.width, bottom.height);
        this.axes = vibrationAxes;

        /* Top: vibration, with the two levels which switch the rotor state. */
        vibrationAxes.grid(ctx, {yLabel: 'Vibration of the accelerometer (m/s²)', xLabels: false, yTicks: 3});
        ctx.save();
        vibrationAxes.clip(ctx);
        for (const level of [tuning.rotors_on_vibration, tuning.rotors_off_vibration]) {
            ctx.strokeStyle = c.ink2;
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.moveTo(vibrationAxes.left, Math.round(vibrationAxes.y(level)) + 0.5);
            ctx.lineTo(vibrationAxes.left + vibrationAxes.width, Math.round(vibrationAxes.y(level)) + 0.5);
            ctx.stroke();
        }
        vibrationAxes.line(ctx, time, vibration, {end: count, color: alpha(c.orange, this.emphasis('vibration'))});
        ctx.restore();
        text(ctx, 'rotors on, above this', vibrationAxes.left + 8, vibrationAxes.y(tuning.rotors_on_vibration) - 5,
            {baseline: 'bottom', color: c.ink2, size: 11.5, halo: true});
        text(ctx, 'rotors off, below this', vibrationAxes.left + 8, vibrationAxes.y(tuning.rotors_off_vibration) + 5,
            {baseline: 'top', color: c.ink2, size: 11.5, halo: true});

        /* Middle: one row per kind of measurement, filled while it is in use. */
        text(ctx, 'Measurements in use', middle.left, middle.top - 9, {baseline: 'bottom', color: c.ink2, size: 12, weight: 500});
        const rowHeight = middle.height / ROWS.length;
        const rowAxes = new Axes(0, WINDOW, 0, 1);
        ROWS.forEach((row, k) => {
            const y = middle.top + k * rowHeight;
            rowAxes.place(middle.left, y + 3, middle.width, rowHeight - 6);
            ctx.fillStyle = c.surface2;
            ctx.fillRect(middle.left, y + 3, middle.width, rowHeight - 6);
            rowAxes.spans(ctx, time, i => this.used(i)[row.key],
                {end: count, color: c.ink2, opacity: 0.38 * this.emphasis(row.key)});
            text(ctx, row.label, middle.left + middle.width - 6, y + rowHeight / 2,
                {align: 'right', baseline: 'middle', color: c.ink, size: 11.5, weight: 500});
        });

        /* Bottom: the speed the filter believes. */
        speedAxes.grid(ctx, {yLabel: 'Speed (m/s)', xFormat: seconds, yTicks: 3});
        ctx.save();
        speedAxes.clip(ctx);
        speedAxes.line(ctx, time, referenceSpeed, {end: count, color: alpha(c.ink, this.emphasis('reference'))});
        speedAxes.line(ctx, time, speed, {end: count, color: alpha(c.blue, this.emphasis('estimate'))});
        ctx.restore();
        legend(ctx, [{label: 'filter', color: c.blue}, {label: 'reference', color: c.ink}],
            speedAxes.left + 150, speedAxes.top - 14);

        /* The time cursor. */
        const t = this.params.t;
        const i = clamp(Math.round(t / recording.duration * (recording.samples - 1)), 0, count - 1);
        for (const axes of [vibrationAxes, speedAxes]) {
            axes.cursor(ctx, t);
        }
        ctx.strokeStyle = c.ink2;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(Math.round(vibrationAxes.x(t)) + 0.5, middle.top);
        ctx.lineTo(Math.round(vibrationAxes.x(t)) + 0.5, middle.top + middle.height);
        ctx.stroke();
        vibrationAxes.dot(ctx, t, Math.min(vibration[i], 2.5), c.orange);
        speedAxes.dot(ctx, t, Math.min(speed[i], 3), c.blue);

        const used = this.used(i);
        const names = ROWS.filter(row => used[row.key]).map(row => row.label);
        const parts = [`${t.toFixed(1)} s`, `vibration ${vibration[i].toFixed(2)} m/s²`, `in use: ${names.join(' and ')}`,
            `filter speed ${speed[i].toFixed(2)} m/s`];
        if (Number.isFinite(referenceSpeed[i])) {
            parts.push(`reference ${referenceSpeed[i].toFixed(2)} m/s`);
        }
        this.readout.textContent = parts.join('  ·  ');
    }
}

registerFigure('modes', ModesFigure);
