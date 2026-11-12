/*
Looking ahead with the real sigma points. The figure takes the filter's
19 sigma points at one moment of the flight and moves each of them
through the process model for up to a second, with no corrections (the
estimator module does this: see ex_predict in explainer/wasm).

Left: how each guess is tilted, compared with the middle guess. Right: the
velocity of each guess across the rotor plane, and the true velocity from
the reference. Guesses which start with a different tilt end with a
different velocity: that is the link which lets a velocity measurement
correct the tilt.
*/

import {registerFigure} from '../core/figure.js';
import {addSlider, addChoice, addReadout} from '../core/controls.js';
import {Axes, text, marker, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {clamp, qMul, qConj, qRotate, qAt, vAt, dot, niceStep} from '../core/math.js';
import {StoryFigure} from './story-kit.js';

const STEPS = 200;              /* one second of samples */
const DEGREES = 180 / Math.PI;

/* The sigma points which differ from the middle one in each kind of state. */
const GROUPS = {
    tilt: [0, 1],
    heading: [2],
    velocity: [3, 4, 5],
    bias: [6, 7, 8],
};

class FanFigure extends StoryFigure {
    build() {
        addSlider(this, 't0', {label: 'Moment in the flight', min: 5, max: 80, value: 30, step: 0.1, format: v => `${v.toFixed(1)} s`});
        addSlider(this, 'ahead', {label: 'Look ahead', min: 0, max: 1, value: 0.6, step: 0.005, format: v => `${v.toFixed(2)} s`});
        addChoice(this, 'group', {
            label: 'Sigma points that differ in',
            options: [
                {value: 'tilt', label: 'tilt'},
                {value: 'heading', label: 'heading'},
                {value: 'velocity', label: 'velocity'},
                {value: 'bias', label: 'gyroscope bias'},
            ],
            value: 'tilt',
        });
        this.readout = addReadout(this);
        this.prediction = null;
        this.predictedIndex = -1;
        this.wanted = -1;
        this.asking = false;
        this.near = null;
        this.onPointer({
            move: at => this.hover(at),
            leave: () => {
                this.near = null;
                this.invalidate();
            },
        });
    }

    ready({recording, run, estimator}) {
        const constants = estimator.constants;
        this.lateral = [
            [constants.lateral_axis_1_x, constants.lateral_axis_1_y, constants.lateral_axis_1_z],
            [constants.lateral_axis_2_x, constants.lateral_axis_2_y, constants.lateral_axis_2_z],
        ];
        /* Keep the moment inside powered flight, with a second of recording after it. */
        let first = 0;
        while (first < run.valid && !run.value('rotor_model_active', first)) {
            first++;
        }
        const low = Math.ceil(recording.time[Math.min(first + 200, run.valid - 1)]);
        const high = Math.floor(recording.duration - 1.5);
        this.inputs.t0.input.min = low;
        this.inputs.t0.input.max = high;
        this.set({t0: clamp(this.params.t0, low, high)}, true);
    }

    changed() {
        if (!this.story) {
            return;
        }
        this.wanted = this.sampleAt(this.params.t0);
        this.ask();
    }

    /* Ask the estimator for the look-ahead from the wanted sample; keep only the newest request. */
    async ask() {
        if (this.asking || this.wanted === this.predictedIndex) {
            return;
        }
        this.asking = true;
        const index = this.wanted;
        const prediction = await this.story.estimator.predict(index, STEPS);
        this.asking = false;
        if (prediction) {
            this.prediction = prediction;
            this.predictedIndex = index;
            this.invalidate();
        }
        if (this.wanted !== index) {
            this.ask();
        }
    }

    /* The state of sigma point k at a snapshot: attitude and velocity. */
    point(snapshot, k) {
        const {data, points, size} = this.prediction;
        const start = (snapshot * points + k) * size;
        return {attitude: Array.from(data.subarray(start, start + 4)), velocity: Array.from(data.subarray(start + 4, start + 7))};
    }

    /* Where a guess sits in each panel: its tilt against the middle guess (degrees), and its velocity (m/s). */
    place(snapshot, k) {
        const guess = this.point(snapshot, k), middle = this.point(snapshot, 0);
        const difference = qMul(guess.attitude, qConj(middle.attitude));
        return {
            tilt: [2 * difference[1] / difference[0] * DEGREES, 2 * difference[2] / difference[0] * DEGREES],
            velocity: [dot(guess.velocity, this.lateral[0]), dot(guess.velocity, this.lateral[1])],
        };
    }

    inGroup(k) {
        if (k === 0) {
            return true;
        }
        const column = (k - 1) % 9;
        return GROUPS[this.params.group].includes(column);
    }

    hover(at) {
        if (!this.layout || !this.prediction) {
            return;
        }
        let best = null, distance = 22;
        const snapshot = this.snapshot;
        for (let k = 0; k < this.prediction.points; k++) {
            if (!this.inGroup(k)) {
                continue;
            }
            const where = this.place(snapshot, k);
            for (const [axes, value] of [[this.layout.tilt, where.tilt], [this.layout.velocity, where.velocity]]) {
                const d = Math.hypot(axes.x(value[0]) - at.x, axes.y(value[1]) - at.y);
                if (d < distance) {
                    distance = d;
                    best = k;
                }
            }
        }
        if (best !== this.near) {
            this.near = best;
            this.invalidate();
        }
    }

    get snapshot() {
        return clamp(Math.round(this.params.ahead / 0.005), 0, this.prediction.snapshots - 1);
    }

    paint(ctx, width, height) {
        const c = this.colors;
        if (!this.prediction) {
            text(ctx, 'Working…', width / 2, height / 2, {align: 'center', baseline: 'middle', color: c.muted, size: 14});
            return;
        }
        const {recording} = this.story;
        const reference = recording.reference;
        const index = this.predictedIndex;
        const snapshots = this.prediction.snapshots;
        const snapshot = this.snapshot;

        /* Two square panels: side by side when there is room, else one above the other. */
        const wide = width >= 560;
        const gap = 46, top = 34, bottom = 44, left = 44, right = 12;
        const side = wide ? Math.min((width - left - right - gap - 34) / 2, height - top - bottom)
            : Math.min(width - left - right, (height - 2 * top - bottom - 10) / 2);
        const tiltBox = {left, top};
        const velocityBox = wide ? {left: left + side + gap + 34, top} : {left, top: top + side + top + 10};

        /* Ranges which hold every guess for the whole second, so the panels do not jump. */
        let tiltReach = 0.2;
        let low = [Infinity, Infinity], high = [-Infinity, -Infinity];
        const grow = value => {
            for (let axis = 0; axis < 2; axis++) {
                low[axis] = Math.min(low[axis], value[axis]);
                high[axis] = Math.max(high[axis], value[axis]);
            }
        };
        for (let s = 0; s < snapshots; s += 10) {
            for (let k = 0; k < this.prediction.points; k++) {
                const where = this.place(s, k);
                tiltReach = Math.max(tiltReach, Math.abs(where.tilt[0]), Math.abs(where.tilt[1]));
                grow(where.velocity);
            }
        }
        const truth = [];
        for (let s = 0; s < snapshots; s++) {
            const i = index + s;
            if (reference.valid[i]) {
                const body = qRotate(qConj(qAt(reference.attitude, i)), vAt(reference.velocity, i));
                truth.push([dot(body, this.lateral[0]), dot(body, this.lateral[1])]);
                grow(truth[truth.length - 1]);
            }
        }
        const span = Math.max(high[0] - low[0], high[1] - low[1], 0.5) * 1.12;
        const centre = [(low[0] + high[0]) / 2, (low[1] + high[1]) / 2];
        const tiltRange = Math.ceil(tiltReach * 1.2 / niceStep(tiltReach, 2)) * niceStep(tiltReach, 2);
        const tiltAxes = new Axes(-tiltRange, tiltRange, -tiltRange, tiltRange).place(tiltBox.left, tiltBox.top, side, side);
        const velocityAxes = new Axes(centre[0] - span / 2, centre[0] + span / 2, centre[1] - span / 2, centre[1] + span / 2)
            .place(velocityBox.left, velocityBox.top, side, side);
        this.layout = {tilt: tiltAxes, velocity: velocityAxes};

        tiltAxes.grid(ctx, {xTicks: 4, yTicks: 4, yLabel: 'Tilt, against the middle sigma point (degrees)'});
        velocityAxes.grid(ctx, {xTicks: 4, yTicks: 4, yLabel: 'Sideways velocity (m/s)', zero: false});

        /* Right: the true velocity over the look-ahead. */
        ctx.save();
        velocityAxes.clip(ctx);
        const shown = Math.min(truth.length, snapshot + 1);
        if (shown > 0) {
            ctx.strokeStyle = alpha(c.ink, this.emphasis('truth'));
            ctx.lineWidth = 2;
            ctx.lineJoin = 'round';
            ctx.beginPath();
            for (let s = 0; s < shown; s++) {
                const px = velocityAxes.x(truth[s][0]), py = velocityAxes.y(truth[s][1]);
                if (s === 0) {
                    ctx.moveTo(px, py);
                } else {
                    ctx.lineTo(px, py);
                }
            }
            ctx.stroke();
        }

        /* The path of each guess, faint for the ones outside the chosen group. */
        for (let k = 0; k < this.prediction.points; k++) {
            const chosen = this.inGroup(k);
            ctx.strokeStyle = alpha(c.green, (chosen ? 0.75 : 0.13) * this.emphasis('guess'));
            ctx.lineWidth = chosen ? 1.5 : 1;
            ctx.beginPath();
            for (let s = 0; s <= snapshot; s += 2) {
                const where = this.place(s, k).velocity;
                if (s === 0) {
                    ctx.moveTo(velocityAxes.x(where[0]), velocityAxes.y(where[1]));
                } else {
                    ctx.lineTo(velocityAxes.x(where[0]), velocityAxes.y(where[1]));
                }
            }
            ctx.stroke();
        }
        ctx.restore();

        /* The guesses now, in both panels. */
        const dots = [];
        for (let k = 0; k < this.prediction.points; k++) {
            dots.push({k, chosen: this.inGroup(k), now: this.place(snapshot, k), start: this.place(0, k)});
        }
        dots.sort((a, b) => Number(a.chosen) - Number(b.chosen));
        for (const [axes, name] of [[tiltAxes, 'tilt'], [velocityAxes, 'velocity']]) {
            ctx.save();
            axes.clip(ctx);
            for (const item of dots) {
                const px = axes.x(item.now[name][0]), py = axes.y(item.now[name][1]);
                if (name === 'velocity' && item.chosen && snapshot > 0) {
                    /* Where this guess started, as a ring in the filter's own colour. */
                    ctx.strokeStyle = alpha(c.blue, 0.9 * this.emphasis('guess'));
                    ctx.lineWidth = 1.5;
                    ctx.beginPath();
                    ctx.arc(axes.x(item.start[name][0]), axes.y(item.start[name][1]), 4, 0, 2 * Math.PI);
                    ctx.stroke();
                }
                const color = name === 'tilt' || snapshot === 0 ? c.blue : c.green;
                if (item.chosen) {
                    marker(ctx, px, py, alpha(color, this.emphasis('guess')), item.k === 0 ? 5.5 : 4.5, item.k === this.near ? 4 : 2);
                } else {
                    marker(ctx, px, py, alpha(color, 0.22), 3, 0);
                }
            }
            ctx.restore();
        }
        if (shown > 0) {
            const end = truth[shown - 1];
            marker(ctx, velocityAxes.x(end[0]), velocityAxes.y(end[1]), alpha(c.ink, this.emphasis('truth')), 4.5);
            text(ctx, 'true', velocityAxes.x(end[0]) + 9, velocityAxes.y(end[1]) - 7, {color: c.ink2, halo: true, size: 12});
        }
        if (this.near !== null) {
            for (const [axes, name] of [[tiltAxes, 'tilt'], [velocityAxes, 'velocity']]) {
                const where = this.place(snapshot, this.near)[name];
                ctx.strokeStyle = c.ink;
                ctx.lineWidth = 1.5;
                ctx.beginPath();
                ctx.arc(axes.x(where[0]), axes.y(where[1]), 9, 0, 2 * Math.PI);
                ctx.stroke();
            }
        }

        legend(ctx, [
            {label: 'now', color: c.blue, ring: true},
            {label: 'later', color: c.green, dot: true},
            {label: 'true', color: c.ink},
        ], velocityAxes.left, velocityAxes.top + side + 32, {size: 11.5});

        /* How far the chosen guesses are apart in velocity, now and at the start. */
        const spreadOf = s => {
            let reach = 0;
            const middle = this.place(s, 0).velocity;
            for (let k = 1; k < this.prediction.points; k++) {
                if (this.inGroup(k)) {
                    const v = this.place(s, k).velocity;
                    reach = Math.max(reach, Math.hypot(v[0] - middle[0], v[1] - middle[1]));
                }
            }
            return reach;
        };
        if (this.near !== null) {
            const where = this.place(snapshot, this.near);
            this.readout.textContent = this.near === 0
                ? `The middle sigma point: velocity ${where.velocity[0].toFixed(2)}, ${where.velocity[1].toFixed(2)} m/s.`
                : `This sigma point: tilted ${Math.hypot(...where.tilt).toFixed(2)}° from the middle one, velocity ${where.velocity[0].toFixed(2)}, ${where.velocity[1].toFixed(2)} m/s.`;
        } else {
            this.readout.textContent = `Farthest of these sigma points from the middle one: ${spreadOf(0).toFixed(2)} m/s now, `
                + `${spreadOf(snapshot).toFixed(2)} m/s after ${this.params.ahead.toFixed(2)} s.`;
        }
    }
}

registerFigure('fan', FanFigure);
