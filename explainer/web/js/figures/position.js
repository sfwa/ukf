/*
Where is the aircraft? The estimator does not say: it estimates attitude,
velocity and gyroscope bias, not position. This figure adds up velocity
over time to get a position (dead reckoning), in two ways, and compares
both with the reference path:

    blue    the filter's velocity, added up
    orange  the raw IMU readings, added up twice with no filter at all

Both start from the true position, velocity and attitude at the chosen
moment. Seen from above.
*/

import {registerFigure} from '../core/figure.js';
import {addSlider, addReadout} from '../core/controls.js';
import {Axes, text, legend, marker} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {clamp, qMul, qRotate, qAt, vAt, qFromRotationVector, add, scale, niceStep} from '../core/math.js';
import {StoryFigure, seconds} from './story-kit.js';

const LONGEST = 20;     /* seconds of dead reckoning at most */
const GRAVITY = [0, 0, -9.80665];
const ERROR_TOP = 10;   /* metres: the top of the error plot */

/*
Both dead-reckoned paths from sample "start", for LONGEST seconds.
"fix" is the fixed turn about the vertical which lines the estimate's
heading up with the reference (see headingFix below).
*/
export function reckon({recording, run}, fix, start) {
    const reference = recording.reference;
    const count = Math.min(Math.round(LONGEST / 0.005) + 1, run.valid - start);
    const paths = {
        start, count,
        truth: new Float32Array(2 * count), filter: new Float32Array(2 * count), raw: new Float32Array(2 * count),
        filterError: new Float32Array(count), rawError: new Float32Array(count),
    };
    let filterPosition = vAt(reference.position, start);
    let rawPosition = filterPosition.slice(), rawVelocity = vAt(reference.velocity, start), rawAttitude = qAt(reference.attitude, start);
    for (let k = 0; k < count; k++) {
        const i = start + k;
        if (k > 0) {
            const dt = recording.time[i] - recording.time[i - 1];
            /* The filter's velocity is in the IMU frame; turn it into the world frame and add it up. */
            const attitude = qMul(fix, run.attitude(i));
            filterPosition = add(filterPosition, scale(qRotate(attitude, run.vector('velocity', i)), dt));
            /* The raw IMU: add up the turn rate, then the acceleration, then the velocity. */
            const rate = scale(add(vAt(recording.gyroscope, i - 1), vAt(recording.gyroscope, i)), 0.5);
            rawAttitude = qMul(rawAttitude, qFromRotationVector(scale(rate, dt)));
            const acceleration = add(qRotate(rawAttitude, vAt(recording.accelerometer, i)), GRAVITY);
            rawPosition = add(rawPosition, add(scale(rawVelocity, dt), scale(acceleration, 0.5 * dt * dt)));
            rawVelocity = add(rawVelocity, scale(acceleration, dt));
        }
        const truth = vAt(reference.position, i);
        paths.truth.set(truth.slice(0, 2), 2 * k);
        paths.filter.set(filterPosition.slice(0, 2), 2 * k);
        paths.raw.set(rawPosition.slice(0, 2), 2 * k);
        paths.filterError[k] = Math.hypot(...filterPosition.map((value, axis) => value - truth[axis]));
        paths.rawError[k] = Math.hypot(...rawPosition.map((value, axis) => value - truth[axis]));
    }
    return paths;
}

/* The fixed turn about the vertical which lines a run's heading up with the reference. */
export function headingFix(run) {
    const offset = run.errors.offset;
    return [Math.cos(offset / 2), 0, 0, -Math.sin(offset / 2)];
}

class PositionFigure extends StoryFigure {
    build() {
        addSlider(this, 't0', {label: 'Start adding up at', min: 5, max: 60, value: 20, step: 0.1, format: v => `${v.toFixed(1)} s`});
        addSlider(this, 'length', {label: 'Carry on for', min: 0, max: LONGEST, value: 10, step: 0.05, format: v => `${v.toFixed(1)} s`});
        this.readout = addReadout(this);
        this.paths = null;
        this.pathsFrom = -1;
    }

    ready({recording, run}) {
        const reference = recording.reference;
        this.headingFix = headingFix(run);

        /* The whole reference path seen from above, and its extent. */
        this.outline = [];
        let low = [Infinity, Infinity], high = [-Infinity, -Infinity];
        for (let i = 0; i < recording.samples; i += 10) {
            if (reference.valid[i]) {
                const x = reference.position[3 * i], y = reference.position[3 * i + 1];
                this.outline.push([x, y]);
                low = [Math.min(low[0], x), Math.min(low[1], y)];
                high = [Math.max(high[0], x), Math.max(high[1], y)];
            }
        }
        this.extent = {low, high};

        let first = 0;
        while (first < run.valid && !run.value('rotor_model_active', first)) {
            first++;
        }
        const from = Math.ceil(recording.time[Math.min(first + 400, run.valid - 1)]);
        const to = Math.floor(recording.duration - LONGEST - 1);
        this.inputs.t0.input.min = from;
        this.inputs.t0.input.max = Math.max(to, from + 1);
        this.set({t0: clamp(this.params.t0, from, Math.max(to, from + 1))}, true);
    }

    paint(ctx, width, height) {
        const c = this.colors;
        const start = this.sampleAt(this.params.t0);
        if (this.pathsFrom !== start) {
            this.paths = reckon(this.story, this.headingFix, start);
            this.pathsFrom = start;
        }
        const paths = this.paths;
        const shown = clamp(Math.round(this.params.length / 0.005) + 1, 1, paths.count);

        /* The map on the left (or on top), the error plot beside (or under) it. */
        const wide = width >= 620;
        const mapBox = wide ? {left: 14, top: 34, width: width * 0.56 - 28, height: height - 62}
            : {left: 14, top: 34, width: width - 28, height: height * 0.6 - 40};
        const errorBox = wide ? {left: width * 0.56 + 46, top: 56, width: width * 0.44 - 62, height: height - 96}
            : {left: 46, top: height * 0.6 + 44, width: width - 62, height: height * 0.4 - 78};

        /* A map with the same scale both ways, which holds the whole reference path. */
        const {low, high} = this.extent;
        const margin = 0.8;
        const spanX = high[0] - low[0] + 2 * margin, spanY = high[1] - low[1] + 2 * margin;
        const metres = Math.max(spanX / mapBox.width, spanY / mapBox.height);
        const cx = (low[0] + high[0]) / 2, cy = (low[1] + high[1]) / 2;
        const map = new Axes(cx - metres * mapBox.width / 2, cx + metres * mapBox.width / 2,
            cy - metres * mapBox.height / 2, cy + metres * mapBox.height / 2).place(mapBox.left, mapBox.top, mapBox.width, mapBox.height);

        ctx.save();
        map.clip(ctx);
        /* A one-metre grid. */
        ctx.strokeStyle = c.grid;
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let x = Math.ceil(map.x0); x <= map.x1; x++) {
            ctx.moveTo(Math.round(map.x(x)) + 0.5, map.top);
            ctx.lineTo(Math.round(map.x(x)) + 0.5, map.top + map.height);
        }
        for (let y = Math.ceil(map.y0); y <= map.y1; y++) {
            ctx.moveTo(map.left, Math.round(map.y(y)) + 0.5);
            ctx.lineTo(map.left + map.width, Math.round(map.y(y)) + 0.5);
        }
        ctx.stroke();
        const trace = (values, count, color, lineWidth) => {
            ctx.strokeStyle = color;
            ctx.lineWidth = lineWidth;
            ctx.lineJoin = 'round';
            ctx.lineCap = 'round';
            ctx.beginPath();
            for (let k = 0; k < count; k++) {
                const px = map.x(values[2 * k]), py = map.y(values[2 * k + 1]);
                if (k === 0) {
                    ctx.moveTo(px, py);
                } else {
                    ctx.lineTo(px, py);
                }
            }
            ctx.stroke();
        };
        trace(this.outline.flat(), this.outline.length, alpha(c.ink2, 0.22), 1.5);
        trace(paths.raw, shown, alpha(c.orange, this.emphasis('raw')), 2);
        trace(paths.truth, shown, alpha(c.ink, this.emphasis('truth')), 2.5);
        trace(paths.filter, shown, alpha(c.blue, this.emphasis('filter')), 2);
        ctx.restore();
        const end = shown - 1;
        marker(ctx, map.x(paths.truth[0]), map.y(paths.truth[1]), c.ink2, 3.5);
        marker(ctx, map.x(paths.truth[2 * end]), map.y(paths.truth[2 * end + 1]), alpha(c.ink, this.emphasis('truth')), 4.5);
        if (map.contains(map.x(paths.filter[2 * end]), map.y(paths.filter[2 * end + 1]))) {
            marker(ctx, map.x(paths.filter[2 * end]), map.y(paths.filter[2 * end + 1]), alpha(c.blue, this.emphasis('filter')), 4.5);
        }
        if (map.contains(map.x(paths.raw[2 * end]), map.y(paths.raw[2 * end + 1]))) {
            marker(ctx, map.x(paths.raw[2 * end]), map.y(paths.raw[2 * end + 1]), alpha(c.orange, this.emphasis('raw')), 4.5);
        }
        /* A scale bar. */
        const bar = niceStep(metres * mapBox.width, 4);
        ctx.strokeStyle = c.ink2;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(map.left + 8, map.top + map.height - 10);
        ctx.lineTo(map.left + 8 + bar / metres, map.top + map.height - 10);
        ctx.stroke();
        text(ctx, `${bar} m`, map.left + 14 + bar / metres, map.top + map.height - 10, {baseline: 'middle', color: c.ink2, size: 11.5, halo: true});
        text(ctx, 'Seen from above', map.left, map.top - 9, {baseline: 'bottom', color: c.ink2, size: 12, weight: 500});
        legend(ctx, [
            {label: 'true path', color: c.ink},
            {label: 'filter', color: c.blue},
            {label: 'raw IMU', color: c.orange},
        ], wide ? errorBox.left : map.left + 130, wide ? 13 : map.top - 15, {size: 12});

        /* How far each is from the true position. */
        const errorAxes = new Axes(0, LONGEST, 0, ERROR_TOP).place(errorBox.left, errorBox.top, errorBox.width, errorBox.height);
        errorAxes.grid(ctx, {yLabel: 'Distance from the true position (m)', xFormat: seconds, xTicks: 4, yTicks: 4});
        const times = Float32Array.from({length: paths.count}, (_, k) => k * 0.005);
        ctx.save();
        errorAxes.clip(ctx);
        errorAxes.line(ctx, times, paths.rawError, {end: shown, color: alpha(c.orange, this.emphasis('raw'))});
        errorAxes.line(ctx, times, paths.filterError, {end: shown, color: alpha(c.blue, this.emphasis('filter'))});
        ctx.restore();
        if (paths.filterError[end] <= ERROR_TOP) {
            errorAxes.dot(ctx, times[end], paths.filterError[end], c.blue);
        }
        if (paths.rawError[end] <= ERROR_TOP) {
            errorAxes.dot(ctx, times[end], paths.rawError[end], c.orange);
        }

        const metresText = value => (value >= 100 ? value.toFixed(0) : value >= 10 ? value.toFixed(1) : value.toFixed(2));
        this.readout.textContent = `After ${(end * 0.005).toFixed(1)} s the filter's velocity, added up, is ${metresText(paths.filterError[end])} m `
            + `from the true position. The raw IMU, added up, is ${metresText(paths.rawError[end])} m away.`;
    }
}

registerFigure('position', PositionFigure);
