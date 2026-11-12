/*
The whole loop of a Kalman filter, on the simplest moving thing: a
quadcopter which flies back and forth along a line. The filter keeps two
numbers, the position and the speed. It predicts all the time: an
accelerometer tells it how the speed changes, the estimate moves on, and
its spread grows. Each position reading then pulls the estimate towards
itself and shrinks the spread.

Lower panel: the spread of the position estimate. It rises between
readings and drops at each one.

The true path and the random part of every sensor value are fixed, so that
the sliders change only what the filter is told.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addToggle, addRow, addReadout} from '../core/controls.js';
import {Axes, text, legend, marker} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {Random, clamp} from '../core/math.js';
import {drawValueBox} from './concept-kit.js';

const DURATION = 20;            /* s */
const STEP = 0.05;              /* the filter predicts this often (s) */
const SAMPLES = Math.round(DURATION / STEP) + 1;
const OUTAGE = [8, 14];         /* the position sensor is blocked between these times (s) */
const SETTLED = 3;              /* errors are measured from this time on (s) */
const ACCELEROMETER_NOISE = 0.3;    /* spread of each accelerometer sample (m/s^2) */

/* Where the quadcopter really is at time t (m), and its speed and acceleration. */
const SLOW = 2 * Math.PI / 16, FAST = 2 * Math.PI / 5.7;
const truePosition = t => 5 + 3.2 * Math.sin(SLOW * t - 0.9) + 0.4 * Math.sin(FAST * t + 1);
const trueSpeed = t => 3.2 * SLOW * Math.cos(SLOW * t - 0.9) + 0.4 * FAST * Math.cos(FAST * t + 1);
const trueAcceleration = t => -3.2 * SLOW * SLOW * Math.sin(SLOW * t - 0.9) - 0.4 * FAST * FAST * Math.sin(FAST * t + 1);

class LoopFigure extends Figure {
    setup() {
        const random = new Random(262);
        this.draws = Array.from({length: 200}, () => random.normal());
        const accelerometer = new Random(17);
        this.accelerometerDraws = Array.from({length: SAMPLES}, () => accelerometer.normal());
        addSlider(this, 'sensorNoise', {label: 'Sensor noise', min: 0.05, max: 2, value: 0.8, step: 0.05,
            format: v => `${v.toFixed(2)} m`});
        addSlider(this, 'interval', {label: 'Time between readings', min: 0.2, max: 3, value: 0.5, step: 0.1,
            format: v => `${v.toFixed(1)} s`});
        addSlider(this, 'modelNoise', {label: 'Model noise', min: 0.01, max: 2, value: 0.07, log: true,
            format: v => v.toFixed(2)});
        const row = addRow(this);
        addToggle(this, 'outage', {label: `Block the position sensor from ${OUTAGE[0]} s to ${OUTAGE[1]} s`, value: false, parent: row});
        addToggle(this, 'accelerometer', {label: 'Use the accelerometer', value: true, parent: row});
        this.readout = addReadout(this);
        this.result = null;
        this.hover = null;
        this.axes = new Axes(0, DURATION, 0, 10);
        this.canvas.classList.add('scrub');
        const point = at => {
            this.hover = this.axes.contains(at.x, at.y, 40) ? clamp(this.axes.fromX(at.x), 0, DURATION) : null;
            this.invalidate();
        };
        this.onPointer({
            down: point,
            drag: point,
            move: point,
            leave: () => {
                this.hover = null;
                this.invalidate();
            },
        });
    }

    changed() {
        this.result = null;
    }

    /* Run the filter over the whole flight with the present settings. */
    simulate() {
        const {sensorNoise, interval, modelNoise, outage, accelerometer} = this.params;
        const time = new Float32Array(SAMPLES), truth = new Float32Array(SAMPLES);
        const estimate = new Float32Array(SAMPLES), spread = new Float32Array(SAMPLES);
        const readings = [];

        /* The estimate (position, speed) and its covariance. The first guess is a poor one. */
        let position = truePosition(0) + 1.5, speed = trueSpeed(0) - 0.4;
        let pp = 1.5 ** 2, pv = 0, vv = 0.5 ** 2;
        const q = modelNoise ** 2;
        let next = 1;       /* number of the next reading; reading k comes at k * interval */

        for (let i = 0; i < SAMPLES; i++) {
            const t = i * STEP;
            if (i > 0) {
                /*
                Predict. The accelerometer says how the speed changes; without it the
                filter must assume that the speed stays as it is. What the model
                does not know adds spread.
                */
                const measured = accelerometer ?
                    trueAcceleration(t - STEP / 2) + ACCELEROMETER_NOISE * this.accelerometerDraws[i] : 0;
                position += speed * STEP + 0.5 * measured * STEP * STEP;
                speed += measured * STEP;
                pp += 2 * STEP * pv + STEP * STEP * vv + q * STEP ** 3 / 3;
                pv += STEP * vv + q * STEP ** 2 / 2;
                vv += q * STEP;
            }
            while (next * interval <= t + 1e-9) {
                const at = next * interval;
                const blocked = outage && at >= OUTAGE[0] && at < OUTAGE[1];
                if (!blocked) {
                    /* Correct: the gain says how far to move towards the reading. */
                    const value = truePosition(at) + sensorNoise * this.draws[next % this.draws.length];
                    const innovation = value - position;
                    const s = pp + sensorNoise ** 2;
                    const gainPosition = pp / s, gainSpeed = pv / s;
                    position += gainPosition * innovation;
                    speed += gainSpeed * innovation;
                    const pp0 = pp, pv0 = pv;
                    pp -= gainPosition * pp0;
                    pv -= gainPosition * pv0;
                    vv -= gainSpeed * pv0;
                    readings.push({time: at, value});
                }
                next++;
            }
            time[i] = t;
            truth[i] = truePosition(t);
            estimate[i] = position;
            spread[i] = Math.sqrt(Math.max(pp, 0));
        }

        /* How wrong the estimate and the raw readings are, once the filter has settled. */
        let estimateSquares = 0, estimateCount = 0, readingSquares = 0, readingCount = 0, widest = 0, inside = 0;
        for (let i = 0; i < SAMPLES; i++) {
            if (time[i] >= SETTLED) {
                estimateSquares += (estimate[i] - truth[i]) ** 2;
                estimateCount++;
                widest = Math.max(widest, spread[i]);
                if (Math.abs(estimate[i] - truth[i]) <= 2 * spread[i]) {
                    inside++;
                }
            }
        }
        for (const reading of readings) {
            if (reading.time >= SETTLED) {
                readingSquares += (reading.value - truePosition(reading.time)) ** 2;
                readingCount++;
            }
        }
        return {
            time, truth, estimate, spread, readings, widest,
            inside: inside / estimateCount,
            estimateError: Math.sqrt(estimateSquares / estimateCount),
            readingError: readingCount > 0 ? Math.sqrt(readingSquares / readingCount) : NaN,
        };
    }

    draw(ctx, width, height) {
        const c = this.colors;
        this.result ??= this.simulate();
        const {time, truth, estimate, spread, readings, widest} = this.result;
        const narrow = width < 560;
        const left = narrow ? 30 : 40, right = width - 18;
        const top = 54, bottom = height - 26;
        const lowerHeight = Math.round((bottom - top) * 0.24);
        const upperBottom = bottom - lowerHeight - 34;
        const axes = this.axes.place(left, top, right - left, upperBottom - top);
        const lower = new Axes(0, DURATION, 0, Math.max(0.4, widest * 1.25)).place(left, bottom - lowerHeight, right - left, lowerHeight);
        const seconds = v => `${v} s`;

        /* Upper panel: position against time. */
        axes.grid(ctx, {xTicks: narrow ? 5 : 10, yTicks: 5, xLabels: false, yLabel: 'position (m)'});
        if (this.params.outage) {
            for (const panel of [axes, lower]) {
                ctx.fillStyle = alpha(c.ink, c.dark ? 0.1 : 0.055);
                ctx.fillRect(panel.x(OUTAGE[0]), panel.top, panel.x(OUTAGE[1]) - panel.x(OUTAGE[0]), panel.height);
            }
            text(ctx, 'sensor blocked', axes.x((OUTAGE[0] + OUTAGE[1]) / 2), axes.top + 14, {align: 'center', color: c.ink2, size: 11.5});
        }
        ctx.save();
        axes.clip(ctx);
        const spreadOn = this.emphasis('spread'), estimateOn = this.emphasis('estimate');
        axes.band(ctx, time, i => estimate[i] - 2 * spread[i], i => estimate[i] + 2 * spread[i],
            {color: c.blue, opacity: 0.16 * Math.max(spreadOn, estimateOn * 0.6)});
        axes.line(ctx, time, truth, {color: alpha(c.ink, 0.9 * this.emphasis('truth')), width: 2});
        axes.line(ctx, time, estimate, {color: alpha(c.blue, estimateOn), width: 2});
        const readingOn = this.emphasis('reading');
        const many = readings.length > 60;
        for (const reading of readings) {
            marker(ctx, axes.x(reading.time), axes.y(reading.value), alpha(c.orange, readingOn), many ? 3 : 4, many ? 0 : 1.5);
        }
        ctx.restore();

        /* Lower panel: the spread of the position estimate. */
        lower.grid(ctx, {xTicks: narrow ? 5 : 10, yTicks: 2, xFormat: seconds, yLabel: 'spread of the estimate (m)'});
        ctx.save();
        lower.clip(ctx);
        lower.line(ctx, time, spread, {color: alpha(c.blue, Math.max(spreadOn, estimateOn * 0.6)), width: 2});
        ctx.restore();

        legend(ctx, [
            {label: narrow ? 'truth' : 'true position', color: c.ink},
            {label: 'reading', color: c.orange, dot: true},
            {label: 'estimate', color: c.blue},
            {label: narrow ? '± 2 spreads' : 'estimate ± 2 spreads', color: c.blue, band: true, width: 0},
        ], left, 14, {size: narrow ? 11.5 : 12.5});

        /* The cursor: the values at one moment. */
        if (this.hover !== null) {
            const i = clamp(Math.round(this.hover / STEP), 0, SAMPLES - 1);
            axes.cursor(ctx, time[i], c.muted);
            lower.cursor(ctx, time[i], c.muted);
            axes.dot(ctx, time[i], truth[i], c.ink, 3.5);
            axes.dot(ctx, time[i], estimate[i], c.blue, 4);
            lower.dot(ctx, time[i], Math.min(spread[i], lower.y1), c.blue, 4);
            drawValueBox(ctx, axes.x(time[i]), axes.top + 6, `${time[i].toFixed(1)} s`, [
                {label: 'true position', value: `${truth[i].toFixed(2)} m`, color: c.ink},
                {label: 'estimate', value: `${estimate[i].toFixed(2)} m`, color: c.blue},
                {label: 'its spread', value: `${spread[i].toFixed(2)} m`},
            ], c, {left: axes.left, top: axes.top, right, bottom: upperBottom});
        }

        const reading = Number.isFinite(this.result.readingError) ? `${this.result.readingError.toFixed(2)} m` : 'none';
        this.readout.textContent = `typical error after ${SETTLED} s:  readings ${reading}  ·  ` +
            `estimate ${this.result.estimateError.toFixed(2)} m  ·  ` +
            `truth inside the band ${Math.round(this.result.inside * 100)}% of the time`;
    }
}

registerFigure('loop', LoopFigure);
