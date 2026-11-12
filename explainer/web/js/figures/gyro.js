/*
The gyroscope on its own: the first seconds of the recording. The aircraft
stands still, but the gyroscope reads a steady turn (its bias). Adding the
readings up makes the attitude lean further and further over. The filter
measures the bias while the aircraft is at rest, and removes it.

Top: tilt from the gyroscope alone, from the filter, and from the
reference. Bottom: the filter's estimate of the bias on one axis.
*/

import {registerFigure} from '../core/figure.js';
import {addReadout} from '../core/controls.js';
import {Axes, text, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {qAt, DEG} from '../core/math.js';
import {StoryFigure, stack, seconds, tiltDegrees, gyroOnlyAttitude} from './story-kit.js';

const WINDOW = 14;      /* seconds shown */
const BIAS_AXIS = 2;    /* the IMU axis with the largest bias */

class GyroFigure extends StoryFigure {
    build() {
        this.params.t = 3;
        this.readout = addReadout(this);
        this.scrub('t', () => this.tiltAxes);
    }

    ready({recording, run, estimator}) {
        const constants = estimator.constants;
        const thrust = [constants.thrust_axis_x, constants.thrust_axis_y, constants.thrust_axis_z];
        const count = Math.min(recording.samples, Math.round(WINDOW / recording.duration * (recording.samples - 1)) + 1);
        this.count = count;
        const gyroOnly = gyroOnlyAttitude(recording, count);
        this.series = {
            gyro: new Float32Array(count),
            filter: new Float32Array(count),
            reference: new Float32Array(count).fill(NaN),
            bias: new Float32Array(count),
            biasSd: new Float32Array(count),
            biasReference: new Float32Array(count).fill(NaN),
        };
        for (let i = 0; i < count; i++) {
            this.series.gyro[i] = tiltDegrees(qAt(gyroOnly, i), thrust);
            this.series.filter[i] = tiltDegrees(run.attitude(i), thrust);
            this.series.bias[i] = run.value('gyro_bias', i, BIAS_AXIS) / DEG;
            this.series.biasSd[i] = run.value('sd', i, 6 + BIAS_AXIS) / DEG;
            if (recording.reference.valid[i]) {
                this.series.reference[i] = tiltDegrees(qAt(recording.reference.attitude, i), thrust);
                this.series.biasReference[i] = recording.reference.gyroBias[3 * i + BIAS_AXIS] / DEG;
            }
        }
    }

    paint(ctx, width, height) {
        const c = this.colors;
        const {recording, run} = this.story;
        const time = recording.time;
        const series = this.series;
        const [top, bottom] = stack(width, height, [1.5, 1], {top: 56, gap: 44});
        const tiltAxes = new Axes(0, WINDOW, 0, 60).place(top.left, top.top, top.width, top.height);
        const biasAxes = new Axes(0, WINDOW, -1, 6).place(bottom.left, bottom.top, bottom.width, bottom.height);
        this.tiltAxes = tiltAxes;

        /* Shade the time when the filter knows the aircraft is at rest. */
        const resting = i => run.value('resting', i) > 0;
        for (const axes of [tiltAxes, biasAxes]) {
            axes.spans(ctx, time, resting, {end: this.count, color: c.ink, opacity: c.dark ? 0.1 : 0.055});
        }
        let firstRest = 0;
        while (firstRest < this.count && !resting(firstRest)) {
            firstRest++;
        }
        if (firstRest < this.count) {
            text(ctx, 'at rest', tiltAxes.x(time[firstRest]) + 6, tiltAxes.top + 6, {baseline: 'top', color: c.muted, size: 11.5});
        }

        tiltAxes.grid(ctx, {yLabel: 'Tilt (degrees)', xLabels: false, yTicks: 4});
        biasAxes.grid(ctx, {yLabel: 'Gyroscope bias, one axis (degrees per second)', xFormat: seconds, yTicks: 3});

        ctx.save();
        tiltAxes.clip(ctx);
        tiltAxes.line(ctx, time, series.reference, {end: this.count, color: alpha(c.ink, this.emphasis('reference'))});
        tiltAxes.line(ctx, time, series.gyro, {end: this.count, color: alpha(c.orange, this.emphasis('gyro'))});
        tiltAxes.line(ctx, time, series.filter, {end: this.count, color: alpha(c.blue, this.emphasis('estimate'))});
        ctx.restore();
        legend(ctx, [
            {label: 'gyroscope alone', color: c.orange},
            {label: 'filter', color: c.blue},
            {label: 'reference', color: c.ink},
        ], tiltAxes.left, 11);

        ctx.save();
        biasAxes.clip(ctx);
        biasAxes.band(ctx, time, i => series.bias[i] - 2 * series.biasSd[i], i => series.bias[i] + 2 * series.biasSd[i],
            {end: this.count, color: c.blue, opacity: 0.18 * this.emphasis('estimate')});
        biasAxes.line(ctx, time, series.biasReference, {end: this.count, color: alpha(c.ink, this.emphasis('reference'))});
        biasAxes.line(ctx, time, series.bias, {end: this.count, color: alpha(c.blue, this.emphasis('estimate'))});
        ctx.restore();

        /* The time cursor and what the lines read there. */
        const t = this.params.t;
        const i = Math.min(this.sampleAt(t), this.count - 1);
        tiltAxes.cursor(ctx, t);
        biasAxes.cursor(ctx, t);
        tiltAxes.dot(ctx, t, Math.min(series.gyro[i], 60), c.orange);
        tiltAxes.dot(ctx, t, series.filter[i], c.blue);
        biasAxes.dot(ctx, t, series.bias[i], c.blue);
        this.readout.innerHTML = '';
        const parts = [
            `${t.toFixed(1)} s`,
            `gyroscope alone ${series.gyro[i].toFixed(1)}°`,
            `filter ${series.filter[i].toFixed(1)}°`,
        ];
        if (Number.isFinite(series.reference[i])) {
            parts.push(`reference ${series.reference[i].toFixed(1)}°`);
        }
        parts.push(`bias estimate ${series.bias[i].toFixed(2)} ± ${(2 * series.biasSd[i]).toFixed(2)} °/s`);
        this.readout.textContent = parts.join('  ·  ');
    }
}

registerFigure('gyro', GyroFigure);
