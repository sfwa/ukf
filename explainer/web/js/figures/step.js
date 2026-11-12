/*
The filter at work, slowed down. In flight it takes 200 steps every
second. At each one it works out what the accelerometer should read across
the rotor plane, compares that with the real reading, and corrects its
state by a small part of the difference.

The plot is the plane of that reading. Green is what the filter expected,
with rings at one and two standard deviations. Orange is what the
accelerometer read, with the readings of the last half second behind it.
The numbers underneath are the real step, taken apart by the estimator
module (ex_explain in explainer/wasm).
*/

import {registerFigure} from '../core/figure.js';
import {addPlayBar, addButton, addRow, addReadout} from '../core/controls.js';
import {Axes, text, marker, arrow, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {clamp, vAt, dot} from '../core/math.js';
import {StoryFigure} from './story-kit.js';

const SLOW = 1 / 20;        /* playback speed: 10 filter steps per second */
const TRAIL = 100;          /* readings kept on screen */
const RANGE = 2.6;          /* m/s^2 from the centre of the plot to its edge */
const DEGREES = 180 / Math.PI;

class StepFigure extends StoryFigure {
    build() {
        addPlayBar(this, 't', {duration: 1, value: 30, format: v => `${v.toFixed(3)} s`});
        const row = addRow(this);
        addButton(this, {label: '‹ One step back', parent: row, onClick: () => this.nudge(-1)});
        addButton(this, {label: 'One step on ›', parent: row, onClick: () => this.nudge(1)});
        this.table = document.createElement('table');
        this.table.className = 'numbers';
        this.controls.append(this.table);
        this.detail = null;         /* the step at this.detailIndex, taken apart by the estimator */
        this.detailIndex = -1;
        this.asking = false;
    }

    ready({recording, run, estimator}) {
        const constants = estimator.constants;
        this.lateral = [
            [constants.lateral_axis_1_x, constants.lateral_axis_1_y, constants.lateral_axis_1_z],
            [constants.lateral_axis_2_x, constants.lateral_axis_2_y, constants.lateral_axis_2_z],
        ];
        this.inputs.t.setDuration(recording.duration);
        this.set({t: Math.min(30, recording.duration / 2)}, true);
    }

    nudge(steps) {
        const recording = this.story.recording;
        const i = clamp(this.sampleAt(this.params.t) + steps, 1, recording.samples - 1);
        this.set({t: recording.time[i], tPlaying: false}, true);
    }

    step(dt) {
        if (!this.story || !this.params.tPlaying) {
            return false;
        }
        const duration = this.story.recording.duration;
        const t = this.params.t + dt * SLOW;
        if (t >= duration) {
            this.set({t: duration, tPlaying: false}, true);
            return false;
        }
        this.set({t}, true);
        return true;
    }

    /* The reading across the rotor plane at sample i, and what the filter expected. */
    reading(i) {
        const {recording, run} = this.story;
        if (!run.value('rotor_model_active', i)) {
            return null;
        }
        const sample = vAt(recording.accelerometer, i);
        const read = [dot(sample, this.lateral[0]), dot(sample, this.lateral[1])];
        const difference = [run.value('innovation', i, 0), run.value('innovation', i, 1)];
        return {
            read,
            difference,
            expected: [read[0] - difference[0], read[1] - difference[1]],
            sd: [run.value('innovation_sd', i, 0), run.value('innovation_sd', i, 1)],
        };
    }

    /* Ask the estimator to take the step at this sample apart; keep only the newest request. */
    async explain(i) {
        if (this.asking || i === this.detailIndex) {
            return;
        }
        this.asking = true;
        const detail = await this.story.estimator.explain(i);
        this.asking = false;
        this.detail = detail;
        this.detailIndex = i;
        this.fillTable();
        const now = Math.max(this.sampleAt(this.params.t), 1);
        if (now !== i) {
            this.explain(now);
        }
    }

    fillTable() {
        const i = Math.max(this.sampleAt(this.params.t), 1);
        const now = this.reading(i);
        const rows = [];
        const pair = values => values.map(v => v.toFixed(2).replace('-', '−')).join(', ');
        if (!now) {
            rows.push(['At this moment the rotor model is not in use', '']);
        } else {
            rows.push(['The accelerometer read', `${pair(now.read)} m/s²`]);
            rows.push(['The filter expected', `${pair(now.expected)} m/s²`]);
            rows.push(['Difference', `${pair(now.difference)} m/s²`]);
            if (this.detail && this.detailIndex === i) {
                const correction = this.detail.column('correction');
                const tilt = Math.hypot(correction[0], correction[1]) * DEGREES;
                const velocity = Math.hypot(correction[3], correction[4], correction[5]);
                const before = Math.sqrt(this.detail.at('predicted_covariance', 0, 0) + this.detail.at('predicted_covariance', 1, 1)) * DEGREES;
                const after = Math.sqrt(this.detail.at('final_covariance', 0, 0) + this.detail.at('final_covariance', 1, 1)) * DEGREES;
                rows.push(['Correction to the tilt at this step', `${tilt.toFixed(3)}°`]);
                rows.push(['Correction to the velocity at this step', `${velocity.toFixed(3)} m/s`]);
                rows.push(['Tilt uncertainty before and after', `${before.toFixed(4)}° and ${after.toFixed(4)}°`]);
            }
        }
        /* Rebuild the table only when its contents change. */
        const contents = JSON.stringify(rows);
        if (contents === this.tableContents) {
            return;
        }
        this.tableContents = contents;
        this.table.replaceChildren(...rows.map(([name, value]) => {
            const row = document.createElement('tr');
            const head = document.createElement('td');
            head.textContent = name;
            const cell = document.createElement('td');
            cell.textContent = value;
            row.append(head, cell);
            return row;
        }));
    }

    paint(ctx, width, height) {
        const c = this.colors;
        const {run} = this.story;
        const i = Math.max(this.sampleAt(this.params.t), 1);
        const side = Math.min(width - 60, height - 78);
        const axes = new Axes(-RANGE, RANGE, -RANGE, RANGE).place((width - side) / 2 + 14, 30, side, side);
        axes.grid(ctx, {xTicks: 4, yTicks: 4, yLabel: 'Sideways accelerometer reading (m/s²)'});

        const now = this.reading(i);
        if (!now) {
            text(ctx, 'The rotor model is not in use here.', axes.left + side / 2, axes.top + side / 2,
                {align: 'center', baseline: 'middle', color: c.muted, size: 13, halo: true});
        } else {
            ctx.save();
            axes.clip(ctx);
            /* What the filter expected, with rings at one and two standard deviations. */
            const cx = axes.x(now.expected[0]), cy = axes.y(now.expected[1]);
            const rx = now.sd[0] / (2 * RANGE) * side, ry = now.sd[1] / (2 * RANGE) * side;
            for (const sigmas of [2, 1]) {
                ctx.fillStyle = alpha(c.green, 0.1 * this.emphasis('expected'));
                ctx.beginPath();
                ctx.ellipse(cx, cy, rx * sigmas, ry * sigmas, 0, 0, 2 * Math.PI);
                ctx.fill();
            }
            /* The readings of the last half second, fading with age. */
            for (let k = TRAIL; k >= 1; k--) {
                const past = this.reading(i - k);
                if (past) {
                    ctx.fillStyle = alpha(c.orange, 0.5 * (1 - k / (TRAIL + 1)) * this.emphasis('reading'));
                    ctx.beginPath();
                    ctx.arc(axes.x(past.read[0]), axes.y(past.read[1]), 2.5, 0, 2 * Math.PI);
                    ctx.fill();
                }
            }
            arrow(ctx, cx, cy, axes.x(now.read[0]), axes.y(now.read[1]), alpha(c.ink2, this.emphasis('difference')), {width: 1.5, head: 8});
            marker(ctx, cx, cy, alpha(c.green, this.emphasis('expected')), 5.5);
            marker(ctx, axes.x(now.read[0]), axes.y(now.read[1]), alpha(c.orange, this.emphasis('reading')), 5.5);
            ctx.restore();
        }
        legend(ctx, [
            {label: 'expected', color: c.green, dot: true},
            {label: 'read', color: c.orange, dot: true},
        ], axes.left, axes.top + side + 34, {size: 11.5});
        text(ctx, `step ${i} of ${run.valid - 1}`, axes.left + side, axes.top + side + 34,
            {align: 'right', baseline: 'middle', color: c.muted, size: 11.5});

        /* The numbers: refresh at once when still, and a few times a second while playing. */
        if (!this.params.tPlaying || Math.abs(i - this.detailIndex) >= 3) {
            this.explain(i);
        }
        this.fillTable();
    }
}

registerFigure('step', StepFigure);
