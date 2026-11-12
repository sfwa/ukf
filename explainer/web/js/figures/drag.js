/*
Rotor drag in the recorded flight. Each dot is half a second of flight: how
fast the aircraft moved across its rotor plane (from the reference system),
and what the accelerometer read in that direction. If the rotor model is
right, the dots lie on a straight line through the centre, and its slope is
the drag rate.
*/

import {registerFigure} from '../core/figure.js';
import {addSlider, addReadout, addButton, addRow} from '../core/controls.js';
import {Axes, text, marker} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {qAt, qConj, qRotate, vAt, dot} from '../core/math.js';
import {StoryFigure} from './story-kit.js';

const BLOCK = 100;      /* samples per dot: half a second */

/*
One dot for each half second of powered flight and each of the two
sideways axes: {x: true speed, y: accelerometer reading, time}. Also the
straight line which fits the dots best.
*/
export function dragPoints({recording, run, estimator}) {
    const constants = estimator.constants;
    const axes = [
        [constants.lateral_axis_1_x, constants.lateral_axis_1_y, constants.lateral_axis_1_z],
        [constants.lateral_axis_2_x, constants.lateral_axis_2_y, constants.lateral_axis_2_z],
    ];
    const reference = recording.reference;
    const points = [];
    for (let start = 0; start + BLOCK <= run.valid; start += BLOCK) {
        const velocity = [0, 0, 0], reading = [0, 0, 0];
        let usable = true;
        for (let i = start; i < start + BLOCK; i++) {
            if (!reference.valid[i] || !run.value('rotor_model_active', i)) {
                usable = false;
                break;
            }
            /* The reference velocity is in the world frame; turn it into the IMU frame. */
            const body = qRotate(qConj(qAt(reference.attitude, i)), vAt(reference.velocity, i));
            for (let k = 0; k < 3; k++) {
                velocity[k] += body[k] / BLOCK;
                reading[k] += recording.accelerometer[3 * i + k] / BLOCK;
            }
        }
        if (usable) {
            for (const axis of axes) {
                points.push({x: dot(velocity, axis), y: dot(reading, axis), time: recording.time[start + BLOCK / 2]});
            }
        }
    }
    const n = points.length;
    const meanX = points.reduce((sum, p) => sum + p.x, 0) / n;
    const meanY = points.reduce((sum, p) => sum + p.y, 0) / n;
    let sxx = 0, sxy = 0;
    for (const p of points) {
        sxx += (p.x - meanX) ** 2;
        sxy += (p.x - meanX) * (p.y - meanY);
    }
    return {points, fit: {slope: sxy / sxx, intercept: meanY - sxy / sxx * meanX}};
}

class DragFigure extends StoryFigure {
    build() {
        addSlider(this, 'rate', {
            label: 'Drag rate of the rule', min: 0, max: 0.6, value: 0.25, step: 0.005,
            format: v => `${v.toFixed(2)} per second`,
        });
        const row = addRow(this);
        addButton(this, {label: "Use the estimator's value", parent: row, onClick: () => {
            this.set({rate: this.story.estimator.constants.drag_rate});
        }});
        this.readout = addReadout(this, {parent: row});
        this.points = [];
        this.near = null;
        this.onPointer({
            move: at => this.hover(at),
            leave: () => {
                this.near = null;
                this.invalidate();
            },
        });
    }

    ready(story) {
        const {points, fit} = dragPoints(story);
        this.points = points;
        this.bestFit = fit;
        this.set({rate: story.estimator.constants.drag_rate}, true);
    }

    hover(at) {
        if (!this.axes) {
            return;
        }
        let best = null, distance = 24;
        for (const point of this.points) {
            const d = Math.hypot(this.axes.x(point.x) - at.x, this.axes.y(point.y) - at.y);
            if (d < distance) {
                distance = d;
                best = point;
            }
        }
        if (best !== this.near) {
            this.near = best;
            this.invalidate();
        }
    }

    paint(ctx, width, height) {
        const c = this.colors;
        const axes = new Axes(-2.6, 2.6, -0.8, 0.8).place(52, 30, width - 70, height - 74);
        this.axes = axes;
        axes.grid(ctx, {
            xTicks: 6, yTicks: 4,
            yLabel: 'Sideways accelerometer reading (m/s²)',
            xLabel: 'sideways speed, from the reference (m/s)',
        });
        ctx.save();
        axes.clip(ctx);
        ctx.fillStyle = alpha(c.orange, 0.5 * this.emphasis('reading'));
        for (const point of this.points) {
            ctx.beginPath();
            ctx.arc(axes.x(point.x), axes.y(point.y), 3, 0, 2 * Math.PI);
            ctx.fill();
        }
        /* The best-fit line, thin and neutral. */
        ctx.strokeStyle = alpha(c.ink2, 0.8 * this.emphasis('fit'));
        ctx.lineWidth = 1.5;
        ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(axes.x(-2.6), axes.y(this.bestFit.intercept - 2.6 * this.bestFit.slope));
        ctx.lineTo(axes.x(2.6), axes.y(this.bestFit.intercept + 2.6 * this.bestFit.slope));
        ctx.stroke();
        ctx.setLineDash([]);
        /* The model's line: reading = -rate x speed. */
        const rate = this.params.rate;
        ctx.strokeStyle = alpha(c.green, this.emphasis('model'));
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(axes.x(-2.6), axes.y(2.6 * rate));
        ctx.lineTo(axes.x(2.6), axes.y(-2.6 * rate));
        ctx.stroke();
        ctx.restore();

        /* Labels at the right-hand ends of the two lines, kept apart. */
        const modelY = axes.y(-2.2 * rate), fitY = axes.y(this.bestFit.intercept + 2.2 * this.bestFit.slope);
        const apart = Math.abs(modelY - fitY) < 16 ? (modelY <= fitY ? -1 : 1) : 0;
        text(ctx, 'rule', axes.x(2.2), modelY - 8 + (apart < 0 ? -6 : apart > 0 ? 22 : 0),
            {align: 'center', baseline: 'bottom', halo: true, color: c.ink2, weight: 500});
        text(ctx, 'best fit', axes.x(2.2), fitY + 8 + (apart < 0 ? 8 : apart > 0 ? -24 : 0),
            {align: 'center', baseline: 'top', halo: true, color: c.ink2});

        if (this.near) {
            marker(ctx, axes.x(this.near.x), axes.y(this.near.y), c.orange, 5);
        }
        const fitText = `Best straight line through the dots: ${(-this.bestFit.slope).toFixed(2)} per second.`;
        this.readout.textContent = this.near
            ? `At ${this.near.time.toFixed(1)} s: speed ${this.near.x.toFixed(2)} m/s, reading ${this.near.y.toFixed(2)} m/s².  ${fitText}`
            : fitText;
    }
}

registerFigure('drag', DragFigure);
