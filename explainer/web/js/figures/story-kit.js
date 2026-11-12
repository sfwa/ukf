/*
Shared parts of the figures which show the recorded flight: loading the
shared recording and run, laying out stacked time plots, and moving a time
cursor with the pointer.
*/

import {Figure} from '../core/figure.js';
import {text} from '../core/plot.js';
import {clamp, qRotate, qMul, qFromRotationVector, qBetween, vAt, scale, add} from '../core/math.js';
import {loadStory} from '../data/story.js';

const DEGREES = 180 / Math.PI;

/* A tick label for a time axis. */
export const seconds = value => `${Number(value.toFixed(3))} s`;

export class StoryFigure extends Figure {
    setup() {
        this.story = null;
        this.failure = null;
        this.build();
        loadStory().then(story => {
            this.story = story;
            this.ready(story);
            this.invalidate();
        }).catch(error => {
            console.error(error);
            this.failure = error.message;
            this.invalidate();
        });
    }

    /* For subclasses. */
    build() {}          /* make the controls */
    ready() {}          /* the recording and the run have arrived */
    paint() {}          /* draw, once they have */

    draw(ctx, width, height) {
        if (!this.story) {
            text(ctx, this.failure ?? 'Loading the recording…', width / 2, height / 2,
                {align: 'center', baseline: 'middle', color: this.colors.muted, size: 14});
            return;
        }
        this.paint(ctx, width, height);
    }

    /* The sample nearest to a time in seconds. */
    sampleAt(t) {
        const recording = this.story.recording;
        return clamp(Math.round(t / recording.duration * (recording.samples - 1)), 0, recording.samples - 1);
    }

    /*
    Let the pointer move the time cursor: pressing or dragging anywhere on
    the canvas sets this.params[name] from the x position on "axes()".
    */
    scrub(name, axes, limits = null) {
        const move = at => {
            const plot = axes();
            if (!plot) {
                return;
            }
            const [low, high] = limits ? limits() : [plot.x0, plot.x1];
            this.set({[name]: clamp(plot.fromX(at.x), low, high)}, true);
        };
        this.canvas.classList.add('scrub');
        this.onPointer({down: move, drag: move});
    }
}

/*
Split a canvas into plot areas stacked top to bottom, which share one time
axis. "weights" gives their relative heights. Returns [{left, top, width,
height}, ...].
*/
export function stack(width, height, weights, margins = {}) {
    const {left = 46, right = 14, top = 26, bottom = 26, gap = 30} = margins;
    const total = weights.reduce((sum, weight) => sum + weight, 0);
    const usable = height - top - bottom - gap * (weights.length - 1);
    const areas = [];
    let y = top;
    for (const weight of weights) {
        const h = usable * weight / total;
        areas.push({left, top: y, width: width - left - right, height: h});
        y += h + gap;
    }
    return areas;
}

/* The angle in degrees between an attitude's thrust axis and straight up. */
export function tiltDegrees(attitude, thrustAxis) {
    return Math.acos(clamp(qRotate(attitude, thrustAxis)[2], -1, 1)) * DEGREES;
}

/*
The attitude which the gyroscope alone would give: start as the estimator
starts (level by the first accelerometer sample) and add up the turn rates,
with no correction of any kind. Returns w x y z per sample, up to "count".
*/
export function gyroOnlyAttitude(recording, count = recording.samples) {
    const out = new Float64Array(4 * count);
    let q = qBetween(vAt(recording.accelerometer, 0), [0, 0, 1]);
    out.set(q, 0);
    for (let i = 1; i < count; i++) {
        const dt = recording.time[i] - recording.time[i - 1];
        const rate = scale(add(vAt(recording.gyroscope, i - 1), vAt(recording.gyroscope, i)), 0.5);
        q = qMul(q, qFromRotationVector(scale(rate, dt)));
        const length = Math.hypot(q[0], q[1], q[2], q[3]);
        q = q.map(value => value / length);
        out.set(q, 4 * i);
    }
    return out;
}
