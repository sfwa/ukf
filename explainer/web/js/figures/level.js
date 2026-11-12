/*
Can the accelerometer be used as a spirit level? Carried by hand, yes: it
feels the hand pushing up against gravity, so it shows which way is up.
In flight, no: the only large force is the rotor thrust, which always
points along the aircraft's own axis, however the aircraft is tilted.

The plot compares the true tilt with the tilt which the accelerometer
would show if it were used as a level.
*/

import {registerFigure} from '../core/figure.js';
import {addChoice, addReadout} from '../core/controls.js';
import {Axes, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {clamp, qAt, vAt, dot, norm} from '../core/math.js';
import {StoryFigure, seconds, tiltDegrees} from './story-kit.js';
import {loadRun} from '../data/story.js';

const HAND_RECORDING = 'MH_01_easy';
const HAND_SECONDS = 19;
const SMOOTH_SECONDS = 0.5;
const DEGREES = 180 / Math.PI;

/* The two tilts for one recording, between two times. */
export function tilts(recording, thrust, from, to) {
    const first = Math.round(from / recording.duration * (recording.samples - 1));
    const last = Math.round(to / recording.duration * (recording.samples - 1));
    const half = Math.round(SMOOTH_SECONDS / 2 / (recording.duration / (recording.samples - 1)));
    const count = last - first + 1;
    const time = new Float32Array(count), truth = new Float32Array(count).fill(NaN), level = new Float32Array(count);

    /* A running mean of the accelerometer takes the rotor vibration out. */
    const sum = [0, 0, 0];
    let low = first - half, high = first - half - 1;
    const push = i => {
        if (i >= 0 && i < recording.samples) {
            for (let k = 0; k < 3; k++) {
                sum[k] += recording.accelerometer[3 * i + k];
            }
        }
    };
    const pop = i => {
        if (i >= 0 && i < recording.samples) {
            for (let k = 0; k < 3; k++) {
                sum[k] -= recording.accelerometer[3 * i + k];
            }
        }
    };
    for (let i = first; i <= last; i++) {
        while (high < i + half) {
            push(++high);
        }
        while (low < i - half) {
            pop(low++);
        }
        const j = i - first;
        time[j] = recording.time[i];
        level[j] = Math.acos(clamp(dot(sum, thrust) / norm(sum), -1, 1)) * DEGREES;
        if (recording.reference.valid[i]) {
            truth[j] = tiltDegrees(qAt(recording.reference.attitude, i), thrust);
        }
    }
    return {time, truth, level, from, to};
}

class LevelFigure extends StoryFigure {
    build() {
        addChoice(this, 'scene', {
            label: '',
            options: [{value: 'flight', label: 'In flight'}, {value: 'hand', label: 'Carried by hand'}],
            value: 'flight',
        });
        this.readout = addReadout(this);
        this.params.t = 30;
        this.scenes = {};
        this.scrub('t', () => this.axes);
    }

    ready({recording, estimator}) {
        const constants = estimator.constants;
        this.thrust = [constants.thrust_axis_x, constants.thrust_axis_y, constants.thrust_axis_z];
        this.scenes.flight = tilts(recording, this.thrust, 0, recording.duration);
    }

    changed() {
        if (this.params.scene === 'hand' && !this.scenes.hand && !this.loadingHand && this.thrust) {
            this.loadingHand = true;
            loadRun(HAND_RECORDING).then(({recording}) => {
                this.scenes.hand = tilts(recording, this.thrust, 0, HAND_SECONDS);
                this.invalidate();
            });
        }
        const scene = this.scenes[this.params.scene];
        if (scene) {
            this.params.t = clamp(this.params.t, scene.from, scene.to);
        }
    }

    paint(ctx, width, height) {
        const c = this.colors;
        const scene = this.scenes[this.params.scene];
        const axes = new Axes(0, 1, 0, 35).place(46, 56, width - 60, height - 84);
        if (!scene) {
            this.axes = null;
            axes.grid(ctx, {xTicks: 0, yLabel: 'Tilt (degrees)'});
            return;
        }
        axes.setRange(scene.from, scene.to, 0, 35);
        this.axes = axes;
        axes.grid(ctx, {yLabel: 'Tilt (degrees)', xFormat: seconds, yTicks: 4});
        ctx.save();
        axes.clip(ctx);
        axes.line(ctx, scene.time, scene.truth, {color: alpha(c.ink, this.emphasis('truth'))});
        axes.line(ctx, scene.time, scene.level, {color: alpha(c.orange, this.emphasis('level'))});
        ctx.restore();
        legend(ctx, [
            {label: 'true tilt (reference)', color: c.ink},
            {label: 'tilt shown by the accelerometer', color: c.orange},
        ], axes.left, 11);

        const t = clamp(this.params.t, scene.from, scene.to);
        const i = clamp(Math.round((t - scene.from) / (scene.to - scene.from) * (scene.time.length - 1)), 0, scene.time.length - 1);
        axes.cursor(ctx, t);
        if (Number.isFinite(scene.truth[i])) {
            axes.dot(ctx, t, Math.min(scene.truth[i], 35), c.ink);
        }
        axes.dot(ctx, t, Math.min(scene.level[i], 35), c.orange);
        const parts = [`${t.toFixed(1)} s`];
        if (Number.isFinite(scene.truth[i])) {
            parts.push(`true tilt ${scene.truth[i].toFixed(1)}°`);
        }
        parts.push(`accelerometer shows ${scene.level[i].toFixed(1)}°`);
        this.readout.textContent = parts.join('  ·  ');
    }
}

registerFigure('level', LevelFigure);
