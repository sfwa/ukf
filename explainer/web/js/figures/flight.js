/*
The opening figure: a recorded flight, replayed. The aircraft's attitude is
the estimator's output, computed in the browser from the IMU samples alone.
Its place in the room is taken from the reference system, because the
estimator does not estimate position.

The small dial shows the tilt: where "up" is as seen from the aircraft,
for the estimate (blue dot) and for the reference (ring).
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addPlayBar, addReadout} from '../core/controls.js';
import {Camera} from '../core/scene.js';
import {Renderer} from '../core/gl.js';
import {Vehicle, VEHICLE_FOOT} from '../core/models.js';
import {alpha} from '../core/theme.js';
import {text, marker} from '../core/plot.js';
import {clamp, qMul, qConj, qRotate, qAt, vAt, dot, DEG} from '../core/math.js';
import {loadStory} from '../data/story.js';

const PATH_EVERY = 8;

class FlightFigure extends Figure {
    setup() {
        this.camera = new Camera({distance: 2.15, azimuth: -125 * DEG, elevation: 20 * DEG, fov: 34 * DEG});
        this.renderer = Renderer.get();
        this.story = null;
        this.failure = null;
        addPlayBar(this, 't', {duration: 1});
        this.readout = addReadout(this);
        this.canvas.classList.add('grab');
        this.onPointer({
            down: () => this.canvas.classList.add('grabbing'),
            drag: (at, delta) => {
                this.camera.orbit(delta.dx, delta.dy);
                this.invalidate();
            },
            up: () => this.canvas.classList.remove('grabbing'),
        });
        loadStory().then(story => this.begin(story)).catch(error => {
            this.failure = error.message;
            this.invalidate();
        });
    }

    begin(story) {
        const {recording, run, estimator} = story;
        this.story = story;
        this.vehicle = this.renderer ? new Vehicle(this.renderer, estimator.constants) : null;
        this.inputs.t.setDuration(recording.duration);

        /* The fixed turn about the vertical which lines the estimate's heading up with the reference. */
        const offset = run.errors.offset;
        this.headingFix = [Math.cos(offset / 2), 0, 0, -Math.sin(offset / 2)];

        const constants = estimator.constants;
        this.right = [constants.lateral_axis_1_x, constants.lateral_axis_1_y, constants.lateral_axis_1_z];
        this.forward = [constants.lateral_axis_2_x, constants.lateral_axis_2_y, constants.lateral_axis_2_z];

        /* The path flown, thinned out for drawing, and the lowest point reached. */
        const reference = recording.reference;
        this.path = [];
        let lowest = Infinity;
        for (let i = 0; i < recording.samples; i += PATH_EVERY) {
            if (reference.valid[i]) {
                this.path.push({index: i, point: vAt(reference.position, i)});
                lowest = Math.min(lowest, reference.position[3 * i + 2]);
            }
        }
        this.floor = lowest - VEHICLE_FOOT;

        this.invalidate();
        if (!window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
            this.set({tPlaying: true}, true);
            this.animate();
        }
    }

    step(dt) {
        if (!this.story || !this.params.tPlaying) {
            return false;
        }
        const duration = this.story.recording.duration;
        const t = this.params.t + dt;
        if (t >= duration) {
            this.set({t: duration, tPlaying: false}, true);
            return false;
        }
        this.set({t}, true);
        return true;
    }

    /* The sample nearest to the figure's time. */
    get sample() {
        const recording = this.story.recording;
        return clamp(Math.round(this.params.t / recording.duration * (recording.samples - 1)), 0, recording.samples - 1);
    }

    /* "Up" as seen from the aircraft, as its rightward and forward parts. */
    tilt(attitude) {
        const up = qRotate(qConj(attitude), [0, 0, 1]);
        return [dot(up, this.right), dot(up, this.forward)];
    }

    draw(ctx, width, height) {
        const c = this.colors;
        if (this.failure || !this.story) {
            text(ctx, this.failure ?? 'Loading the recording…', width / 2, height / 2,
                {align: 'center', baseline: 'middle', color: c.muted, size: 14});
            return;
        }
        const {recording, run} = this.story;
        const reference = recording.reference;
        const i = this.sample;
        const position = vAt(reference.position, i);
        const estimate = qMul(this.headingFix, run.attitude(Math.min(i, run.valid - 1)));
        const flying = run.value('rotor_vibration', i) > 0.65;

        const camera = this.camera;
        camera.target = [position[0], position[1], position[2] - 0.03];
        camera.update(width, height);

        /* The floor, the path and the shadow go under the aircraft. */
        camera.grid(ctx, {centre: position, z: this.floor, extent: 4.5, step: 0.5, color: c.grid});
        ctx.save();
        ctx.lineWidth = 1.5;
        ctx.lineJoin = 'round';
        ctx.strokeStyle = alpha(c.ink2, 0.2);
        ctx.beginPath();
        camera.path(ctx, this.path.length, k => this.path[k].point);
        ctx.stroke();
        ctx.strokeStyle = alpha(c.ink2, 0.85);
        ctx.lineWidth = 2;
        ctx.beginPath();
        let flown = 0;
        while (flown < this.path.length && this.path[flown].index <= i) {
            flown++;
        }
        camera.path(ctx, flown + 1, k => (k < flown ? this.path[k].point : position));
        ctx.stroke();
        ctx.restore();
        const lift = position[2] - VEHICLE_FOOT - this.floor;
        camera.shadow(ctx, position, this.floor, 0.36 + 0.1 * Math.min(lift, 2), c.dark ? [0, 0, 0] : c.inkRgb,
            (c.dark ? 0.55 : 0.26) / (1 + 0.7 * lift));

        if (this.vehicle) {
            this.renderer.begin(width, height, this.dpr, camera.viewProjection);
            this.vehicle.draw(this.renderer, c, position, estimate, {spinning: flying});
            this.renderer.finish(ctx, width, height);
        }

        this.drawDial(ctx, width, height, estimate, reference.valid[i] ? qAt(reference.attitude, i) : null);

        const error = run.errors.inclination[i];
        const parts = [flying ? 'rotors on' : 'rotors off'];
        if (Number.isFinite(error)) {
            parts.push(`tilt error ${error.toFixed(1)}°`);
        }
        this.readout.textContent = parts.join('  ·  ');
    }

    /* The tilt dial in the corner. */
    drawDial(ctx, width, height, estimate, reference) {
        const c = this.colors;
        const radius = clamp(width * 0.085, 44, 70);
        const cx = width - radius - 18, cy = radius + 30;
        const range = 25 * DEG;     /* the edge of the dial */
        ctx.save();
        ctx.fillStyle = alpha(c.surface, 0.82);
        ctx.beginPath();
        ctx.arc(cx, cy, radius + 6, 0, 2 * Math.PI);
        ctx.fill();
        ctx.strokeStyle = c.grid;
        ctx.lineWidth = 1;
        for (const ring of [10, 20]) {
            ctx.beginPath();
            ctx.arc(cx, cy, radius * ring * DEG / range, 0, 2 * Math.PI);
            ctx.stroke();
        }
        ctx.strokeStyle = c.axis;
        ctx.beginPath();
        ctx.arc(cx, cy, radius, 0, 2 * Math.PI);
        ctx.moveTo(cx - radius, cy);
        ctx.lineTo(cx + radius, cy);
        ctx.moveTo(cx, cy - radius);
        ctx.lineTo(cx, cy + radius);
        ctx.stroke();
        ctx.restore();

        const place = attitude => {
            const [right, forward] = this.tilt(attitude);
            const x = Math.asin(clamp(right, -1, 1)) / range, y = Math.asin(clamp(forward, -1, 1)) / range;
            const length = Math.hypot(x, y);
            const limit = length > 1 ? 1 / length : 1;
            return [cx + x * limit * radius, cy - y * limit * radius];
        };
        if (reference) {
            const [x, y] = place(reference);
            ctx.save();
            ctx.strokeStyle = alpha(c.ink, this.emphasis('reference'));
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.arc(x, y, 8, 0, 2 * Math.PI);
            ctx.stroke();
            ctx.restore();
        }
        const [x, y] = place(estimate);
        marker(ctx, x, y, alpha(c.blue, this.emphasis('estimate')), 4.5, 0);
        text(ctx, 'tilt', cx, cy - radius - 12, {align: 'center', baseline: 'bottom', color: c.ink2, size: 12, weight: 500});
        text(ctx, '20°', cx + radius * 20 * DEG / range + 3, cy + 3, {baseline: 'top', color: c.muted, size: 10.5});
    }
}

registerFigure('flight', FlightFigure);
