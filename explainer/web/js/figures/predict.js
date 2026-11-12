/*
Waiting makes a belief stale. At the start the filter knows quite well
where the quadcopter is. The quadcopter then flies on at about one metre
per second, pushed about by the air, and no new reading comes. The model
can say where the quadcopter probably is now, but it is less and less
sure: the bell moves along and spreads out.

The faint quadcopters are the possibilities which the model allows. The
solid one is where the quadcopter really is; the filter does not know which
of the possibilities that is.

Every path is fixed in advance, so the model-noise slider stretches the
same paths instead of drawing new ones.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addPlayBar, addReadout} from '../core/controls.js';
import {Axes, text, legend, marker} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {Random, clamp, DEG} from '../core/math.js';
import {drawQuadSide, drawBell, drawBracket, reducedMotion} from './concept-kit.js';

const TRACK = 12;               /* length of the position axis (m) */
const START = 2;                /* centre of the belief at the start (m) */
const START_SPREAD = 0.2;       /* its spread (m) */
const SPEED = 1;                /* m/s */
const DURATION = 6;             /* s */
const PATH_STEP = 0.05;         /* s */
const GHOSTS = 30;

/* A random walk which has a spread of sqrt(t) after t seconds. */
function randomWalk(random) {
    const steps = Math.round(DURATION / PATH_STEP);
    const walk = new Float32Array(steps + 1);
    for (let i = 1; i <= steps; i++) {
        walk[i] = walk[i - 1] + Math.sqrt(PATH_STEP) * random.normal();
    }
    return walk;
}

function walkAt(walk, t) {
    const position = clamp(t / PATH_STEP, 0, walk.length - 1);
    const left = Math.min(Math.floor(position), walk.length - 2);
    return walk[left] + (walk[left + 1] - walk[left]) * (position - left);
}

class PredictFigure extends Figure {
    setup() {
        const random = new Random(77);
        this.ghosts = [];
        for (let i = 0; i < GHOSTS; i++) {
            this.ghosts.push({offset: random.normal(), walk: randomWalk(random), lift: random.uniform() * 2 - 1});
        }
        /* The real quadcopter is one more possibility, a little ahead of the centre. */
        this.truth = {offset: 0.5, walk: randomWalk(new Random(5))};
        addPlayBar(this, 't', {duration: DURATION, value: reducedMotion() ? 4 : 0});
        addSlider(this, 'noise', {label: 'Model noise', min: 0.05, max: 1, value: 0.35, step: 0.01,
            format: v => `${v.toFixed(2)} m/√s`});
        this.readout = addReadout(this);
        if (!reducedMotion()) {
            this.set({tPlaying: true}, true);
            this.animate();
        }
    }

    step(dt) {
        if (!this.params.tPlaying) {
            return false;
        }
        const t = this.params.t + dt;
        if (t >= DURATION) {
            this.set({t: DURATION, tPlaying: false}, true);
            return false;
        }
        this.set({t}, true);
        return true;
    }

    /* Where one possibility is after t seconds. */
    place(possibility, t) {
        return START + START_SPREAD * possibility.offset + SPEED * t + this.params.noise * walkAt(possibility.walk, t);
    }

    draw(ctx, width, height) {
        const c = this.colors;
        const t = clamp(this.params.t, 0, DURATION);
        const centre = START + SPEED * t;
        const spread = Math.sqrt(START_SPREAD ** 2 + this.params.noise ** 2 * t);
        const narrow = width < 560;
        const left = 16, right = width - 16;
        const sceneTop = 36, sceneBottom = Math.round(height * 0.4);
        const base = height - (narrow ? 62 : 58);
        const axes = new Axes(0, TRACK, 0, 1).place(left, sceneBottom + 22, right - left, base - sceneBottom - 22);
        const x = value => axes.x(value);

        axes.grid(ctx, {xTicks: narrow ? 6 : 12, yTicks: 0, xFormat: v => `${v} m`});
        ctx.strokeStyle = c.axis;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(left, base + 0.5);
        ctx.lineTo(right, base + 0.5);
        ctx.stroke();

        /*
        The scene, above the same position axis. Upper row: the possibilities
        which the model allows, drawn over each other, so that the picture
        blurs as they spread. Lower row: the real quadcopter.
        */
        const scene = sceneBottom - sceneTop;
        const size = narrow ? 44 : 68;
        const predictionOn = this.emphasis('prediction');
        ctx.save();
        ctx.beginPath();
        ctx.rect(left - 6, 0, right - left + 12, height);
        ctx.clip();
        for (const ghost of this.ghosts) {
            drawQuadSide(ctx, x(this.place(ghost, t)), sceneTop + scene * (0.3 + 0.07 * ghost.lift), size,
                5 * DEG, c, {tint: c.green, opacity: 0.13 * predictionOn});
        }
        const truthAt = this.place(this.truth, t);
        drawQuadSide(ctx, x(truthAt), sceneTop + scene * 0.8, size, 5 * DEG, c, {opacity: this.emphasis('truth')});
        ctx.restore();

        /* The belief at the start, and the prediction now. Both bells have the same area. */
        const unit = (base - axes.top) * 0.96 * START_SPREAD;
        drawBell(ctx, {mean: START, sd: START_SPREAD, along: x, base, peak: unit / START_SPREAD,
            color: c.blue, fill: 0.1, opacity: this.emphasis('start') * (t > 0.02 ? 0.7 : 1), from: 0, to: TRACK});
        if (t > 0.02) {
            drawBell(ctx, {mean: centre, sd: spread, along: x, base, peak: unit / spread,
                color: c.green, fill: 0.16, opacity: predictionOn, from: 0, to: TRACK});
        }

        /* The true position on the axis, and the spread of the prediction as a length. */
        marker(ctx, x(truthAt), base, alpha(c.ink, this.emphasis('truth')), 4.5, 2);
        const mark = base + 30;
        const end = Math.min(x(centre + spread), right);
        drawBracket(ctx, x(centre), mark, end, mark, alpha(c.ink2, predictionOn));
        const label = `spread ${spread.toFixed(2)} m`;
        const roomRight = end + 8 + 86 < right;
        text(ctx, label, roomRight ? end + 8 : x(centre) - 8, mark, {baseline: 'middle', align: roomRight ? 'left' : 'right',
            color: alpha(c.ink2, predictionOn), size: 12});

        legend(ctx, [
            {label: narrow ? 'start' : 'belief at the start', color: c.blue},
            {label: narrow ? 'prediction' : 'prediction now', color: c.green},
            {label: narrow ? 'truth' : 'where it really is', color: c.ink2, dot: true},
        ], left, 16);

        this.readout.textContent =
            `after ${t.toFixed(1)} s  ·  predicted position ${centre.toFixed(1)} m  ·  spread ${spread.toFixed(2)} m`;
    }
}

registerFigure('predict', PredictFigure);
