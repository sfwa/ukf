/*
How a filter learns a thing which it never measures. The belief here has
two parts: where the quadcopter is, and how fast it moves. Each dot is one
possibility, and the ellipse holds most of them.

Waiting moves every dot sideways by its own speed, so the fast ones get
ahead and the cloud leans over: position and speed are now tied together.
A position reading then keeps the dots near the reading. Because the cloud
leans, the dots which are left also agree on the speed. The reading
measured only position, but the speed is better known too.

The dots, the true state and the random part of the reading are fixed, so
that the controls change only the time waited and the sensor.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addToggle, addReadout} from '../core/controls.js';
import {Axes, text, legend, marker, ellipse} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {Random, covarianceEllipse} from '../core/math.js';
import {drawBracket} from './concept-kit.js';

/* The belief at the start: position and speed, each with its own spread, not yet tied together. */
const START = {position: 2, speed: 1, positionSpread: 0.4, speedSpread: 0.45};
/* The real quadcopter is one possibility among the others. */
const TRUTH = {position: 2.25, speed: 1.4};
const READING_DRAW = 0.4;       /* the random part of the reading, in units of the sensor's spread */
const DOTS = 260;
const SIGMAS = 2;               /* the ellipses hold the belief out to two spreads */

class HiddenFigure extends Figure {
    setup() {
        const random = new Random(31);
        this.dots = Array.from({length: DOTS}, () => ({position: random.normal(), speed: random.normal()}));
        addSlider(this, 'wait', {label: 'Time waited', min: 0, max: 3, value: 0, step: 0.02,
            format: v => `${v.toFixed(1)} s`});
        addToggle(this, 'reading', {label: 'Take a position reading', value: false});
        addSlider(this, 'sensorNoise', {label: 'Sensor noise', min: 0.05, max: 1, value: 0.2, step: 0.01,
            format: v => `${v.toFixed(2)} m`});
        this.readout = addReadout(this);
    }

    /* The belief after waiting, and after the reading: centres, covariances and the reading itself. */
    get beliefs() {
        const t = this.params.wait, noise = this.params.sensorNoise;
        const vv = START.speedSpread ** 2;
        const predicted = {
            position: START.position + START.speed * t,
            speed: START.speed,
            pp: START.positionSpread ** 2 + vv * t * t,
            pv: vv * t,
            vv,
        };
        const truth = {position: TRUTH.position + TRUTH.speed * t, speed: TRUTH.speed};
        const value = truth.position + READING_DRAW * noise;
        const s = predicted.pp + noise * noise;
        const gainPosition = predicted.pp / s, gainSpeed = predicted.pv / s;
        const corrected = {
            position: predicted.position + gainPosition * (value - predicted.position),
            speed: predicted.speed + gainSpeed * (value - predicted.position),
            pp: predicted.pp - gainPosition * predicted.pp,
            pv: predicted.pv - gainPosition * predicted.pv,
            vv: predicted.vv - gainSpeed * predicted.pv,
        };
        return {predicted, corrected, truth, value};
    }

    draw(ctx, width, height) {
        const c = this.colors;
        const {wait, reading, sensorNoise} = this.params;
        const {predicted, corrected, truth, value} = this.beliefs;
        const narrow = width < 560;
        const left = narrow ? 38 : 50, right = width - 14, top = 44, bottom = height - (narrow ? 40 : 44);
        const axes = new Axes(0, 8.5, -0.4, 2.6).place(left, top, right - left, bottom - top);
        axes.grid(ctx, {xTicks: narrow ? 4 : 8, yTicks: 5, xFormat: v => `${v} m`, yLabel: 'speed (m/s)',
            xLabel: 'position'});

        const outline = (belief, color, key, options = {}) => {
            const shape = covarianceEllipse(belief.pp, belief.pv, belief.vv, SIGMAS);
            /* The ellipse is worked out in metres and m/s; the two axes have different scales on screen. */
            const sx = axes.width / (axes.x1 - axes.x0), sy = axes.height / (axes.y1 - axes.y0);
            ctx.save();
            ctx.translate(axes.x(belief.position), axes.y(belief.speed));
            ctx.scale(sx, -sy);
            ctx.beginPath();
            ctx.ellipse(0, 0, Math.max(shape.major, 1e-4), Math.max(shape.minor, 1e-4), shape.angle, 0, 2 * Math.PI);
            ctx.restore();
            const on = this.emphasis(key);
            if (options.fill) {
                ctx.fillStyle = alpha(color, options.fill * on);
                ctx.fill();
            }
            ctx.strokeStyle = alpha(color, on * (options.opacity ?? 1));
            ctx.lineWidth = 2;
            if (options.dash) {
                ctx.setLineDash(options.dash);
            }
            ctx.stroke();
            ctx.setLineDash([]);
        };

        ctx.save();
        axes.clip(ctx);

        /* The reading: a band of positions which the sensor allows. */
        if (reading) {
            const on = this.emphasis('reading');
            ctx.fillStyle = alpha(c.orange, 0.13 * on);
            ctx.fillRect(axes.x(value - 2 * sensorNoise), axes.top, axes.x(value + 2 * sensorNoise) - axes.x(value - 2 * sensorNoise), axes.height);
            ctx.strokeStyle = alpha(c.orange, on);
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(axes.x(value), axes.top);
            ctx.lineTo(axes.x(value), axes.top + axes.height);
            ctx.stroke();
        }

        /* The possibilities. With a reading, those which the reading makes unlikely fade away. */
        const predictionOn = this.emphasis('prediction');
        for (const dot of this.dots) {
            const speed = START.speed + START.speedSpread * dot.speed;
            const position = START.position + START.positionSpread * dot.position + speed * wait;
            const kept = reading ? Math.exp(-0.5 * ((position - value) / sensorNoise) ** 2) : 1;
            const shade = reading ? c.blue : (wait > 0.01 ? c.green : c.blue);
            ctx.fillStyle = alpha(reading && kept < 0.5 ? c.green : shade, (reading ? 0.12 + 0.6 * kept : 0.5) * predictionOn);
            ctx.beginPath();
            ctx.arc(axes.x(position), axes.y(speed), narrow ? 2 : 2.5, 0, 2 * Math.PI);
            ctx.fill();
        }

        /* The belief at the start, after waiting, and after the reading. */
        const start = {position: START.position, speed: START.speed, pp: START.positionSpread ** 2, pv: 0, vv: START.speedSpread ** 2};
        if (wait > 0.01) {
            outline(start, c.blue, 'start', {opacity: 0.55, dash: [5, 4]});
            outline(predicted, c.green, 'prediction', {fill: reading ? 0 : 0.08});
        } else {
            outline(start, c.blue, 'start', {fill: 0.08});
        }
        if (reading) {
            outline(corrected, c.blue, 'estimate', {fill: 0.16});
        }

        /* Where the quadcopter really is, and how fast it really moves. */
        const tx = axes.x(truth.position), ty = axes.y(truth.speed);
        ctx.strokeStyle = alpha(c.surface, this.emphasis('truth'));
        ctx.lineWidth = 5;
        ctx.lineCap = 'round';
        for (const stroke of [c.surface, c.ink]) {
            ctx.strokeStyle = alpha(stroke, this.emphasis('truth'));
            ctx.lineWidth = stroke === c.ink ? 2 : 5;
            ctx.beginPath();
            ctx.moveTo(tx - 5, ty - 5);
            ctx.lineTo(tx + 5, ty + 5);
            ctx.moveTo(tx - 5, ty + 5);
            ctx.lineTo(tx + 5, ty - 5);
            ctx.stroke();
        }
        ctx.restore();

        /* How well the speed is known, as a bracket beside the speed axis. */
        const now = reading ? corrected : predicted;
        const speedSpread = Math.sqrt(now.vv);
        const bracketX = right - 10;
        drawBracket(ctx, bracketX, axes.y(now.speed - speedSpread), bracketX, axes.y(now.speed + speedSpread),
            alpha(reading ? c.blue : (wait > 0.01 ? c.green : c.blue), 1), {width: 2});
        text(ctx, `speed ± ${speedSpread.toFixed(2)} m/s`, bracketX - 9, axes.y(now.speed + speedSpread) - 7,
            {align: 'right', color: c.ink2, size: 12, halo: true});
        if (reading) {
            text(ctx, 'reading', axes.x(value) + 6, axes.top + 13, {color: c.ink2, size: 12, halo: true});
        }

        const entries = [{label: narrow ? 'start' : 'belief at the start', color: c.blue, dash: wait > 0.01 ? [4, 3] : null}];
        if (wait > 0.01) {
            entries.push({label: narrow ? 'waited' : 'after waiting', color: c.green});
        }
        if (reading) {
            entries.push({label: narrow ? 'reading' : 'position reading', color: c.orange});
            entries.push({label: narrow ? 'estimate' : 'after the reading', color: c.blue});
        }
        legend(ctx, entries, left, 14, {size: narrow ? 11.5 : 12.5});
        text(ctx, '×', narrow ? right - 44 : right - 118, 15, {color: c.ink, size: 15, weight: 700, baseline: 'middle'});
        text(ctx, narrow ? 'truth' : 'the true state', narrow ? right - 32 : right - 104, 14, {color: c.ink2, size: narrow ? 11.5 : 12.5, baseline: 'middle'});

        const before = `position ± ${Math.sqrt(predicted.pp).toFixed(2)} m, speed ± ${Math.sqrt(predicted.vv).toFixed(2)} m/s`;
        this.readout.textContent = reading ?
            `before the reading: ${before}  ·  after: position ± ${Math.sqrt(corrected.pp).toFixed(2)} m, ` +
            `speed ± ${Math.sqrt(corrected.vv).toFixed(2)} m/s` :
            `the belief now: ${before}`;
    }
}

registerFigure('hidden', HiddenFigure);
