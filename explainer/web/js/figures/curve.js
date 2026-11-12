/*
A belief pushed through a curve. A quadcopter 2 m above the floor looks at
a marker on the floor; its camera gives the angle between straight down and
the marker. The distance to the marker along the floor is then
height x tan(angle), which is a curve, not a straight line.

The angle is not known exactly (blue bell). What do we then know about the
distance? The figure compares three answers:

    grey bars     the truth: many angles drawn from the bell, each one
                  turned into a distance
    dashed line   pretend the curve is a straight line (its tangent)
    green         the unscented transform: three sigma points

The sigma points use the library's default settings (kappa = 3), which is
what the multirotor filter uses too.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addToggle, addReadout} from '../core/controls.js';
import {Axes, text, marker, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {DEG} from '../core/math.js';
import {drawBell} from './concept-kit.js';

const HEIGHT = 2;               /* metres above the floor */
const KAPPA = 3;                /* the library's default */
const ANGLE_LIMIT = 85;         /* degrees: the camera cannot see further out than this */
const SAMPLES = 6000;
const BINS = 60;
const DISTANCE_TOP = 9;         /* metres at the top of the plot */

const distance = angle => HEIGHT * Math.tan(angle * DEG);

/* Evenly spread points of a standard bell curve (its quantiles), so the "truth" does not flicker. */
const QUANTILES = (() => {
    /* The inverse of the normal distribution function (Acklam's rational approximation). */
    const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239];
    const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
    const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
    const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
    const inverse = p => {
        if (p < 0.02425) {
            const q = Math.sqrt(-2 * Math.log(p));
            return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
        }
        if (p > 1 - 0.02425) {
            return -inverse(1 - p);
        }
        const q = p - 0.5, r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
    };
    return Float64Array.from({length: SAMPLES}, (_, i) => inverse((i + 0.5) / SAMPLES));
})();

class CurveFigure extends Figure {
    setup() {
        addSlider(this, 'angle', {label: 'Angle to the marker', min: 10, max: 60, value: 50, step: 0.5, format: v => `${v.toFixed(0)}°`});
        addSlider(this, 'spread', {label: 'Its uncertainty', min: 1, max: 10, value: 8, step: 0.1, format: v => `${v.toFixed(1)}°`});
        addToggle(this, 'line', {label: 'Show the straight-line answer', value: true});
        this.table = document.createElement('table');
        this.table.className = 'numbers';
        this.controls.append(this.table);
        this.readout = addReadout(this);
    }

    /* The three answers for the current belief. */
    solve() {
        const {angle, spread} = this.params;
        /* The truth: push every sample through the curve. */
        const values = new Float64Array(SAMPLES);
        let sum = 0;
        for (let i = 0; i < SAMPLES; i++) {
            values[i] = distance(Math.min(angle + spread * QUANTILES[i], ANGLE_LIMIT));
            sum += values[i];
        }
        const mean = sum / SAMPLES;
        let squares = 0;
        for (let i = 0; i < SAMPLES; i++) {
            squares += (values[i] - mean) ** 2;
        }
        /* The straight line: the curve's own value and slope at the middle angle. */
        const slope = HEIGHT / Math.cos(angle * DEG) ** 2 * DEG;      /* metres per degree */
        /* Sigma points: the middle, and one each side at sqrt(1 + kappa) spreads. */
        const reach = Math.sqrt(1 + KAPPA) * spread;
        const points = [angle, angle + reach, angle - reach].map(a => Math.min(a, ANGLE_LIMIT));
        const weights = [KAPPA / (1 + KAPPA), 0.5 / (1 + KAPPA), 0.5 / (1 + KAPPA)];
        const moved = points.map(distance);
        const sigmaMean = moved.reduce((total, value, i) => total + weights[i] * value, 0);
        const sigmaVariance = moved.reduce((total, value, i) => total + weights[i] * (value - sigmaMean) ** 2, 0);
        return {
            values,
            truth: {mean, sd: Math.sqrt(squares / SAMPLES)},
            line: {mean: distance(angle), sd: slope * spread},
            sigma: {mean: sigmaMean, sd: Math.sqrt(sigmaVariance), points, moved},
        };
    }

    draw(ctx, width, height) {
        const c = this.colors;
        const {angle, spread} = this.params;
        const result = this.solve();
        const showLine = this.params.line;

        /* The plot, with room under it for the belief and beside it for the answers. */
        const side = Math.max(96, Math.min(170, width * 0.2));
        const under = Math.min(74, height * 0.18);
        const axes = new Axes(0, 90, 0, DISTANCE_TOP).place(46, 50, width - 46 - side - 16, height - 50 - under - 30);
        axes.grid(ctx, {xTicks: 6, yTicks: 4, xFormat: v => `${v}°`, yLabel: 'Distance to the marker along the floor (m)'});
        legend(ctx, [
            {label: 'belief about the angle', color: c.blue},
            {label: 'sigma points', color: c.green, dot: true},
            {label: 'truth', color: c.muted, band: true},
        ].concat(showLine ? [{label: 'straight line', color: c.ink2, dash: [5, 4], width: 1.5}] : []), axes.left, 11);

        /* The curve, and its tangent at the middle angle. */
        ctx.save();
        axes.clip(ctx);
        ctx.strokeStyle = alpha(c.ink, this.emphasis('curve'));
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let a = 0; a <= 88; a += 0.5) {
            const px = axes.x(a), py = axes.y(distance(a));
            if (a === 0) {
                ctx.moveTo(px, py);
            } else {
                ctx.lineTo(px, py);
            }
        }
        ctx.stroke();
        if (showLine) {
            const slope = result.line.sd / spread;
            ctx.strokeStyle = alpha(c.ink2, 0.9 * this.emphasis('line'));
            ctx.lineWidth = 1.5;
            ctx.setLineDash([5, 4]);
            ctx.beginPath();
            ctx.moveTo(axes.x(angle - 40), axes.y(result.line.mean - 40 * slope));
            ctx.lineTo(axes.x(angle + 40), axes.y(result.line.mean + 40 * slope));
            ctx.stroke();
            ctx.setLineDash([]);
        }
        /* Each sigma point: up from the angle to the curve, then across to the distance. */
        ctx.strokeStyle = alpha(c.green, 0.55 * this.emphasis('sigma'));
        ctx.lineWidth = 1;
        result.sigma.points.forEach((point, i) => {
            const px = axes.x(point), py = axes.y(Math.min(result.sigma.moved[i], DISTANCE_TOP * 2));
            ctx.beginPath();
            ctx.moveTo(px, axes.top + axes.height);
            ctx.lineTo(px, py);
            ctx.lineTo(axes.left + axes.width, py);
            ctx.stroke();
        });
        ctx.restore();

        /* Under the plot: the belief about the angle, and the sigma points on it. */
        const base = axes.top + axes.height + 24 + under;
        drawBell(ctx, {mean: angle, sd: spread, along: v => axes.x(v), base, peak: under - 6, color: c.blue,
            opacity: this.emphasis('belief'), from: 0, to: 90});
        ctx.strokeStyle = c.axis;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(axes.left, base + 0.5);
        ctx.lineTo(axes.left + axes.width, base + 0.5);
        ctx.stroke();
        result.sigma.points.forEach((point, i) => {
            marker(ctx, axes.x(point), base, alpha(c.blue, this.emphasis('sigma')), i === 0 ? 5 : 4.5);
            marker(ctx, axes.x(point), axes.y(Math.min(result.sigma.moved[i], DISTANCE_TOP)), alpha(c.green, this.emphasis('sigma')), i === 0 ? 5 : 4.5);
        });

        /* Beside the plot: the three answers about the distance. */
        const edge = axes.left + axes.width + 14;
        const counts = new Float64Array(BINS);
        for (const value of result.values) {
            const bin = Math.floor(value / DISTANCE_TOP * BINS);
            if (bin >= 0 && bin < BINS) {
                counts[bin]++;
            }
        }
        const tallest = Math.max(...counts);
        const barHeight = axes.height / BINS;
        ctx.fillStyle = alpha(c.muted, 0.42 * this.emphasis('truth'));
        for (let bin = 0; bin < BINS; bin++) {
            const length = counts[bin] / tallest * (side - 24);
            if (length > 0.3) {
                ctx.fillRect(edge, axes.y((bin + 1) / BINS * DISTANCE_TOP) + 0.5, length, Math.max(barHeight - 1, 1));
            }
        }
        /* Bells with the same area as the bars, so that their heights can be compared. */
        const binWidth = DISTANCE_TOP / BINS;
        const peakOf = sd => (side - 24) * (SAMPLES * binWidth / (sd * Math.sqrt(2 * Math.PI))) / tallest;
        ctx.save();
        ctx.beginPath();
        ctx.rect(edge - 2, axes.top, side, axes.height);
        ctx.clip();
        if (showLine && result.line.sd > 1e-6) {
            drawBell(ctx, {mean: result.line.mean, sd: result.line.sd, along: v => axes.y(v), base: edge, peak: Math.min(peakOf(result.line.sd), side * 3),
                vertical: true, color: c.ink2, fill: 0, width: 1.5, dash: [5, 4], opacity: this.emphasis('line')});
        }
        drawBell(ctx, {mean: result.sigma.mean, sd: result.sigma.sd, along: v => axes.y(v), base: edge, peak: Math.min(peakOf(result.sigma.sd), side * 3),
            vertical: true, color: c.green, fill: 0.12, opacity: this.emphasis('sigma')});
        ctx.restore();
        ctx.strokeStyle = c.axis;
        ctx.beginPath();
        ctx.moveTo(edge - 0.5, axes.top);
        ctx.lineTo(edge - 0.5, axes.top + axes.height);
        ctx.stroke();
        text(ctx, 'how likely', edge + 2, axes.top + axes.height + 6, {baseline: 'top', color: c.muted, size: 11.5});

        /* The numbers. */
        const cell = value => `${value.mean.toFixed(2)} m ± ${value.sd.toFixed(2)} m`;
        const rows = [['', 'Distance'], ['Truth', cell(result.truth)], ['Sigma points', cell(result.sigma)]];
        if (showLine) {
            rows.push(['Straight line', cell(result.line)]);
        }
        this.table.replaceChildren(...rows.map((cells, r) => {
            const row = document.createElement('tr');
            cells.forEach((value, k) => {
                const node = document.createElement(r === 0 || k === 0 ? 'th' : 'td');
                node.textContent = value;
                row.append(node);
            });
            return row;
        }));
    }
}

registerFigure('curve', CurveFigure);
