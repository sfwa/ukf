/*
Sigma points in two dimensions. A ground station measures how far away a
quadcopter is and in which direction. Both numbers are a little uncertain
(left). Where is the quadcopter on the map (right)?

The map position is a curved function of distance and direction, so a
round belief becomes a bent cloud (grey dots). The figure shows two ways to
sum that cloud up as a centre and an ellipse:

    green   the unscented transform: five sigma points, each put on the
            map, then their weighted mean and spread
    dashed  the straight-line answer: the curve replaced by its tangent

The sigma points use the library's default settings (kappa = 3).
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addToggle, addReadout} from '../core/controls.js';
import {Axes, text, marker, ellipse, legend} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {DEG, Random, covarianceEllipse} from '../core/math.js';

const RANGE = 8;                /* metres from the station to the quadcopter */
const BEARING = 55;             /* degrees from the map's x axis */
const KAPPA = 3;
const CLOUD = 1400;

/* The same standard-normal pairs every time, so the cloud only changes when a slider moves. */
const NORMALS = (() => {
    const random = new Random(7);
    return Float64Array.from({length: 2 * CLOUD}, () => random.normal());
})();

const toMap = (range, bearing) => [range * Math.cos(bearing * DEG), range * Math.sin(bearing * DEG)];

/* The weighted mean and covariance of points on the map. */
function summarise(points, weights) {
    const mean = [0, 0];
    points.forEach((point, i) => {
        mean[0] += weights[i] * point[0];
        mean[1] += weights[i] * point[1];
    });
    let xx = 0, xy = 0, yy = 0;
    points.forEach((point, i) => {
        const dx = point[0] - mean[0], dy = point[1] - mean[1];
        xx += weights[i] * dx * dx;
        xy += weights[i] * dx * dy;
        yy += weights[i] * dy * dy;
    });
    return {mean, xx, xy, yy};
}

class Sigma2dFigure extends Figure {
    setup() {
        addSlider(this, 'bearingSpread', {label: 'Uncertainty of the direction', min: 1, max: 30, value: 15, step: 0.25, format: v => `${v.toFixed(0)}°`});
        addSlider(this, 'rangeSpread', {label: 'Uncertainty of the distance', min: 0.05, max: 1.5, value: 0.3, step: 0.01, format: v => `${v.toFixed(2)} m`});
        addToggle(this, 'line', {label: 'Show the straight-line answer', value: true});
        this.readout = addReadout(this);
    }

    solve() {
        const {bearingSpread, rangeSpread} = this.params;
        /* The cloud: the truth, as near as many samples can show it. */
        const cloud = [];
        for (let i = 0; i < CLOUD; i++) {
            cloud.push(toMap(RANGE + rangeSpread * NORMALS[2 * i], BEARING + bearingSpread * NORMALS[2 * i + 1]));
        }
        const truth = summarise(cloud, cloud.map(() => 1 / CLOUD));

        /* Five sigma points: the middle, and one each way along each axis, sqrt(2 + kappa) spreads out. */
        const reach = Math.sqrt(2 + KAPPA);
        const before = [
            [RANGE, BEARING],
            [RANGE + reach * rangeSpread, BEARING], [RANGE, BEARING + reach * bearingSpread],
            [RANGE - reach * rangeSpread, BEARING], [RANGE, BEARING - reach * bearingSpread],
        ];
        const weights = [KAPPA / (2 + KAPPA), ...Array(4).fill(0.5 / (2 + KAPPA))];
        const after = before.map(([range, bearing]) => toMap(range, bearing));
        const sigma = summarise(after, weights);

        /* The straight-line answer: the tangent of the curve at the middle. */
        const cosine = Math.cos(BEARING * DEG), sine = Math.sin(BEARING * DEG);
        const a = rangeSpread, b = RANGE * bearingSpread * DEG;      /* spreads along and across the line of sight */
        const line = {
            mean: toMap(RANGE, BEARING),
            xx: (cosine * a) ** 2 + (sine * b) ** 2,
            xy: cosine * sine * (a * a - b * b),
            yy: (sine * a) ** 2 + (cosine * b) ** 2,
        };
        return {cloud, truth, before, after, sigma, line};
    }

    draw(ctx, width, height) {
        const c = this.colors;
        const {bearingSpread, rangeSpread} = this.params;
        const result = this.solve();
        const showLine = this.params.line;

        /* Two panels: side by side when there is room, else one above the other. */
        const wide = width >= 600;
        const top = 54, bottom = 34;
        const left = wide ? {left: 46, top, width: width * 0.36 - 46, height: height - top - bottom}
            : {left: 46, top, width: width - 60, height: height * 0.36 - top};
        const mapSide = wide ? Math.min(width * 0.64 - 60, height - top - bottom) : Math.min(width - 60, height * 0.64 - 64);
        const right = wide ? {left: width * 0.36 + 44, top} : {left: 46, top: height * 0.36 + 50};

        /* Left: the belief as measured. */
        const measured = new Axes(BEARING - 75, BEARING + 75, RANGE - 4, RANGE + 4).place(left.left, left.top, left.width, left.height);
        measured.grid(ctx, {xTicks: 4, yTicks: 4, xFormat: v => `${v}°`, yLabel: 'Distance (m), against direction', zero: false});
        ctx.save();
        measured.clip(ctx);
        const rx = bearingSpread / 150 * left.width, ry = rangeSpread / 8 * left.height;
        for (const sigmas of [2, 1]) {
            ellipse(ctx, measured.x(BEARING), measured.y(RANGE), rx * sigmas, ry * sigmas, 0,
                {fill: alpha(c.blue, 0.1 * this.emphasis('belief')), stroke: sigmas === 1 ? alpha(c.blue, this.emphasis('belief')) : null});
        }
        result.before.forEach(([range, bearing], i) => {
            marker(ctx, measured.x(bearing), measured.y(range), alpha(c.blue, this.emphasis('sigma')), i === 0 ? 5 : 4.5);
        });
        ctx.restore();

        /* Right: the map. The station is at the origin. */
        const map = new Axes(-3, 9, -1, 11).place(right.left, right.top, mapSide, mapSide);
        map.grid(ctx, {xTicks: 4, yTicks: 4, yLabel: 'Map (m)'});
        ctx.save();
        map.clip(ctx);
        ctx.fillStyle = alpha(c.muted, 0.42 * this.emphasis('cloud'));
        for (const [x, y] of result.cloud) {
            ctx.beginPath();
            ctx.arc(map.x(x), map.y(y), 1.6, 0, 2 * Math.PI);
            ctx.fill();
        }
        const scale = mapSide / 12;         /* pixels per metre */
        const outline = (summary, options) => {
            const shape = covarianceEllipse(summary.xx, summary.xy, summary.yy, 2);
            ellipse(ctx, map.x(summary.mean[0]), map.y(summary.mean[1]), shape.major * scale, shape.minor * scale, shape.angle, options);
        };
        if (showLine) {
            outline(result.line, {stroke: alpha(c.ink2, 0.9 * this.emphasis('line')), width: 1.5, dash: [5, 4]});
            marker(ctx, map.x(result.line.mean[0]), map.y(result.line.mean[1]), alpha(c.ink2, this.emphasis('line')), 3, 0);
        }
        outline(result.sigma, {stroke: alpha(c.green, this.emphasis('sigma')), fill: alpha(c.green, 0.08 * this.emphasis('sigma'))});
        result.after.forEach(([x, y], i) => {
            marker(ctx, map.x(x), map.y(y), alpha(c.green, this.emphasis('sigma')), i === 0 ? 5 : 4.5);
        });
        /* The weighted mean of the sigma points, and the true centre of the cloud. */
        ctx.strokeStyle = alpha(c.green, this.emphasis('sigma'));
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(map.x(result.sigma.mean[0]), map.y(result.sigma.mean[1]), 8, 0, 2 * Math.PI);
        ctx.stroke();
        const tx = map.x(result.truth.mean[0]), ty = map.y(result.truth.mean[1]);
        ctx.strokeStyle = alpha(c.ink, this.emphasis('cloud'));
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(tx - 5, ty - 5);
        ctx.lineTo(tx + 5, ty + 5);
        ctx.moveTo(tx - 5, ty + 5);
        ctx.lineTo(tx + 5, ty - 5);
        ctx.stroke();
        /* The ground station. */
        ctx.fillStyle = c.ink;
        ctx.fillRect(map.x(0) - 5, map.y(0) - 5, 10, 10);
        ctx.restore();
        text(ctx, 'station', map.x(0) + 10, map.y(0) + 4, {baseline: 'top', color: c.ink2, halo: true, size: 12});

        legend(ctx, [
            {label: 'belief', color: c.blue},
            {label: 'sigma points', color: c.green, dot: true},
            {label: 'what really happens', color: c.muted, dot: true},
        ].concat(showLine ? [{label: 'straight line', color: c.ink2, dash: [5, 4], width: 1.5}] : []), 46, 11);

        const miss = summary => Math.hypot(summary.mean[0] - result.truth.mean[0], summary.mean[1] - result.truth.mean[1]);
        this.readout.textContent = `Centre of the cloud (×): the sigma points miss it by ${miss(result.sigma).toFixed(2)} m`
            + (showLine ? `, the straight line by ${miss(result.line).toFixed(2)} m.` : '.');
    }
}

registerFigure('sigma2d', Sigma2dFigure);
