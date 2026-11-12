/*
A sensor reading is a guess. A quadcopter hovers at a fixed height, and a
range sensor under it reports the height again and again. No two readings
agree. Each new reading enters on the left, moves to the right as it ages,
and is counted in the pile on the right, which grows into a bell curve.

The random part of every reading is fixed in advance, so that moving the
noise slider stretches the same readings instead of drawing new ones.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addButton, addRow, addReadout} from '../core/controls.js';
import {Axes, text, legend, marker, arrow} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {Random, clamp} from '../core/math.js';
import {drawQuadSide, quadSensorPoint, drawBell, reducedMotion, meanAndSpread} from './concept-kit.js';

const TRUE_HEIGHT = 1;          /* m */
const MOST_READINGS = 1000;
const STILL_READINGS = 400;     /* shown at once when the reader wants less motion */
const SHOWN_IN_STRIP = 44;      /* how many of the newest readings are drawn as dots */
const BINS_PER_SPREAD = 2.5;    /* the pile's bins are spread / 2.5 wide */

class NoiseFigure extends Figure {
    setup() {
        const random = new Random(20241);
        this.draws = new Float32Array(MOST_READINGS);
        for (let i = 0; i < MOST_READINGS; i++) {
            this.draws[i] = random.normal();
        }
        addSlider(this, 'noise', {label: 'Sensor noise', min: 0.03, max: 0.3, value: 0.15, step: 0.005,
            format: v => `${v.toFixed(2)} m`});
        const row = addRow(this);
        addButton(this, {label: 'Start again', parent: row, onClick: () => this.action('restart')});
        this.readout = addReadout(this, {parent: row});
        /* How many readings have been taken so far; it grows while the figure runs. */
        this.params.count = reducedMotion() ? STILL_READINGS : 0;
        if (!reducedMotion()) {
            this.animate();
        }
    }

    action(name) {
        if (name === 'restart') {
            this.transitions.delete('count');
            this.params.count = 0;
            this.animate();
            this.invalidate();
        }
    }

    step(dt) {
        const count = this.params.count;
        if (count >= MOST_READINGS) {
            return false;
        }
        /* Slow at first, so that single readings can be seen, then faster. */
        const rate = Math.min(5 + 0.3 * count, 80);
        this.params.count = Math.min(MOST_READINGS, count + rate * dt);
        return true;
    }

    draw(ctx, width, height) {
        const c = this.colors;
        const sigma = this.params.noise;
        const count = clamp(this.params.count, 0, MOST_READINGS);
        const taken = Math.floor(count);
        const narrow = width < 560;

        /* Three parts side by side: the scene, the newest readings, the pile. */
        const left = 40, top = 34, bottom = height - 14;
        const sceneWidth = narrow ? 78 : Math.min(200, width * 0.21);
        const pileWidth = narrow ? width * 0.25 : width * 0.27;
        const stripLeft = left + sceneWidth, stripRight = width - pileWidth - (narrow ? 8 : 16);
        const pileLeft = stripRight + (narrow ? 8 : 16), pileRight = width - 6;
        const axes = new Axes(0, 1, -0.22, 2.2).place(left, top, width - left - 6, bottom - top);
        const y = value => axes.y(value);

        /* Height grid and the ground. */
        axes.grid(ctx, {xTicks: 0, yTicks: 4, zero: false, yFormat: v => `${v.toFixed(1)} m`});
        ctx.fillStyle = alpha(c.ink, c.dark ? 0.1 : 0.05);
        ctx.fillRect(left, y(0), width - left - 6, bottom - y(0));
        ctx.strokeStyle = c.axis;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(left, y(0));
        ctx.lineTo(width - 6, y(0));
        ctx.stroke();

        /* The band which holds about two readings in three: one spread each side of the truth. */
        const spreadOn = this.emphasis('spread');
        ctx.fillStyle = alpha(c.orange, 0.1 * spreadOn);
        ctx.fillRect(stripLeft, y(TRUE_HEIGHT + sigma), pileRight - stripLeft, y(TRUE_HEIGHT - sigma) - y(TRUE_HEIGHT + sigma));

        /* The true height. */
        ctx.strokeStyle = alpha(c.ink, 0.85 * this.emphasis('truth'));
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(left, y(TRUE_HEIGHT));
        ctx.lineTo(pileRight, y(TRUE_HEIGHT));
        ctx.stroke();

        /* The quadcopter and the beam of its range sensor. */
        const quadX = left + sceneWidth / 2, quadSize = Math.min(sceneWidth * 0.82, 120);
        const [sensorX, sensorY] = quadSensorPoint(quadX, y(TRUE_HEIGHT), quadSize, 0);
        const pulse = count >= MOST_READINGS || taken === 0 ? 0 : 1 - (count - taken);
        ctx.strokeStyle = alpha(c.orange, (0.3 + 0.5 * pulse) * this.emphasis('reading'));
        ctx.lineWidth = 2;
        ctx.setLineDash([2, 5]);
        ctx.lineCap = 'round';
        ctx.beginPath();
        ctx.moveTo(sensorX, sensorY + 6);
        ctx.lineTo(sensorX, y(0) - 2);
        ctx.stroke();
        ctx.setLineDash([]);
        drawQuadSide(ctx, quadX, y(TRUE_HEIGHT) - quadSize * 0.03, quadSize, 0, c, {sensor: true});

        /* The newest readings: each enters on the left and ages to the right. */
        const spacing = (stripRight - stripLeft - 14) / SHOWN_IN_STRIP;
        const readingOn = this.emphasis('reading');
        for (let k = Math.max(0, taken - SHOWN_IN_STRIP - 1); k < taken; k++) {
            const age = count - 1 - k;
            const px = stripLeft + 10 + age * spacing;
            if (px > stripRight - 2) {
                continue;
            }
            const fade = clamp(1 - age / SHOWN_IN_STRIP, 0, 1);
            const value = TRUE_HEIGHT + sigma * this.draws[k];
            /* A new reading grows to full size as it arrives. */
            const grow = k === taken - 1 && count > taken ? clamp((count - taken) * 4, 0.3, 1) : 1;
            marker(ctx, px, y(value), alpha(c.orange, (0.3 + 0.7 * fade) * readingOn), 4 * grow, k === taken - 1 ? 2 : 0);
        }

        /* The pile: how many readings fell at each height, and the bell curve they tend to. */
        const binWidth = sigma / BINS_PER_SPREAD;
        const bins = new Map();
        for (let k = 0; k < taken; k++) {
            const bin = Math.round(this.draws[k] * BINS_PER_SPREAD);
            bins.set(bin, (bins.get(bin) ?? 0) + 1);
        }
        /* A bell's middle bin holds this share of all readings; scale so that it fills the panel. */
        const middleShare = 1 / (BINS_PER_SPREAD * Math.sqrt(2 * Math.PI));
        const fullLength = (pileRight - pileLeft) * 0.9;
        const perReading = fullLength / (Math.max(taken, 40) * middleShare);
        const thickness = Math.max(1, Math.abs(y(binWidth) - y(0)) - 1);
        ctx.fillStyle = alpha(c.orange, 0.5 * readingOn);
        for (const [bin, number] of bins) {
            const centre = y(TRUE_HEIGHT + bin * binWidth);
            ctx.fillRect(pileLeft, centre - thickness / 2, Math.min(number * perReading, pileRight - pileLeft), thickness);
        }
        drawBell(ctx, {mean: TRUE_HEIGHT, sd: sigma, along: y, base: pileLeft, peak: fullLength, vertical: true,
            color: c.orange, fill: 0, opacity: this.emphasis('bell'), sigmas: 3.6});
        ctx.strokeStyle = c.axis;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(pileLeft + 0.5, top);
        ctx.lineTo(pileLeft + 0.5, y(0));
        ctx.stroke();

        /* The spread, marked as a distance from the truth to the edge of the band. */
        const markX = stripLeft + (stripRight - stripLeft) * (narrow ? 0.5 : 0.62);
        if (Math.abs(y(TRUE_HEIGHT + sigma) - y(TRUE_HEIGHT)) > 9) {
            arrow(ctx, markX, y(TRUE_HEIGHT) - 1, markX, y(TRUE_HEIGHT + sigma) + 1, alpha(c.ink2, spreadOn), {width: 1.5, head: 6});
        }
        text(ctx, `spread ${sigma.toFixed(2)} m`, markX + 7, y(TRUE_HEIGHT + sigma) - 6,
            {color: alpha(c.ink2, spreadOn), size: 12, halo: true});

        /* Labels. */
        legend(ctx, [
            {label: 'reading', color: c.orange, dot: true},
            {label: 'true height', color: c.ink, width: 1.5},
        ], left + 2, 14);
        if (!narrow) {
            text(ctx, 'newest', stripLeft + 6, bottom - 6, {color: c.muted, size: 11.5});
            text(ctx, 'older', stripRight - 4, bottom - 6, {color: c.muted, size: 11.5, align: 'right'});
        }
        text(ctx, narrow ? 'count' : 'how often', pileLeft + 4, bottom - 6, {color: c.muted, size: 11.5});

        const heights = new Float32Array(taken);
        for (let k = 0; k < taken; k++) {
            heights[k] = TRUE_HEIGHT + sigma * this.draws[k];
        }
        const {mean, spread} = meanAndSpread(heights);
        this.readout.textContent = taken < 2 ? `${taken} reading${taken === 1 ? '' : 's'}` :
            `${taken} readings  ·  average ${mean.toFixed(3)} m  ·  spread ${spread.toFixed(3)} m`;
    }
}

registerFigure('noise', NoiseFigure);
