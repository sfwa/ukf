/*
Two guesses about one position, and how to join them. The model's
prediction and the sensor's reading are each a bell curve: a centre and a
spread. Multiplying the two curves gives the estimate: a third bell, which
lies between them, nearer to the one with the smaller spread, and is
narrower than either.

The share of the way from the prediction to the reading is the gain.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addRow, addReadout} from '../core/controls.js';
import {Axes, text, legend, marker} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {clamp} from '../core/math.js';
import {drawBell} from './concept-kit.js';

const TRACK = 10;       /* the position axis runs from 0 to this many metres */
const GRAB = 30;        /* a handle can be taken from this far away (pixels) */

class FuseFigure extends Figure {
    setup() {
        const metres = v => `${v.toFixed(2)} m`;
        const first = addRow(this), second = addRow(this);
        addSlider(this, 'prediction', {label: 'Prediction', min: 0.5, max: TRACK - 0.5, value: 3.5, step: 0.05,
            format: metres, parent: first});
        addSlider(this, 'predictionSpread', {label: 'its spread', min: 0.2, max: 2.5, value: 1.2, step: 0.05,
            format: metres, parent: first});
        addSlider(this, 'reading', {label: 'Reading', min: 0.5, max: TRACK - 0.5, value: 6.5, step: 0.05,
            format: metres, parent: second});
        addSlider(this, 'sensorSpread', {label: 'its spread', min: 0.2, max: 2.5, value: 0.6, step: 0.05,
            format: metres, parent: second});
        /* Two sliders share a row where there is room for both. */
        for (const name of ['prediction', 'predictionSpread', 'reading', 'sensorSpread']) {
            this.inputs[name].row.style.flex = '1 1 17rem';
            this.inputs[name].row.style.gridTemplateColumns = 'minmax(4.6rem, auto) 1fr minmax(3.6rem, auto)';
        }
        this.readout = addReadout(this);
        this.axes = new Axes(0, TRACK, 0, 1);
        this.dragged = false;
        this.held = null;

        const nearest = at => {
            let name = null, distance = GRAB;
            for (const candidate of ['prediction', 'reading']) {
                const gap = Math.abs(at.x - this.axes.x(this.params[candidate]));
                if (gap < distance) {
                    distance = gap;
                    name = candidate;
                }
            }
            return name;
        };
        this.onPointer({
            down: at => {
                this.held = nearest(at);
                if (this.held === null) {
                    return false;
                }
                this.canvas.classList.add('grabbing');
                return true;
            },
            drag: at => {
                this.dragged = true;
                this.transitions.delete(this.held);
                this.set({[this.held]: clamp(this.axes.fromX(at.x), 0.5, TRACK - 0.5)}, true);
            },
            up: () => {
                this.held = null;
                this.canvas.classList.remove('grabbing');
            },
            move: at => this.canvas.classList.toggle('grab', nearest(at) !== null),
            leave: () => this.canvas.classList.remove('grab'),
        });
    }

    /* The estimate which the prediction and the reading together support. */
    get estimate() {
        const {prediction, reading, predictionSpread, sensorSpread} = this.params;
        const gain = predictionSpread ** 2 / (predictionSpread ** 2 + sensorSpread ** 2);
        return {
            gain,
            centre: prediction + gain * (reading - prediction),
            spread: Math.sqrt((1 - gain) * predictionSpread ** 2),
        };
    }

    draw(ctx, width, height) {
        const c = this.colors;
        const {prediction, reading, predictionSpread, sensorSpread} = this.params;
        const {gain, centre, spread} = this.estimate;
        const narrow = width < 560;
        const left = 16, right = width - 16, top = 40;
        const rail = height - (narrow ? 44 : 48);       /* the line on which the handles sit */
        const base = rail - 46;                         /* the position axis */
        const axes = this.axes.place(left, top, right - left, base - top);
        const x = value => axes.x(value);

        axes.grid(ctx, {xTicks: narrow ? 5 : 10, yTicks: 0, xFormat: v => `${v} m`});
        ctx.strokeStyle = c.axis;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(left, base + 0.5);
        ctx.lineTo(right, base + 0.5);
        ctx.stroke();

        /* The tallest bell (always the estimate) fills the plot; the others keep their true proportions. */
        const unit = (base - top) * 0.94 * spread;
        const bells = [
            {key: 'prediction', centre: prediction, spread: predictionSpread, color: c.green, fill: 0.13},
            {key: 'reading', centre: reading, spread: sensorSpread, color: c.orange, fill: 0.13},
            {key: 'estimate', centre, spread, color: c.blue, fill: 0.2},
        ];
        for (const bell of bells) {
            const on = this.emphasis(bell.key);
            drawBell(ctx, {mean: bell.centre, sd: bell.spread, along: x, base, peak: unit / bell.spread,
                color: bell.color, fill: bell.fill, opacity: on, from: 0, to: TRACK});
            /* A thin line ties each bell to its mark on the rail. */
            ctx.strokeStyle = alpha(bell.color, 0.45 * on);
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.moveTo(x(bell.centre), base - unit / bell.spread);
            ctx.lineTo(x(bell.centre), rail);
            ctx.stroke();
        }

        /* The rail: from the prediction to the reading, with the estimate part of the way along. */
        const gainOn = this.emphasis('gain');
        ctx.lineCap = 'round';
        ctx.strokeStyle = alpha(c.axis, gainOn);
        ctx.lineWidth = 4;
        ctx.beginPath();
        ctx.moveTo(x(prediction), rail);
        ctx.lineTo(x(reading), rail);
        ctx.stroke();
        ctx.strokeStyle = alpha(c.ink2, gainOn);
        ctx.beginPath();
        ctx.moveTo(x(prediction), rail);
        ctx.lineTo(x(centre), rail);
        ctx.stroke();
        marker(ctx, x(prediction), rail, alpha(c.green, this.emphasis('prediction')), 8, 2);
        marker(ctx, x(reading), rail, alpha(c.orange, this.emphasis('reading')), 8, 2);
        marker(ctx, x(centre), rail, alpha(c.blue, this.emphasis('estimate')), 5.5, 2);

        const percent = Math.round(gain * 100);
        const sentence = narrow ? `gain ${gain.toFixed(2)}` :
            `gain ${gain.toFixed(2)}: the estimate is ${percent}% of the way from the prediction to the reading`;
        const middle = clamp((x(prediction) + x(reading)) / 2, left + (narrow ? 40 : 250), right - (narrow ? 40 : 250));
        text(ctx, sentence, middle, rail + 26, {align: 'center', color: alpha(c.ink2, gainOn), size: 12.5});
        if (!this.dragged) {
            const hintAt = prediction < reading ? x(prediction) - 14 : x(prediction) + 14;
            text(ctx, 'drag', hintAt, rail + 4, {align: prediction < reading ? 'right' : 'left', color: c.muted, size: 11.5});
        }

        legend(ctx, [
            {label: narrow ? 'prediction' : 'prediction (model)', color: c.green},
            {label: narrow ? 'reading' : 'reading (sensor)', color: c.orange},
            {label: 'estimate', color: c.blue},
        ], left, 16);

        this.readout.textContent =
            `gain ${gain.toFixed(2)}  ·  estimate ${centre.toFixed(2)} m  ·  spread ${spread.toFixed(2)} m`;
    }
}

registerFigure('fuse', FuseFigure);
