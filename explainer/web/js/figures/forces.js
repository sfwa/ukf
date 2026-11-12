/*
What an accelerometer feels. It does not feel gravity. It feels every other
force on the aircraft: the hand or the ground that holds it up, or in
flight the rotor thrust and the rotor drag.

Side view. Gravity is drawn in ink because the sensor does not feel it;
the forces it does feel are orange; the velocity is blue.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addChoice, addReadout} from '../core/controls.js';
import {text, arrow} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {DEG} from '../core/math.js';

const GRAVITY = 9.80665;
const DRAG_RATE = 0.25;         /* per second, as in the estimator */
const DRAG_DRAWN = 5;           /* the drag arrow is drawn this many times longer, to be visible */

/* A multirotor seen from the side, centred on the origin of the current transform. */
function drawVehicle(ctx, size, c, spinning) {
    const arm = size * 0.5;
    ctx.save();
    ctx.lineCap = 'round';
    ctx.strokeStyle = c.body;
    ctx.fillStyle = c.body;
    ctx.lineWidth = size * 0.045;
    ctx.beginPath();
    ctx.moveTo(-arm, 0);
    ctx.lineTo(arm, 0);
    ctx.stroke();
    /* Body and legs. */
    ctx.beginPath();
    ctx.roundRect(-size * 0.16, -size * 0.07, size * 0.32, size * 0.14, size * 0.04);
    ctx.fill();
    ctx.lineWidth = size * 0.025;
    ctx.beginPath();
    for (const side of [-1, 1]) {
        ctx.moveTo(side * size * 0.1, size * 0.05);
        ctx.lineTo(side * size * 0.2, size * 0.2);
    }
    ctx.stroke();
    /* Motors and rotors. */
    for (const side of [-1, 1]) {
        ctx.fillStyle = c.body2;
        ctx.fillRect(side * arm - size * 0.03, -size * 0.07, size * 0.06, size * 0.07);
        ctx.fillStyle = alpha(c.body3, spinning ? 0.45 : 0.9);
        ctx.beginPath();
        ctx.ellipse(side * arm, -size * 0.085, size * (spinning ? 0.22 : 0.2), size * (spinning ? 0.022 : 0.012), 0, 0, 2 * Math.PI);
        ctx.fill();
    }
    ctx.restore();
}

class ForcesFigure extends Figure {
    setup() {
        addChoice(this, 'scene', {
            label: '',
            options: [{value: 'flight', label: 'In flight'}, {value: 'hand', label: 'Held in the hand'}],
            value: 'flight',
        });
        addSlider(this, 'tilt', {label: 'Tilt', min: -30, max: 30, value: 15, step: 0.5, format: v => `${v.toFixed(0)}°`});
        addSlider(this, 'speed', {label: 'Speed', min: -4, max: 4, value: 2, step: 0.05, format: v => `${v.toFixed(1)} m/s`});
        this.readout = addReadout(this);
    }

    changed() {
        /* A vehicle held in the hand does not move in this picture. */
        this.inputs.speed.input.disabled = this.params.scene === 'hand';
        this.inputs.speed.row.style.opacity = this.params.scene === 'hand' ? 0.4 : 1;
    }

    draw(ctx, width, height) {
        const c = this.colors;
        const flying = this.params.scene === 'flight';
        const tilt = this.params.tilt * DEG;
        const speed = flying ? this.params.speed : 0;
        const cx = width / 2, cy = height * 0.52;
        const size = Math.min(width * 0.36, height * 0.62, 230);
        const unit = Math.min(height * 0.36, 140) / GRAVITY;     /* pixels for 1 m/s^2 */

        /* Body axes on screen: "up" along the thrust, "across" along the rotor plane. */
        const up = [Math.sin(tilt), -Math.cos(tilt)], across = [Math.cos(tilt), Math.sin(tilt)];

        /* The forces on each kilogram of the aircraft, as accelerations. */
        const thrust = flying ? GRAVITY / Math.cos(tilt) : 0;
        const acrossSpeed = speed * Math.cos(tilt);
        const drag = flying ? -DRAG_RATE * acrossSpeed : 0;
        const felt = flying ? [up[0] * thrust + across[0] * drag, up[1] * thrust + across[1] * drag] : [0, -GRAVITY];

        ctx.save();
        ctx.translate(cx, cy);
        ctx.rotate(tilt);
        drawVehicle(ctx, size, c, flying);
        ctx.restore();
        if (!flying) {
            /* A hand under the vehicle: a simple cupped shape. */
            ctx.save();
            ctx.translate(cx, cy);
            ctx.rotate(tilt);
            ctx.strokeStyle = c.muted;
            ctx.lineWidth = 3;
            ctx.lineCap = 'round';
            ctx.beginPath();
            ctx.moveTo(-size * 0.3, size * 0.24);
            ctx.quadraticCurveTo(0, size * 0.32, size * 0.3, size * 0.24);
            ctx.stroke();
            ctx.restore();
        }

        const label = (string, x, y, align = 'left') => text(ctx, string, x, y, {align, baseline: 'middle', halo: true, color: c.ink, size: 12.5});

        /* Gravity: not felt by the sensor. */
        arrow(ctx, cx, cy, cx, cy + GRAVITY * unit, alpha(c.ink, this.emphasis('gravity')), {width: 2.5});
        label('gravity', cx + 8, cy + GRAVITY * unit - 4);

        if (flying) {
            arrow(ctx, cx, cy, cx + up[0] * thrust * unit, cy + up[1] * thrust * unit, alpha(c.orange, this.emphasis('felt')), {width: 3});
            label('thrust', cx + up[0] * thrust * unit + 8, cy + up[1] * thrust * unit + 4);
            if (Math.abs(drag) > 0.01) {
                const length = drag * unit * DRAG_DRAWN;
                const ox = cx + up[0] * size * 0.16, oy = cy + up[1] * size * 0.16;
                arrow(ctx, ox, oy, ox + across[0] * length, oy + across[1] * length, alpha(c.orange, this.emphasis('felt')), {width: 3});
                label(`drag (drawn ${DRAG_DRAWN} times longer)`, ox + across[0] * length + (drag > 0 ? 10 : -10), oy + across[1] * length - 2,
                    drag > 0 ? 'left' : 'right');
            }
            if (Math.abs(speed) > 0.05) {
                const vy = cy + size * 0.42;
                arrow(ctx, cx, vy, cx + speed * unit * 2.2, vy, alpha(c.blue, this.emphasis('velocity')), {width: 3});
                label('velocity', cx + speed * unit * 2.2 + (speed > 0 ? 8 : -8), vy, speed > 0 ? 'left' : 'right');
            }
        } else {
            arrow(ctx, cx, cy, cx, cy - GRAVITY * unit, alpha(c.orange, this.emphasis('felt')), {width: 3});
            label('the hand holds it up', cx + 8, cy - GRAVITY * unit + 4);
        }

        /* What the sensor reads, in the aircraft's own axes, and what a level would make of it. */
        const alongUp = felt[0] * up[0] + felt[1] * up[1], alongAcross = felt[0] * across[0] + felt[1] * across[1];
        const levelTilt = Math.atan2(-alongAcross, alongUp) / DEG;
        const shown = `${Math.abs(levelTilt).toFixed(1)}°`;
        this.readout.textContent = `The accelerometer reads ${alongUp.toFixed(1)} m/s² along the thrust axis and `
            + `${alongAcross.toFixed(2).replace('-', '−')} m/s² across it. Used as a level it would show a tilt of ${shown}. `
            + `The true tilt is ${Math.abs(this.params.tilt).toFixed(0)}°.`;
    }
}

registerFigure('forces', ForcesFigure);
