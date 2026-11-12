/*
Uncertainty about attitude. The aircraft is drawn at the filter's best
guess. The blue arrows and the blue wedge show the other attitudes it must
allow for: tilted a little to each side (the arrows are their "up"), or
turned a little about the vertical (the wedge spans their headings). They
are placed as the library places its sigma points: each is the best guess
turned about one world axis.

Tilt can be measured against gravity. Heading cannot, with these sensors:
turning about the vertical changes nothing that an accelerometer or a
gyroscope at rest can see.
*/

import {Figure, registerFigure} from '../core/figure.js';
import {addSlider, addReadout} from '../core/controls.js';
import {Camera} from '../core/scene.js';
import {Renderer} from '../core/gl.js';
import {Vehicle, VEHICLE_FOOT, drawArrow} from '../core/models.js';
import {text} from '../core/plot.js';
import {alpha} from '../core/theme.js';
import {qRotate, qFromAxisAngle, qBetween, DEG} from '../core/math.js';
import {loadStory} from '../data/story.js';

class AttitudeFigure extends Figure {
    setup() {
        this.camera = new Camera({distance: 3.0, azimuth: -118 * DEG, elevation: 28 * DEG, fov: 30 * DEG, target: [0, 0, 0.1]});
        this.renderer = Renderer.get();
        this.vehicle = null;
        this.spread = Math.sqrt(12);
        addSlider(this, 'tilt', {label: 'Tilt uncertainty', min: 0, max: 6, value: 3, step: 0.05, format: v => `${v.toFixed(1)}°`});
        addSlider(this, 'heading', {label: 'Heading uncertainty', min: 0, max: 30, value: 12, step: 0.25, format: v => `${v.toFixed(0)}°`});
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
        loadStory().then(({estimator}) => {
            this.constants = estimator.constants;
            this.spread = estimator.constants.sigma_spread;
            if (this.renderer) {
                this.vehicle = new Vehicle(this.renderer, estimator.constants);
            }
            /* The attitude which stands the vehicle level: its thrust axis straight up. */
            const thrust = [this.constants.thrust_axis_x, this.constants.thrust_axis_y, this.constants.thrust_axis_z];
            this.level = qBetween(thrust, [0, 0, 1]);
            this.forward = [this.constants.lateral_axis_2_x, this.constants.lateral_axis_2_y, this.constants.lateral_axis_2_z];
            /* Look at the vehicle from in front and to one side, so the heading wedge opens towards the reader. */
            const nose = qRotate(this.level, this.forward);
            this.camera.azimuth = Math.atan2(nose[1], nose[0]) - 38 * DEG;
            this.camera.target = [0.28 * nose[0], 0.28 * nose[1], 0.1];
            this.invalidate();
        }).catch(error => console.error(error));
    }

    draw(ctx, width, height) {
        const c = this.colors;
        if (!this.vehicle) {
            text(ctx, this.renderer ? 'Loading…' : 'This figure needs WebGL.', width / 2, height / 2,
                {align: 'center', baseline: 'middle', color: c.muted, size: 14});
            return;
        }
        const camera = this.camera;
        camera.update(width, height);
        const floor = -VEHICLE_FOOT;
        camera.grid(ctx, {centre: [0, 0, 0], z: floor, extent: 1.6, step: 0.25, color: c.grid});

        /* How far the other sigma points are from the estimate: the library's spread times the uncertainty. */
        const tilt = this.params.tilt * DEG * this.spread, heading = this.params.heading * DEG * this.spread;
        const nose = qRotate(this.level, this.forward);
        const bearing = Math.atan2(nose[1], nose[0]);

        /* Heading: a wedge on the floor, from one heading guess to the other. */
        const reach = 0.95;
        const onFloor = angle => [reach * Math.cos(angle), reach * Math.sin(angle), floor];
        if (heading > 1e-4) {
            const wedge = [[0, 0, floor]];
            for (let k = 0; k <= 32; k++) {
                wedge.push(onFloor(bearing - heading + 2 * heading * k / 32));
            }
            wedge.push([0, 0, floor]);
            ctx.save();
            ctx.fillStyle = alpha(c.blue, 0.16 * this.emphasis('heading'));
            ctx.beginPath();
            camera.path(ctx, wedge.length, k => wedge[k]);
            ctx.fill();
            ctx.restore();
            for (const sign of [-1, 1]) {
                camera.line(ctx, [[0, 0, floor], onFloor(bearing + sign * heading)], alpha(c.blue, this.emphasis('heading')), 2);
            }
        }
        camera.line(ctx, [[0, 0, floor], onFloor(bearing)], alpha(c.ink2, 0.9), 1.5);
        camera.shadow(ctx, [0, 0, 0], floor, 0.34, c.dark ? [0, 0, 0] : c.inkRgb, c.dark ? 0.5 : 0.2);

        /* Tilt: the "up" direction of each tilt guess, as thin arrows around the best guess's. */
        const height3 = 0.5;
        this.renderer.begin(width, height, this.dpr, camera.viewProjection);
        this.vehicle.draw(this.renderer, c, [0, 0, 0], this.level, {});
        drawArrow(this.renderer, [0, 0, 0.05], [0, 0, height3], c.ink2Rgb, {thickness: 0.5});
        const tips = [];
        if (tilt > 1e-4) {
            for (const [axis, sign] of [[[1, 0, 0], 1], [[0, 1, 0], 1], [[1, 0, 0], -1], [[0, 1, 0], -1]]) {
                const up = qRotate(qFromAxisAngle(axis, sign * tilt), [0, 0, height3]);
                drawArrow(this.renderer, [0, 0, 0.05], up, c.blueRgb, {thickness: 0.4, opacity: this.emphasis('tilt')});
                tips.push([up[0], up[1], up[2] + 0.05]);
            }
        }
        this.renderer.finish(ctx, width, height);

        /* A ring through the tips of the tilt arrows. */
        if (tips.length > 0) {
            const radius = Math.hypot(tips[0][0], tips[0][1]), z = tips[0][2];
            const ring = Array.from({length: 49}, (_, k) => [radius * Math.cos(2 * Math.PI * k / 48), radius * Math.sin(2 * Math.PI * k / 48), z]);
            camera.line(ctx, ring, alpha(c.blue, 0.55 * this.emphasis('tilt')), 1.5);
        }
        const label = (string, point, dy) => {
            const at = camera.project(point);
            if (at) {
                text(ctx, string, at[0], at[1] + dy, {align: 'center', baseline: 'middle', color: c.ink, halo: true, size: 12.5});
            }
        };
        label('up', [0, 0, height3 + 0.05], -12);
        label('nose', onFloor(bearing), 14);
        this.readout.textContent = `The other sigma points are ${this.spread.toFixed(2)} spreads from the estimate: `
            + `${(this.params.tilt * this.spread).toFixed(1)}° of tilt, ${(this.params.heading * this.spread).toFixed(0)}° of heading.`;
    }
}

registerFigure('attitude', AttitudeFigure);
