/*
The camera of a three-dimensional figure, and line drawing in its space.

Solid shapes are drawn with WebGL (gl.js). Lines, grids and labels are
drawn on the figure's ordinary canvas by projecting their points with the
same camera, so they stay sharp at any size.

The world's +Z axis points up.
*/

import {add, scale, clamp, mLookAt, mPerspective, mMul, mTransformPoint, DEG} from './math.js';

const NEAR = 0.05;

export class Camera {
    constructor(options = {}) {
        this.target = options.target ?? [0, 0, 0];
        this.distance = options.distance ?? 3;
        this.azimuth = options.azimuth ?? -60 * DEG;
        this.elevation = options.elevation ?? 25 * DEG;
        this.fov = options.fov ?? 32 * DEG;
        this.minElevation = options.minElevation ?? -5 * DEG;
        this.maxElevation = options.maxElevation ?? 89 * DEG;
        this.update(1, 1);
    }

    /* Turn the view, as a drag of (dx, dy) pixels does. */
    orbit(dx, dy) {
        this.azimuth -= dx * 0.008;
        this.elevation = clamp(this.elevation + dy * 0.008, this.minElevation, this.maxElevation);
    }

    get eye() {
        const c = Math.cos(this.elevation);
        return add(this.target, scale([c * Math.cos(this.azimuth), c * Math.sin(this.azimuth), Math.sin(this.elevation)], this.distance));
    }

    /* Work out the matrices for a picture of the given size in CSS pixels. */
    update(width, height) {
        this.width = width;
        this.height = height;
        const view = mLookAt(this.eye, this.target, [0, 0, 1]);
        const projection = mPerspective(this.fov, width / height, NEAR, 200);
        this.viewProjection = mMul(projection, view);
        return this.viewProjection;
    }

    /* A world point in clip coordinates [x, y, z, w]. */
    clip(point) {
        return mTransformPoint(this.viewProjection, point);
    }

    /* Clip coordinates to CSS pixels. */
    pixel(clip) {
        return [(clip[0] / clip[3] * 0.5 + 0.5) * this.width, (0.5 - clip[1] / clip[3] * 0.5) * this.height];
    }

    /* A world point in CSS pixels, or null if it is behind the camera. */
    project(point) {
        const clip = this.clip(point);
        return clip[3] > NEAR ? this.pixel(clip) : null;
    }

    /* How many CSS pixels one metre covers at the given world point. */
    pixelsPerMetre(point) {
        const clip = this.clip(point);
        return this.height / (2 * Math.tan(this.fov / 2) * Math.max(clip[3], NEAR));
    }

    /*
    Add a line through world points to the current canvas path. Parts
    behind the camera are cut off. "read(i)" returns point i.
    */
    path(ctx, count, read) {
        let previous = null;
        let pen = false;
        for (let i = 0; i < count; i++) {
            const point = read(i);
            const current = point ? this.clip(point) : null;
            if (previous && current) {
                let a = previous, b = current;
                const aIn = a[3] > NEAR, bIn = b[3] > NEAR;
                if (aIn || bIn) {
                    if (!aIn || !bIn) {
                        /* Cut the segment where it passes the camera plane. */
                        const t = (NEAR - a[3]) / (b[3] - a[3]);
                        const cut = [0, 1, 2, 3].map(k => a[k] + (b[k] - a[k]) * t);
                        cut[3] = NEAR;
                        if (aIn) {
                            b = cut;
                        } else {
                            a = cut;
                            pen = false;
                        }
                    }
                    const pa = this.pixel(a), pb = this.pixel(b);
                    if (!pen) {
                        ctx.moveTo(pa[0], pa[1]);
                    }
                    ctx.lineTo(pb[0], pb[1]);
                    pen = bIn;
                } else {
                    pen = false;
                }
            } else if (!current) {
                pen = false;
            }
            previous = current;
        }
    }

    /* Stroke one line through an array of world points. */
    line(ctx, points, color, width = 1.5) {
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = width;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';
        ctx.beginPath();
        this.path(ctx, points.length, i => points[i]);
        ctx.stroke();
        ctx.restore();
    }

    /*
    A square grid on the horizontal plane at height z. The lines fade out
    with distance from the centre, so the grid has no hard edge.
    */
    grid(ctx, {centre, z, extent, step, color, width = 1}) {
        const x0 = Math.floor((centre[0] - extent) / step) * step, x1 = Math.ceil((centre[0] + extent) / step) * step;
        const y0 = Math.floor((centre[1] - extent) / step) * step, y1 = Math.ceil((centre[1] + extent) / step) * step;
        /* Collect short pieces of line in a few groups of equal strength. */
        const LEVELS = 5;
        const groups = Array.from({length: LEVELS}, () => []);
        const piece = (ax, ay, bx, by) => {
            const mx = (ax + bx) / 2 - centre[0], my = (ay + by) / 2 - centre[1];
            const strength = 1 - (mx * mx + my * my) / (extent * extent);
            if (strength > 0) {
                groups[Math.min(LEVELS - 1, Math.floor(strength * LEVELS))].push([[ax, ay, z], [bx, by, z]]);
            }
        };
        for (let x = x0; x <= x1 + 1e-9; x += step) {
            for (let y = y0; y < y1 - 1e-9; y += step) {
                piece(x, y, x, y + step);
            }
        }
        for (let y = y0; y <= y1 + 1e-9; y += step) {
            for (let x = x0; x < x1 - 1e-9; x += step) {
                piece(x, y, x + step, y);
            }
        }
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = width;
        groups.forEach((pieces, level) => {
            ctx.globalAlpha = (level + 1) / LEVELS;
            ctx.beginPath();
            for (const ends of pieces) {
                this.path(ctx, 2, i => ends[i]);
            }
            ctx.stroke();
        });
        ctx.restore();
    }

    /*
    A soft round shadow on the plane at height z, under a world point.
    "color" is an [r, g, b] triple from 0 to 1 and "strength" its opacity
    at the centre.
    */
    shadow(ctx, point, z, radius, color, strength) {
        const centre = this.project([point[0], point[1], z]);
        const alongX = this.project([point[0] + radius, point[1], z]);
        const alongY = this.project([point[0], point[1] + radius, z]);
        if (!centre || !alongX || !alongY) {
            return;
        }
        /* Map a unit circle onto the ellipse which the round shadow becomes on screen. */
        const rgb = color.map(value => Math.round(value * 255)).join(', ');
        ctx.save();
        ctx.transform(alongX[0] - centre[0], alongX[1] - centre[1], alongY[0] - centre[0], alongY[1] - centre[1],
            centre[0], centre[1]);
        const gradient = ctx.createRadialGradient(0, 0, 0, 0, 0, 1);
        gradient.addColorStop(0, `rgba(${rgb}, ${strength})`);
        gradient.addColorStop(0.55, `rgba(${rgb}, ${strength * 0.55})`);
        gradient.addColorStop(1, `rgba(${rgb}, 0)`);
        ctx.fillStyle = gradient;
        ctx.beginPath();
        ctx.arc(0, 0, 1, 0, 2 * Math.PI);
        ctx.fill();
        ctx.restore();
    }
}
