/*
Drawing charts on a canvas: a plot area with axes, lines, bands and labels.

Marks follow one set of rules in every figure: data lines are 2 px, grid
lines are solid hairlines, uncertainty is a light wash of the line's own
colour, and labels are written in neutral ink beside a coloured mark.
*/

import {colors, alpha} from './theme.js';
import {niceStep} from './math.js';

export const FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';

/* Write text. Options: align, baseline, color, size, weight, halo. */
export function text(ctx, string, x, y, options = {}) {
    const {align = 'left', baseline = 'alphabetic', size = 12.5, weight = 400, halo = false} = options;
    ctx.save();
    ctx.font = `${weight} ${size}px ${FONT}`;
    ctx.textAlign = align;
    ctx.textBaseline = baseline;
    if (halo) {
        ctx.lineJoin = 'round';
        ctx.lineWidth = 3.5;
        ctx.strokeStyle = colors().surface;
        ctx.strokeText(string, x, y);
    }
    ctx.fillStyle = options.color ?? colors().ink2;
    ctx.fillText(string, x, y);
    ctx.restore();
}

/* A legend entry: a short stroke (or a dot, a ring or a band) of the mark's colour, then its name in ink. */
export function key(ctx, label, x, y, color, options = {}) {
    const {dot = false, band = false, align = 'left', size = 12.5} = options;
    ctx.save();
    ctx.font = `400 ${size}px ${FONT}`;
    const width = ctx.measureText(label).width + 22;
    const left = align === 'right' ? x - width : x;
    if (band) {
        ctx.fillStyle = alpha(color, 0.3);
        ctx.fillRect(left, y - 5, 16, 10);
    }
    if (options.ring) {
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.arc(left + 8, y, 4, 0, 2 * Math.PI);
        ctx.stroke();
    } else if (dot) {
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.arc(left + 8, y, 4, 0, 2 * Math.PI);
        ctx.fill();
    } else if (!band) {
        ctx.strokeStyle = color;
        ctx.lineWidth = options.width ?? 2;
        ctx.lineCap = 'round';
        if (options.dash) {
            ctx.setLineDash(options.dash);
        }
        ctx.beginPath();
        ctx.moveTo(left + 1, y);
        ctx.lineTo(left + 15, y);
        ctx.stroke();
    }
    ctx.restore();
    text(ctx, label, left + 22, y, {baseline: 'middle', size, color: colors().ink2});
    return width;
}

/*
A row of legend entries: [{label, color, dot, band, dash}, ...]. On a narrow
canvas the row wraps onto further lines. Returns the number of lines.
*/
export function legend(ctx, entries, x, y, options = {}) {
    const right = options.maxRight ?? (ctx.canvas.clientWidth || Infinity) - 4;
    const size = options.size ?? 12.5;
    let at = x, line = y, lines = 1;
    for (const entry of entries) {
        ctx.save();
        ctx.font = `400 ${size}px ${FONT}`;
        const width = ctx.measureText(entry.label).width + 22;
        ctx.restore();
        if (at > x && at + width > right) {
            at = x;
            line += 15;
            lines++;
        }
        at += key(ctx, entry.label, at, line, entry.color, {...options, ...entry}) + 14;
    }
    return lines;
}

export function formatNumber(value, step) {
    const decimals = Math.max(0, Math.min(6, -Math.floor(Math.log10(step) + 1e-9)));
    const string = value.toFixed(decimals);
    return string === '-0' ? '0' : string.replace('-', '−');
}

export class Axes {
    constructor(x0, x1, y0, y1) {
        this.setRange(x0, x1, y0, y1);
        this.place(0, 0, 1, 1);
    }

    setRange(x0, x1, y0, y1) {
        this.x0 = x0;
        this.x1 = x1;
        this.y0 = y0;
        this.y1 = y1;
        return this;
    }

    /* The plot area, in CSS pixels. */
    place(left, top, width, height) {
        this.left = left;
        this.top = top;
        this.width = width;
        this.height = height;
        return this;
    }

    x(value) {
        return this.left + (value - this.x0) / (this.x1 - this.x0) * this.width;
    }

    y(value) {
        return this.top + this.height - (value - this.y0) / (this.y1 - this.y0) * this.height;
    }

    fromX(pixel) {
        return this.x0 + (pixel - this.left) / this.width * (this.x1 - this.x0);
    }

    fromY(pixel) {
        return this.y0 + (this.top + this.height - pixel) / this.height * (this.y1 - this.y0);
    }

    contains(px, py, margin = 0) {
        return px >= this.left - margin && px <= this.left + this.width + margin &&
            py >= this.top - margin && py <= this.top + this.height + margin;
    }

    clip(ctx) {
        ctx.beginPath();
        ctx.rect(this.left, this.top, this.width, this.height);
        ctx.clip();
    }

    /*
    Grid lines and tick labels. Options: xTicks and yTicks (how many to aim
    for; 0 for none), xFormat and yFormat, xLabel and yLabel, zero (draw the
    y = 0 line a little stronger), frame.
    */
    grid(ctx, options = {}) {
        const c = colors();
        const {xTicks = 6, yTicks = 4, zero = true, labels = true} = options;
        ctx.save();
        ctx.lineWidth = 1;
        const yStep = yTicks > 0 ? niceStep(this.y1 - this.y0, yTicks) : 0;
        const xStep = xTicks > 0 ? niceStep(this.x1 - this.x0, xTicks) : 0;
        if (yStep > 0) {
            for (let v = Math.ceil(this.y0 / yStep - 1e-9) * yStep; v <= this.y1 + yStep * 1e-9; v += yStep) {
                const py = Math.round(this.y(v)) + 0.5;
                ctx.strokeStyle = zero && Math.abs(v) < yStep * 1e-6 ? c.axis : c.grid;
                ctx.beginPath();
                ctx.moveTo(this.left, py);
                ctx.lineTo(this.left + this.width, py);
                ctx.stroke();
                if (labels) {
                    const string = options.yFormat ? options.yFormat(v) : formatNumber(v, yStep);
                    text(ctx, string, this.left - 7, py, {align: 'right', baseline: 'middle', color: c.muted, size: 11.5});
                }
            }
        }
        if (xStep > 0) {
            for (let v = Math.ceil(this.x0 / xStep - 1e-9) * xStep; v <= this.x1 + xStep * 1e-9; v += xStep) {
                const px = Math.round(this.x(v)) + 0.5;
                if (options.xGrid !== false) {
                    ctx.strokeStyle = c.grid;
                    ctx.beginPath();
                    ctx.moveTo(px, this.top);
                    ctx.lineTo(px, this.top + this.height);
                    ctx.stroke();
                }
                if (labels && options.xLabels !== false) {
                    const string = options.xFormat ? options.xFormat(v) : formatNumber(v, xStep);
                    text(ctx, string, px, this.top + this.height + 6, {align: 'center', baseline: 'top', color: c.muted, size: 11.5});
                }
            }
        }
        if (options.xLabel) {
            text(ctx, options.xLabel, this.left + this.width, this.top + this.height + 22,
                {align: 'right', baseline: 'top', color: c.muted, size: 11.5});
        }
        if (options.yLabel) {
            text(ctx, options.yLabel, this.left, this.top - 9, {align: 'left', baseline: 'bottom', color: c.ink2, size: 12, weight: 500});
        }
        ctx.restore();
    }

    /*
    A line through (xs[i], ys[i]). Gaps (NaN) break the line. Options:
    color, width, start, end, dash, and stride/offset to read one column of
    a flat array of rows.
    */
    line(ctx, xs, ys, options = {}) {
        const {color = colors().ink, width = 2, start = 0, end = xs.length, stride = 1, offset = 0, scale = 1} = options;
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = width;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';
        if (options.dash) {
            ctx.setLineDash(options.dash);
        }
        ctx.beginPath();
        const count = end - start;
        const columns = Math.max(1, Math.floor(this.width));
        let pen = false;
        if (count <= columns * 64) {
            for (let i = start; i < end; i++) {
                const value = ys[i * stride + offset] * scale;
                if (!Number.isFinite(value)) {
                    pen = false;
                    continue;
                }
                const px = this.x(xs[i]), py = this.y(value);
                if (pen) {
                    ctx.lineTo(px, py);
                } else {
                    ctx.moveTo(px, py);
                    pen = true;
                }
            }
        } else {
            /* A very long series: draw each pixel column's first, lowest, highest and last value. */
            let column = null, first = 0, last = 0, low = 0, high = 0;
            const flush = () => {
                if (column === null) {
                    return;
                }
                const px = column + 0.5;
                if (pen) {
                    ctx.lineTo(px, first);
                } else {
                    ctx.moveTo(px, first);
                    pen = true;
                }
                if (high - low > 0.5) {
                    ctx.lineTo(px, low);
                    ctx.lineTo(px, high);
                }
                ctx.lineTo(px, last);
                column = null;
            };
            for (let i = start; i < end; i++) {
                const value = ys[i * stride + offset] * scale;
                if (!Number.isFinite(value)) {
                    flush();
                    pen = false;
                    continue;
                }
                const at = Math.floor(this.x(xs[i]));
                const py = this.y(value);
                if (at !== column) {
                    flush();
                    column = at;
                    first = low = high = py;
                }
                last = py;
                low = Math.min(low, py);
                high = Math.max(high, py);
            }
            flush();
        }
        ctx.stroke();
        ctx.restore();
    }

    /*
    A filled band between centre - half and centre + half (or between two
    arrays when "low" and "high" are given instead).
    */
    band(ctx, xs, low, high, options = {}) {
        const {color = colors().blue, opacity = 0.16, start = 0, end = xs.length} = options;
        const columns = Math.max(1, Math.floor(this.width));
        const every = Math.max(1, Math.floor((end - start) / (columns * 2)));
        const top = [], bottom = [];
        for (let i = start; i < end; i += every) {
            /* Keep the extremes inside each skipped group, so thin spikes survive. */
            let lo = Infinity, hi = -Infinity;
            const stop = Math.min(i + every, end);
            for (let j = i; j < stop; j++) {
                const a = low(j), b = high(j);
                if (Number.isFinite(a) && Number.isFinite(b)) {
                    lo = Math.min(lo, a);
                    hi = Math.max(hi, b);
                }
            }
            if (hi >= lo) {
                const px = this.x(xs[i]);
                top.push(px, this.y(hi));
                bottom.push(px, this.y(lo));
            } else if (top.length > 0) {
                this.#fill(ctx, top, bottom, color, opacity);
                top.length = 0;
                bottom.length = 0;
            }
        }
        this.#fill(ctx, top, bottom, color, opacity);
    }

    #fill(ctx, top, bottom, color, opacity) {
        if (top.length < 4) {
            return;
        }
        ctx.save();
        ctx.fillStyle = alpha(color, opacity);
        ctx.beginPath();
        ctx.moveTo(top[0], top[1]);
        for (let i = 2; i < top.length; i += 2) {
            ctx.lineTo(top[i], top[i + 1]);
        }
        for (let i = bottom.length - 2; i >= 0; i -= 2) {
            ctx.lineTo(bottom[i], bottom[i + 1]);
        }
        ctx.closePath();
        ctx.fill();
        ctx.restore();
    }

    /* Shade the time spans where flag(i) is true, as a quiet full-height block. */
    spans(ctx, xs, flag, options = {}) {
        const {color = colors().ink, opacity = 0.06, start = 0, end = xs.length, top = this.top, height = this.height} = options;
        ctx.save();
        ctx.fillStyle = alpha(color, opacity);
        let from = null;
        let drawnTo = -Infinity;
        for (let i = start; i <= end; i++) {
            const on = i < end && flag(i);
            if (on && from === null) {
                from = i;
            } else if (!on && from !== null) {
                /* Start no further left than the last block ended, so blocks never overlap and darken. */
                const x0 = Math.max(this.x(xs[from]), drawnTo);
                const x1 = Math.max(this.x(xs[Math.min(i, xs.length - 1)]), x0 + 1);
                ctx.fillRect(x0, top, x1 - x0, height);
                drawnTo = x1;
                from = null;
            }
        }
        ctx.restore();
    }

    /* A dot with a ring of the surface colour, so that it reads over lines. */
    dot(ctx, x, y, color, radius = 4.5) {
        const px = this.x(x), py = this.y(y);
        ctx.save();
        ctx.fillStyle = colors().surface;
        ctx.beginPath();
        ctx.arc(px, py, radius + 2, 0, 2 * Math.PI);
        ctx.fill();
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.arc(px, py, radius, 0, 2 * Math.PI);
        ctx.fill();
        ctx.restore();
    }

    /* A vertical hairline at x, the cursor of a time plot. */
    cursor(ctx, x, color = colors().ink2) {
        const px = Math.round(this.x(x)) + 0.5;
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(px, this.top);
        ctx.lineTo(px, this.top + this.height);
        ctx.stroke();
        ctx.restore();
    }
}

/* A filled circle with a surface-coloured ring, in pixel coordinates. */
export function marker(ctx, px, py, color, radius = 4.5, ring = 2) {
    ctx.save();
    if (ring > 0) {
        ctx.fillStyle = colors().surface;
        ctx.beginPath();
        ctx.arc(px, py, radius + ring, 0, 2 * Math.PI);
        ctx.fill();
    }
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(px, py, radius, 0, 2 * Math.PI);
    ctx.fill();
    ctx.restore();
}

/* An arrow from (x0, y0) to (x1, y1) in pixels. */
export function arrow(ctx, x0, y0, x1, y1, color, options = {}) {
    const {width = 2.5, head = 9} = options;
    const length = Math.hypot(x1 - x0, y1 - y0);
    if (length < 0.5) {
        return;
    }
    const ux = (x1 - x0) / length, uy = (y1 - y0) / length;
    const size = Math.min(head, length * 0.6);
    ctx.save();
    ctx.strokeStyle = color;
    ctx.fillStyle = color;
    ctx.lineWidth = width;
    ctx.lineCap = 'round';
    ctx.beginPath();
    ctx.moveTo(x0, y0);
    ctx.lineTo(x1 - ux * size * 0.8, y1 - uy * size * 0.8);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x1 - ux * size - uy * size * 0.5, y1 - uy * size + ux * size * 0.5);
    ctx.lineTo(x1 - ux * size + uy * size * 0.5, y1 - uy * size - ux * size * 0.5);
    ctx.closePath();
    ctx.fill();
    ctx.restore();
}

/* The outline of an ellipse given by covarianceEllipse(), centred on (px, py). */
export function ellipse(ctx, px, py, major, minor, angle, options = {}) {
    const {stroke = null, fill = null, width = 2, dash = null} = options;
    ctx.save();
    ctx.beginPath();
    /* The canvas's y axis points down, so the angle is mirrored. */
    ctx.ellipse(px, py, Math.max(major, 0.01), Math.max(minor, 0.01), -angle, 0, 2 * Math.PI);
    if (fill) {
        ctx.fillStyle = fill;
        ctx.fill();
    }
    if (stroke) {
        ctx.strokeStyle = stroke;
        ctx.lineWidth = width;
        if (dash) {
            ctx.setLineDash(dash);
        }
        ctx.stroke();
    }
    ctx.restore();
}
