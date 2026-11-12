/*
Shared pieces of the concept figures (the figures at the start of the
article, which use made-up numbers instead of a recording): a quadcopter
seen from the side, bell curves, brackets and a small value box.
*/

import {alpha} from '../core/theme.js';
import {text, FONT} from '../core/plot.js';

/* True if the reader has asked the browser for less motion. */
export function reducedMotion() {
    return window.matchMedia('(prefers-reduced-motion: reduce)').matches;
}

function roundedRect(ctx, x, y, width, height, radius) {
    const r = Math.min(radius, width / 2, height / 2);
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.lineTo(x + width - r, y);
    ctx.arcTo(x + width, y, x + width, y + r, r);
    ctx.lineTo(x + width, y + height - r);
    ctx.arcTo(x + width, y + height, x + width - r, y + height, r);
    ctx.lineTo(x + r, y + height);
    ctx.arcTo(x, y + height, x, y + height - r, r);
    ctx.lineTo(x, y + r);
    ctx.arcTo(x, y, x + r, y, r);
    ctx.closePath();
}

/*
A quadcopter seen from the side, centred on (x, y) in pixels. "size" is its
width from rotor tip to rotor tip. A positive tilt leans it to the right
(clockwise on the screen), in radians.

Options: tint (draw every part in this one colour), opacity, spinning
(false draws the rotors as still blades), sensor (true draws the small
range sensor under the body).
*/
export function drawQuadSide(ctx, x, y, size, tilt, colors, options = {}) {
    const {tint = null, opacity = 1, spinning = true, sensor = false} = options;
    const s = size;
    const frame = tint ?? colors.body, dark = tint ?? colors.body2, rotor = tint ?? colors.body3;
    ctx.save();
    ctx.translate(x, y);
    ctx.rotate(tilt);
    ctx.globalAlpha *= opacity;
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';

    /* Legs and feet. */
    ctx.strokeStyle = frame;
    ctx.lineWidth = Math.max(1, s * 0.022);
    ctx.beginPath();
    for (const side of [-1, 1]) {
        ctx.moveTo(side * s * 0.07, s * 0.05);
        ctx.lineTo(side * s * 0.17, s * 0.17);
        ctx.moveTo(side * s * 0.12, s * 0.17);
        ctx.lineTo(side * s * 0.22, s * 0.17);
    }
    ctx.stroke();

    /* The arm from rotor to rotor. */
    ctx.fillStyle = frame;
    roundedRect(ctx, -s * 0.38, -s * 0.018, s * 0.76, s * 0.036, s * 0.018);
    ctx.fill();

    /* Motors, then rotors: a blurred disc while turning, a blade when not. */
    for (const side of [-1, 1]) {
        ctx.fillStyle = dark;
        roundedRect(ctx, side * s * 0.36 - s * 0.03, -s * 0.085, s * 0.06, s * 0.085, s * 0.012);
        ctx.fill();
        ctx.fillStyle = rotor;
        ctx.globalAlpha *= spinning ? 0.55 : 1;
        ctx.beginPath();
        ctx.ellipse(side * s * 0.36, -s * 0.1, s * 0.14, spinning ? s * 0.02 : s * 0.011, 0, 0, 2 * Math.PI);
        ctx.fill();
        ctx.globalAlpha /= spinning ? 0.55 : 1;
    }

    /* The body. */
    ctx.fillStyle = dark;
    roundedRect(ctx, -s * 0.12, -s * 0.07, s * 0.24, s * 0.13, s * 0.045);
    ctx.fill();
    if (sensor) {
        ctx.fillStyle = frame;
        roundedRect(ctx, -s * 0.03, s * 0.055, s * 0.06, s * 0.035, s * 0.01);
        ctx.fill();
    }
    ctx.restore();
}

/* Where the range sensor of a quadcopter drawn by drawQuadSide() is, in pixels. */
export function quadSensorPoint(x, y, size, tilt) {
    const down = size * 0.09;
    return [x - Math.sin(tilt) * down, y + Math.cos(tilt) * down];
}

/*
A bell curve along an axis. "along(value)" gives the pixel position of a
value on the axis and "base" the pixel position of the axis itself; the
bell rises from the axis by "peak" pixels at its centre. With vertical set,
the axis runs up the screen and the bell rises sideways (a negative peak
rises to the left; on a level axis a positive peak rises upwards).

Options: color, fill (opacity of the wash under the curve; 0 for none),
width, dash, opacity, sigmas (how far out to draw), from and to (clip the
curve to this range of values).
*/
export function drawBell(ctx, options) {
    const {mean, sd, along, base, peak, vertical = false, color, fill = 0.14, width = 2, dash = null,
        opacity = 1, sigmas = 4} = options;
    const low = Math.max(mean - sigmas * sd, options.from ?? -Infinity);
    const high = Math.min(mean + sigmas * sd, options.to ?? Infinity);
    if (!(high > low)) {
        return;
    }
    const steps = 96;
    const points = [];
    for (let i = 0; i <= steps; i++) {
        const value = low + (high - low) * i / steps;
        const rise = peak * Math.exp(-0.5 * ((value - mean) / sd) ** 2);
        points.push(vertical ? [base + rise, along(value)] : [along(value), base - rise]);
    }
    const first = vertical ? [base, along(low)] : [along(low), base];
    const last = vertical ? [base, along(high)] : [along(high), base];
    ctx.save();
    if (fill > 0) {
        ctx.fillStyle = alpha(color, fill * opacity);
        ctx.beginPath();
        ctx.moveTo(first[0], first[1]);
        for (const [px, py] of points) {
            ctx.lineTo(px, py);
        }
        ctx.lineTo(last[0], last[1]);
        ctx.closePath();
        ctx.fill();
    }
    if (width > 0) {
        ctx.strokeStyle = alpha(color, opacity);
        ctx.lineWidth = width;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';
        if (dash) {
            ctx.setLineDash(dash);
        }
        ctx.beginPath();
        points.forEach(([px, py], i) => (i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py)));
        ctx.stroke();
    }
    ctx.restore();
}

/*
A measuring bracket between two pixel positions: a line with a short bar
at each end. Used to mark a spread.
*/
export function drawBracket(ctx, x0, y0, x1, y1, color, options = {}) {
    const {width = 1.5, bar = 5} = options;
    const length = Math.hypot(x1 - x0, y1 - y0);
    if (length < 0.5) {
        return;
    }
    const nx = -(y1 - y0) / length * bar, ny = (x1 - x0) / length * bar;
    ctx.save();
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.lineCap = 'round';
    ctx.beginPath();
    ctx.moveTo(x0, y0);
    ctx.lineTo(x1, y1);
    ctx.moveTo(x0 - nx, y0 - ny);
    ctx.lineTo(x0 + nx, y0 + ny);
    ctx.moveTo(x1 - nx, y1 - ny);
    ctx.lineTo(x1 + nx, y1 + ny);
    ctx.stroke();
    ctx.restore();
}

/*
A small box of values beside a cursor. "rows" is a list of {label, value,
color}; a row with a colour gets a short stroke of that colour as its key.
The box is kept inside the rectangle "bounds" = {left, top, right, bottom}.
*/
export function drawValueBox(ctx, x, y, title, rows, colors, bounds) {
    const size = 12, lineHeight = 17, pad = 8;
    ctx.save();
    ctx.font = `400 ${size}px ${FONT}`;
    let inner = title ? ctx.measureText(title).width : 0;
    for (const row of rows) {
        inner = Math.max(inner, 20 + ctx.measureText(`${row.label}  ${row.value}`).width);
    }
    const width = inner + 2 * pad, height = (rows.length + (title ? 1 : 0)) * lineHeight + 2 * pad - 3;
    let left = x + 12;
    if (left + width > bounds.right) {
        left = x - 12 - width;
    }
    left = Math.max(bounds.left, left);
    const top = Math.min(Math.max(y, bounds.top), bounds.bottom - height);
    ctx.fillStyle = alpha(colors.surface, 0.94);
    ctx.strokeStyle = colors.axis;
    ctx.lineWidth = 1;
    roundedRect(ctx, left + 0.5, top + 0.5, width, height, 6);
    ctx.fill();
    ctx.stroke();
    ctx.restore();
    let line = top + pad + size - 2;
    if (title) {
        text(ctx, title, left + pad, line, {color: colors.muted, size});
        line += lineHeight;
    }
    for (const row of rows) {
        if (row.color) {
            ctx.save();
            ctx.strokeStyle = row.color;
            ctx.lineWidth = 2;
            ctx.lineCap = 'round';
            ctx.beginPath();
            ctx.moveTo(left + pad + 1, line - 4);
            ctx.lineTo(left + pad + 13, line - 4);
            ctx.stroke();
            ctx.restore();
        }
        text(ctx, row.label, left + pad + 20, line, {color: colors.ink2, size});
        text(ctx, row.value, left + width - pad, line, {color: colors.ink, size, weight: 600, align: 'right'});
        line += lineHeight;
    }
}

/* The mean and the spread (standard deviation) of the first "count" values. */
export function meanAndSpread(values, count = values.length) {
    let sum = 0;
    for (let i = 0; i < count; i++) {
        sum += values[i];
    }
    const mean = count > 0 ? sum / count : NaN;
    let squares = 0;
    for (let i = 0; i < count; i++) {
        squares += (values[i] - mean) ** 2;
    }
    return {mean, spread: count > 1 ? Math.sqrt(squares / (count - 1)) : NaN};
}

/* A minus sign which lines up with digits, for numbers shown to the reader. */
export function signed(value, digits) {
    const string = value.toFixed(digits);
    return (Number(string) === 0 ? string.replace('-', '') : string).replace('-', '−');
}
