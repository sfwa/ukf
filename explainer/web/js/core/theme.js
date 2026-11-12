/*
Colours for the canvases. The values live in css/style.css as custom
properties, so that the light and dark themes are defined in one place;
this module reads them and tells the figures when they change.

The meaning of the three hues is the same in every figure:
    blue    what the filter believes (its state)
    orange  what a sensor reads
    green   what the model predicts
Everything else is drawn in neutral ink.
*/

const NAMES = ['surface', 'surface-2', 'ink', 'ink-2', 'muted', 'grid', 'axis', 'blue', 'orange', 'green',
    'body', 'body-2', 'body-3', 'floor'];

let cache = null;
const listeners = new Set();

function parse(value) {
    const hex = value.trim().replace('#', '');
    const full = hex.length === 3 ? [...hex].map(c => c + c).join('') : hex;
    return [0, 2, 4].map(i => parseInt(full.slice(i, i + 2), 16) / 255);
}

function read() {
    const style = getComputedStyle(document.documentElement);
    const colors = {};
    for (const name of NAMES) {
        const value = style.getPropertyValue(`--${name}`).trim();
        const key = name.replace(/-(\w)/g, (_, c) => c.toUpperCase());
        colors[key] = value;
        colors[`${key}Rgb`] = parse(value);
    }
    colors.dark = style.getPropertyValue('--is-dark').trim() === '1';
    return colors;
}

/* The current colours: for example colors().blue is a CSS colour string and
colors().blueRgb is the same colour as three numbers from 0 to 1. */
export function colors() {
    if (cache === null) {
        cache = read();
    }
    return cache;
}

/* A colour string with the given opacity. */
export function alpha(color, opacity) {
    const [r, g, b] = typeof color === 'string' ? parse(color) : color;
    return `rgba(${Math.round(r * 255)}, ${Math.round(g * 255)}, ${Math.round(b * 255)}, ${opacity})`;
}

/* A colour between two others; t = 0 gives the first. */
export function blend(a, b, t) {
    const from = typeof a === 'string' ? parse(a) : a;
    const to = typeof b === 'string' ? parse(b) : b;
    return from.map((value, i) => value + (to[i] - value) * t);
}

export function onThemeChange(listener) {
    listeners.add(listener);
}

function changed() {
    cache = null;
    for (const listener of listeners) {
        listener();
    }
}

if (typeof window !== 'undefined') {
    window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', changed);
    new MutationObserver(changed).observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
}
