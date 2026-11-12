/*
The base of every interactive figure, and the loop that draws them.

A figure is a <figure data-figure="name"> element in the article. Its
script supplies a class which extends Figure. The base class gives it a
canvas which fits the page and the screen's pixel density, controls under
the canvas, pointer input, and values which the text can change.

Figures are set up when they first come near the screen, and are drawn only
while they are on it and something has changed.
*/

import {colors, onThemeChange} from './theme.js';
import {clamp, smooth} from './math.js';

const figures = new Map();
const classes = new Map();
const active = new Set();
let frameRequested = false;
let lastFrame = 0;

export class Figure {
    constructor(root) {
        this.root = root;
        this.id = root.id;
        this.stage = root.querySelector('.fig-stage');
        this.controls = root.querySelector('.fig-controls');
        this.canvas = document.createElement('canvas');
        this.canvas.setAttribute('role', 'img');
        const caption = root.querySelector('figcaption');
        if (caption) {
            this.canvas.setAttribute('aria-label', caption.textContent.trim().replace(/\s+/g, ' '));
        }
        this.stage.append(this.canvas);
        this.ctx = this.canvas.getContext('2d');
        this.width = 0;
        this.height = 0;
        this.dpr = 1;
        this.visible = false;
        this.dirty = true;
        this.animating = false;
        this.started = false;
        this.highlighted = null;
        this.params = {};
        this.inputs = {};
        this.transitions = new Map();
        this.pointer = null;
        this.time = 0;
    }

    /* For subclasses. */
    setup() {}                  /* called once, before the first draw */
    draw() {}                   /* paint this.ctx; sizes are in CSS pixels */
    step() { return false; }    /* advance by dt seconds; return true to keep animating */
    changed() {}                /* called after any value in this.params changes */

    /* Ask for a redraw. */
    invalidate() {
        this.dirty = true;
        schedule();
    }

    /* Start calling step() every frame until it returns false. */
    animate() {
        this.animating = true;
        schedule();
    }

    /*
    Change named values, as a control or a link in the text does. Numbers
    glide to their new value unless "instant" is set.
    */
    set(values, instant = false) {
        for (const [name, value] of Object.entries(values)) {
            if (typeof value === 'number' && typeof this.params[name] === 'number' && !instant) {
                this.transitions.set(name, {from: this.params[name], to: value, progress: 0});
            } else {
                this.params[name] = value;
                this.inputs[name]?.show(value);
            }
        }
        this.changed();
        this.invalidate();
        if (this.transitions.size > 0) {
            this.animate();
        }
    }

    /* Emphasise the part of the figure which a word in the text names. */
    highlight(key) {
        if (this.highlighted !== key) {
            this.highlighted = key;
            this.invalidate();
        }
    }

    /* The opacity of a part: full unless another part is emphasised. */
    emphasis(key) {
        return this.highlighted === null || this.highlighted === key ? 1 : 0.22;
    }

    get colors() {
        return colors();
    }

    /*
    Pointer input on the canvas. The handlers receive positions in CSS
    pixels from the canvas's top left corner. "drag" is called only while a
    button or finger is down; return false from "down" to ignore a press.
    */
    onPointer({down, drag, up, move, leave}) {
        const position = event => {
            const rect = this.canvas.getBoundingClientRect();
            return {x: event.clientX - rect.left, y: event.clientY - rect.top};
        };
        this.canvas.addEventListener('pointerdown', event => {
            const at = position(event);
            if (down && down(at, event) === false) {
                return;
            }
            this.pointer = {...at, id: event.pointerId};
            this.canvas.setPointerCapture(event.pointerId);
            event.preventDefault();
        });
        this.canvas.addEventListener('pointermove', event => {
            const at = position(event);
            if (this.pointer && this.pointer.id === event.pointerId) {
                drag?.(at, {dx: at.x - this.pointer.x, dy: at.y - this.pointer.y}, event);
                this.pointer.x = at.x;
                this.pointer.y = at.y;
            } else {
                move?.(at, event);
            }
        });
        const release = event => {
            if (this.pointer && this.pointer.id === event.pointerId) {
                this.pointer = null;
                up?.(position(event), event);
            }
        };
        this.canvas.addEventListener('pointerup', release);
        this.canvas.addEventListener('pointercancel', release);
        this.canvas.addEventListener('pointerleave', event => leave?.(position(event), event));
    }

    /* Internal: make the canvas match its box and the screen. */
    fit() {
        const rect = this.stage.getBoundingClientRect();
        const dpr = Math.min(window.devicePixelRatio || 1, 2);
        const width = Math.round(rect.width), height = Math.round(rect.height);
        if (width === this.width && height === this.height && dpr === this.dpr) {
            return;
        }
        this.width = width;
        this.height = height;
        this.dpr = dpr;
        this.canvas.width = Math.round(width * dpr);
        this.canvas.height = Math.round(height * dpr);
        this.dirty = true;
    }

    /* Internal: one frame. Returns true if the figure wants another. */
    frame(dt) {
        if (!this.started) {
            this.started = true;
            this.setup();
        }
        this.time += dt;
        if (this.transitions.size > 0) {
            for (const [name, transition] of this.transitions) {
                transition.progress = clamp(transition.progress + dt / 0.45, 0, 1);
                const value = transition.from + (transition.to - transition.from) * smooth(transition.progress);
                this.params[name] = value;
                this.inputs[name]?.show(value);
                if (transition.progress >= 1) {
                    this.transitions.delete(name);
                }
            }
            this.changed();
            this.dirty = true;
        }
        if (this.animating) {
            const more = this.step(dt) === true;
            this.animating = more || this.transitions.size > 0;
            this.dirty = true;
        }
        if (this.dirty && this.width > 0 && this.height > 0) {
            this.dirty = false;
            this.ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
            this.ctx.clearRect(0, 0, this.width, this.height);
            this.draw(this.ctx, this.width, this.height);
        }
        return this.animating || this.dirty;
    }
}

function schedule() {
    if (!frameRequested) {
        frameRequested = true;
        requestAnimationFrame(loop);
    }
}

function loop(now) {
    frameRequested = false;
    /* A long pause (a hidden tab) must not become one huge step. */
    const dt = lastFrame === 0 ? 0 : Math.min((now - lastFrame) / 1000, 0.05);
    lastFrame = now;
    let more = false;
    for (const figure of active) {
        more = figure.frame(dt) || more;
    }
    if (more) {
        schedule();
    } else {
        lastFrame = 0;
    }
}

/* Tell the page which class draws the figure called "name". */
export function registerFigure(name, figureClass) {
    classes.set(name, figureClass);
}

export function getFigure(id) {
    return figures.get(id);
}

export function allFigures() {
    return [...figures.values()];
}

/* Find every figure in the page and start watching it. */
export function startFigures() {
    const resize = new ResizeObserver(entries => {
        for (const entry of entries) {
            const figure = figures.get(entry.target.closest('figure').id);
            if (figure) {
                figure.fit();
                if (figure.visible) {
                    schedule();
                }
            }
        }
    });
    const visibility = new IntersectionObserver(entries => {
        for (const entry of entries) {
            const figure = figures.get(entry.target.id);
            figure.visible = entry.isIntersecting;
            if (figure.visible) {
                figure.fit();
                active.add(figure);
                figure.dirty = true;
                schedule();
            } else {
                active.delete(figure);
            }
        }
    }, {rootMargin: '300px 0px'});

    for (const root of document.querySelectorAll('figure[data-figure]')) {
        const figureClass = classes.get(root.dataset.figure);
        if (!figureClass) {
            console.warn(`No figure is registered as "${root.dataset.figure}"`);
            continue;
        }
        const figure = new figureClass(root);
        figures.set(root.id, figure);
        resize.observe(figure.stage);
        visibility.observe(root);
    }

    onThemeChange(() => {
        for (const figure of figures.values()) {
            figure.themeChanged?.();
            figure.dirty = true;
        }
        schedule();
    });
}
