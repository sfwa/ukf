/*
The links between the text and the figures.

    <span class="k k-blue" data-hl="estimate">estimate</span>
        Pointing at this word emphasises that part of the nearest figure.

    <button class="act" data-fig="fig-fuse" data-set='{"noise": 0.1}'>…</button>
        Pressing this phrase changes values of the named figure.

    <span class="live" data-live="tilt-rms"></span>
        A number which a figure fills in, with setLive().
*/

import {getFigure, allFigures} from './figure.js';

/* The figure a word refers to: the named one, or else the nearest on the page. */
function figureFor(node) {
    if (node.dataset.fig) {
        return getFigure(node.dataset.fig);
    }
    const box = node.getBoundingClientRect();
    const middle = (box.top + box.bottom) / 2;
    let nearest = null, distance = Infinity;
    for (const figure of allFigures()) {
        const rect = figure.root.getBoundingClientRect();
        const gap = middle < rect.top ? rect.top - middle : middle > rect.bottom ? middle - rect.bottom : 0;
        if (gap < distance) {
            distance = gap;
            nearest = figure;
        }
    }
    return nearest;
}

export function startLinks() {
    for (const node of document.querySelectorAll('.k[data-hl]')) {
        node.tabIndex = 0;
        let target = null;
        const on = () => {
            target = figureFor(node);
            target?.highlight(node.dataset.hl);
        };
        const off = () => {
            target?.highlight(null);
            target = null;
        };
        node.addEventListener('pointerenter', on);
        node.addEventListener('pointerleave', off);
        node.addEventListener('focus', on);
        node.addEventListener('blur', off);
    }

    for (const button of document.querySelectorAll('button.act')) {
        button.type = 'button';
        button.addEventListener('click', () => {
            const figure = figureFor(button);
            if (!figure) {
                return;
            }
            if (button.dataset.set) {
                figure.set(JSON.parse(button.dataset.set));
            }
            if (button.dataset.do) {
                figure.action?.(button.dataset.do);
            }
            /* Bring the figure into view if the reader cannot see it change. */
            const rect = figure.root.getBoundingClientRect();
            if (rect.bottom < 80 || rect.top > window.innerHeight - 80) {
                figure.root.scrollIntoView({block: 'center'});
            }
        });
    }
}

/* Fill in every number in the text which is marked data-live="name". */
export function setLive(name, value) {
    for (const node of document.querySelectorAll(`.live[data-live="${name}"]`)) {
        node.textContent = value;
    }
}
