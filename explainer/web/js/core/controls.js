/*
The controls under a figure: sliders, switches, choices, buttons and a
play bar. Each one is an ordinary HTML form control, so it works with a
keyboard and a screen reader.

A control is tied to one named value of its figure. Moving it calls
figure.set(); when the value changes for another reason (a link in the
text, an animation) the control follows.
*/

function element(tag, className, text) {
    const node = document.createElement(tag);
    if (className) {
        node.className = className;
    }
    if (text !== undefined) {
        node.textContent = text;
    }
    return node;
}

/*
A slider. Options: label, min, max, step, value, format (value to text),
and log (true to spread a wide range evenly by its powers of ten).
*/
export function addSlider(figure, name, options) {
    const {label, min, max, value, format = v => String(v), log = false} = options;
    const toSlider = v => (log ? Math.log10(v) : v);
    const fromSlider = v => (log ? 10 ** v : v);
    const row = element('label', 'ctl ctl-slider');
    const input = element('input');
    const output = element('output', 'ctl-value');
    input.type = 'range';
    input.min = toSlider(min);
    input.max = toSlider(max);
    input.step = options.step ?? (toSlider(max) - toSlider(min)) / 200;
    row.append(element('span', 'ctl-label', label), input, output);
    (options.parent ?? figure.controls).append(row);

    const show = v => {
        input.value = toSlider(v);
        output.textContent = format(v);
        /* The filled part of the track. */
        const fraction = (toSlider(v) - toSlider(min)) / (toSlider(max) - toSlider(min));
        input.style.setProperty('--fill', `${Math.min(Math.max(fraction, 0), 1) * 100}%`);
    };
    input.addEventListener('input', () => {
        figure.transitions.delete(name);
        figure.set({[name]: fromSlider(Number(input.value))}, true);
    });
    figure.params[name] = value;
    figure.inputs[name] = {show, input, row};
    show(value);
    return figure.inputs[name];
}

/* An on/off switch. */
export function addToggle(figure, name, {label, value, parent}) {
    const row = element('label', 'ctl ctl-toggle');
    const input = element('input');
    input.type = 'checkbox';
    row.append(input, element('span', 'ctl-switch'), element('span', 'ctl-label', label));
    (parent ?? figure.controls).append(row);
    const show = v => {
        input.checked = Boolean(v);
    };
    input.addEventListener('change', () => figure.set({[name]: input.checked}, true));
    figure.params[name] = value;
    figure.inputs[name] = {show, input, row};
    show(value);
    return figure.inputs[name];
}

/* A choice between a few named options: [{value, label}, ...]. */
export function addChoice(figure, name, {label, options, value, parent}) {
    const row = element('div', 'ctl ctl-choice');
    const group = element('div', 'seg');
    group.setAttribute('role', 'group');
    group.setAttribute('aria-label', label);
    const buttons = options.map(option => {
        const button = element('button', '', option.label);
        button.type = 'button';
        button.addEventListener('click', () => figure.set({[name]: option.value}, true));
        group.append(button);
        return button;
    });
    if (label) {
        row.append(element('span', 'ctl-label', label));
    }
    row.append(group);
    (parent ?? figure.controls).append(row);
    const show = v => {
        buttons.forEach((button, i) => button.setAttribute('aria-pressed', String(options[i].value === v)));
    };
    figure.params[name] = value;
    figure.inputs[name] = {show, row, buttons};
    show(value);
    return figure.inputs[name];
}

/* A plain button. */
export function addButton(figure, {label, onClick, parent, kind = ''}) {
    const button = element('button', `ctl-button ${kind}`.trim(), label);
    button.type = 'button';
    button.addEventListener('click', onClick);
    (parent ?? figure.controls).append(button);
    return button;
}

/* A row which holds several controls side by side. */
export function addRow(figure, className = '') {
    const row = element('div', `ctl-row ${className}`.trim());
    figure.controls.append(row);
    return row;
}

/* A line of text among the controls, for a live result. */
export function addReadout(figure, {parent, className = ''} = {}) {
    const node = element('div', `ctl-readout ${className}`.trim());
    node.setAttribute('aria-live', 'polite');
    (parent ?? figure.controls).append(node);
    return node;
}

const PLAY_ICON = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M4 2.5v11l9-5.5z"/></svg>';
const PAUSE_ICON = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M3.5 2.5h3v11h-3zM9.5 2.5h3v11h-3z"/></svg>';

/*
A play button with a time slider. The figure's value "name" is the time in
seconds; "playing" is kept in figure.params[`${name}Playing`]. The figure's
own step() moves the time while playing.
*/
export function addPlayBar(figure, name, {duration, value = 0, format = v => `${v.toFixed(1)} s`, parent}) {
    const row = element('div', 'ctl ctl-time');
    const button = element('button', 'ctl-play');
    button.type = 'button';
    const input = element('input');
    input.type = 'range';
    input.min = 0;
    input.max = duration;
    input.step = 'any';
    input.setAttribute('aria-label', 'Time');
    const output = element('output', 'ctl-value');
    row.append(button, input, output);
    (parent ?? figure.controls).append(row);
    const playingName = `${name}Playing`;

    const showPlaying = playing => {
        button.innerHTML = playing ? PAUSE_ICON : PLAY_ICON;
        button.setAttribute('aria-label', playing ? 'Pause' : 'Play');
    };
    const show = v => {
        input.value = v;
        output.textContent = format(v);
        input.style.setProperty('--fill', `${Math.min(Math.max(v / Number(input.max), 0), 1) * 100}%`);
    };
    button.addEventListener('click', () => {
        const playing = !figure.params[playingName];
        /* Pressing play at the end starts again from the beginning. */
        if (playing && figure.params[name] >= Number(input.max) - 1e-6) {
            figure.set({[name]: 0}, true);
        }
        figure.set({[playingName]: playing}, true);
        if (playing) {
            figure.animate();
        }
    });
    input.addEventListener('input', () => {
        figure.set({[name]: Number(input.value), [playingName]: false}, true);
    });
    figure.params[name] = value;
    figure.params[playingName] = false;
    figure.inputs[name] = {show, input, row, setDuration: d => {
        input.max = d;
        show(figure.params[name]);
    }};
    figure.inputs[playingName] = {show: showPlaying};
    show(value);
    showPlaying(false);
    return figure.inputs[name];
}
