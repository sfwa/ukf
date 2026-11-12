/*
What the figures with real data share: the list of recordings, the first
recording, and one run of the estimator over it with the benchmark's own
tuning. They are loaded once, the first time a figure asks.
*/

import {fetchRecording} from './pack.js';
import {Estimator} from './estimator.js';

const WASM_URL = new URL('../../ukf_explainer.wasm', import.meta.url);
const DATA_URL = new URL('../../data/', import.meta.url);

let index = null;
let story = null;
const recordings = new Map();

/* The list of recordings on this site: [{name, file, samples, duration_s, split, reference, bytes}, ...]. */
export function loadIndex() {
    index ??= fetch(new URL('recordings.json', DATA_URL)).then(response => {
        if (!response.ok) {
            throw new Error(`Cannot load the list of recordings (${response.status})`);
        }
        return response.json();
    }).then(data => data.recordings);
    return index;
}

/* A recording by name; each is fetched once. */
export async function loadRecording(name) {
    if (!recordings.has(name)) {
        const list = await loadIndex();
        const entry = list.find(item => item.name === name);
        if (!entry) {
            throw new Error(`There is no recording called ${name}`);
        }
        recordings.set(name, fetchRecording(new URL(entry.file, DATA_URL)));
    }
    return recordings.get(name);
}

/* A new estimator of its own, for a figure which changes tuning or sensors. */
export function createEstimator() {
    return Estimator.create(WASM_URL);
}

/*
The shared story: {recording, estimator, run}. The estimator here keeps
the first recording and the benchmark's tuning, so figures may ask it to
take any step of that run apart.
*/
export function loadStory() {
    story ??= (async () => {
        const [list, estimator] = await Promise.all([loadIndex(), createEstimator()]);
        const recording = await loadRecording(list[0].name);
        await estimator.load(recording);
        const run = await estimator.run();
        return {recording, estimator, run};
    })();
    return story;
}

let side = null;
let queue = Promise.resolve();
const runs = new Map();

/*
A run of the estimator over any recording, with the benchmark's tuning:
{recording, run}. One more estimator is shared by all such runs, so they
are made one after another.
*/
export function loadRun(name) {
    if (!runs.has(name)) {
        side ??= createEstimator();
        const result = queue.then(async () => {
            const [estimator, recording] = await Promise.all([side, loadRecording(name)]);
            await estimator.load(recording);
            return {recording, run: await estimator.run(), constants: estimator.constants, tuning: estimator.defaultTuning};
        });
        queue = result.catch(() => {});
        runs.set(name, result);
    }
    return runs.get(name);
}
