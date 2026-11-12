/*
Numbers in the text which come from the estimator itself. Each is written
into the page as a default, and replaced here by the value from the
estimator module and from the run over the first recording, so the text
cannot drift away from what the code does.
*/

import {loadStory} from '../data/story.js';
import {setLive} from '../core/links.js';
import {qAt} from '../core/math.js';
import {gyroOnlyAttitude, tiltDegrees} from './story-kit.js';
import {tilts} from './level.js';
import {dragPoints} from './drag.js';
import {reckon, headingFix} from './position.js';

const DEGREES = 180 / Math.PI;

/* A weight such as 0.041666 written as the fraction 1/24. */
function fraction(weight) {
    const parts = Math.round(1 / weight);
    return Math.abs(1 / parts - weight) < 1e-9 ? `1/${parts}` : weight.toFixed(3);
}

export async function startFacts() {
    const {recording, run, estimator} = await loadStory();
    const constants = estimator.constants, tuning = estimator.defaultTuning;
    const step = recording.duration / (recording.samples - 1);
    setLive('tilt-rms', run.errors.rms.inclination.toFixed(2));
    setLive('heading-rms', run.errors.rms.heading.toFixed(1));
    setLive('error-size', String(constants.error_size));
    setLive('sigma-count', String(constants.sigma_count));
    setLive('sigma-spread', constants.sigma_spread.toFixed(2));
    setLive('weight-centre', fraction(constants.weight_mean_centre));
    setLive('weight-outer', fraction(constants.weight_mean_outer));
    setLive('heading-limit', (tuning.heading_sigma_limit * DEGREES).toFixed(0));
    setLive('sample-rate', (1 / step).toFixed(0));
    setLive('step-ms', (step * 1000).toFixed(0));
    setLive('tilt-leak', (constants.gravity * Math.sin(1 / DEGREES)).toFixed(2));

    const sampleAt = t => Math.round(t / step);
    const thrust = [constants.thrust_axis_x, constants.thrust_axis_y, constants.thrust_axis_z];

    /* The filter's own uncertainty halfway through the flight: tilt (per axis) and heading. */
    const middle = Math.round(recording.samples / 2);
    setLive('sd-tilt', (Math.hypot(run.value('sd', middle, 0), run.value('sd', middle, 1)) / Math.SQRT2 * DEGREES).toFixed(1));
    setLive('sd-heading', (run.value('sd', middle, 2) * DEGREES).toFixed(0));

    /* The gyroscope: what it reads while the aircraft is at rest, and where adding it up leads. */
    const window = sampleAt(14);
    const mean = [0, 0, 0];
    let rested = 0;
    for (let i = 0; i < window; i++) {
        if (run.value('resting', i)) {
            for (let axis = 0; axis < 3; axis++) {
                mean[axis] += recording.gyroscope[3 * i + axis];
            }
            rested++;
        }
    }
    if (rested > 0) {
        setLive('gyro-rest-rate', (Math.hypot(...mean) / rested * DEGREES).toFixed(1));
    }
    const gyroOnly = gyroOnlyAttitude(recording, sampleAt(10) + 1);
    setLive('gyro-only-tilt-10s', tiltDegrees(qAt(gyroOnly, sampleAt(10)), thrust).toFixed(0));

    /* The accelerometer as a spirit level, while the rotor model is in use. */
    const level = tilts(recording, thrust, 0, recording.duration);
    let largestTruth = 0, largestLevel = 0;
    for (let i = 0; i < level.time.length; i++) {
        if (run.value('rotor_model_active', i) && Number.isFinite(level.truth[i])) {
            largestTruth = Math.max(largestTruth, level.truth[i]);
            largestLevel = Math.max(largestLevel, level.level[i]);
        }
    }
    setLive('level-true-max', largestTruth.toFixed(0));
    setLive('level-accel-max', largestLevel.toFixed(0));

    /* Rotor drag: the estimator's rate, and the slope which fits this recording best. */
    setLive('drag-rate', constants.drag_rate.toFixed(2));
    setLive('drag-fit', (-dragPoints({recording, run, estimator}).fit.slope).toFixed(2));

    /* Dead reckoning from 20 s into the recording. */
    const paths = reckon({recording, run}, headingFix(run), sampleAt(20));
    const after = seconds => Math.min(Math.round(seconds / step), paths.count - 1);
    setLive('reckon-raw-2s', paths.rawError[after(2)].toFixed(1));
    setLive('reckon-raw-10s', paths.rawError[after(10)].toFixed(0));
    setLive('reckon-filter-10s', paths.filterError[after(10)].toFixed(0));
}
