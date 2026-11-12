/*
Start of the page: load every figure's script, then find the figures in
the article and bring them to life.
*/

import {startFigures} from './core/figure.js';
import {startLinks} from './core/links.js';
import {startFacts} from './figures/facts.js';
import './figures/flight.js';
import './figures/noise.js';
import './figures/fuse.js';
import './figures/predict.js';
import './figures/loop.js';
import './figures/hidden.js';
import './figures/curve.js';
import './figures/sigma2d.js';
import './figures/attitude.js';
import './figures/gyro.js';
import './figures/forces.js';
import './figures/level.js';
import './figures/drag.js';
import './figures/fan.js';
import './figures/modes.js';
import './figures/step.js';
import './figures/playground.js';
import './figures/position.js';

startFigures();
startLinks();
startFacts().catch(error => console.error(error));
