/* Test page for the concept figures (not part of the article). */
import {startFigures} from './core/figure.js';
import {startLinks} from './core/links.js';

/* A figure which is not written yet must not stop the others from loading. */
const names = ['noise', 'fuse', 'predict', 'loop', 'hidden', 'curve', 'sigma2d'];
await Promise.all(names.map(name => import(`./figures/${name}.js`).catch(error => {
    console.warn(`figure "${name}" did not load: ${error.message}`);
})));

startFigures();
startLinks();
