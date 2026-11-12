/*
The estimator runs here, off the page's main thread, so that the page
stays smooth while a recording is processed. The page talks to this worker
through estimator.js.
*/

import {EstimatorModule} from './wasm.js';

let estimator = null;

self.onmessage = async ({data}) => {
    const {id, type} = data;
    try {
        if (type === 'init') {
            estimator = await EstimatorModule.create(data.module);
            self.postMessage({id, result: estimator.describe()});
        } else if (type === 'load') {
            estimator.load(data.time, data.gyroscope, data.accelerometer);
            self.postMessage({id, result: true});
        } else if (type === 'run') {
            const started = performance.now();
            const result = estimator.run(data.tuning);
            result.elapsed = performance.now() - started;
            self.postMessage({id, result}, [result.trace.buffer]);
        } else if (type === 'predict') {
            const prediction = estimator.predict(data.index, data.steps);
            self.postMessage({id, result: prediction}, prediction ? [prediction.buffer] : []);
        } else if (type === 'explain') {
            const step = estimator.explain(data.index);
            self.postMessage({id, result: step}, step ? [step.buffer] : []);
        }
    } catch (error) {
        self.postMessage({id, error: String(error?.message ?? error)});
    }
};
