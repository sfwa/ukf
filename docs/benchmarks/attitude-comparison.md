# Recorded attitude comparison

Recorded 2026-11-11. Full suite: **912 evaluations**, 39 BROAD and 18 Sassari recordings, both axis configurations, all 8 algorithms/variants.

Values are mean per-recording RMSE in degrees over author-marked movement intervals. Heading/total use fixed initial yaw alignment; six-axis heading is drift, not absolute heading accuracy.

## Validation split

New-estimator development uses BROAD 01, 08, 10, 24, 28 and Sassari slow (all six sensors). The table below evaluates the remaining 34 BROAD and 12 Sassari recordings. These labels govern the new implementation's iteration; baseline results on all recordings were already known.

| Dataset | Axes | Algorithm | Trials | Inclination | Heading | Total |
|---|---:|---|---:|---:|---:|---:|
| broad | 6 | ukf | 34 | 10.430 | 75.850 | 76.659 |
| broad | 6 | srukf | 34 | 10.402 | 75.952 | 76.760 |
| broad | 6 | robust-rest | 34 | 4.234 | 9.180 | 10.778 |
| broad | 6 | robust-acc | 34 | 0.861 | 6.061 | 6.205 |
| broad | 6 | robust | 34 | 0.876 | 6.061 | 6.208 |
| broad | 6 | madgwick | 34 | 2.207 | 23.516 | 23.794 |
| broad | 6 | mahony | 34 | 11.056 | 71.510 | 72.620 |
| broad | 6 | vqf | 34 | 0.827 | 6.137 | 6.277 |
| broad | 9 | ukf | 34 | 7.341 | 12.992 | 15.078 |
| broad | 9 | srukf | 34 | 7.341 | 12.992 | 15.078 |
| broad | 9 | robust-rest | 34 | 4.234 | 6.125 | 7.571 |
| broad | 9 | robust-acc | 34 | 0.861 | 3.520 | 3.663 |
| broad | 9 | robust | 34 | 0.876 | 2.930 | 3.087 |
| broad | 9 | madgwick | 34 | 2.491 | 4.873 | 5.591 |
| broad | 9 | mahony | 34 | 9.603 | 28.782 | 30.444 |
| broad | 9 | vqf | 34 | 0.827 | 2.403 | 2.589 |
| sassari | 6 | ukf | 12 | 6.173 | 97.271 | 97.424 |
| sassari | 6 | srukf | 12 | 6.173 | 97.271 | 97.424 |
| sassari | 6 | robust-rest | 12 | 4.313 | 3.288 | 5.684 |
| sassari | 6 | robust-acc | 12 | 2.407 | 3.507 | 4.398 |
| sassari | 6 | robust | 12 | 2.457 | 3.507 | 4.436 |
| sassari | 6 | madgwick | 12 | 4.040 | 27.492 | 28.335 |
| sassari | 6 | mahony | 12 | 10.595 | 66.774 | 67.828 |
| sassari | 6 | vqf | 12 | 2.326 | 3.500 | 4.343 |
| sassari | 9 | ukf | 12 | 4.869 | 5.956 | 7.724 |
| sassari | 9 | srukf | 12 | 4.869 | 5.956 | 7.724 |
| sassari | 9 | robust-rest | 12 | 4.313 | 5.081 | 6.818 |
| sassari | 9 | robust-acc | 12 | 2.407 | 4.061 | 4.823 |
| sassari | 9 | robust | 12 | 2.457 | 4.168 | 4.952 |
| sassari | 9 | madgwick | 12 | 3.941 | 4.311 | 5.917 |
| sassari | 9 | mahony | 12 | 9.965 | 15.701 | 18.702 |
| sassari | 9 | vqf | 12 | 2.326 | 3.845 | 4.580 |

## Full dataset comparison

| Dataset | Axes | Algorithm | Trials | Inclination | Heading | Total |
|---|---:|---|---:|---:|---:|---:|
| broad | 6 | ukf | 39 | 9.866 | 76.368 | 77.103 |
| broad | 6 | srukf | 39 | 9.843 | 76.504 | 77.237 |
| broad | 6 | robust-rest | 39 | 4.014 | 9.517 | 10.966 |
| broad | 6 | robust-acc | 39 | 0.892 | 6.232 | 6.374 |
| broad | 6 | robust | 39 | 0.906 | 6.232 | 6.376 |
| broad | 6 | madgwick | 39 | 2.158 | 22.987 | 23.257 |
| broad | 6 | mahony | 39 | 10.326 | 68.586 | 69.618 |
| broad | 6 | vqf | 39 | 0.857 | 6.311 | 6.448 |
| broad | 9 | ukf | 39 | 6.978 | 12.048 | 14.077 |
| broad | 9 | srukf | 39 | 6.978 | 12.048 | 14.077 |
| broad | 9 | robust-rest | 39 | 4.014 | 5.780 | 7.148 |
| broad | 9 | robust-acc | 39 | 0.892 | 3.487 | 3.637 |
| broad | 9 | robust | 39 | 0.906 | 2.881 | 3.051 |
| broad | 9 | madgwick | 39 | 2.411 | 4.605 | 5.316 |
| broad | 9 | mahony | 39 | 8.934 | 26.454 | 28.014 |
| broad | 9 | vqf | 39 | 0.857 | 2.397 | 2.593 |
| sassari | 6 | ukf | 18 | 4.988 | 88.484 | 88.599 |
| sassari | 6 | srukf | 18 | 4.988 | 88.484 | 88.599 |
| sassari | 6 | robust-rest | 18 | 3.486 | 3.547 | 5.301 |
| sassari | 6 | robust-acc | 18 | 2.163 | 3.694 | 4.415 |
| sassari | 6 | robust | 18 | 2.196 | 3.695 | 4.440 |
| sassari | 6 | madgwick | 18 | 3.377 | 31.782 | 32.441 |
| sassari | 6 | mahony | 18 | 8.099 | 69.223 | 69.950 |
| sassari | 6 | vqf | 18 | 2.095 | 3.637 | 4.323 |
| sassari | 9 | ukf | 18 | 4.000 | 5.082 | 6.536 |
| sassari | 9 | srukf | 18 | 4.000 | 5.082 | 6.536 |
| sassari | 9 | robust-rest | 18 | 3.486 | 4.306 | 5.693 |
| sassari | 9 | robust-acc | 18 | 2.163 | 3.621 | 4.329 |
| sassari | 9 | robust | 18 | 2.196 | 3.692 | 4.415 |
| sassari | 9 | madgwick | 18 | 3.296 | 3.848 | 5.143 |
| sassari | 9 | mahony | 18 | 7.548 | 11.583 | 13.942 |
| sassari | 9 | vqf | 18 | 2.095 | 3.573 | 4.235 |

## Nine-axis conditions

Mean per-recording total RMSE in degrees. The final stage adds both magnetic rejection and acceleration guards; protection can trade accuracy on clean data for disturbance tolerance.

| Dataset | Condition | Recordings | robust-acc | robust | VQF |
|---|---|---:|---:|---:|---:|
| broad | disturbed | 16 | 4.574 | 3.146 | 3.185 |
| broad | undisturbed | 23 | 2.986 | 2.985 | 2.181 |
| sassari | fast | 6 | 5.670 | 5.932 | 5.184 |
| sassari | medium | 6 | 3.977 | 3.972 | 3.977 |
| sassari | slow | 6 | 3.341 | 3.341 | 3.545 |

The baseline `ukf`/`srukf` use the original fixed-noise joint attitude/bias model. `robust-rest`, `robust-acc`, and `robust` are cumulative estimator variants; the first stage also changes the bias/heading architecture, so its difference is not a rest-only ablation. These measurements are not an optimized ranking of filter families.

[Methodology and reproduction](../attitude-benchmark.md) · [Per-trial metrics](attitude-per-trial.csv) · [Source hashes, dependencies and provenance](attitude-provenance.json)

Linux x86_64, GCC 14.2.0, CMake 3.31.10, Release, double precision; Python 3.12.14 and pinned reference packages. The provenance revision identifies the checked-out commit; the dirty flag and source hashes identify local changes. CI regenerates reports for the checked-out commit. Timing is diagnostic and not a speed ranking.
