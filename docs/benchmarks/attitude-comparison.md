# Recorded attitude comparison

Recorded 2026-11-12. Full suite: **570 evaluations**, 39 BROAD and 18 Sassari recordings, both axis configurations, all 5 algorithms/variants.

Values are mean per-recording RMSE in degrees over author-marked movement intervals. Heading/total use fixed initial yaw alignment; six-axis heading is drift, not absolute heading accuracy.

## Validation split

The UKF estimator is tuned on BROAD 01, 08, 10, 24, 28 and Sassari slow (all six sensors) only. The table below evaluates the remaining 34 BROAD and 12 Sassari recordings. Reference-algorithm results on all recordings were already known.

| Dataset | Axes | Algorithm | Trials | Inclination | Heading | Total |
|---|---:|---|---:|---:|---:|---:|
| broad | 6 | ukf | 34 | 0.920 | 6.189 | 6.374 |
| broad | 6 | srukf | 34 | 0.920 | 6.189 | 6.374 |
| broad | 6 | madgwick | 34 | 2.207 | 23.516 | 23.794 |
| broad | 6 | mahony | 34 | 11.056 | 71.510 | 72.620 |
| broad | 6 | vqf | 34 | 0.827 | 6.137 | 6.277 |
| broad | 9 | ukf | 34 | 0.920 | 2.266 | 2.486 |
| broad | 9 | srukf | 34 | 0.920 | 2.266 | 2.486 |
| broad | 9 | madgwick | 34 | 2.491 | 4.873 | 5.591 |
| broad | 9 | mahony | 34 | 9.603 | 28.782 | 30.444 |
| broad | 9 | vqf | 34 | 0.827 | 2.403 | 2.589 |
| sassari | 6 | ukf | 12 | 1.979 | 3.867 | 4.514 |
| sassari | 6 | srukf | 12 | 1.979 | 3.867 | 4.514 |
| sassari | 6 | madgwick | 12 | 4.040 | 27.492 | 28.335 |
| sassari | 6 | mahony | 12 | 10.595 | 66.774 | 67.828 |
| sassari | 6 | vqf | 12 | 2.326 | 3.500 | 4.343 |
| sassari | 9 | ukf | 12 | 1.979 | 3.990 | 4.511 |
| sassari | 9 | srukf | 12 | 1.979 | 3.990 | 4.511 |
| sassari | 9 | madgwick | 12 | 3.941 | 4.311 | 5.917 |
| sassari | 9 | mahony | 12 | 9.965 | 15.701 | 18.702 |
| sassari | 9 | vqf | 12 | 2.326 | 3.845 | 4.580 |

## Full dataset comparison

| Dataset | Axes | Algorithm | Trials | Inclination | Heading | Total |
|---|---:|---|---:|---:|---:|---:|
| broad | 6 | ukf | 39 | 0.936 | 6.218 | 6.404 |
| broad | 6 | srukf | 39 | 0.936 | 6.218 | 6.404 |
| broad | 6 | madgwick | 39 | 2.158 | 22.987 | 23.257 |
| broad | 6 | mahony | 39 | 10.326 | 68.586 | 69.618 |
| broad | 6 | vqf | 39 | 0.857 | 6.311 | 6.448 |
| broad | 9 | ukf | 39 | 0.936 | 2.244 | 2.472 |
| broad | 9 | srukf | 39 | 0.936 | 2.244 | 2.472 |
| broad | 9 | madgwick | 39 | 2.411 | 4.605 | 5.316 |
| broad | 9 | mahony | 39 | 8.934 | 26.454 | 28.014 |
| broad | 9 | vqf | 39 | 0.857 | 2.397 | 2.593 |
| sassari | 6 | ukf | 18 | 1.861 | 3.862 | 4.423 |
| sassari | 6 | srukf | 18 | 1.861 | 3.862 | 4.423 |
| sassari | 6 | madgwick | 18 | 3.377 | 31.782 | 32.441 |
| sassari | 6 | mahony | 18 | 8.099 | 69.223 | 69.950 |
| sassari | 6 | vqf | 18 | 2.095 | 3.637 | 4.323 |
| sassari | 9 | ukf | 18 | 1.861 | 3.689 | 4.199 |
| sassari | 9 | srukf | 18 | 1.861 | 3.689 | 4.199 |
| sassari | 9 | madgwick | 18 | 3.296 | 3.848 | 5.143 |
| sassari | 9 | mahony | 18 | 7.548 | 11.583 | 13.942 |
| sassari | 9 | vqf | 18 | 2.095 | 3.573 | 4.235 |

## Conditions

Mean per-recording total RMSE in degrees. BROAD recordings 24 to 39 contain deliberate magnetic or mechanical disturbances; Sassari recordings differ in rotation speed.

| Dataset | Axes | Condition | Recordings | UKF | SR-UKF | Madgwick | Mahony | VQF |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| broad | 6 | disturbed | 16 | 5.996 | 5.996 | 16.290 | 75.443 | 6.038 |
| broad | 9 | disturbed | 16 | 2.689 | 2.689 | 5.677 | 36.618 | 3.185 |
| broad | 6 | undisturbed | 23 | 6.689 | 6.689 | 28.104 | 65.566 | 6.733 |
| broad | 9 | undisturbed | 23 | 2.322 | 2.322 | 5.065 | 22.028 | 2.181 |
| sassari | 6 | fast | 6 | 5.380 | 5.380 | 28.199 | 53.491 | 4.959 |
| sassari | 9 | fast | 6 | 5.177 | 5.177 | 7.054 | 14.191 | 5.184 |
| sassari | 6 | medium | 6 | 3.648 | 3.648 | 28.471 | 82.166 | 3.727 |
| sassari | 9 | medium | 6 | 3.845 | 3.845 | 4.780 | 23.213 | 3.977 |
| sassari | 6 | slow | 6 | 4.241 | 4.241 | 40.653 | 74.192 | 4.282 |
| sassari | 9 | slow | 6 | 3.574 | 3.574 | 3.596 | 4.422 | 3.545 |

`ukf` and `srukf` run the same estimator model on `UKF::Core` and `UKF::SquareRootCore`. These measurements are not an optimized ranking of filter families.

[Methodology and reproduction](../attitude-benchmark.md) · [Per-trial metrics](attitude-per-trial.csv) · [Source hashes, dependencies and provenance](attitude-provenance.json)

macOS 27.2 arm64, Apple clang 21.0.0, CMake 4.4.4, Release, patched Eigen 3.2.8, double precision; Python 3.12.16 and pinned reference packages. The provenance revision identifies the checked-out commit; the dirty flag and source hashes identify local changes. CI regenerates reports for the checked-out commit. Timing is diagnostic and not a speed ranking.
