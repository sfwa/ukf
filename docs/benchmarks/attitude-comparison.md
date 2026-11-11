# Recorded attitude comparison

Recorded 2026-11-11. Full suite: **570 evaluations**, 39 BROAD and 18 Sassari recordings, both axis configurations, all five algorithms.

Values are mean per-recording RMSE in degrees over author-marked movement intervals. Heading/total use fixed initial yaw alignment; six-axis heading is drift, not absolute heading accuracy.

| Dataset | Axes | Algorithm | Inclination | Heading | Total |
|---|---:|---|---:|---:|---:|
| broad | 6 | ukf | 9.866 | 76.368 | 77.103 |
| broad | 6 | srukf | 9.843 | 76.504 | 77.237 |
| broad | 6 | madgwick | 2.158 | 22.987 | 23.257 |
| broad | 6 | mahony | 10.326 | 68.586 | 69.618 |
| broad | 6 | vqf | 0.857 | 6.311 | 6.448 |
| broad | 9 | ukf | 6.978 | 12.048 | 14.077 |
| broad | 9 | srukf | 6.978 | 12.048 | 14.077 |
| broad | 9 | madgwick | 2.411 | 4.605 | 5.316 |
| broad | 9 | mahony | 8.934 | 26.454 | 28.014 |
| broad | 9 | vqf | 0.857 | 2.397 | 2.593 |
| sassari | 6 | ukf | 4.988 | 88.484 | 88.599 |
| sassari | 6 | srukf | 4.988 | 88.484 | 88.599 |
| sassari | 6 | madgwick | 3.377 | 31.782 | 32.441 |
| sassari | 6 | mahony | 8.099 | 69.223 | 69.950 |
| sassari | 6 | vqf | 2.095 | 3.637 | 4.323 |
| sassari | 9 | ukf | 4.000 | 5.082 | 6.536 |
| sassari | 9 | srukf | 4.000 | 5.082 | 6.536 |
| sassari | 9 | madgwick | 3.296 | 3.848 | 5.143 |
| sassari | 9 | mahony | 7.548 | 11.583 | 13.942 |
| sassari | 9 | vqf | 2.095 | 3.573 | 4.235 |

These results describe the supplied fixed-noise UKF model and selected reference settings. They are not an optimized ranking of filter families.

The fixed-noise UKFs do not reject disturbed measurements. Large six-axis heading errors are retained: heading is unobservable and gyro bias/model errors can accumulate during motion.

[Methodology and reproduction](../attitude-benchmark.md) · [Per-trial metrics](attitude-per-trial.csv) · [Source hashes, dependencies and provenance](attitude-provenance.json)

Measured locally with GCC 14.2.0, CMake 3.31.10, Release, double precision. The provenance revision identifies the checked-out commit; the dirty flag and source hashes identify local changes. CI regenerates reports for the checked-out commit. Timing is diagnostic and not a speed ranking.
