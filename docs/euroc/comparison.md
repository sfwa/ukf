# Recorded EuRoC MAV attitude comparison

Recorded 2026-11-11. **121 evaluations**: all eleven sequences, eleven algorithms/variants, six-axis IMU input only.

Values are arithmetic means of per-sequence RMSEs in degrees. Heading and total error use one fixed initial yaw alignment; absolute heading is unobservable. The primary process-model target is inclination accuracy.

## Primary validation: externally constrained Vicon attitude

V1_03_difficult and V2_01_easy, V2_02_medium, V2_03_difficult were reserved from model iteration. The reference combines Vicon poses and IMU in an offline batch estimate; it is not wholly independent of the input IMU.

| Algorithm | Sequences | Inclination | Heading | Total |
|---|---:|---:|---:|---:|
| ukf | 4 | 4.801 | 69.437 | 69.638 |
| srukf | 4 | 4.801 | 69.437 | 69.638 |
| robust-rest | 4 | 5.221 | 102.869 | 102.962 |
| robust-acc | 4 | 13.468 | 102.018 | 102.644 |
| robust | 4 | 13.517 | 102.018 | 102.648 |
| mav-srukf | 4 | 1.130 | 2.844 | 3.116 |
| mav-srukf-gravity | 4 | 1.631 | 3.311 | 3.779 |
| madgwick | 4 | 8.812 | 100.709 | 100.995 |
| mahony | 4 | 3.295 | 92.538 | 92.587 |
| vqf | 4 | 5.483 | 111.540 | 111.639 |
| vqf-mav | 4 | 1.615 | 88.528 | 88.785 |

## Primary Vicon subset, including development

The five Vicon sequences excluding V1_01_easy. V1_02_medium is a development recording.

| Algorithm | Sequences | Inclination | Heading | Total |
|---|---:|---:|---:|---:|
| ukf | 5 | 5.004 | 60.366 | 60.663 |
| srukf | 5 | 5.004 | 60.366 | 60.663 |
| robust-rest | 5 | 5.313 | 97.684 | 97.794 |
| robust-acc | 5 | 13.200 | 96.480 | 97.148 |
| robust | 5 | 13.244 | 96.480 | 97.152 |
| mav-srukf | 5 | 1.056 | 2.531 | 2.791 |
| mav-srukf-gravity | 5 | 1.540 | 2.891 | 3.360 |
| madgwick | 5 | 8.589 | 95.678 | 95.973 |
| mahony | 5 | 3.564 | 87.548 | 87.615 |
| vqf | 5 | 5.464 | 107.252 | 107.357 |
| vqf-mav | 5 | 1.570 | 75.119 | 75.333 |

## All eleven recordings (descriptive)

Includes five Machine Hall recordings with no independent optical orientation observations and V1_01_easy with its documented canonical-reference orientation issue.

| Algorithm | Sequences | Inclination | Heading | Total |
|---|---:|---:|---:|---:|
| ukf | 11 | 4.000 | 83.389 | 83.550 |
| srukf | 11 | 3.999 | 82.754 | 82.915 |
| robust-rest | 11 | 4.783 | 105.663 | 105.742 |
| robust-acc | 11 | 13.035 | 105.176 | 105.741 |
| robust | 11 | 13.063 | 105.176 | 105.743 |
| mav-srukf | 11 | 1.393 | 8.523 | 8.858 |
| mav-srukf-gravity | 11 | 1.881 | 10.532 | 10.946 |
| madgwick | 11 | 6.985 | 104.757 | 104.939 |
| mahony | 11 | 3.184 | 100.243 | 100.289 |
| vqf | 11 | 5.106 | 109.725 | 109.812 |
| vqf-mav | 11 | 1.946 | 60.936 | 61.127 |

## Per-sequence model comparison

Inclination RMSE in degrees. The gravity control shares the MAV estimator configuration; the default and adapted VQF configurations remain separately visible.

| Sequence | Split | Gravity control | MAV SR-UKF | VQF | Adapted VQF |
|---|---|---:|---:|---:|---:|
| MH_01_easy | development | 1.277 | 1.194 | 5.219 | 2.127 |
| MH_02_easy | validation | 2.749 | 2.739 | 5.213 | 2.249 |
| MH_03_medium | validation | 1.944 | 1.090 | 5.304 | 2.183 |
| MH_04_difficult | validation | 1.671 | 1.118 | 5.241 | 2.681 |
| MH_05_difficult | validation | 1.390 | 0.942 | 5.077 | 1.467 |
| V1_01_easy | validation | 3.961 | 2.964 | 2.795 | 2.847 |
| V1_02_medium | development | 1.176 | 0.759 | 5.384 | 1.391 |
| V1_03_difficult | validation | 1.385 | 0.966 | 4.554 | 1.417 |
| V2_01_easy | validation | 1.328 | 1.145 | 5.650 | 1.649 |
| V2_02_medium | validation | 2.035 | 1.285 | 5.724 | 1.555 |
| V2_03_difficult | validation | 1.775 | 1.123 | 6.005 | 1.838 |

Development used MH_01_easy and V1_02_medium. All estimators process complete sensor recordings; no ground-truth position, velocity or bias enters an estimator, and no reference-derived flight mask selects favorable intervals. Configuration changes and model limitations are described in the methodology.

[Methodology](../euroc-benchmark.md) · [Per-sequence metrics](per-sequence.csv) · [Configuration and provenance](provenance.json)

Linux x86_64, GCC 14.2.0, CMake 3.31.10, Release, double precision; Python 3.12.14 and pinned reference packages. The recorded source and executable hashes identify the replayed implementation. Timing is diagnostic, not a cross-language speed ranking.
