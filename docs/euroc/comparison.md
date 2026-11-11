# Recorded EuRoC MAV attitude comparison

Recorded 2026-11-12. **88 evaluations**: all eleven sequences, 8 algorithms/variants, six-axis IMU input only.

Values are arithmetic means of per-sequence RMSEs in degrees. Heading and total error use one fixed initial yaw alignment; absolute heading is unobservable. The primary process-model target is inclination accuracy.

## Primary validation: externally constrained Vicon attitude

V1_03_difficult and V2_01_easy, V2_02_medium, V2_03_difficult were not used for tuning; the methodology records two structural choices made after validation results had been seen. The reference combines Vicon poses and IMU in an offline batch estimate; it is not wholly independent of the input IMU.

| Algorithm | Sequences | Inclination | Heading | Total |
|---|---:|---:|---:|---:|
| ukf | 4 | 1.511 | 76.698 | 76.718 |
| srukf | 4 | 1.511 | 76.698 | 76.718 |
| mav-srukf | 4 | 1.172 | 2.901 | 3.157 |
| mav-srukf-gravity | 4 | 1.543 | 5.129 | 5.400 |
| madgwick | 4 | 8.812 | 100.709 | 100.995 |
| mahony | 4 | 3.295 | 92.538 | 92.587 |
| vqf | 4 | 5.483 | 111.540 | 111.639 |
| vqf-mav | 4 | 1.615 | 88.528 | 88.785 |

## Primary Vicon subset, including development

The five Vicon sequences excluding V1_01_easy. V1_02_medium is a development recording.

| Algorithm | Sequences | Inclination | Heading | Total |
|---|---:|---:|---:|---:|
| ukf | 5 | 1.586 | 67.501 | 67.528 |
| srukf | 5 | 1.586 | 67.501 | 67.528 |
| mav-srukf | 5 | 1.109 | 2.646 | 2.894 |
| mav-srukf-gravity | 5 | 1.609 | 5.671 | 5.932 |
| madgwick | 5 | 8.589 | 95.678 | 95.973 |
| mahony | 5 | 3.564 | 87.548 | 87.615 |
| vqf | 5 | 5.464 | 107.252 | 107.357 |
| vqf-mav | 5 | 1.570 | 75.119 | 75.333 |

## All eleven recordings (descriptive)

Includes five Machine Hall recordings with no independent optical orientation observations and V1_01_easy with its documented canonical-reference orientation issue.

| Algorithm | Sequences | Inclination | Heading | Total |
|---|---:|---:|---:|---:|
| ukf | 11 | 1.717 | 66.861 | 66.886 |
| srukf | 11 | 1.717 | 66.861 | 66.886 |
| mav-srukf | 11 | 1.228 | 3.461 | 3.710 |
| mav-srukf-gravity | 11 | 1.801 | 8.116 | 8.433 |
| madgwick | 11 | 6.985 | 104.757 | 104.939 |
| mahony | 11 | 3.184 | 100.243 | 100.289 |
| vqf | 11 | 5.106 | 109.725 | 109.812 |
| vqf-mav | 11 | 1.946 | 60.936 | 61.127 |

## Per-sequence model comparison

Inclination RMSE in degrees. The gravity control is the MAV estimator with its rotor model switched off; the default and adapted VQF configurations remain separately visible.

| Sequence | Split | Gravity control | MAV SR-UKF | VQF | Adapted VQF |
|---|---|---:|---:|---:|---:|
| MH_01_easy | development | 1.231 | 1.084 | 5.219 | 2.127 |
| MH_02_easy | validation | 1.236 | 1.095 | 5.213 | 2.249 |
| MH_03_medium | validation | 2.607 | 0.932 | 5.304 | 2.183 |
| MH_04_difficult | validation | 1.809 | 0.914 | 5.241 | 2.681 |
| MH_05_difficult | validation | 1.641 | 0.911 | 5.077 | 1.467 |
| V1_01_easy | validation | 3.239 | 3.025 | 2.795 | 2.847 |
| V1_02_medium | development | 1.875 | 0.856 | 5.384 | 1.391 |
| V1_03_difficult | validation | 1.658 | 1.067 | 4.554 | 1.417 |
| V2_01_easy | validation | 1.132 | 1.255 | 5.650 | 1.649 |
| V2_02_medium | validation | 1.634 | 1.299 | 5.724 | 1.555 |
| V2_03_difficult | validation | 1.749 | 1.066 | 6.005 | 1.838 |

Development used MH_01_easy and V1_02_medium. All estimators process complete sensor recordings; no ground-truth position, velocity or bias enters an estimator, and no reference-derived flight mask selects favorable intervals. Configuration changes and model limitations are described in the methodology.

[Methodology](../euroc-benchmark.md) · [Per-sequence metrics](per-sequence.csv) · [Configuration and provenance](provenance.json)

Release build on macOS 27.2 arm64 with Apple clang 21.0.0, C++17, warnings treated as errors, patched Eigen 3.2.8, UKF double precision and Python 3.12.16. The recorded source and executable hashes identify the replayed implementation. Timing is diagnostic, not a cross-language speed ranking.
