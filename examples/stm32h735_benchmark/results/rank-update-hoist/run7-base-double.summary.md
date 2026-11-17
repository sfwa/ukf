STM32H735 on-device **compute-only** measurements.

Input: synthetic_500hz_rest_motion_magnetic_disturbance; 6001 samples, mean input interval 2 ms; float64 filters. The comparison budget is 2 ms per update; input timestamps are not retimed.

| Filter | Pass | Mean µs | P95 µs | P99 µs | Max µs | Max / 2 ms | Updates > 2 ms | Replay valid |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| UKF | 1 | 312.418 | 364.000 | 364.000 | 365.000 | 18.25% | 0 | yes |
| SR-UKF | 1 | 330.806 | 386.000 | 386.000 | 389.000 | 19.45% | 0 | yes |
| UKF | 2 | 312.429 | 364.000 | 364.000 | 365.000 | 18.25% | 0 | yes |
| SR-UKF | 2 | 330.807 | 385.000 | 386.000 | 389.000 | 19.45% | 0 | yes |

Each row is one replay; percentiles are not combined across passes. Sensor I/O, output logging, sample loading and task-release scheduling are outside the timed region.

All four replays completed within the observed compute budget. This does not establish a worst-case end-to-end deadline guarantee.
