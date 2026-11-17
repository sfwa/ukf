STM32H735 on-device **compute-only** measurements.

Input: synthetic_500hz_rest_motion_magnetic_disturbance; 6001 samples, mean input interval 2 ms; float64 filters. The comparison budget is 2 ms per update; input timestamps are not retimed.

| Filter | Pass | Mean µs | P95 µs | P99 µs | Max µs | Max / 2 ms | Updates > 2 ms | Replay valid |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| UKF | 1 | 312.309 | 364.000 | 364.000 | 365.000 | 18.25% | 0 | yes |
| SR-UKF | 1 | 326.439 | 378.000 | 379.000 | 381.000 | 19.05% | 0 | yes |
| UKF | 2 | 312.307 | 364.000 | 364.000 | 365.000 | 18.25% | 0 | yes |
| SR-UKF | 2 | 326.433 | 378.000 | 379.000 | 382.000 | 19.10% | 0 | yes |

Each row is one replay; percentiles are not combined across passes. Sensor I/O, output logging, sample loading and task-release scheduling are outside the timed region.

All four replays completed within the observed compute budget. This does not establish a worst-case end-to-end deadline guarantee.
