STM32H735 on-device **compute-only** measurements.

Input: synthetic_500hz_rest_motion_magnetic_disturbance; 6001 samples, mean input interval 2 ms; float32 filters. The comparison budget is 2 ms per update; input timestamps are not retimed.

| Filter | Pass | Mean µs | P95 µs | P99 µs | Max µs | Max / 2 ms | Updates > 2 ms | Replay valid |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| UKF | 1 | 166.672 | 191.000 | 192.000 | 193.000 | 9.65% | 0 | yes |
| SR-UKF | 1 | 163.567 | 191.000 | 191.000 | 192.000 | 9.60% | 0 | yes |
| UKF | 2 | 166.675 | 191.000 | 192.000 | 193.000 | 9.65% | 0 | yes |
| SR-UKF | 2 | 163.564 | 191.000 | 191.000 | 192.000 | 9.60% | 0 | yes |

Each row is one replay; percentiles are not combined across passes. Sensor I/O, output logging, sample loading and task-release scheduling are outside the timed region.

All four replays completed within the observed compute budget. This does not establish a worst-case end-to-end deadline guarantee.
