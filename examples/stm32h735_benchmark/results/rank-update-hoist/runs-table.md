| Run | Build | Filter | Pass | Mean µs | P95 | P99 | Max | Replay valid |
|---:|---|---|---:|---:|---:|---:|---:|---|
| 1 | base-float | UKF | 1 | 166.267 | 191 | 191 | 193 | True |
| 1 | base-float | SR-UKF | 1 | 165.525 | 194 | 194 | 195 | True |
| 1 | base-float | SR-UKF | 2 | 165.529 | 194 | 194 | 195 | True |
| 1 | base-float | UKF | 2 | 166.267 | 191 | 192 | 193 | True |
| 2 | hoist-float | UKF | 1 | 166.662 | 191 | 192 | 193 | True |
| 2 | hoist-float | SR-UKF | 1 | 163.564 | 191 | 191 | 192 | True |
| 2 | hoist-float | SR-UKF | 2 | 163.562 | 191 | 191 | 192 | True |
| 2 | hoist-float | UKF | 2 | 166.670 | 191 | 192 | 193 | True |
| 3 | base-float | UKF | 1 | 166.276 | 191 | 192 | 192 | True |
| 3 | base-float | SR-UKF | 1 | 165.521 | 194 | 194 | 195 | True |
| 3 | base-float | SR-UKF | 2 | 165.540 | 194 | 194 | 195 | True |
| 3 | base-float | UKF | 2 | 166.268 | 191 | 192 | 192 | True |
| 4 | hoist-float | UKF | 1 | 166.672 | 191 | 192 | 193 | True |
| 4 | hoist-float | SR-UKF | 1 | 163.567 | 191 | 191 | 192 | True |
| 4 | hoist-float | SR-UKF | 2 | 163.564 | 191 | 191 | 192 | True |
| 4 | hoist-float | UKF | 2 | 166.675 | 191 | 192 | 193 | True |
| 5 | base-double | UKF | 1 | 312.407 | 364 | 364 | 365 | True |
| 5 | base-double | SR-UKF | 1 | 330.819 | 386 | 386 | 389 | True |
| 5 | base-double | SR-UKF | 2 | 330.811 | 386 | 386 | 389 | True |
| 5 | base-double | UKF | 2 | 312.411 | 364 | 364 | 365 | True |
| 6 | hoist-double | UKF | 1 | 312.309 | 364 | 364 | 365 | True |
| 6 | hoist-double | SR-UKF | 1 | 326.439 | 378 | 379 | 381 | True |
| 6 | hoist-double | SR-UKF | 2 | 326.433 | 378 | 379 | 382 | True |
| 6 | hoist-double | UKF | 2 | 312.307 | 364 | 364 | 365 | True |
| 7 | base-double | UKF | 1 | 312.418 | 364 | 364 | 365 | True |
| 7 | base-double | SR-UKF | 1 | 330.806 | 386 | 386 | 389 | True |
| 7 | base-double | SR-UKF | 2 | 330.807 | 385 | 386 | 389 | True |
| 7 | base-double | UKF | 2 | 312.429 | 364 | 364 | 365 | True |
| 8 | hoist-double | UKF | 1 | 312.309 | 364 | 364 | 365 | True |
| 8 | hoist-double | SR-UKF | 1 | 326.441 | 378 | 379 | 382 | True |
| 8 | hoist-double | SR-UKF | 2 | 326.443 | 378 | 379 | 382 | True |
| 8 | hoist-double | UKF | 2 | 312.310 | 364 | 364 | 365 | True |

| Run | Build | suite_ok | stack_min_free_bytes | log_overflow | optimization | summarize exit |
|---:|---|---|---:|---|---|---:|
| 1 | base-float | True | 55992 | False | O2 | 0 |
| 2 | hoist-float | True | 55992 | False | O2 | 0 |
| 3 | base-float | True | 55992 | False | O2 | 0 |
| 4 | hoist-float | True | 55992 | False | O2 | 0 |
| 5 | base-double | True | 47448 | False | O2 | 0 |
| 6 | hoist-double | True | 47448 | False | O2 | 0 |
| 7 | base-double | True | 47448 | False | O2 | 0 |
| 8 | hoist-double | True | 47448 | False | O2 | 0 |

| Build | Filter | final_attitude_wxyz | checksum | distinct outcomes |
|---|---|---|---:|---:|
| base-double | SR-UKF | [-0.509771947, 0.271514054, -0.360215446, 0.732569119] | 4000.11481903 | 1 |
| base-double | UKF | [-0.509771947, 0.271514054, -0.360215446, 0.732569119] | 4000.11481903 | 1 |
| base-float | SR-UKF | [-0.509767592, 0.27151382, -0.360224187, 0.732567966] | 4000.07817621 | 1 |
| base-float | UKF | [-0.509771526, 0.271512389, -0.360221654, 0.732567012] | 4000.07878503 | 1 |
| hoist-double | SR-UKF | [-0.509771947, 0.271514054, -0.360215446, 0.732569119] | 4000.11481903 | 1 |
| hoist-double | UKF | [-0.509771947, 0.271514054, -0.360215446, 0.732569119] | 4000.11481903 | 1 |
| hoist-float | SR-UKF | [-0.509768188, 0.27151373, -0.360224158, 0.732567608] | 4000.07729478 | 1 |
| hoist-float | UKF | [-0.509771526, 0.271512389, -0.360221654, 0.732567012] | 4000.07878503 | 1 |

| Precision | Filter | Base mean µs (4 replays) | Hoist mean µs (4 replays) | Change |
|---|---|---:|---:|---:|
| float | UKF | 166.269 | 166.670 | +0.241% |
| float | SR-UKF | 165.529 | 163.564 | -1.187% |
| double | UKF | 312.416 | 312.309 | -0.034% |
| double | SR-UKF | 330.811 | 326.439 | -1.322% |

Notes:
UART capture equals the DTCM log byte for byte in every run; all summarize.py exits 0.
