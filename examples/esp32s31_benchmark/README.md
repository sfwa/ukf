# ESP32-S31 attitude-estimator benchmark

This ESP-IDF application compares the recorded-data example's **UKF and
SR-UKF in single precision** against a **2,000 µs computation budget**.
It runs the same [`Attitude_Estimator`](../attitude_benchmark/attitude_estimator.h),
including its RK4 process model, rest detection, gyro-bias estimation,
magnetic disturbance checks and separate heading filter. No smaller substitute
filter is used. The main filter has nine covariance dimensions and 19 sigma
points; the heading filter has one dimension.

The firmware measures computation on the board or in an emulator. The host
executable checks correctness and reporting; **host timing is not an ESP32-S31
estimate**. [Simulated cloud results](#simulated-cloud-results-2026-11-14) are
recorded below; no physical ESP32-S31 timing measurements are available.

## Build and run on ESP32-S31

Use an ESP-IDF checkout/toolchain with ESP32-S31 support; consult Espressif's
[support status](https://developer.espressif.com/hardware/esp32s31/).
The application reuses the repository's patched Eigen dependency, downloaded
by CMake on the first build. From the repository root:

```sh
. /path/to/esp-idf/export.sh
cd examples/esp32s31_benchmark
idf.py set-target esp32s31
idf.py build
idf.py -p /dev/cu.usbmodemYOUR_BOARD flash monitor
```

The benchmark starts automatically at boot and stops after four measured
replays. Reset the board to run it again. To save a run, start monitoring
before resetting the board:

```sh
idf.py -p /dev/cu.usbmodemYOUR_BOARD monitor | tee benchmark.log
# Exit the monitor with Ctrl-]. Then:
python summarize.py benchmark.log
```

The summary prints UKF/SR-UKF mean, p95, p99, maximum, maximum budget usage,
and overruns for each pass. Capture one complete boot per log; the summary
rejects mixed or incomplete suites. Its exit status is nonzero if a complete
hardware run fails numerically or exceeds the computation budget.

Defaults select 320 MHz, `-O2`, two CPU cores, internal SRAM, no dynamic
frequency scaling and no PSRAM. The benchmark runs on core 1 in a task with
a 64 KiB stack. Compile-time checks reject incompatible settings. Wi-Fi,
Bluetooth and sensor drivers are not started. The platform record includes
the chip revision, IDF/compiler versions, configured CPU frequency, task core,
and heap/stack information.

Use `idf.py size` to inspect the image and memory footprint. Large embedded
recordings may require a larger app partition; ESP-IDF's
`CONFIG_PARTITION_TABLE_SINGLE_APP_LARGE=y` selects the 1.5 MiB single-app
layout when the attached flash can accommodate it.

## SIMULATED cloud results (2026-11-14)

**These are emulator timings, not measurements of physical ESP32-S31 hardware.
They do not establish that a board can meet a 2 ms deadline or sustain 500 Hz.**
The cross-compiled firmware ran in official Espressif `esp-emu` 1.0.0 on
Linux x86_64 with normal ROM/bootloader startup. It used the recorded-data
estimator with the default **synthetic** 6,001-sample, 500 Hz fixture, not a
sensor recording.

All four replays completed: 6,000 measured updates each, 24,000 total, with
`platform_complete.suite_ok=true` and emulator exit code 0. `suite_ok` confirms
numerical replay success, not deadline compliance. Each row below is a separate
replay; percentiles are not combined across passes. All timings are **simulated
microseconds**.

| Filter | Pass | Mean µs | P95 µs | P99 µs | Max µs | Updates > 2 ms |
|---|---:|---:|---:|---:|---:|---:|
| UKF | 1 | 248.573 | 290 | 290 | 290 | 0 |
| UKF | 2 | 248.573 | 290 | 290 | 290 | 0 |
| SR-UKF | 1 | 280.511 | 331 | 331 | 333 | 0 |
| SR-UKF | 2 | 280.507 | 331 | 331 | 333 | 0 |

SR-UKF used about 13% more simulated time per update than UKF on this workload.
The largest simulated interval was 333 µs, or 16.65% of the 2,000 µs comparison
budget. This comparison does not assume a cycle-accurate model of silicon;
physical performance and end-to-end scheduling still require board measurements.

The source was base commit `e826de4dc65212035f1c6a938f03c6afc95f74cd` plus the
benchmark snapshot, whose gzip SHA-256 was
`9079714c13585a277dbc6b3a8055105db48102447cd7e25dc30b441c1984cf75`.
The verified build configuration was ESP-IDF **v6.1** at commit
`fff9895c82d744c7237be8847347bdd1b07c6643`, GCC **15.2.0** from
`esp-15.2.0_20251204`, `UKF_SINGLE_PRECISION`, `EIGEN_NO_MALLOC`, `-mabi=ilp32f`
and `-O2`; UART confirmed float32, 320 MHz and execution on core 1.
The application binary was 449,776 bytes (449,408 bytes before padding), and
the merged flash image was 515,312 bytes. All 64 relevant host tests passed
in the cloud: 46 attitude-benchmark tests and 18 ESP32-S31 benchmark tests.

The seven complete [UART JSON records](results/esp-emu-v1.0.0.jsonl) are
preserved in original execution order. Boot chatter is omitted; JSON text is
unchanged, with LF line endings and a final newline. The file's SHA-256 is
`f9a6856a5d9f4d775748bb3d077e12280cbc8789c7a28b478bbac07cae32d075`.
From this directory, reproduce the report with:

```sh
python summarize.py --emulated results/esp-emu-v1.0.0.jsonl
```

Always use `--emulated` for emulator UART output. The firmware still identifies
its target as `esp32s31` under emulation, so the platform field alone cannot
distinguish an emulator from a board. `--emulated`
labels the report as simulated and excludes simulated deadline results from
the hardware pass/fail decision; numerical replay failures still fail.

No filter, workload, C++ or CMake correction was needed for the cloud run.
The compiler was killed during the parallel build; serial compilation
with `IDF_PY_BUILD_JOBS=1` succeeded. Cloud host validation used
`-DUKF_WARNINGS_AS_ERRORS=OFF` for an Eigen dependency warning; this was specific
to that environment, not a change to the repository's default warning policy.
With this ESP-IDF version, `idf.py --preview merge-bin -o merged_flash.bin`
writes `build/merged_flash.bin`, because the output path is relative to the
build directory.

## What is measured

Each sample is prepared before timing. The timed interval contains:

```cpp
estimator.update(previous, current);
estimator.attitude();
```

`esp_timer_get_time()` measures elapsed microseconds, including interrupts
that occur inside this interval. Timer overhead is retained, not subtracted.
A compiler barrier keeps the output calculation inside the interval.
Initialization is measured separately. Input generation/decoding, output
validation, statistics, sorting, logging and task delays are outside it.
The filter math uses `float`; input timestamps stay `double` to preserve time
resolution, with elapsed time converted to `real_t` by the existing estimator.

Each replay warms the code with up to 256 updates on a disposable estimator,
then creates a fresh estimator and measures every update from the start.
The order is UKF, SR-UKF, SR-UKF, UKF, so both filters run first in one pass.
State is never carried across replays. A short task delay every 32 samples
lets the idle task/watchdog run. Playback is **not paced at 500 Hz**: this
measures computation needed per sample, not a live acquisition/scheduling loop.

JSON lines report:

- Mean, minimum, nearest-rank p50/p95/p99, maximum, and the number of updates
  taking **more than 2,000 µs**. Exactly 2,000 µs is within the stated budget.
- Mean/maximum budget usage and minimum observed headroom, which is negative
  if the maximum exceeds 2 ms.
- Separate counts, means and maxima for rest/motion and accepted/rejected
  magnetic readings. Rest adds a three-axis bias measurement and can cost more.
- Numerical validity, warm-up/initialization status, failing sample indices,
  and a consumed output checksum. Failed/partial replays cannot pass the budget.

`platform_complete.suite_ok` means the replay succeeded numerically. The
per-run `observed_compute_within_budget` fields give the timing verdict.
A pass establishes only that the observed updates fit; it is not a worst-case
execution-time proof. Rare covariance fallback paths are not deliberately
forced. Repeat on the intended board and with the application's competing
workloads. Sensor transfer time, scheduling latency and output delivery must
also fit in a real system's 2 ms period.

## Default 500 Hz workload

Without an input header, both filters receive the same deterministic synthetic
sequence: 6,001 samples at 2 ms intervals (12 s of sensor time), containing
tilted rest with gyro bias, multi-axis rotation, a changed magnetic field,
a magnetic outage, and recovery. This exercises the recorded-data estimator
without requiring a sensor or redistributing a dataset. It is a workload
fixture, not an accuracy benchmark or an actual recording.

The shared timing buffer uses 48,000 bytes of internal SRAM for this sequence.
The estimator objects and temporaries use the task stack. `EIGEN_NO_MALLOC`
is defined; state and bounded measurement matrices have fixed-capacity storage.

## Replay a recorded CSV

Use the sensor CSV format accepted by the existing recorded-data example:

```text
time,gx,gy,gz,ax,ay,az,mx,my,mz
```

Time is in seconds, gyro in rad/s, acceleration in m/s², and magnetic field
in any consistent unit. All sensor axes must already share the same body
frame. Export a contiguous window from the repository root:

```sh
python examples/esp32s31_benchmark/export_recording.py sensors.csv /tmp/ukf-recording.h --start 0 --count 6001
cd examples/esp32s31_benchmark
idf.py -DUKF_BENCHMARK_RECORDING_HEADER=/tmp/ukf-recording.h build
idf.py -p /dev/cu.usbmodemYOUR_BOARD flash monitor
```

The exporter validates the CSV, rounds sensor readings to float, rebases
timestamps, and embeds a source SHA-256 and starting row. It preserves elapsed
times: a 100 Hz recording stays 100 Hz and is **not** silently retimed to
500 Hz. Its update computation is still compared with 2 ms, but rate-dependent
filter behaviour may differ. For a direct 500 Hz workload, use a recording
acquired at that rate or the default fixture.

Select 2 to 8,193 samples. `--start` is zero-based; `--count` must be available
in full. Omitting `--count` selects the remaining rows and fails if the cap
is exceeded. The cap bounds the timing buffer at 64 KiB. The recording itself
is const data in flash; it is not loaded wholesale into RAM. Sample conversion
is outside timing. The starting window must permit initialization from the
accelerometer; choose a representative lead-in for the filter to settle.

Clear the cached header selection to restore the default fixture:

```sh
idf.py -DUKF_BENCHMARK_RECORDING_HEADER= build
```

## Host validation

Normal repository builds require no ESP-IDF installation:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target esp32s31_benchmark_host --parallel 2
python -m unittest discover -s examples/esp32s31_benchmark -p 'test_*.py' -v
build/examples/esp32s31_benchmark/esp32s31_benchmark_host > /tmp/ukf-host.jsonl
python examples/esp32s31_benchmark/summarize.py /tmp/ukf-host.jsonl
```

The host target also accepts CMake's `-DUKF_BENCHMARK_RECORDING_HEADER=...`
option. Tests cover successful replay by both cores, branch coverage on the
default fixture, repeated-pass determinism, numerical agreement, deadline
boundary accounting, percentile selection, and recording export/reporting.
They do not assert desktop timing thresholds. The existing attitude benchmark
CI runs these host checks; it does not build, flash or time an ESP32-S31.

Espressif references: [timing and optimization](https://docs.espressif.com/projects/esp-idf/en/latest/esp32s31/api-guides/performance/speed.html),
[ESP Timer](https://docs.espressif.com/projects/esp-idf/en/latest/esp32s31/api-reference/system/esp_timer.html).
