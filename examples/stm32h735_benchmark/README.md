# STM32H735 attitude-estimator benchmark

This bare-metal application runs the **same float32 UKF and SR-UKF estimator
and workload** as the [ESP32-S31 benchmark](../esp32s31_benchmark/), by including
its shared `benchmark.h` directly. It compares `update()` plus `attitude()`
against a **2,000 µs computation budget** on the STM32H735G-DK Discovery kit.

## Measured STM32H735 hardware results

On 2026-11-14, the cross-compiled application completed on a physical
**STM32H735IGK6, revision Z, at 550 MHz**. These are **hardware measurements**,
using the shared synthetic sensor fixture; they are not simulated timings or
a replay of a real sensor recording. The estimator uses float32, hardware
floating point, `-O2`, and enabled instruction/data caches.

The measured build used library headers from `b75b6a6` (source base
`878c1dc` plus this platform port), before the SR-UKF optimizations in
[#92](https://github.com/sfwa/ukf/pull/92). These recorded timings do not
describe that newer core revision; rebuild and rerun to measure it.

| Filter | Pass | Updates | Mean (µs) | p95 (µs) | p99 (µs) | Maximum (µs) | Updates >2 ms |
|---|---:|---:|---:|---:|---:|---:|---:|
| UKF | 1 | 6,000 | 161.621 | 187 | 188 | 188 | 0 |
| SR-UKF | 1 | 6,000 | 190.616 | 226 | 227 | 228 | 0 |
| SR-UKF | 2 | 6,000 | 190.624 | 226 | 227 | 228 | 0 |
| UKF | 2 | 6,000 | 161.625 | 187 | 188 | 189 | 0 |

Across the two equal-length replays, mean time was **161.623 µs for UKF**
and **190.620 µs for SR-UKF**: SR-UKF took approximately **17.9% longer**.
Percentiles above are per replay, not combined. Both completed all 12,000
measured updates per filter without exceeding 2 ms. The largest observed
UKF update consumed **9.45%** of the budget, leaving **1,811 µs**; SR-UKF
consumed **11.40%**, leaving **1,772 µs**. This supports a 500 Hz estimator
computation rate under this workload; it is not a complete system deadline
or worst-case timing guarantee.

All four replays passed numerical checks. Minimum free stack was 53,664 of
65,536 bytes, with no capture overflow or recorded CPU fault. The 5,005-byte
UART capture matched the retained DTCM output byte for byte after completion.
The original 1 MiB internal flash was backed up twice with matching hashes
before programming; the backup is under `SFWA/toolchains/board-backups`.
The board is left with the benchmark application in internal flash.

Raw evidence: [UART JSONL](results/stm32h735-550mhz.jsonl) and
[build/run manifest](results/stm32h735-550mhz-manifest.json), including firmware
and source hashes. Reproduce the report from the repository root with:

```sh
python examples/esp32s31_benchmark/summarize.py \
  examples/stm32h735_benchmark/results/stm32h735-550mhz.jsonl
```

## Toolchain and build

The tools are installed in `SFWA/toolchains`, alongside the `SFWA/ukf` checkout
and outside the UKF repository:

| Component | Installation directory / revision |
|---|---|
| Arm GNU Toolchain 14.3.Rel1, GCC 14.3.1 | `arm-gnu-toolchain-14.3.rel1-darwin-arm64-arm-none-eabi` |
| Official STM32CubeProgrammer 2.22.0 | `programmer-2.22.0+st.1-aarch64-darwin` |
| xPack OpenOCD 0.12.0-7 | `xpack-openocd-0.12.0-7` |
| STM32CubeH7 v1.13.0 | `STM32CubeH7`, commit `5abb9764b32e11a6557b90bf39531528019b5761` |
| STM32H7 HAL submodule | `a1996eed9172b59887bafaaa0ea1816ea14d48b5` |
| STM32H7 CMSIS device submodule | `8f922cdc7cc6de2344e75ddd657889f4ff761790` |

The workspace's `toolchains/arm-tools-manifest.json` and
`toolchains/st-vendor-tools-manifest.json` record download URLs, package
versions and archive hashes. The compiler and debugger are
native macOS arm64 installations; no global package or shell-profile changes
are required. On another host, install its Arm GNU Toolchain package and pass
`-DARM_TOOLCHAIN_ROOT=/path/to/toolchain` when configuring.

To reproduce the minimal SDK checkout from the SFWA directory:

```sh
git clone --depth 1 --branch v1.13.0 --filter=blob:none --sparse \
  https://github.com/STMicroelectronics/STM32CubeH7.git toolchains/STM32CubeH7
git -C toolchains/STM32CubeH7 sparse-checkout set \
  Drivers/CMSIS Drivers/STM32H7xx_HAL_Driver Projects/STM32H735G-DK
git -C toolchains/STM32CubeH7 submodule update --init --depth 1 \
  Drivers/CMSIS/Device/ST/STM32H7xx Drivers/STM32H7xx_HAL_Driver
git -C toolchains/STM32CubeH7 rev-parse HEAD
```

From the UKF repository root, prepare its patched Eigen 3.2.8 dependency and
cross-compile. In this workspace, `.venv/bin/cmake` is the installed CMake:

```sh
.venv/bin/cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
.venv/bin/cmake --build build --target eigen3
.venv/bin/cmake -S examples/stm32h735_benchmark -B build-stm32h735 \
  -DCMAKE_TOOLCHAIN_FILE=toolchain.cmake -DCMAKE_BUILD_TYPE=Release
.venv/bin/cmake --build build-stm32h735 --parallel 2
```

The build emits `ukf_stm32h735_benchmark.elf`, `.bin`, `.map`, and
`compile_commands.json` in `build-stm32h735`. Override `STM32CUBE_ROOT` or
`UKF_EIGEN_DIR` to use equivalent pinned dependencies elsewhere. The ordinary
repository host build skips this cross-compilation target.

## Flash, capture, and summarize

Power the board off, set **SW1 (BOOT0) to 0**, and power it back on through its
**ST-LINK USB** connector before running the benchmark. SW1=1 selects the ROM
bootloader instead of the programmed application. This board has bootloader
V9.2, which can lock `PWR_CR3` until power-off; directly starting flash code
after visiting that bootloader is insufficient to clear the lock. See
[AN2606, STM32H72xxx/73xxx bootloader versions](https://www.st.com/resource/en/application_note/an2606-stm32-microcontroller-system-memory-boot-mode-stmicroelectronics.pdf),
table 144 in revision 71. The firmware expects the board's default direct-SMPS
hardware configuration.

On this Mac, connecting through a USB hub resolved the original USB transport
failure. The embedded ST-LINK-V3E was then updated to **V3J17M11**, and official
STM32CubeProgrammer 2.22.0 successfully identified the MCU and programmed and
verified internal flash. OpenOCD is installed, but the vendor programmer is
the verified flashing route for this setup.

Unmount the ST-LINK programming volume before a reset or firmware update if
the host has mounted it. The following command replaces
the resident application in internal flash with this benchmark; it does not
program external Octo-SPI flash or change option bytes. From the repository
root, after the ST-LINK is accessible:

```sh
../toolchains/programmer-2.22.0+st.1-aarch64-darwin/bin/STM32_Programmer_CLI \
  --connect port=SWD mode=UR reset=HWrst freq=1000 \
  --download build-stm32h735/ukf_stm32h735_benchmark.elf --verify -hardRst
```

Capture the ST-LINK virtual COM port at **115200 baud, 8N1**, starting the
capture before the programmer's `-hardRst` or a subsequent board reset. The
device appears as `/dev/cu.usbmodem...` on macOS; select the connected board's
actual port. Firmware output uses USART3 on PD8/PD9,
alternate function 7. Each boot produces one complete JSON suite and then
stops. Summarize a captured log with:

```sh
python examples/esp32s31_benchmark/summarize.py stm32h735.jsonl
```

The shared summarizer recognizes `platform:"stm32h735"` as physical hardware
and requires all four replay records plus `platform_complete`. Its exit status
is nonzero for numerical failure, incomplete replay, or an observed update
over 2 ms. Do not pass `--emulated` for this board run.

The same output is retained in a 16 KiB DTCM buffer for SWD capture when the
virtual COM port is unavailable. Obtain current addresses from the ELF's
symbols: `benchmark_log`, `benchmark_log_length`, `benchmark_log_overflow`,
and `benchmark_complete`. DTCM is uncached, so debugger reads do not require
cache maintenance. Completion is 0 while running, 1 for a valid completed
suite, 2 for replay failure, or 3 for setup/fault failure. A value of 1 is
**not** a deadline verdict; summarize the recorded durations. Read exactly
`benchmark_log_length` bytes and require zero overflow. Let the MCU run
uninterrupted during measurement; halt only after completion to collect data.

## Configuration and interpretation

- Cortex-M7 at 550 MHz: the board's 25 MHz HSE, PLL M=5/N=110/P=1, voltage
  scale 0, AHB at 275 MHz, APBs at 137.5 MHz, and flash latency 3.
- **Direct SMPS supply**, matching the board's default solder-bridge layout.
  Vendor startup selects this before RAM initialization; the application
  verifies the resulting CPU/bus clock configuration. `USER_VECT_TAB_ADDRESS`
  explicitly places the interrupt vector table at `0x08000000`.
- `-O2`, `-mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard`, float32 estimator
  math, `EIGEN_NO_MALLOC`, and no fast-math option. Full newlib preserves the
  shared harness's size and 64-bit integer JSON formatting.
- Instruction/data caches enabled; code in internal flash, data and the
  48,000-byte timing buffer in AXI SRAM, and a 64 KiB stack in DTCM. No external
  RAM, display, networking, sensor drivers, or RTOS is started. The 1 kHz
  SysTick interrupt remains enabled and is included when it interrupts timing.

The timer extends the 32-bit DWT cycle counter before converting it to integer
microseconds using the verified CPU clock. Reads occur at each update and
yield, within its approximately 7.8-second wrap interval. Resolution is 1 µs;
timer overhead is retained. Input generation, logging, validation and sorting
remain outside the measured region. Playback is not paced at 500 Hz.

The default input is the same **synthetic** 6,001-sample, 500 Hz fixture:
rest, rotation, magnetic disturbance, outage and recovery. It exercises the
recorded-data estimator but is not an actual recording or an accuracy test.
Each fresh replay follows 256 disposable warm-up updates. The four measured
replays run UKF, SR-UKF, SR-UKF, UKF, with 6,000 updates each.

For a recording, use the existing `export_recording.py` and pass
`-DUKF_BENCHMARK_RECORDING_HEADER=/absolute/path/recording.h` to CMake. Its
timestamps are preserved. See the shared benchmark's
[recorded CSV instructions](../esp32s31_benchmark/README.md#replay-a-recorded-csv).

`platform_complete.suite_ok` establishes numerical completion and intact
capture/stack checks. Budget status comes from each replay's timing results.
Passing this compute-only workload does not establish a worst-case or complete
500 Hz system deadline: sensor transfers, scheduling and other application
work must also fit. ESP32-S31 emulator timings are not hardware measurements
and should not be used to infer a physical cross-MCU speed ratio.

References: [ST board and user documentation](https://www.st.com/en/evaluation-tools/stm32h735g-dk.html),
[STM32H735 550 MHz specification](https://www.st.com/en/microcontrollers-microprocessors/stm32h725-735.html),
[pinned ST clock example](https://github.com/STMicroelectronics/STM32CubeH7/blob/5abb9764b32e11a6557b90bf39531528019b5761/Projects/STM32H735G-DK/Examples/UART/UART_Printf/Src/main.c),
and [CMSIS DWT register definitions](https://github.com/STMicroelectronics/STM32CubeH7/blob/5abb9764b32e11a6557b90bf39531528019b5761/Drivers/CMSIS/Include/core_cm7.h).
The ST example uses 520 MHz (N=104); this benchmark uses N=110 for 550 MHz.
