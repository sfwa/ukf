# Rank-update reciprocal hoist on the STM32H735 (2026-11-17)

Base is `d0e10a9`; "hoist" is the same tree with `rank_update_lower` in
`include/UKF/RankUpdate.h` using one reciprocal and multiplies instead of the
two divisions in the update branch, and a reciprocal computed once per column
instead of a per-row division in the downdate branch. The change is
adopted in the same commit as these results;
[manifest.json](manifest.json) records both file hashes, the toolchain, flags,
board and firmware hashes.

All 8 runs: `-O2`, STM32H735IGK6 rev Z at 550 MHz, programmed and run in the
order base, hoist, base, hoist for float, then the same for double.
Each mean below combines 2 runs × 2 replays of 6,000 updates.

| Precision | Filter | Base mean µs | Hoist mean µs | Change |
|---|---|---:|---:|---:|
| float32 | UKF | 166.269 | 166.670 | +0.241% |
| float32 | SR-UKF | 165.529 | 163.564 | −1.187% |
| float64 | UKF | 312.416 | 312.309 | −0.034% |
| float64 | SR-UKF | 330.811 | 326.439 | −1.322% |

Repeat runs of the same firmware agreed within 0.025 µs on means; SR-UKF
p95/p99/max dropped by 3 µs (float) and 7–8 µs (double). Only the SR-UKF
calls `rank_update_lower` (through `Detail::update_root` for the central
sigma point and the posterior downdates, and `Detail::whiten`). The standard UKF's
output is bit-identical between base and hoist, so its small shifts are a code
placement effect, not arithmetic. The float SR-UKF final quaternion and checksum
changed in the 7th significant figure (rounding of reciprocal multiplication);
the double results printed identically.

Code: `vdiv` count falls by 2 in each precision (float 143 → 141 `vdiv.f32`;
double 151 → 149 `vdiv.f64`). Out-of-line `rank_update_lower` grows by 8 bytes
(float) and 16 bytes (double), and the copy inlined into `whiten` by 8 bytes.
Total `.text` is unchanged (123,140 float, 128,316 double) because alignment
padding absorbs the growth.

Files: `runN-<build>.jsonl` are the UART logs (byte-identical to the DTCM copy),
`runN-<build>.summary.md` the `summarize.py` output, `<build>.size-vdiv.txt` the
`arm-none-eabi-size` output and whole-ELF `vdiv` counts, and
[runs-table.md](runs-table.md) all per-replay rows, stack, final attitudes and
checksums.
