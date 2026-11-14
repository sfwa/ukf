# Embedded attitude-estimator code size

These cross-linked probes measure the **specific float32 attitude estimator
used by the STM32H735 and ESP32-S31 benchmarks**, separately for UKF and
SR-UKF. A header-only template library has no single universal code size:
the state, measurement types, model, compiler and calling context matter.

Measured results and compiler versions are in the
[STM32H735 README](../stm32h735_benchmark/README.md#estimator-code-size) and
[ESP32-S31 README](../esp32s31_benchmark/README.md#estimator-code-size).
Both size measurements include the core optimizations from #92. The earlier
ESP32 emulator timing capture uses the older core and remains separate.

## Measurement boundary

`probe.cpp` exposes three C functions: construct/initialize in caller-owned
storage, update from runtime sensor inputs, and return the attitude quaternion.
The magnetometer is enabled. It uses the actual `Attitude_Estimator` template,
including its attitude and heading filters, rest/magnetic checks and Eigen
operations. There is no synthetic trajectory for the compiler to fold away.

Each target/optimization builds an empty adapter, UKF and SR-UKF independently.
The linker explicitly retains all three entry points and garbage-collects
unreachable sections. No timing firmware, input generator, recording, sorting,
logging, HAL, RTOS or MCU startup code is included.

The reports distinguish:

- **Code:** retained executable sections from `probe.o`, including inlined
  estimator/Eigen code and the three thin API wrappers. Embedded constants
  such as ARM literal pools also count here when stored in executable sections.
- **Read-only data:** retained constants from the same object. These bytes are
  separate from code, but also occupy flash. Neither column includes external
  math/compiler/C/C++ runtime archive code.
- **Added linked flash:** all allocated, file-backed section payloads in the
  final estimator ELF, minus the empty adapter ELF. This also includes needed
  runtime helpers, unwind/initialization tables and initialized-data load
  bytes. Section-internal padding is included; gaps between sections, ELF
  headers, debug information, bootloaders, image headers and partition padding
  are excluded. Empty baselines are 32 bytes for ARM and 44 bytes for ESP32,
  including 8 bytes of fixed metadata.

The runtime contribution includes allocation/error support reachable from
generic Eigen paths. ESP32 also links C++ allocation/exception support from
the shipped Picolibc-compatible `libstdc++`, even though the probe is compiled
with exceptions disabled. This is a **standalone link cost**, not the marginal
addition to an existing application that may already share helpers or use
ESP-IDF ROM bindings. Do not add UKF and SR-UKF totals to estimate an application
containing both: it can share instantiated functions and runtime routines.

The linker script uses artificial addresses and reserves no stack or heap.
Its zero-sized heap symbols satisfy runtime link requirements; **these ELFs
are not bootable or runnable firmware**. Stock runtime stubs can produce link
warnings, saved in each build's `build.log`. No unresolved symbols are allowed
or suppressed. This does not alter or flash either timing application.

Both estimator instances require 3,116 bytes of caller-owned storage, aligned
to 4 bytes, in the recorded builds. That is not peak RAM: stack, runtime data,
TLS and any allocation paths are separate. JSON `writable_allocated_payload_bytes`
also includes writable-flagged initialization arrays in flash and is not a
RAM figure; `zero_fill_bytes_including_tls` includes any `.tbss`.

## Optimization settings

The two configurations differ only in **`-O2` versus `-Os`**. Both use:

```text
-std=c++17 -DNDEBUG -DUKF_SINGLE_PRECISION -DEIGEN_NO_MALLOC
-ffunction-sections -fdata-sections -fno-lto
-fno-exceptions -fno-rtti -fno-threadsafe-statics -fno-use-cxa-atexit
link: -nostdlib --gc-sections --no-undefined
```

There is no fast-math flag. Assertions are disabled (`-DNDEBUG`); therefore
`EIGEN_NO_MALLOC` assertions do not prove absence of allocation. The probes use
the repository-patched Eigen 3.2.8. Compiler, source, Eigen-header and selected
runtime archive hashes, full commands, allocated sections and ELF hashes are
saved in each JSON report. Library paths come from the actual linker map,
because Picolibc specs redirect searches that `-print-file-name` may not reflect.

STM32H735 uses Arm GNU Toolchain 14.3.Rel1 / GCC 14.3.1 and:

```text
-mcpu=cortex-m7 -mthumb -mfpu=fpv5-d16 -mfloat-abi=hard
libraries: -lm -lc -lgcc -lnosys
```

ESP32-S31 uses Espressif GCC 15.2.0 (`esp-15.2.0_20251204`) and:

```text
-march=rv32imafc_zicsr_zifencei_zaamo_zalrsc_xesploop_xespv
-mabi=ilp32f -mtune=esp-base --specs=picolibc.specs
-D__STDC_WANT_LIB_EXT1__=0 -fstrict-volatile-bitfields
-fno-jump-tables -fno-tree-switch-conversion
libraries: -lstdc++ -lm -lc -lgcc
```

These ESP32 architecture/runtime settings come from the pinned
[ESP-IDF target rules](https://github.com/espressif/esp-idf/blob/fff9895c82d744c7237be8847347bdd1b07c6643/components/soc/project_include.cmake),
[S31 capabilities](https://github.com/espressif/esp-idf/blob/fff9895c82d744c7237be8847347bdd1b07c6643/components/soc/esp32s31/include/soc/Kconfig.soc_caps.in),
[compiler options](https://github.com/espressif/esp-idf/blob/fff9895c82d744c7237be8847347bdd1b07c6643/CMakeLists.txt),
and [Picolibc defaults](https://github.com/espressif/esp-idf/blob/fff9895c82d744c7237be8847347bdd1b07c6643/components/esp_libc/Kconfig).
The probe deliberately normalizes both targets to C++17 with assertions
disabled. It does not inherit ESP-IDF's sdkconfig: that IDF revision defaults
to GNU C++26 and enabled assertions, unlike these size probes. Consequently
these results are not exact symbol attribution inside the ESP timing ELF.
The standalone API boundary can also change inlining relative to either
combined timing firmware. **No `-Os` timings have been measured.**

## Reproduce

From the repository root, with patched Eigen already prepared by the normal
build and toolchains installed under the sibling `SFWA/toolchains` directory:

```sh
python examples/embedded_code_size/measure.py --target stm32h735 --optimization O2
python examples/embedded_code_size/measure.py --target stm32h735 --optimization Os
python examples/embedded_code_size/measure.py --target esp32s31 --optimization O2
python examples/embedded_code_size/measure.py --target esp32s31 --optimization Os
```

Use `--compiler /path/to/TARGET-g++` on other installations, `--eigen` to select
the patched headers, or `--output` to change the build directory. Default
outputs are `build-code-size/TARGET-OPT/{empty,ukf,srukf}/probe.{o,elf,map}`
and `build-code-size/TARGET-OPT/results.json`. Python 3.9 or newer suffices;
no Python dependencies, ESP-IDF checkout or connected board are required.

The script validates ELF32 little-endian target identity, retention of all
three APIs, estimator state metadata and absence of unresolved symbols.
Tests exercise code/constants/runtime/data separation and exclusion of BSS
and debug information from flash totals:

```sh
python -m unittest discover -s examples/embedded_code_size -p 'test_*.py' -v
```
