#!/usr/bin/env python3
"""Cross-link isolated float32 attitude estimators and report retained bytes.

The generated ELFs are size probes, not bootable MCU applications. No flashing
or execution is performed. Only the Python standard library is required.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGETS = {
    "stm32h735": {
        "machine": 40,
        "compiler": "arm-gnu-toolchain-14.3.rel1-darwin-arm64-arm-none-eabi/bin/arm-none-eabi-g++",
        "flags": ["-mcpu=cortex-m7", "-mthumb", "-mfpu=fpv5-d16", "-mfloat-abi=hard"],
        "compile_flags": [],
        "link_flags": [],
        "libraries": ["-lm", "-lc", "-lgcc", "-lnosys"],
    },
    "esp32s31": {
        "machine": 243,
        "compiler": "riscv32-esp-elf-15.2.0_20251204-aarch64-apple-darwin/bin/riscv32-esp-elf-g++",
        "flags": ["-march=rv32imafc_zicsr_zifencei_zaamo_zalrsc_xesploop_xespv",
                  "-mabi=ilp32f", "-mtune=esp-base"],
        "compile_flags": ["--specs=picolibc.specs", "-D__STDC_WANT_LIB_EXT1__=0",
                          "-fstrict-volatile-bitfields", "-fno-jump-tables", "-fno-tree-switch-conversion"],
        "link_flags": ["--specs=picolibc.specs", "-fno-rtti"],
        "libraries": ["-lstdc++", "-lm", "-lc", "-lgcc"],
    },
}
ROOT_SYMBOLS = ["ukf_size_initialise", "ukf_size_update", "ukf_size_attitude"]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def elf_sections(data):
    """Read section payloads/flags from an ELF32 little-endian executable."""
    if len(data) < 52 or data[:7] != b"\x7fELF\x01\x01\x01":
        raise ValueError("Expected ELF32 little-endian version 1")
    header = struct.unpack_from("<16sHHIIIIIHHHHHH", data)
    _, kind, machine, _, _, _, offset, _, _, _, _, stride, count, names_index = header
    if kind != 2 or stride != 40 or not count or names_index >= count:
        raise ValueError("Unsupported ELF executable section table")
    if offset + count * stride > len(data):
        raise ValueError("Truncated ELF section table")
    entries = [struct.unpack_from("<IIIIIIIIII", data, offset + i * stride) for i in range(count)]
    names = entries[names_index]
    strings = data[names[4]:names[4] + names[5]]
    sections = []
    for name, kind, flags, address, start, size, _, _, _, _ in entries:
        if kind != 8 and start + size > len(data):
            raise ValueError("Truncated ELF section payload")
        end = strings.find(b"\0", name)
        if name >= len(strings) or end < 0:
            raise ValueError("Invalid ELF section name")
        sections.append({"name": strings[name:end].decode("ascii"), "type": kind,
                         "flags": flags, "address": address, "bytes": size,
                         "payload": b"" if kind == 8 else data[start:start + size]})
    return machine, sections


def section_totals(sections):
    allocated = [s for s in sections if s["flags"] & 2]
    def named(name):
        return sum(s["bytes"] for s in allocated if s["name"] == name)
    return {
        "estimator_code_bytes": named(".ukf_text"),
        "estimator_rodata_bytes": named(".ukf_rodata"),
        "estimator_code_and_rodata_bytes": named(".ukf_text") + named(".ukf_rodata"),
        "runtime_code_bytes": named(".runtime_text"),
        "runtime_rodata_bytes": named(".runtime_rodata"),
        "total_flash_payload_bytes": sum(s["bytes"] for s in allocated if s["type"] != 8),
        # SHF_WRITE also covers init/fini arrays in flash, so this is not a
        # RAM-use figure. NOBITS includes thread-local .tbss, where present.
        "writable_allocated_payload_bytes": sum(s["bytes"] for s in allocated if s["flags"] & 1 and s["type"] != 8),
        "zero_fill_bytes_including_tls": sum(s["bytes"] for s in allocated if s["type"] == 8),
    }


def measure(args):
    target = TARGETS[args.target]
    compiler = (args.compiler or ROOT.parent / "toolchains" / target["compiler"]).resolve()
    eigen = args.eigen.resolve()
    output = (args.output or ROOT / "build-code-size" / f"{args.target}-{args.optimization}").resolve()
    if not compiler.is_file() or not (eigen / "Eigen/Core").is_file():
        raise ValueError("Compiler or Eigen headers missing; pass --compiler/--eigen")
    output.mkdir(parents=True, exist_ok=True)
    nm = compiler.with_name(compiler.name.removesuffix("g++") + "nm")
    commands = []
    def normalized(text):
        for path, replacement in ((compiler.parent.parent, "{toolchain}"), (output, "{output}"),
                                  (eigen, "{eigen}"), (ROOT, "{repo}")):
            text = text.replace(str(path), replacement)
        return text
    def run(command, directory=output, record=True):
        command = [str(x) for x in command]
        if record:
            commands.append({"cwd": normalized(str(directory)), "argv": [normalized(x) for x in command]})
        result = subprocess.run(command, cwd=directory, capture_output=True, text=True, check=True)
        if result.stderr:
            with (directory / "build.log").open("a") as log:
                log.write(result.stderr)
        return result.stdout.strip()
    arch = target["flags"]
    common = [*arch, *target["compile_flags"], f"-{args.optimization}", "-std=c++17", "-DNDEBUG", "-DUKF_SINGLE_PRECISION",
              "-DEIGEN_NO_MALLOC", "-ffunction-sections", "-fdata-sections", "-fno-lto",
              "-fno-exceptions", "-fno-rtti", "-fno-threadsafe-statics", "-fno-use-cxa-atexit",
              "-I", str(ROOT / "include"), "-isystem", str(eigen),
              "-include", str(ROOT / "cmake/EigenGccWarnings.h"), "-Wall", "-Wextra"]
    version = run([compiler, "--version"], record=False).splitlines()[0]
    records = []
    library_paths = set()
    for variant, name in enumerate(("empty", "ukf", "srukf")):
        directory = output / name
        directory.mkdir(exist_ok=True)
        (directory / "build.log").write_text("")
        print(f"{args.target} -{args.optimization}: building {name}", flush=True)
        run([compiler, *common, f"-DUKF_SIZE_VARIANT={variant}", "-c", HERE / "probe.cpp", "-o", "probe.o"], directory)
        run([compiler, *arch, *target["link_flags"], "-fno-lto", "-nostdlib", "-T", HERE / "probe.ld", "probe.o",
             "-Wl,--gc-sections", "-Wl,--no-undefined", "-Wl,-Map=probe.map",
             *[f"-Wl,--undefined={symbol}" for symbol in ROOT_SYMBOLS],
             "-Wl,--start-group", *target["libraries"], "-Wl,--end-group", "-o", "probe.elf"], directory)
        # Specs may redirect libc searches: GCC -print-file-name can report a
        # different archive. Record the archives the linker actually loaded.
        for line in (directory / "probe.map").read_text().splitlines():
            match = re.fullmatch(r"LOAD (.+\.a)", line)
            if match:
                library_paths.add(Path(match[1]).resolve())
        if run([nm, "-u", "probe.elf"], directory):
            raise ValueError(f"Unresolved symbols in {name}")
        symbols = run([nm, "--defined-only", "probe.elf"], directory)
        if any(not any(line.endswith(" " + symbol) for line in symbols.splitlines()) for symbol in ROOT_SYMBOLS):
            raise ValueError(f"Missing retained API entry in {name}")
        machine, sections = elf_sections((directory / "probe.elf").read_bytes())
        if machine != target["machine"]:
            raise ValueError("Wrong target architecture in size probe")
        metadata = next(s["payload"] for s in sections if s["name"] == ".probe_metadata")
        state_bytes, alignment = struct.unpack("<II", metadata)
        totals = section_totals(sections)
        if bool(state_bytes) != bool(variant) or not totals["estimator_code_bytes"]:
            raise ValueError("Estimator or baseline was unexpectedly eliminated")
        record = {"variant": name, **totals, "instance_storage_bytes": state_bytes,
                  "instance_alignment_bytes": alignment, "elf_sha256": sha256(directory / "probe.elf"),
                  "object_sha256": sha256(directory / "probe.o"),
                  "allocated_sections": [{k: v for k, v in s.items() if k != "payload"}
                                         for s in sections if s["flags"] & 2]}
        records.append(record)
    baseline = records[0]
    for record in records:
        record["incremental_flash_bytes"] = record["total_flash_payload_bytes"] - baseline["total_flash_payload_bytes"]
        record["estimator_incremental_bytes"] = record["estimator_code_and_rodata_bytes"] - baseline["estimator_code_and_rodata_bytes"]
    sources = [HERE / "probe.cpp", HERE / "probe.ld", Path(__file__), ROOT / "cmake/EigenGccWarnings.h",
               *sorted((ROOT / "include/UKF").glob("*.h")),
               *[ROOT / "examples/attitude_benchmark" / name for name in
                 ("attitude_estimator.h", "imu_sample.h", "sensor_checks.h")]]
    libraries = [{"path": normalized(str(path)), "sha256": sha256(path)} for path in sorted(library_paths)]
    if not libraries:
        raise ValueError("No runtime archive provenance found in linker maps")
    eigen_hash = hashlib.sha256()
    for path in sorted((eigen / "Eigen").rglob("*")):
        if path.is_file():
            eigen_hash.update(str(path.relative_to(eigen)).encode() + b"\0" + path.read_bytes())
    report = {"schema_version": 1, "measured_at_utc": datetime.now(timezone.utc).isoformat(),
              "target": args.target, "compiler": version, "compiler_sha256": sha256(compiler),
              "source_commit": run(["git", "rev-parse", "HEAD"], ROOT, False),
              "source_tree_dirty": bool(run(["git", "status", "--porcelain"], ROOT, False)),
              "source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in sources},
              "eigen": {"version": "3.2.8 with repository patches", "headers_sha256": eigen_hash.hexdigest()},
              "configuration": {"optimization": args.optimization, "precision": "float32", "standard": "c++17",
                                "assertions": False, "lto": False, "fast_math": False,
                                "section_garbage_collection": True, "target_flags": arch},
              "libraries": libraries, "commands": commands, "measurements": records,
              "scope": "Isolated attitude estimator plus Eigen and thin API wrappers; no benchmark, startup, HAL, RTOS, logging or sensor I/O. Incremental flash includes required linked math/runtime helpers. Not a bootable firmware or a timing measurement."}
    (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print("\n| Filter | Code B | Read-only B | Estimator + read-only B | Added linked flash B | State B |")
    print("|---|---:|---:|---:|---:|---:|")
    for record in records[1:]:
        print(f"| {record['variant']} | {record['estimator_code_bytes']} | {record['estimator_rodata_bytes']} | "
              f"{record['estimator_code_and_rodata_bytes']} | {record['incremental_flash_bytes']} | {record['instance_storage_bytes']} |")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", required=True, choices=TARGETS)
    parser.add_argument("--compiler", type=Path, help="Target g++ executable; defaults to the SFWA toolchain installation")
    parser.add_argument("--eigen", type=Path, default=ROOT / "build/src/eigen3", help="Patched Eigen 3.2.8 headers")
    parser.add_argument("--optimization", choices=("O2", "Os"), default="O2")
    parser.add_argument("--output", type=Path, help="Build/report directory; default: build-code-size/TARGET-OPT")
    try:
        measure(parser.parse_args())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n{getattr(error, 'stdout', '') or ''}{getattr(error, 'stderr', '') or ''}")
