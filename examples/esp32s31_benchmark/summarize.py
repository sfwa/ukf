#!/usr/bin/env python3
"""Summarize one complete ESP32-S31, STM32H735, or host benchmark log.

Exit 0: complete valid replay (and within the compute budget on hardware).
Exit 1: numerical replay failure or hardware compute-budget failure.
With --emulated, timings are SIMULATED and never establish silicon performance.
Exit 2: malformed, incomplete, or mixed benchmark suites.
"""

import argparse
import json
import math
from pathlib import Path
import re
import sys


RECORD_TYPES = {"benchmark_input", "attitude_benchmark", "platform", "platform_complete", "benchmark_error"}
MARKER = re.compile(r'"type"\s*:\s*"(' + "|".join(sorted(RECORD_TYPES)) + r')"')
ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
EXPECTED = {("ukf", 1), ("srukf", 1), ("ukf", 2), ("srukf", 2)}
FIRMWARE_PLATFORMS = {"esp32s31": "ESP32-S31", "stm32h735": "STM32H735"}


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def number(value, name):
    if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be a finite non-negative number")
    return value


def boolean(value, name):
    if type(value) is not bool:
        raise ValueError(f"{name} must be a boolean")
    return value


def reject_constant(value):
    raise ValueError(f"Non-finite JSON constant: {value}")


def parse_log(text):
    records = {kind: [] for kind in RECORD_TYPES}
    for line_number, line in enumerate(text.splitlines(), 1):
        line = ANSI.sub("", line)
        # Boot/monitor messages and unrelated JSON are not benchmark data.
        if not MARKER.search(line):
            continue
        start = line.find("{")
        try:
            record = json.loads(line[start:], parse_constant=reject_constant)
        except (ValueError, json.JSONDecodeError) as error:
            raise ValueError(f"Malformed benchmark JSON at line {line_number}: {error}") from None
        if not isinstance(record, dict) or record.get("type") not in RECORD_TYPES:
            raise ValueError(f"Invalid benchmark record at line {line_number}")
        records[record["type"]].append(record)

    if records["benchmark_error"]:
        raise ValueError("Benchmark reported an input/setup error: " +
                         str(records["benchmark_error"][0].get("reason", "unknown")))
    if len(records["benchmark_input"]) != 1:
        raise ValueError("Expected exactly one benchmark_input record; capture one complete suite per log")
    source = records["benchmark_input"][0]
    platform = source.get("platform")
    if platform != "host" and platform not in FIRMWARE_PLATFORMS:
        raise ValueError("Unsupported benchmark platform")
    if len(records["platform"]) > 1 or len(records["platform_complete"]) > 1:
        raise ValueError("Multiple platform records; capture one complete suite per log")
    for record in records["platform"] + records["attitude_benchmark"]:
        if record.get("platform") != platform:
            raise ValueError("Mixed benchmark platforms")
    if platform in FIRMWARE_PLATFORMS and (len(records["platform"]) != 1 or
                                           len(records["platform_complete"]) != 1):
        raise ValueError(f"{FIRMWARE_PLATFORMS[platform]} log requires platform and platform_complete records")

    samples = integer(source.get("samples"), "samples", 2)
    if source.get("passes") != 2:
        raise ValueError("Expected two benchmark passes")
    if not isinstance(source.get("source"), str):
        raise ValueError("Input source name is missing")
    for key in ("dt_min_s", "dt_max_s", "dt_mean_s"):
        if number(source.get(key), key) <= 0:
            raise ValueError("Input intervals must be positive")
    if not source["dt_min_s"] <= source["dt_mean_s"] <= source["dt_max_s"]:
        raise ValueError("Inconsistent input timing metadata")

    runs = records["attitude_benchmark"]
    keys = [(run.get("filter"), run.get("pass")) for run in runs]
    if len(keys) != 4 or set(keys) != EXPECTED:
        raise ValueError("Expected exactly UKF and SR-UKF passes 1 and 2; missing or duplicate results")
    for run in runs:
        if (run.get("precision") != "float32" or run.get("axes") != 9 or
                run.get("timed_scope") != "update_and_attitude" or run.get("budget_us") != 2000 or
                run.get("expected_updates") != samples - 1):
            raise ValueError("Run metadata does not match the float32, nine-axis, 2000 us benchmark input")
        boolean(run.get("replay_valid"), "replay_valid")
        boolean(run.get("observed_compute_within_budget"), "observed_compute_within_budget")
        timing = run.get("timing")
        if not isinstance(timing, dict):
            raise ValueError("Missing timing summary")
        count = integer(timing.get("count"), "timing.count")
        over = integer(timing.get("over_budget"), "timing.over_budget")
        mean = number(timing.get("mean_us"), "timing.mean_us")
        maximum = number(timing.get("max_us"), "timing.max_us")
        p95 = number(run.get("p95_us"), "p95_us")
        p99 = number(run.get("p99_us"), "p99_us")
        if count > samples - 1 or over > count or mean > maximum or not p95 <= p99 <= maximum:
            raise ValueError("Inconsistent timing summary")
    complete = records["platform_complete"]
    if complete:
        boolean(complete[0].get("suite_ok"), "platform_complete.suite_ok")
    return source, runs, complete


def markdown(value):
    return str(value).replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def summarize(text, *, emulated=False):
    source, runs, complete = parse_log(text)
    if emulated and source["platform"] != "esp32s31":
        raise ValueError("--emulated requires an ESP32-S31 firmware log")
    hardware = source["platform"] in FIRMWARE_PLATFORMS and not emulated
    lines = ["**SIMULATED ESP32-S31 firmware timings — NOT physical silicon measurements.**" if emulated else
             f"{FIRMWARE_PLATFORMS[source['platform']]} on-device **compute-only** measurements." if hardware else
             "**HOST functional check — NOT ESP32-S31 target performance.**", "",
             f"Input: {markdown(source['source'])}; {source['samples']} samples, "
             f"mean input interval {source['dt_mean_s'] * 1000:.6g} ms. "
             "The comparison budget is 2 ms per update; input timestamps are not retimed.", "",
             "| Filter | Pass | Mean µs | P95 µs | P99 µs | Max µs | Max / 2 ms | Updates > 2 ms | Replay valid |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    all_valid = True
    all_within_budget = True
    for run in sorted(runs, key=lambda row: (row["pass"], 0 if row["filter"] == "ukf" else 1)):
        timing = run["timing"]
        valid = run["replay_valid"] and timing["count"] == run["expected_updates"]
        within_budget = (run["observed_compute_within_budget"] and timing["over_budget"] == 0 and
                         timing["max_us"] <= 2000)
        all_valid = all_valid and valid
        all_within_budget = all_within_budget and within_budget
        label = "UKF" if run["filter"] == "ukf" else "SR-UKF"
        lines.append(f"| {label} | {run['pass']} | {timing['mean_us']:.3f} | {run['p95_us']:.3f} | "
                     f"{run['p99_us']:.3f} | {timing['max_us']:.3f} | {timing['max_us'] / 20:.2f}% | "
                     f"{timing['over_budget']} | {'yes' if valid else 'NO'} |")
    if complete:
        all_valid = all_valid and complete[0]["suite_ok"]
    lines += ["", "Each row is one replay; percentiles are not combined across passes. "
              "Sensor I/O, output logging, sample loading and task-release scheduling are outside the timed region."]
    if not all_valid:
        lines += ["", "**FAIL: numerical replay failed or an update sequence did not complete.**"]
    elif hardware and not all_within_budget:
        lines += ["", "**FAIL: at least one replay did not meet the observed 2 ms compute budget.**"]
    elif hardware:
        lines += ["", "All four replays completed within the observed compute budget. "
                  "This does not establish a worst-case end-to-end deadline guarantee."]
    elif emulated:
        lines += ["", "Functional replay passed. All timings and overrun counts above are SIMULATED. "
                  "Emulation cannot establish whether physical ESP32-S31 silicon meets 2 ms; "
                  "simulated budget results are not used as a pass/fail gate."]
    else:
        lines += ["", "Functional replay passed. Host timing does not establish ESP32-S31 performance; "
                  "the host budget result is deliberately not used as a pass/fail gate."]
    return "\n".join(lines) + "\n", 0 if all_valid and (not hardware or all_within_budget) else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path, help="Serial log containing exactly one complete suite")
    parser.add_argument("--emulated", action="store_true",
                        help="Label ESP32-S31 firmware timings as SIMULATED; gate only functional success")
    args = parser.parse_args(argv)
    try:
        report, status = summarize(args.log.read_text(encoding="utf-8", errors="replace"), emulated=args.emulated)
    except (OSError, ValueError, TypeError) as error:
        print(f"summarize.py: {error}", file=sys.stderr)
        return 2
    print(report, end="")
    return status


if __name__ == "__main__":
    sys.exit(main())
