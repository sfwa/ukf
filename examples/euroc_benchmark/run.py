#!/usr/bin/env python3
"""Causal six-axis EuRoC attitude comparison with fixed recording-group splits."""
from __future__ import annotations

import argparse
import csv
from fnmatch import fnmatchcase
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time as clock

import numpy as np
from vqf import VQF

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ATTITUDE = HERE.parent / "attitude_benchmark"


def load_module(name, path):
    """Use explicit paths: both benchmark directories contain a datasets.py."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


datasets = load_module("euroc_datasets", HERE / "datasets.py")
metrics = load_module("euroc_common_metrics", ATTITUDE / "metrics.py")
references = load_module("euroc_common_references", ATTITUDE / "references.py")

ALGORITHMS = ("ukf", "srukf", "robust-rest", "robust-acc", "robust", "mav-srukf", "mav-srukf-gravity",
              "madgwick", "mahony", "vqf", "vqf-mav")
BASELINE_CPP = frozenset(("ukf", "srukf", "robust-rest", "robust-acc", "robust"))
ROBUST_CPP = frozenset(("robust-rest", "robust-acc", "robust"))
MAV_CPP = frozenset(("mav-srukf", "mav-srukf-gravity"))
CPP_ALGORITHMS = BASELINE_CPP | MAV_CPP
DIAGNOSTIC_ALGORITHMS = ROBUST_CPP | MAV_CPP
DEVELOPMENT = frozenset(("MH_01_easy", "V1_02_medium"))
SPLIT_PROTOCOL = {
    "name": "euroc-mav-development-v1",
    "development": sorted(DEVELOPMENT),
    "validation": "The other nine complete EuRoC sequences.",
    "grouping": "Whole recordings; no sample-level random split.",
    "qualification": "Freeze the MAV estimator before its first validation evaluation.",
}
VQF_INTERVAL = 1 / 200
CLOCK_QUANTIZATION_TOLERANCE = 1e-6
VQF_MAV_PARAMETERS = {"biasClip": 10.0, "biasSigmaInit": 5.0}
PRIMARY_PROTOCOL = {
    "selection": "Vicon recordings except V1_01_easy, whose canonical orientation is known to be inaccurate.",
    "reason": "Machine Hall orientation was reconstructed using IMU data and is not an independent attitude reference.",
    "coverage": "Five Vicon recordings overall; four in the validation split.",
    "primary_metric": "inclination_rmse_deg",
    "metric_reason": "The vehicle model targets inclination; six-axis yaw is unobservable and heading/total can be dominated by bias drift.",
    "reference_limit": "All published state references include offline IMU fusion; Vicon additionally observes optical orientation.",
}


def source_hashes():
    paths = list(HERE.iterdir()) + list(ATTITUDE.iterdir())
    paths += list((ROOT / "include").rglob("*.h"))
    paths += [ROOT / "CMakeLists.txt", ROOT / "examples/CMakeLists.txt"]
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(set(paths)) if path.is_file()
            and path.suffix in (".py", ".cpp", ".h", ".json", ".txt")
            and not path.name.startswith("test_")}


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sequence_split(name):
    return "development" if name in DEVELOPMENT else "validation"


def independent_attitude(name):
    return name.startswith(("V1_", "V2_")) and name != "V1_01_easy"


def select_sequences(names, split="all", patterns=(), smoke=False):
    """Select before loading, so development runs do not fetch validation data."""
    if split not in ("all", "development", "validation"):
        raise ValueError(f"Unknown split {split}")
    if smoke and (split != "all" or patterns):
        raise ValueError("--smoke cannot be combined with --split or --sequences")
    if smoke:
        if "MH_01_easy" not in names:
            raise ValueError("Smoke sequence missing from manifest")
        return ["MH_01_easy"]
    selected, matched = [], set()
    for name in names:
        if split != "all" and sequence_split(name) != split:
            continue
        matching = {pattern for pattern in patterns if fnmatchcase(name, pattern)}
        if patterns and not matching:
            continue
        selected.append(name)
        matched.update(matching)
    if set(patterns) - matched:
        raise ValueError(f"Sequence patterns matched nothing in the selected split: {sorted(set(patterns)-matched)}")
    if not selected:
        raise ValueError("No recordings selected")
    return selected


def describe(executable, algorithm):
    result = subprocess.run([str(executable), "--describe-filter", algorithm],
                            check=True, capture_output=True, text=True, timeout=30)
    description = json.loads(result.stdout)
    if not isinstance(description, dict) or description.get("stage") != algorithm:
        raise ValueError(f"Invalid parameter description for {algorithm}")
    return description


def cpp_estimate(executable, trial, algorithm, directory, diagnostics=None):
    source, destination = directory / "input.csv", directory / "output.csv"
    # Only measured IMU data enters the executable; no reference state, optical
    # velocity, bias, or position is serialized. EuRoC has no published mag data.
    values = np.column_stack((trial.time, trial.gyr, trial.acc, np.zeros_like(trial.acc)))
    np.savetxt(source, values, delimiter=",", comments="", fmt="%.17g",
               header="time,gx,gy,gz,ax,ay,az,mx,my,mz")
    command = [str(executable), "--input", str(source), "--output", str(destination),
               "--filter", algorithm, "--axes", "6"]
    destination.unlink(missing_ok=True)
    if diagnostics is not None:
        diagnostics.parent.mkdir(parents=True, exist_ok=True)
        diagnostics.unlink(missing_ok=True)
        command.extend(("--diagnostics", str(diagnostics)))
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=600)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"{algorithm} failed on {trial.name}: {error.stderr.strip()}") from error
    if diagnostics is not None and not diagnostics.is_file():
        raise ValueError(f"{algorithm} did not write requested diagnostics")
    output = np.loadtxt(destination, delimiter=",", skiprows=1, ndmin=2)
    if output.shape != (len(trial.time), 5) or not np.allclose(
            output[:, 0], trial.time, atol=1e-8, rtol=0):
        raise ValueError("C++ output does not match input timestamps")
    return output[:, 1:], float(json.loads(result.stderr)["elapsed_seconds"])


def vqf_clock(time):
    """Remove published timestamp quantization, not physical sampling jitter.

    EuRoC IMU timestamps are quantized at 256 ns, although the sensor runs at
    200 Hz. VQF requires one fixed update period. Accept only a clock lying
    within 1 microsecond of that exact grid; reject dropped samples or actual
    variable-rate recordings instead of silently resampling sensor values.
    """
    time = np.asarray(time, dtype=float)
    if time.ndim != 1 or len(time) < 2 or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError("Expected finite increasing IMU timestamps")
    nominal = time[0] + np.arange(len(time)) * VQF_INTERVAL
    deviation = float(np.max(np.abs(time-nominal)))
    interval_deviation = float(np.max(np.abs(np.diff(time)-VQF_INTERVAL)))
    if deviation > CLOCK_QUANTIZATION_TOLERANCE or interval_deviation > CLOCK_QUANTIZATION_TOLERANCE:
        raise ValueError("EuRoC VQF adapter only accepts sub-microsecond 200 Hz clock quantization; "
                         "a timing gap or variable-rate recording requires a different protocol")
    return nominal, {
        "policy": "Nominal 200 Hz VQF clock; unchanged IMU samples, no interpolation or resampling.",
        "nominal_interval_seconds": VQF_INTERVAL,
        "max_timestamp_deviation_seconds": deviation,
        "max_interval_deviation_seconds": interval_deviation,
        "allowed_quantization_seconds": CLOCK_QUANTIZATION_TOLERANCE,
    }


def vqf_mav_estimate(time, gyro, acc):
    """Same online VQF adapter with the MAV example's bias range and prior SD."""
    gyro, acc = np.ascontiguousarray(gyro), np.ascontiguousarray(acc)
    filt = VQF(float(time[1]-time[0]), **VQF_MAV_PARAMETERS)
    quat = np.empty((len(time), 4))
    filt.updateAcc(acc[0])
    quat[0] = filt.getQuat6D()
    for index in range(1, len(time)):
        filt.updateGyr((gyro[index-1]+gyro[index])*0.5)
        filt.updateAcc(acc[index])
        quat[index] = filt.getQuat6D()
    return quat


def aggregate(rows):
    result = []
    for algorithm in ALGORITHMS:
        selected = [row for row in rows if row["algorithm"] == algorithm]
        if not selected:
            continue
        item = dict(dataset="euroc", axes=6, algorithm=algorithm, sequences=len(selected))
        for metric in ("inclination", "heading", "total"):
            item[metric + "_rmse_deg"] = float(np.mean([row[metric + "_rmse_deg"] for row in selected]))
        result.append(item)
    return result


def grouped(rows, field):
    return {name: aggregate([row for row in rows if row[field] == name])
            for name in sorted({row[field] for row in rows})}


def paired_comparisons(rows, candidate="mav-srukf", reference="vqf"):
    """Paired recording-level differences; no must-win or significance claim."""
    pairs = {}
    for row in rows:
        if row["algorithm"] in (candidate, reference):
            pair = pairs.setdefault(row["sequence"], {})
            if row["algorithm"] in pair:
                raise ValueError("Duplicate algorithm/recording in paired comparison")
            pair[row["algorithm"]] = row
    present = {algorithm for pair in pairs.values() for algorithm in pair}
    if candidate not in present or reference not in present:
        return None
    complete = [pair for pair in pairs.values() if candidate in pair and reference in pair]
    if len(complete) != len(pairs):
        raise ValueError("Cannot compare algorithms with different recording coverage")
    result = dict(candidate=candidate, reference=reference, sequences=len(complete))
    for metric in ("inclination", "heading", "total"):
        key = metric + "_rmse_deg"
        candidate_values = np.array([pair[candidate][key] for pair in complete])
        reference_values = np.array([pair[reference][key] for pair in complete])
        result[metric] = {
            "mean_difference_deg": float(np.mean(candidate_values-reference_values)),
            "ratio_of_means": float(np.mean(candidate_values)/np.mean(reference_values))
            if np.mean(reference_values) else None,
            "candidate_wins": int(np.count_nonzero(candidate_values < reference_values)),
        }
    return result


def table(summary):
    if not summary:
        return ["No selected recordings in this subset."]
    lines = ["| Algorithm | Sequences | Inclination | Heading | Total |",
             "|---|---:|---:|---:|---:|"]
    for row in summary:
        lines.append(f"| {row['algorithm']} | {row['sequences']} | "
                     f"{row['inclination_rmse_deg']:.3f} | {row['heading_rmse_deg']:.3f} | "
                     f"{row['total_rmse_deg']:.3f} |")
    return lines


def markdown(report):
    lines = ["# EuRoC MAV attitude comparison", "",
             f"Scope: **{report['scope']}**. Revision: `{report['revision']}`" +
             (" (working tree modified)." if report["dirty"] else "."), "",
             "Mean per-sequence RMSE in degrees; each complete recording has equal weight. "
             "The published data provides gyro and accelerometer measurements, so this is six-axis only. "
             "Heading/total use one fixed initial yaw alignment; heading measures subsequent drift.", "",
             "## Primary Vicon-attitude reference subset", "",
             "Vicon sequences except V1_01_easy, whose canonical orientation has a known quality issue. "
             "Machine Hall has no independent optical orientation observation and is reported separately. "
             "All published state references include offline IMU fusion, including Vicon. "
             "The primary vehicle-model target is inclination RMSE; heading and total are retained "
             "but can be dominated by unobservable-yaw drift and gyro-bias handling.", ""]
    lines += table(report["summary_primary"])
    if report["summary_primary_validation"]:
        lines += ["", "### Primary validation sequences", ""]
        lines += table(report["summary_primary_validation"])
    lines += ["", "## All selected sequences (descriptive)", ""]
    lines += table(report["summary"])
    if "validation" in report["summary_by_split"]:
        lines += ["", "## Validation sequences", ""]
        lines += table(report["summary_by_split"]["validation"])
    lines += ["", f"Scoring starts after {report['burn_in_seconds']:g} s and the initial "
              "one-second ground-truth alignment window. All valid reference samples are scored; "
              "there is no optical-truth-derived flight selection or full-trajectory alignment. "
              "Every filter processes the complete IMU sequence, including gaps in reference coverage.", "",
              "Development uses MH_01_easy and V1_02_medium; the other nine sequences are validation. "
              "No reference orientation, position, velocity or bias is supplied to any estimator.", "",
              "The original VQF uses online defaults. `vqf-mav` changes `biasClip` to 10°/s "
              "and `biasSigmaInit` to 5°/s, widening its bias allowance and matching the MAV initial prior SD. "
              "The native 10°/s rest threshold is not a hard bias-state clip. "
              "This controls for a sensor-bias specification advantage over the default VQF. "
              "VQF receives a nominal 200 Hz clock only after the source timestamps are verified "
              "within 1 µs of that grid; all IMU values are unchanged. Other estimators use measured intervals.", "",
              "C++ update-loop and Python wrapper times include different overheads and are not a "
              "cross-language speed ranking. Optional diagnostic writes occur inside the native timed loop.", ""]
    warnings = [(name, metadata.get("ground_truth_warning") or metadata.get("quality_warning"))
                for name, metadata in report["sequence_metadata"].items()]
    if any(warning for _, warning in warnings):
        lines += ["## Reference limitations", ""]
        lines += [f"- {name}: {warning}" for name, warning in warnings if warning]
        lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=ROOT / "build/examples/attitude_benchmark/attitude_benchmark")
    parser.add_argument("--mav-executable", type=Path,
                        default=ROOT / "build/examples/euroc_benchmark/euroc_benchmark")
    parser.add_argument("--cache", type=Path, default=ROOT / ".benchmark-cache/euroc")
    parser.add_argument("--output", type=Path, default=ROOT / "benchmark-results/euroc")
    parser.add_argument("--algorithms", nargs="+", choices=ALGORITHMS, default=list(ALGORITHMS))
    parser.add_argument("--split", choices=("all", "development", "validation"), default="all")
    parser.add_argument("--sequences", nargs="+", default=[], help="Quoted whole-recording name globs")
    parser.add_argument("--smoke", action="store_true", help="Complete MH_01_easy recording only")
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument("--burn-in", type=float, default=5.0)
    args = parser.parse_args()
    if len(set(args.algorithms)) != len(args.algorithms):
        parser.error("--algorithms must not contain duplicates")
    if not np.isfinite(args.burn_in) or args.burn_in < 1:
        parser.error("--burn-in must be finite and >= 1 second")
    try:
        selected = select_sequences(datasets.SEQUENCES, args.split, args.sequences, args.smoke)
    except ValueError as error:
        parser.error(str(error))
    executables = {}
    if BASELINE_CPP.intersection(args.algorithms):
        executables["attitude_benchmark"] = args.executable.resolve()
    if MAV_CPP.intersection(args.algorithms):
        executables["euroc_benchmark"] = args.mav_executable.resolve()
    for name, path in executables.items():
        if not path.is_file():
            parser.error(f"Build {name} first: {path}")
    args.output.mkdir(parents=True, exist_ok=True)
    sources = source_hashes()
    binaries = {name: file_hash(path) for name, path in executables.items()}
    parameters = {name: describe(args.executable.resolve(), name)
                  for name in sorted(ROBUST_CPP.intersection(args.algorithms))}
    for name in sorted(MAV_CPP.intersection(args.algorithms)):
        parameters[name] = describe(args.mav_executable.resolve(), name)
    rows, sequence_metadata = [], {}
    with tempfile.TemporaryDirectory(prefix="euroc-benchmark-") as temporary:
        for trial in datasets.load_trials(args.cache, selected):
            sequence_metadata[trial.name] = dict(trial.metadata)
            reference_clock = trial.time
            if {"vqf", "vqf-mav"}.intersection(args.algorithms):
                reference_clock, timing_metadata = vqf_clock(trial.time)
                sequence_metadata[trial.name]["vqf_clock"] = timing_metadata
            for algorithm in args.algorithms:
                print(f"euroc/{trial.name} 6-axis {algorithm}", flush=True)
                start = clock.perf_counter()
                if algorithm in CPP_ALGORITHMS:
                    executable = args.mav_executable if algorithm in MAV_CPP else args.executable
                    diagnostics = (args.output / "diagnostics" / f"{trial.name}-{algorithm}.csv") \
                        if args.diagnostics and algorithm in DIAGNOSTIC_ALGORITHMS else None
                    quat, elapsed = cpp_estimate(executable.resolve(), trial, algorithm,
                                                 Path(temporary), diagnostics)
                    timing_scope = "cpp_update_loop_with_diagnostics_io" if diagnostics else "cpp_update_loop"
                else:
                    quat = (vqf_mav_estimate(reference_clock, trial.gyr, trial.acc) if algorithm == "vqf-mav" else
                            references.run_reference(algorithm, reference_clock if algorithm == "vqf" else trial.time,
                                                     trial.gyr, trial.acc, np.zeros_like(trial.acc), 6))
                    elapsed = clock.perf_counter()-start
                    timing_scope = "python_reference_call"
                if quat.shape != (len(trial.time), 4) or not np.isfinite(quat).all():
                    raise ValueError(f"Invalid {algorithm} output on {trial.name}")
                if not np.allclose(np.linalg.norm(quat, axis=1), 1, atol=1e-5, rtol=0):
                    raise ValueError(f"Non-unit {algorithm} output on {trial.name}")
                summary, _ = metrics.compute_errors(quat, trial.quat_gt, trial.valid, trial.time,
                                                     burn_in=args.burn_in)
                rows.append(dict(dataset="euroc", sequence=trial.name, axes=6, algorithm=algorithm,
                                 split=sequence_split(trial.name), environment=trial.name.split("_")[0],
                                 difficulty=trial.name.rsplit("_", 1)[1],
                                 primary_attitude_reference=independent_attitude(trial.name), **summary,
                                 elapsed_seconds=elapsed, timing_scope=timing_scope))
    wanted = {(name, algorithm) for name in selected for algorithm in args.algorithms}
    actual = {(row["sequence"], row["algorithm"]) for row in rows}
    if wanted != actual or len(rows) != len(wanted):
        raise ValueError("Incomplete or duplicate benchmark coverage")
    if sources != source_hashes():
        raise RuntimeError("Benchmark sources changed during replay; rerun for consistent provenance")
    if binaries != {name: file_hash(path) for name, path in executables.items()}:
        raise RuntimeError("Estimator binary changed during replay; rerun for consistent provenance")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip())
    scope = ("smoke (not the full benchmark)" if args.smoke else "custom recording selection" if args.sequences
             else f"{args.split} split" if args.split != "all" else "full EuRoC dataset")
    reference_parameters = dict(references.PARAMETERS)
    reference_parameters["vqf-mav"] = {
        "implementation": "vqf 2.0.1", "parameters": VQF_MAV_PARAMETERS,
        "other_parameters": "online defaults", "adapter": "same midpoint gyro and initial accelerometer update as vqf",
    }
    primary = [row for row in rows if row["primary_attitude_reference"]]
    report = dict(scope=scope, revision=revision, dirty=dirty, source_sha256=sources,
                  binary_sha256=binaries, platform=platform.platform(), machine=platform.machine(),
                  python=sys.version,
                  packages={name: importlib.metadata.version(name) for name in ("numpy", "scipy", "h5py", "ahrs", "vqf")},
                  reference_parameters=reference_parameters, cpp_parameters=parameters,
                  algorithms=args.algorithms, axes=[6], selected_sequences=selected,
                  split=args.split, split_protocol=SPLIT_PROTOCOL, sequence_patterns=args.sequences,
                  sequence_metadata=sequence_metadata, burn_in_seconds=args.burn_in,
                  heading_alignment_seconds=1.0, diagnostics_enabled=args.diagnostics,
                  summary=aggregate(rows), summary_by_split=grouped(rows, "split"),
                  primary_attitude_protocol=PRIMARY_PROTOCOL, summary_primary=aggregate(primary),
                  summary_primary_validation=aggregate([row for row in primary if row["split"] == "validation"]),
                  summary_by_environment=grouped(rows, "environment"),
                  summary_by_difficulty=grouped(rows, "difficulty"),
                  paired_comparison_to_vqf=paired_comparisons(rows),
                  paired_comparison_primary={name: paired_comparisons(primary, reference=name)
                                              for name in ("vqf", "vqf-mav")},
                  paired_comparison_primary_validation={name: paired_comparisons(
                      [row for row in primary if row["split"] == "validation"], reference=name)
                      for name in ("vqf", "vqf-mav")},
                  paired_model_effect_primary=paired_comparisons(primary, reference="mav-srukf-gravity"),
                  paired_model_effect_primary_validation=paired_comparisons(
                      [row for row in primary if row["split"] == "validation"], reference="mav-srukf-gravity"),
                  paired_comparison_by_split={name: paired_comparisons([row for row in rows if row["split"] == name])
                                              for name in sorted({row["split"] for row in rows})},
                  per_sequence=rows)
    (args.output / "results.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (args.output / "per_sequence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    comparison = markdown(report)
    (args.output / "comparison.md").write_text(comparison)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as stream:
            stream.write(comparison)
    print(f"Wrote {len(rows)} results to {args.output}")


if __name__ == "__main__":
    main()
