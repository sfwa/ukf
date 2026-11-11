#!/usr/bin/env python3
"""Replay pinned datasets and produce a reproducible accuracy comparison."""
import argparse
import csv
from fnmatch import fnmatchcase
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time as clock

import numpy as np

from datasets import load_trials
from metrics import compute_errors
from references import run_reference, PARAMETERS

ALGORITHMS = ("ukf", "srukf", "madgwick", "mahony", "vqf")
CPP_ALGORITHMS = frozenset(("ukf", "srukf"))
DEVELOPMENT_BROAD = frozenset(("01", "08", "10", "24", "28"))
DEVELOPMENT_SASSARI = frozenset(("slow",))
SPLIT_PROTOCOL = {
    "name": "attitude-development-v1",
    "development_broad": sorted(DEVELOPMENT_BROAD),
    "development_sassari": sorted(DEVELOPMENT_SASSARI),
    "grouping": "All Sassari sensors from one physical motion stay in the same split.",
    "qualification": "The estimator is tuned on the development recordings only; reference results were already known on all recordings.",
}
ROOT = Path(__file__).resolve().parents[2]


def source_hashes():
    paths = list(Path(__file__).parent.iterdir()) + list((ROOT / "include").rglob("*.h"))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths) if p.suffix in (".py", ".cpp", ".h", ".json", ".txt")
            and not p.name.startswith("test_")}


def describe_cpp_algorithms(executable, algorithms):
    """Capture the executable's fixed parameters once, before any replay."""
    descriptions = {}
    for algorithm in sorted(CPP_ALGORITHMS.intersection(algorithms)):
        result = subprocess.run([str(executable), "--describe-filter", algorithm],
                                check=True, capture_output=True, text=True, timeout=30)
        description = json.loads(result.stdout)
        if not isinstance(description, dict) or description.get("stage") != algorithm:
            raise ValueError(f"Invalid parameter description for {algorithm}")
        descriptions[algorithm] = description
    return descriptions


def trial_split(dataset, name):
    if dataset == "broad":
        development = name.split("_", 1)[0] in DEVELOPMENT_BROAD
    elif dataset == "sassari":
        development = name.split("_", 1)[0] in DEVELOPMENT_SASSARI
    else:
        raise ValueError(f"Unknown dataset {dataset}")
    return "development" if development else "validation"


def trial_condition(dataset, name):
    if dataset == "broad":
        return "undisturbed" if int(name.split("_", 1)[0]) <= 23 else "disturbed"
    return name.split("_", 1)[0]


def select_trials(trials, split="all", patterns=()):
    """Select whole recording names; explicit patterns are marked custom in reports."""
    if split not in ("all", "development", "validation"):
        raise ValueError(f"Unknown split {split}")
    matched = set()
    for trial in trials:
        if split != "all" and trial_split(trial.dataset, trial.name) != split:
            continue
        names = (trial.name, f"{trial.dataset}/{trial.name}")
        matching = {pattern for pattern in patterns if any(fnmatchcase(name, pattern) for name in names)}
        if patterns and not matching:
            continue
        matched.update(matching)
        yield trial
    if set(patterns) - matched:
        raise ValueError(f"Trial patterns matched no recording in the selected split: {sorted(set(patterns)-matched)}")


def cpp_estimate(executable, trial, algorithm, axes, directory, diagnostics=None):
    source, dest = directory / "input.csv", directory / "output.csv"
    values = np.column_stack((trial.time, trial.gyr, trial.acc, trial.mag))
    np.savetxt(source, values, delimiter=",", comments="",
               header="time,gx,gy,gz,ax,ay,az,mx,my,mz", fmt="%.17g")
    command = [str(executable), "--input", str(source), "--output", str(dest),
               "--filter", algorithm, "--axes", str(axes)]
    if diagnostics is not None:
        diagnostics.parent.mkdir(parents=True, exist_ok=True)
        command.extend(("--diagnostics", str(diagnostics)))
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=600)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"{algorithm} failed on {trial.dataset}/{trial.name}: {error.stderr.strip()}") from error
    if diagnostics is not None and not diagnostics.is_file():
        raise ValueError(f"{algorithm} did not write requested diagnostics")
    output = np.loadtxt(dest, delimiter=",", skiprows=1, ndmin=2)
    if output.shape != (len(trial.time), 5) or not np.allclose(
            output[:, 0], trial.time, rtol=0, atol=1e-8):
        raise ValueError("C++ output does not match input timestamps")
    return output[:, 1:], json.loads(result.stderr)["elapsed_seconds"]


def aggregate(rows):
    """Equal weight per recording/sensor, not per sample or dataset duration."""
    output = []
    for dataset in sorted({r["dataset"] for r in rows}):
        for axes in sorted({r["axes"] for r in rows}):
            for algorithm in ALGORITHMS:
                subset = [r for r in rows if (r["dataset"], r["axes"], r["algorithm"])
                          == (dataset, axes, algorithm)]
                if not subset:
                    continue
                item = dict(dataset=dataset, axes=axes, algorithm=algorithm, trials=len(subset))
                for key in ("inclination_rmse_deg", "heading_rmse_deg", "total_rmse_deg"):
                    item[key] = float(np.mean([r[key] for r in subset]))
                output.append(item)
    return output


def grouped_summary(rows, field):
    return {group: aggregate([row for row in rows if row[field] == group])
            for group in sorted({row[field] for row in rows})}


def relative_comparisons(summary):
    """Accuracy ratios against VQF only; no speed comparison or implied acceptance threshold."""
    output = []
    for row in summary:
        match = next((r for r in summary if r["dataset"] == row["dataset"] and
                      r["axes"] == row["axes"] and r["algorithm"] == "vqf"), None)
        if row["algorithm"] not in CPP_ALGORITHMS or match is None:
            continue
        result = {key: row[key] for key in ("dataset", "axes", "algorithm")}
        result["reference"] = "vqf"
        for metric in ("inclination", "heading", "total"):
            denominator = match[f"{metric}_rmse_deg"]
            result[f"{metric}_rmse_ratio"] = row[f"{metric}_rmse_deg"] / denominator if denominator else None
        output.append(result)
    return output


def markdown(report):
    lines = ["# Attitude benchmark comparison", "",
             f"Scope: **{report['scope']}**. Commit: `{report['revision']}`" +
             (" (working tree modified)." if report["dirty"] else "."), "",
             "Mean per-trial RMSE in degrees (each sensor recording has equal weight).",
             "Heading and total error use one fixed yaw alignment in the first second with optical truth; inclination does not.",
             "Six-axis heading measures drift from that alignment, not absolute heading accuracy.", "",
             "| Dataset | Axes | Algorithm | Trials | Inclination | Heading | Total |",
             "|---|---:|---|---:|---:|---:|---:|"]
    for r in report["summary"]:
        lines.append(f"| {r['dataset']} | {r['axes']} | {r['algorithm']} | {r['trials']} | "
                     f"{r['inclination_rmse_deg']:.3f} | {r['heading_rmse_deg']:.3f} | "
                     f"{r['total_rmse_deg']:.3f} |")
    lines += ["", f"Selection split: **{report['split']}**. Explicit selectors are custom subsets. "
              "The UKF estimator is tuned on the development recordings only; reference results were already known on all recordings.",
              "", f"Scoring starts after both {report['burn_in_seconds']:g} s and the optical heading-alignment window. "
              "Only author-marked movement intervals with valid optical references are scored; "
              "filters still process the entire recording, including rest. "
              "Input sensors are never replaced by ground truth.",
              "", "Timing in results.json is diagnostic: C++ update-loop and Python reference-call "
              "times include different overheads and are not a cross-language speed ranking. "
              "Rows marked cpp_update_loop_with_diagnostics_io include in-loop diagnostic formatting and writes.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=ROOT / "build/examples/attitude_benchmark/attitude_benchmark")
    parser.add_argument("--cache", type=Path, default=ROOT / ".benchmark-cache")
    parser.add_argument("--output", type=Path, default=ROOT / "benchmark-results")
    parser.add_argument("--datasets", nargs="+", choices=["broad", "sassari"],
                        default=["broad", "sassari"])
    parser.add_argument("--algorithms", nargs="+", choices=ALGORITHMS, default=list(ALGORITHMS))
    parser.add_argument("--axes", nargs="+", choices=[6, 9], type=int, default=[6, 9])
    parser.add_argument("--smoke", action="store_true", help="One recording from each dataset")
    parser.add_argument("--split", choices=("all", "development", "validation"), default="all",
                        help="Recording-group split for staged estimator development")
    parser.add_argument("--trials", nargs="+", default=[],
                        help="Quoted name globs, e.g. '01_*' or 'sassari/slow_*'; marks the report as custom")
    parser.add_argument("--diagnostics", action="store_true",
                        help="Write UKF estimator state traces below OUTPUT/diagnostics")
    parser.add_argument("--burn-in", type=float, default=5.0)
    args = parser.parse_args()
    if args.smoke and (args.split != "all" or args.trials):
        parser.error("--smoke cannot be combined with --split or --trials")
    for name in ("datasets", "algorithms", "axes"):
        if len(set(getattr(args, name))) != len(getattr(args, name)):
            parser.error(f"--{name} must not contain duplicates")
    if not np.isfinite(args.burn_in) or args.burn_in < 1:
        parser.error("--burn-in must be finite and >= 1 s (the heading alignment window)")
    executable = args.executable.resolve()
    if any(a in CPP_ALGORITHMS for a in args.algorithms) and not executable.is_file():
        parser.error(f"Build the attitude_benchmark target first: {executable}")
    args.output.mkdir(parents=True, exist_ok=True)
    sources = source_hashes()
    binary_hash = hashlib.sha256(executable.read_bytes()).hexdigest() if executable.exists() else None
    cpp_parameters = describe_cpp_algorithms(executable, args.algorithms)
    rows = []
    selected = []
    with tempfile.TemporaryDirectory(prefix="attitude-benchmark-") as temporary:
        trials = select_trials(load_trials(args.cache, args.datasets, smoke=args.smoke), args.split, args.trials)
        for trial in trials:
            selected.append(f"{trial.dataset}/{trial.name}")
            for axes in args.axes:
                for algorithm in args.algorithms:
                    print(f"{trial.dataset}/{trial.name} {axes}-axis {algorithm}", flush=True)
                    start = clock.perf_counter()
                    if algorithm in CPP_ALGORITHMS:
                        diagnostics = (args.output / "diagnostics" / trial.dataset /
                                       f"{trial.name}-{axes}-{algorithm}.csv") if args.diagnostics else None
                        quat, elapsed = cpp_estimate(executable, trial, algorithm, axes, Path(temporary), diagnostics)
                        timing_scope = "cpp_update_loop_with_diagnostics_io" if diagnostics is not None else "cpp_update_loop"
                    else:
                        quat = run_reference(algorithm, trial.time, trial.gyr, trial.acc, trial.mag, axes)
                        elapsed = clock.perf_counter() - start
                        timing_scope = "python_reference_call"
                    # A failed filter is a failed benchmark, never silently dropped from the average.
                    if quat.shape != (len(trial.time), 4) or not np.isfinite(quat).all():
                        raise ValueError(f"Invalid output from {algorithm} on {trial.name}")
                    if not np.allclose(np.linalg.norm(quat, axis=1), 1, atol=1e-5, rtol=0):
                        raise ValueError(f"Non-unit quaternion from {algorithm} on {trial.name}")
                    summary, _ = compute_errors(quat, trial.quat_gt, trial.valid, trial.time,
                                                burn_in=args.burn_in, score_mask=trial.movement)
                    rows.append(dict(dataset=trial.dataset, trial=trial.name, axes=axes,
                                     algorithm=algorithm, split=trial_split(trial.dataset, trial.name),
                                     condition=trial_condition(trial.dataset, trial.name),
                                     **summary, elapsed_seconds=elapsed,
                                     timing_scope=timing_scope))
    if not rows:
        raise ValueError("No recordings evaluated")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip())
    if sources != source_hashes():
        raise RuntimeError("Benchmark sources changed during this run; rerun for consistent provenance")
    if executable.exists() and hashlib.sha256(executable.read_bytes()).hexdigest() != binary_hash:
        raise RuntimeError("Estimator binary changed during this run; rerun for consistent provenance")
    scope = ("smoke (not the full benchmark)" if args.smoke else
             "custom recording selection" if args.trials else
             f"{args.split} split" if args.split != "all" else "full selected datasets")
    summary = aggregate(rows)
    report = dict(scope=scope,
                  revision=revision, dirty=dirty, source_sha256=sources,
                  binary_sha256=binary_hash,
                  platform=platform.platform(), python=sys.version, machine=platform.machine(),
                  packages={p: importlib.metadata.version(p) for p in ("numpy", "scipy", "h5py", "ahrs", "vqf")},
                  reference_parameters=PARAMETERS,
                  cpp_parameters=cpp_parameters,
                  cpp_configuration="Fixed tuning in attitude_estimator.h; ukf and srukf run the same model on UKF::Core and UKF::SquareRootCore.",
                  datasets=args.datasets, axes=args.axes, algorithms=args.algorithms,
                  split=args.split, split_protocol=SPLIT_PROTOCOL,
                  trial_patterns=args.trials, selected_trials=selected,
                  diagnostics_enabled=args.diagnostics,
                  burn_in_seconds=args.burn_in, heading_alignment_seconds=1.0,
                  summary=summary, summary_by_split=grouped_summary(rows, "split"),
                  summary_by_condition=grouped_summary(rows, "condition"),
                  comparisons=relative_comparisons(summary), per_trial=rows)
    (args.output / "results.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (args.output / "per_trial.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary_path = args.output / "comparison.md"
    summary_path.write_text(markdown(report))
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as stream:
            stream.write(summary_path.read_text())
    print(f"Wrote {len(rows)} results to {args.output}")


if __name__ == "__main__":
    main()
