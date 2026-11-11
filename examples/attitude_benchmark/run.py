#!/usr/bin/env python3
"""Replay pinned datasets and produce a reproducible accuracy comparison."""
import argparse
import csv
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
ROOT = Path(__file__).resolve().parents[2]


def source_hashes():
    paths = list(Path(__file__).parent.iterdir()) + list((ROOT / "include").rglob("*.h"))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths) if p.suffix in (".py", ".cpp", ".h", ".json", ".txt")
            and not p.name.startswith("test_")}


def cpp_estimate(executable, trial, algorithm, axes, directory):
    source, dest = directory / "input.csv", directory / "output.csv"
    values = np.column_stack((trial.time, trial.gyr, trial.acc, trial.mag))
    np.savetxt(source, values, delimiter=",", comments="",
               header="time,gx,gy,gz,ax,ay,az,mx,my,mz", fmt="%.17g")
    result = subprocess.run([str(executable), "--input", str(source), "--output", str(dest),
                             "--filter", algorithm, "--axes", str(axes)],
                            check=True, capture_output=True, text=True, timeout=600)
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
    lines += ["", f"Scoring starts after {report['burn_in_seconds']:g} s. "
              "Only author-marked movement intervals with valid optical references are scored; "
              "filters still process the entire recording, including rest. "
              "Input sensors are never replaced by ground truth.",
              "", "Timing in results.json is diagnostic: C++ update-loop and Python reference-call "
              "times include different overheads and are not a cross-language speed ranking.", ""]
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
    parser.add_argument("--burn-in", type=float, default=5.0)
    args = parser.parse_args()
    if not np.isfinite(args.burn_in) or args.burn_in < 1:
        parser.error("--burn-in must be finite and >= 1 s (the heading alignment window)")
    executable = args.executable.resolve()
    if any(a in args.algorithms for a in ("ukf", "srukf")) and not executable.is_file():
        parser.error(f"Build the attitude_benchmark target first: {executable}")
    args.output.mkdir(parents=True, exist_ok=True)
    sources = source_hashes()
    rows = []
    with tempfile.TemporaryDirectory(prefix="attitude-benchmark-") as temporary:
        for trial in load_trials(args.cache, args.datasets, smoke=args.smoke):
            for axes in args.axes:
                for algorithm in args.algorithms:
                    print(f"{trial.dataset}/{trial.name} {axes}-axis {algorithm}", flush=True)
                    start = clock.perf_counter()
                    if algorithm in ("ukf", "srukf"):
                        quat, elapsed = cpp_estimate(executable, trial, algorithm, axes, Path(temporary))
                        timing_scope = "cpp_update_loop"
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
                                     algorithm=algorithm, **summary, elapsed_seconds=elapsed,
                                     timing_scope=timing_scope))
    if not rows:
        raise ValueError("No recordings evaluated")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip())
    if sources != source_hashes():
        raise RuntimeError("Benchmark sources changed during this run; rerun for consistent provenance")
    report = dict(scope="smoke (not the full benchmark)" if args.smoke else "full selected datasets",
                  revision=revision, dirty=dirty, source_sha256=sources,
                  binary_sha256=hashlib.sha256(executable.read_bytes()).hexdigest() if executable.exists() else None,
                  platform=platform.platform(), python=sys.version, machine=platform.machine(),
                  packages={p: importlib.metadata.version(p) for p in ("numpy", "scipy", "h5py", "ahrs", "vqf")},
                  reference_parameters=PARAMETERS,
                  datasets=args.datasets, axes=args.axes, algorithms=args.algorithms,
                  burn_in_seconds=args.burn_in, heading_alignment_seconds=1.0,
                  summary=aggregate(rows), per_trial=rows)
    (args.output / "results.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (args.output / "per_trial.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
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
