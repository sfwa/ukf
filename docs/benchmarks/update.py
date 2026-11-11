#!/usr/bin/env python3
"""Refresh the checked-in comparison from two complete run.py result files.

python docs/benchmarks/update.py benchmark-results/broad/results.json \
    benchmark-results/sassari/results.json
"""
import argparse
import csv
import datetime
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs=2, type=Path)
    parser.add_argument("--build-description", default="Compiler/build details were not supplied.")
    args = parser.parse_args()
    reports = [json.loads(path.read_text()) for path in args.reports]
    by_dataset = {}
    for report in reports:
        if len(report["datasets"]) != 1:
            raise ValueError("Supply one full report per dataset")
        dataset = report["datasets"][0]
        expected = {"broad": 39, "sassari": 18}[dataset]
        algorithms = {"ukf", "srukf", "madgwick", "mahony", "vqf"}
        if "smoke" in report["scope"] or set(report["axes"]) != {6, 9} or set(report["algorithms"]) != algorithms:
            raise ValueError("Only full, five-algorithm, six/nine-axis runs can update the recorded comparison")
        rows = report["per_trial"]
        trials = {r["trial"] for r in rows}
        wanted = {(trial, axes, algorithm) for trial in trials for axes in (6, 9) for algorithm in algorithms}
        actual = {(r["trial"], r["axes"], r["algorithm"]) for r in rows}
        if len(trials) != expected or actual != wanted or len(rows) != len(wanted):
            raise ValueError(f"Incomplete or duplicate {dataset} results")
        by_dataset[dataset] = report
    if set(by_dataset) != {"broad", "sassari"}:
        raise ValueError("Need both BROAD and Sassari")
    for key in ("source_sha256", "binary_sha256", "packages", "reference_parameters",
                "burn_in_seconds", "heading_alignment_seconds"):
        if reports[0][key] != reports[1][key]:
            raise ValueError(f"Reports have different {key}; rerun with identical sources and settings")
    destination = Path(__file__).resolve().parent
    metadata = {name: {k: v for k, v in report.items() if k != "per_trial"}
                for name, report in by_dataset.items()}
    for details in metadata.values():
        details["build_description"] = args.build_description
    (destination / "attitude-provenance.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    rows = [row for dataset in ("broad", "sassari") for row in by_dataset[dataset]["per_trial"]]
    with (destination / "attitude-per-trial.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lines = ["# Recorded attitude comparison", "",
             f"Recorded {datetime.date.today().isoformat()}. Full suite: **570 evaluations**, "
             "39 BROAD and 18 Sassari recordings, both axis configurations, all five algorithms.", "",
             "Values are mean per-recording RMSE in degrees over author-marked movement intervals. "
             "Heading/total use fixed initial yaw alignment; six-axis heading is drift, not absolute heading accuracy.", "",
             "| Dataset | Axes | Algorithm | Inclination | Heading | Total |",
             "|---|---:|---|---:|---:|---:|"]
    for dataset in ("broad", "sassari"):
        for row in by_dataset[dataset]["summary"]:
            lines.append(f"| {dataset} | {row['axes']} | {row['algorithm']} | "
                         f"{row['inclination_rmse_deg']:.3f} | {row['heading_rmse_deg']:.3f} | "
                         f"{row['total_rmse_deg']:.3f} |")
    lines += ["", "These results describe the supplied fixed-noise UKF model and selected reference settings. "
              "They are not an optimized ranking of filter families.", "",
              "The fixed-noise UKFs do not reject disturbed measurements. Large six-axis heading errors "
              "are retained: heading is unobservable and gyro bias/model errors can accumulate during motion.", "",
              "[Methodology and reproduction](../attitude-benchmark.md) · "
              "[Per-trial metrics](attitude-per-trial.csv) · "
              "[Source hashes, dependencies and provenance](attitude-provenance.json)", "",
              args.build_description + " "
              "The provenance revision identifies the checked-out commit; the dirty flag and source hashes identify local changes. "
              "CI regenerates reports for the checked-out commit. Timing is diagnostic and not a speed ranking.", ""]
    (destination / "attitude-comparison.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
