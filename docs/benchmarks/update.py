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

ALGORITHMS = {"ukf", "srukf", "robust-rest", "robust-acc", "robust", "madgwick", "mahony", "vqf"}
ROOT = Path(__file__).resolve().parents[2]


def expected_trials():
    manifest = json.loads((ROOT / "examples/attitude_benchmark/datasets.lock.json").read_text())
    names = {"broad": set(), "sassari": set()}
    for entry in manifest["files"]:
        if entry["dataset"] == "broad":
            names["broad"].add(entry["name"])
        else:
            names["sassari"].update(f"{entry['name']}_{sensor}"
                                   for sensor in ("XS1", "XS2", "AP1", "AP2", "SH1", "SH2"))
    return names


def validate_reports(reports):
    by_dataset = {}
    expected = expected_trials()
    for report in reports:
        if len(report["datasets"]) != 1:
            raise ValueError("Supply one full report per dataset")
        dataset = report["datasets"][0]
        if dataset not in expected or dataset in by_dataset:
            raise ValueError("Need exactly one BROAD and one Sassari report")
        if (report["scope"] != "full selected datasets" or report.get("split") != "all"
                or report.get("trial_patterns") or set(report["axes"]) != {6, 9}
                or set(report["algorithms"]) != ALGORITHMS):
            raise ValueError("Only full, eight-algorithm, six/nine-axis runs can update the recorded comparison")
        rows = report["per_trial"]
        wanted = {(trial, axes, algorithm) for trial in expected[dataset]
                  for axes in (6, 9) for algorithm in ALGORITHMS}
        actual = {(r["trial"], r["axes"], r["algorithm"]) for r in rows}
        if actual != wanted or len(rows) != len(wanted) or any(r["dataset"] != dataset for r in rows):
            raise ValueError(f"Incomplete or duplicate {dataset} results")
        by_dataset[dataset] = report
    if set(by_dataset) != {"broad", "sassari"}:
        raise ValueError("Need both BROAD and Sassari")
    for key in ("source_sha256", "binary_sha256", "packages", "reference_parameters",
                "burn_in_seconds", "heading_alignment_seconds", "split_protocol", "cpp_configuration", "cpp_parameters"):
        if reports[0][key] != reports[1][key]:
            raise ValueError(f"Reports have different {key}; rerun with identical sources and settings")
    return by_dataset


def table(summary):
    lines = ["| Dataset | Axes | Algorithm | Trials | Inclination | Heading | Total |",
             "|---|---:|---|---:|---:|---:|---:|"]
    for row in summary:
        lines.append(f"| {row['dataset']} | {row['axes']} | {row['algorithm']} | {row['trials']} | "
                     f"{row['inclination_rmse_deg']:.3f} | {row['heading_rmse_deg']:.3f} | "
                     f"{row['total_rmse_deg']:.3f} |")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs=2, type=Path)
    parser.add_argument("--build-description", default="Compiler/build details were not supplied.")
    args = parser.parse_args()
    reports = [json.loads(path.read_text()) for path in args.reports]
    by_dataset = validate_reports(reports)
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
             f"Recorded {datetime.date.today().isoformat()}. Full suite: **{len(rows)} evaluations**, "
             f"39 BROAD and 18 Sassari recordings, both axis configurations, all {len(ALGORITHMS)} algorithms/variants.", "",
             "Values are mean per-recording RMSE in degrees over author-marked movement intervals. "
             "Heading/total use fixed initial yaw alignment; six-axis heading is drift, not absolute heading accuracy.", "",
             "## Validation split", "",
             "New-estimator development uses BROAD 01, 08, 10, 24, 28 and Sassari slow (all six sensors). "
             "The table below evaluates the remaining 34 BROAD and 12 Sassari recordings. "
             "These labels govern the new implementation's iteration; baseline results on all recordings were already known.", ""]
    lines += table([row for dataset in ("broad", "sassari")
                    for row in by_dataset[dataset]["summary_by_split"]["validation"]])
    lines += ["", "## Full dataset comparison", ""]
    lines += table([row for dataset in ("broad", "sassari") for row in by_dataset[dataset]["summary"]])
    lines += ["", "## Nine-axis conditions", "",
              "Mean per-recording total RMSE in degrees. The final stage adds both magnetic rejection "
              "and acceleration guards; protection can trade accuracy on clean data for disturbance tolerance.", "",
              "| Dataset | Condition | Recordings | robust-acc | robust | VQF |",
              "|---|---|---:|---:|---:|---:|"]
    for dataset in ("broad", "sassari"):
        for condition, summary in by_dataset[dataset]["summary_by_condition"].items():
            selected = {row["algorithm"]: row for row in summary if row["axes"] == 9}
            lines.append(f"| {dataset} | {condition} | {selected['robust']['trials']} | "
                         f"{selected['robust-acc']['total_rmse_deg']:.3f} | "
                         f"{selected['robust']['total_rmse_deg']:.3f} | "
                         f"{selected['vqf']['total_rmse_deg']:.3f} |")
    lines += ["", "The baseline `ukf`/`srukf` use the original fixed-noise joint attitude/bias model. "
              "`robust-rest`, `robust-acc`, and `robust` are cumulative estimator variants; "
              "the first stage also changes the bias/heading architecture, so its difference is not a rest-only ablation. "
              "These measurements are not an optimized ranking of filter families.", "",
              "[Methodology and reproduction](../attitude-benchmark.md) · "
              "[Per-trial metrics](attitude-per-trial.csv) · "
              "[Source hashes, dependencies and provenance](attitude-provenance.json)", "",
              args.build_description + " "
              "The provenance revision identifies the checked-out commit; the dirty flag and source hashes identify local changes. "
              "CI regenerates reports for the checked-out commit. Timing is diagnostic and not a speed ranking.", ""]
    (destination / "attitude-comparison.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
