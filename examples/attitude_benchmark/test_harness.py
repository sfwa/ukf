"""Regression checks for recording-group splits and complete snapshot provenance."""
from copy import deepcopy
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest

from datasets import manifest, SASSARI_SENSORS
import run


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("attitude_snapshot", ROOT / "docs/benchmarks/update.py")
snapshot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(snapshot)


def recording_names():
    trials = []
    for entry in manifest()["files"]:
        names = ([entry["name"]] if entry["dataset"] == "broad" else
                 [f"{entry['name']}_{sensor}" for sensor in SASSARI_SENSORS])
        trials.extend(SimpleNamespace(dataset=entry["dataset"], name=name) for name in names)
    return trials


class SelectionTests(unittest.TestCase):
    def test_development_and_validation_are_disjoint_complete_motion_groups(self):
        all_trials = recording_names()
        development = list(run.select_trials(all_trials, "development"))
        validation = list(run.select_trials(all_trials, "validation"))
        ids = lambda trials: {(trial.dataset, trial.name) for trial in trials}
        self.assertEqual(len(development), 11)
        self.assertEqual(len(validation), 46)
        self.assertFalse(ids(development) & ids(validation))
        self.assertEqual(ids(development) | ids(validation), ids(all_trials))
        self.assertEqual({trial.name for trial in development if trial.dataset == "sassari"},
                         {f"slow_{sensor}" for sensor in SASSARI_SENSORS})
        for speed in ("slow", "medium", "fast"):
            self.assertEqual(len({run.trial_split("sassari", f"{speed}_{s}") for s in SASSARI_SENSORS}), 1)

    def test_name_and_qualified_patterns_select_whole_recordings(self):
        selected = list(run.select_trials(recording_names(), patterns=("01_*", "sassari/slow_*")))
        self.assertEqual(len(selected), 7)
        self.assertEqual(sum(trial.dataset == "sassari" for trial in selected), 6)

    def test_unmatched_pattern_or_pattern_outside_split_fails(self):
        for split, pattern in (("all", "typo*"), ("validation", "01_*")):
            with self.subTest(split=split, pattern=pattern):
                with self.assertRaisesRegex(ValueError, "matched no recording"):
                    list(run.select_trials(recording_names(), split, (pattern,)))


class SnapshotTests(unittest.TestCase):
    def reports(self):
        reports = []
        for dataset, names in snapshot.expected_trials().items():
            report = dict(datasets=[dataset], scope="full selected datasets", split="all",
                          trial_patterns=[], axes=[6, 9], algorithms=list(run.ALGORITHMS))
            report["per_trial"] = [dict(dataset=dataset, trial=name, axes=axes, algorithm=algorithm)
                                   for name in names for axes in (6, 9) for algorithm in run.ALGORITHMS]
            for key in ("source_sha256", "binary_sha256", "packages", "reference_parameters",
                        "burn_in_seconds", "heading_alignment_seconds", "split_protocol", "cpp_configuration", "cpp_parameters"):
                report[key] = {"fixture": "same"}
            reports.append(report)
        return reports

    def test_complete_staged_reports_are_accepted(self):
        self.assertEqual(set(run.ALGORITHMS), snapshot.ALGORITHMS)
        reports = self.reports()
        validated = snapshot.validate_reports(reports)
        self.assertEqual(sum(len(r["per_trial"]) for r in validated.values()), 912)

    def test_repeated_row_cannot_hide_a_missing_recording(self):
        reports = self.reports()
        reports[0]["per_trial"][-1] = deepcopy(reports[0]["per_trial"][0])
        with self.assertRaisesRegex(ValueError, "Incomplete or duplicate"):
            snapshot.validate_reports(reports)

    def test_same_count_but_wrong_recording_names_fails(self):
        reports = self.reports()
        old_name = reports[0]["per_trial"][0]["trial"]
        for row in reports[0]["per_trial"]:
            if row["trial"] == old_name:
                row["trial"] = "not_in_the_pinned_dataset"
        with self.assertRaisesRegex(ValueError, "Incomplete or duplicate"):
            snapshot.validate_reports(reports)

    def test_partial_splits_and_different_builds_cannot_replace_snapshot(self):
        reports = self.reports()
        reports[0]["split"] = "development"
        with self.assertRaisesRegex(ValueError, "Only full"):
            snapshot.validate_reports(reports)
        reports = self.reports()
        reports[0]["binary_sha256"] = "different"
        with self.assertRaisesRegex(ValueError, "different binary_sha256"):
            snapshot.validate_reports(reports)

    def test_different_executable_parameter_descriptions_fail(self):
        reports = self.reports()
        reports[0]["cpp_parameters"] = {"fixture": "different settings"}
        with self.assertRaisesRegex(ValueError, "different cpp_parameters"):
            snapshot.validate_reports(reports)


if __name__ == "__main__":
    unittest.main()
