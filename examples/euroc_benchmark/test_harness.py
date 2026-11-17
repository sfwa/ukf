"""Checks of recording selection, paired summaries and sensor-only replay input."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("euroc_runner", HERE / "run.py")
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def row(sequence, algorithm, error, samples=100):
    return dict(sequence=sequence, algorithm=algorithm, samples=samples,
                inclination_rmse_deg=error, heading_rmse_deg=2*error,
                total_rmse_deg=3*error)


class RecordingSelectionTests(unittest.TestCase):
    def test_whole_recording_splits_are_disjoint_and_complete(self):
        names = runner.datasets.SEQUENCES
        development = runner.select_sequences(names, split="development")
        validation = runner.select_sequences(names, split="validation")
        self.assertEqual(set(development), {"MH_01_easy", "V1_02_medium"})
        self.assertEqual(len(names), 11)
        self.assertEqual(len(validation), 9)
        self.assertFalse(set(development) & set(validation))
        self.assertEqual(set(development) | set(validation), set(names))
        self.assertEqual(development, [name for name in names if name in development])

    def test_smoke_selects_one_complete_development_recording(self):
        self.assertEqual(runner.select_sequences(runner.datasets.SEQUENCES, smoke=True), ["MH_01_easy"])
        for options in ({"split": "development"}, {"patterns": ("MH_*",)}):
            with self.subTest(options=options), self.assertRaisesRegex(ValueError, "cannot be combined"):
                runner.select_sequences(runner.datasets.SEQUENCES, smoke=True, **options)
        with self.assertRaisesRegex(ValueError, "missing from manifest"):
            runner.select_sequences(["V1_02_medium"], smoke=True)

    def test_patterns_select_named_recordings_before_loading(self):
        self.assertEqual(runner.select_sequences(runner.datasets.SEQUENCES, patterns=("V2_*",)),
                         ["V2_01_easy", "V2_02_medium", "V2_03_difficult"])
        self.assertEqual(runner.select_sequences(runner.datasets.SEQUENCES, split="development",
                                                 patterns=("MH_*",)), ["MH_01_easy"])
        for split, patterns in (("all", ("misspelled*",)), ("validation", ("MH_01_easy",)),
                                ("all", ("V2_*", "missing*"))):
            with self.subTest(split=split, patterns=patterns):
                with self.assertRaisesRegex(ValueError, "matched nothing"):
                    runner.select_sequences(runner.datasets.SEQUENCES, split=split, patterns=patterns)

    def test_primary_attitude_comparison_excludes_dependent_or_known_bad_reference(self):
        primary = {name for name in runner.datasets.SEQUENCES if runner.independent_attitude(name)}
        self.assertEqual(primary, {"V1_02_medium", "V1_03_difficult", "V2_01_easy",
                                   "V2_02_medium", "V2_03_difficult"})
        self.assertEqual(len(primary - runner.DEVELOPMENT), 4)


class SummaryTests(unittest.TestCase):
    def test_summary_weights_recordings_equally_instead_of_samples(self):
        rows = [row("short", "mav-srukf", 1, samples=10),
                row("long", "mav-srukf", 9, samples=100000), row("short", "vqf", 3)]
        result = {item["algorithm"]: item for item in runner.aggregate(rows)}
        self.assertEqual(result["mav-srukf"]["sequences"], 2)
        self.assertEqual(result["mav-srukf"]["inclination_rmse_deg"], 5)
        self.assertEqual(result["mav-srukf"]["total_rmse_deg"], 15)
        self.assertEqual(result["vqf"]["inclination_rmse_deg"], 3)

    def test_paired_comparison_matches_names_and_uses_ratio_of_means(self):
        rows = [row("A", "mav-srukf", 1), row("B", "vqf", 4),
                row("B", "mav-srukf", 8), row("A", "vqf", 2), row("unrelated", "ukf", 100)]
        comparison = runner.paired_comparisons(rows)
        self.assertEqual(comparison["sequences"], 2)
        self.assertEqual(comparison["inclination"]["mean_difference_deg"], 1.5)
        self.assertEqual(comparison["inclination"]["ratio_of_means"], 1.5)
        self.assertEqual(comparison["inclination"]["candidate_wins"], 1)
        self.assertEqual(comparison["total"]["mean_difference_deg"], 4.5)

    def test_incomplete_pairs_cannot_disappear_from_comparison(self):
        with self.assertRaisesRegex(ValueError, "different recording coverage"):
            runner.paired_comparisons([row("A", "mav-srukf", 1), row("A", "vqf", 2),
                                       row("B", "mav-srukf", 3)])
        self.assertIsNone(runner.paired_comparisons([row("A", "vqf", 2)]))

    def test_disjoint_or_duplicate_recording_results_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "different recording coverage"):
            runner.paired_comparisons([row("A", "mav-srukf", 1), row("B", "vqf", 2)])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            runner.paired_comparisons([row("A", "mav-srukf", 1), row("A", "vqf", 2),
                                       row("A", "mav-srukf", 9)])


class ReferenceClockTests(unittest.TestCase):
    def test_known_timestamp_quantization_uses_nominal_clock_without_interpolation(self):
        # Construct integer epoch timestamps on the source's 256 ns grid before
        # converting to relative seconds; no large epoch float is needed.
        ideal_ns = np.int64(1403636579763555584) + np.arange(1000, dtype=np.int64)*5000000
        quantized_ns = ((ideal_ns+128)//256)*256
        source_time = (quantized_ns-quantized_ns[0]).astype(float)*1e-9
        original = source_time.copy()
        nominal, metadata = runner.vqf_clock(source_time)
        np.testing.assert_allclose(nominal, np.arange(1000)*.005, rtol=0, atol=1e-14)
        np.testing.assert_array_equal(source_time, original)
        self.assertGreater(metadata["max_timestamp_deviation_seconds"], 0)
        self.assertLess(metadata["max_timestamp_deviation_seconds"], 1e-6)
        self.assertEqual(metadata["nominal_interval_seconds"], .005)

    def test_real_sampling_gap_or_jitter_cannot_be_hidden_by_nominal_clock(self):
        nominal = np.arange(1000)*.005
        jittered = nominal.copy()
        jittered[30] += .0001
        for source in (np.delete(nominal, 30), np.arange(1000)*.006, jittered):
            with self.subTest(intervals=np.diff(source)[:3]):
                with self.assertRaisesRegex(ValueError, "timing gap or variable-rate"):
                    runner.vqf_clock(source)


class NativeReplayTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.trial = SimpleNamespace(
            name="sensor_only_fixture", time=np.array([0., .005, .010]),
            gyr=np.array([[.01, .02, .03], [.04, .05, .06], [.07, .08, .09]]),
            acc=np.array([[1., 2., 9.], [2., 3., 8.], [3., 4., 7.]]),
            # These deliberately conspicuous fields must not reach the native
            # estimator, even if a caller provides them on a Trial object.
            mag=np.full((3, 3), 12345.), quat_gt=np.full((3, 4), 67890.),
            reference_velocity=np.full((3, 3), 99999.),
        )

    def fake_process(self, command, **kwargs):
        self.assertTrue(kwargs["check"])
        self.assertEqual(command[command.index("--axes")+1], "6")
        source = Path(command[command.index("--input")+1])
        self.assertEqual(source.read_text().splitlines()[0], "time,gx,gy,gz,ax,ay,az,mx,my,mz")
        actual = np.loadtxt(source, delimiter=",", skiprows=1)
        expected = np.column_stack((self.trial.time, self.trial.gyr, self.trial.acc, np.zeros((3, 3))))
        np.testing.assert_array_equal(actual, expected)
        output = Path(command[command.index("--output")+1])
        quaternions = np.tile([1., 0., 0., 0.], (3, 1))
        np.savetxt(output, np.column_stack((self.trial.time, quaternions)), delimiter=",",
                   header="time,qw,qx,qy,qz", comments="", fmt="%.17g")
        return SimpleNamespace(stderr=json.dumps({"elapsed_seconds": .125}))

    def test_native_input_contains_only_measured_imu_and_unused_zero_mag(self):
        with patch.object(runner.subprocess, "run", side_effect=self.fake_process) as process:
            quaternions, elapsed = runner.cpp_estimate(Path("estimator"), self.trial,
                                                       "mav-srukf", self.directory)
        process.assert_called_once()
        np.testing.assert_array_equal(quaternions, np.tile([1., 0., 0., 0.], (3, 1)))
        self.assertEqual(elapsed, .125)

    def test_requested_diagnostics_must_exist(self):
        with patch.object(runner.subprocess, "run", side_effect=self.fake_process):
            with self.assertRaisesRegex(ValueError, "did not write requested diagnostics"):
                runner.cpp_estimate(Path("estimator"), self.trial, "mav-srukf", self.directory,
                                    self.directory / "diagnostics" / "states.csv")

    def test_prior_output_cannot_masquerade_as_current_estimate(self):
        destination = self.directory / "output.csv"
        np.savetxt(destination, np.column_stack((self.trial.time, np.tile([1., 0., 0., 0.], (3, 1)))),
                   delimiter=",", header="time,qw,qx,qy,qz", comments="")
        diagnostics = self.directory / "states.csv"
        diagnostics.write_text("stale diagnostic from previous execution\n")

        def no_artifacts(command, **kwargs):
            self.assertFalse(destination.exists(), "Old attitude output must be removed before execution")
            self.assertFalse(diagnostics.exists(), "Old diagnostics must be removed before execution")
            return SimpleNamespace(stderr=json.dumps({"elapsed_seconds": .1}))

        with patch.object(runner.subprocess, "run", side_effect=no_artifacts):
            with self.assertRaisesRegex(ValueError, "did not write requested diagnostics"):
                runner.cpp_estimate(Path("estimator"), self.trial, "mav-srukf", self.directory, diagnostics)

    def test_native_failures_are_not_silently_dropped(self):
        failure = subprocess.CalledProcessError(1, ["estimator"], stderr="invalid process state")
        with patch.object(runner.subprocess, "run", side_effect=failure):
            with self.assertRaisesRegex(RuntimeError, "mav-srukf failed on sensor_only_fixture: invalid process state"):
                runner.cpp_estimate(Path("estimator"), self.trial, "mav-srukf", self.directory)


if __name__ == "__main__":
    unittest.main()
