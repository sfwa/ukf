"""Validate replay and deadline accounting; never assert host performance."""

import json
import math
import os
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = Path(os.environ.get(
    "ESP32S31_BENCHMARK_EXECUTABLE",
    ROOT / "build/examples/esp32s31_benchmark/esp32s31_benchmark_host",
))


class BenchmarkTests(unittest.TestCase):
    def test_deadline_boundary_and_percentiles(self):
        subprocess.run([str(EXECUTABLE), "--check-timing"], check=True, timeout=30)

    def test_both_filters_replay_all_samples_and_report_consistent_results(self):
        result = subprocess.run([str(EXECUTABLE)], capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("not ESP32-S31 performance measurements", result.stderr)
        records = [json.loads(line) for line in result.stdout.splitlines()]
        metadata, *runs = records
        self.assertEqual(metadata["type"], "benchmark_input")
        self.assertEqual(metadata["platform"], "host")
        self.assertEqual(metadata["sigma_points"], 19)
        self.assertEqual(metadata["covariance_dimensions"], 9)
        self.assertEqual([(run["filter"], run["pass"]) for run in runs],
                         [("ukf", 1), ("srukf", 1), ("srukf", 2), ("ukf", 2)])
        for run in runs:
            self.assertEqual(run["precision"], "float32")
            self.assertEqual(run["budget_us"], 2000)
            self.assertTrue(run["replay_valid"])
            self.assertEqual(run["failed_sample"], 0)
            timing = run["timing"]
            self.assertEqual(timing["count"], metadata["samples"] - 1)
            self.assertEqual(timing["count"], run["expected_updates"])
            self.assertEqual(run["rest"]["count"] + run["moving"]["count"], timing["count"])
            self.assertEqual(run["magnetic_accepted"]["count"] + run["magnetic_rejected"]["count"],
                             timing["count"])
            self.assertEqual(run["observed_compute_within_budget"], timing["over_budget"] == 0)
            self.assertEqual(run["min_headroom_us"], 2000 - timing["max_us"])
            self.assertAlmostEqual(run["max_budget_percent"], timing["max_us"] / 20, places=3)
            quantiles = [timing["min_us"], run["p50_us"], run["p95_us"], run["p99_us"], timing["max_us"]]
            self.assertEqual(quantiles, sorted(quantiles))
            self.assertAlmostEqual(sum(x*x for x in run["final_attitude_wxyz"]), 1, places=5)
            self.assertTrue(math.isfinite(run["checksum"]))
        if metadata["source"].startswith("synthetic_"):
            self.assertAlmostEqual(metadata["dt_mean_s"], 0.002)
            for run in runs:
                for branch in ("rest", "moving", "magnetic_accepted", "magnetic_rejected"):
                    self.assertGreater(run[branch]["count"], 100)
            a, b = runs[0]["final_attitude_wxyz"], runs[1]["final_attitude_wxyz"]
            # Compare equivalent quaternion signs without the ill-conditioned acos near 1.
            distance = min(sum((x-y)**2 for x, y in zip(a, b)),
                           sum((x+y)**2 for x, y in zip(a, b)))
            self.assertLess(distance, 1e-6)
        for first, second in ((runs[0], runs[3]), (runs[1], runs[2])):
            self.assertEqual(first["final_attitude_wxyz"], second["final_attitude_wxyz"])
            self.assertEqual(first["checksum"], second["checksum"])


if __name__ == "__main__":
    unittest.main()
