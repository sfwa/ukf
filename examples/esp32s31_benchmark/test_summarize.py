"""Interpret logged results without relying on host performance or hardware."""

import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("summarize", HERE / "summarize.py")
SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY)


def suite(platform="esp32s31"):
    records = [{"type": "benchmark_input", "platform": platform, "source": "test.csv",
                "samples": 101, "passes": 2, "dt_min_s": 0.002, "dt_max_s": 0.002,
                "dt_mean_s": 0.002}]
    for algorithm, iteration in (("ukf", 1), ("srukf", 1), ("srukf", 2), ("ukf", 2)):
        records.append({"type": "attitude_benchmark", "platform": platform,
                        "filter": algorithm, "pass": iteration, "precision": "float32",
                        "axes": 9, "timed_scope": "update_and_attitude", "budget_us": 2000,
                        "expected_updates": 100, "replay_valid": True,
                        "observed_compute_within_budget": True,
                        "timing": {"count": 100, "mean_us": 500.5, "max_us": 2000, "over_budget": 0},
                        "p95_us": 700, "p99_us": 900})
    if platform in ("esp32s31", "stm32h735"):
        records.insert(0, {"type": "platform", "platform": platform})
        records.append({"type": "platform_complete", "suite_ok": True})
    return records


def log(records):
    return "\n".join(json.dumps(record) for record in records)


class SummarizeTests(unittest.TestCase):
    def test_hardware_success_and_exact_budget_boundary(self):
        for platform, label in (("esp32s31", "ESP32-S31"), ("stm32h735", "STM32H735")):
            with self.subTest(platform=platform):
                report, status = SUMMARY.summarize(log(suite(platform)))
                self.assertEqual(status, 0)
                self.assertIn(f"{label} on-device **compute-only**", report)
                self.assertIn("100.00%", report)
                self.assertEqual(report.count("| UKF |"), 2)
                self.assertEqual(report.count("| SR-UKF |"), 2)
                self.assertIn("percentiles are not combined", report)

    def test_emulated_overrun_is_simulated_and_not_a_silicon_failure(self):
        records = suite()
        records[2]["timing"].update(max_us=12000, over_budget=1)
        records[2]["observed_compute_within_budget"] = False
        report, status = SUMMARY.summarize(log(records), emulated=True)
        self.assertEqual(status, 0)
        self.assertIn("SIMULATED", report)
        self.assertNotIn("on-device", report)
        self.assertIn("cannot establish", report)
        self.assertIn("| 1 | yes |", report)

    def test_emulated_requires_all_four_replays_and_final_success(self):
        records = suite()
        with self.assertRaises(ValueError):
            SUMMARY.summarize(log(records[:-1]), emulated=True)
        with self.assertRaises(ValueError):
            SUMMARY.summarize(log(records[:3] + records[-1:]), emulated=True)
        records[-1]["suite_ok"] = False
        report, status = SUMMARY.summarize(log(records), emulated=True)
        self.assertEqual(status, 1)
        self.assertIn("FAIL: numerical", report)
        records = suite()
        records[2]["timing"]["count"] = 99
        self.assertEqual(SUMMARY.summarize(log(records), emulated=True)[1], 1)

    def test_emulated_rejects_host_and_cli_labels_firmware(self):
        for platform in ("host", "stm32h735"):
            with self.subTest(platform=platform), self.assertRaisesRegex(ValueError, "requires an ESP32-S31"):
                SUMMARY.summarize(log(suite(platform)), emulated=True)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "serial.log"
            path.write_text(log(suite()))
            result = subprocess.run([sys.executable, str(HERE / "summarize.py"),
                                     "--emulated", str(path)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0)
            self.assertIn("SIMULATED", result.stdout)
            self.assertNotIn("on-device", result.stdout)

    def test_hardware_overrun_fails(self):
        for platform in ("esp32s31", "stm32h735"):
            with self.subTest(platform=platform):
                records = suite(platform)
                run = records[2]
                run["timing"]["max_us"] = 2001
                run["timing"]["over_budget"] = 1
                run["observed_compute_within_budget"] = False
                report, status = SUMMARY.summarize(log(records))
                self.assertEqual(status, 1)
                self.assertIn("100.05%", report)
                self.assertIn("did not meet", report)

    def test_stm32_cli_labels_hardware_and_rejects_failed_completion(self):
        records = suite("stm32h735")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "serial.log"
            path.write_text(log(records))
            result = subprocess.run([sys.executable, str(HERE / "summarize.py"), str(path)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0)
            self.assertIn("STM32H735 on-device", result.stdout)
            self.assertNotIn("ESP32", result.stdout)
            self.assertNotIn("SIMULATED", result.stdout)
            records[-1]["suite_ok"] = False
            path.write_text(log(records))
            result = subprocess.run([sys.executable, str(HERE / "summarize.py"), str(path)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 1)
            self.assertIn("FAIL: numerical", result.stdout)

    def test_stm32_requires_complete_unmixed_firmware_suite(self):
        records = suite("stm32h735")
        candidates = [records[1:], records[:-1], records[:3] + records[-1:], records + records]
        mixed = copy.deepcopy(records)
        mixed[2]["platform"] = "esp32s31"
        candidates.append(mixed)
        for candidate in candidates:
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                SUMMARY.summarize(log(candidate))

    def test_host_overrun_is_not_a_target_performance_failure(self):
        records = suite("host")
        records[1]["timing"]["max_us"] = 12000
        records[1]["timing"]["over_budget"] = 1
        records[1]["observed_compute_within_budget"] = False
        report, status = SUMMARY.summarize(log(records))
        self.assertEqual(status, 0)
        self.assertIn("NOT ESP32-S31 target performance", report)
        self.assertIn("not used as a pass/fail gate", report)

    def test_numerical_failure_and_partial_replay_fail(self):
        for platform in ("host", "esp32s31", "stm32h735"):
            for invalid in ({"replay_valid": False}, {"timing": {"count": 50, "mean_us": 200,
                                                                             "max_us": 2000, "over_budget": 0}}):
                records = suite(platform)
                next(row for row in records if row["type"] == "attitude_benchmark").update(invalid)
                report, status = SUMMARY.summarize(log(records))
                self.assertEqual(status, 1)
                self.assertIn("| NO |", report)

    def test_missing_duplicate_mixed_and_conflicting_records_rejected(self):
        records = suite()
        bad_logs = [[], records[:-1], records[:3], records + records,
                    records + [copy.deepcopy(records[2])]]
        mixed = copy.deepcopy(records)
        mixed[2]["platform"] = "host"
        bad_logs.append(mixed)
        mismatch = copy.deepcopy(records)
        mismatch[2]["expected_updates"] = 99
        bad_logs.append(mismatch)
        for candidate in bad_logs:
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                SUMMARY.summarize(log(candidate))

    def test_monitor_noise_ignored_but_malformed_result_rejected(self):
        text = 'I (100) boot: starting\n{"unrelated": true}\n' + log(suite()) + '\n--- monitor stopped ---'
        _, status = SUMMARY.summarize(text)
        self.assertEqual(status, 0)
        with self.assertRaisesRegex(ValueError, "Malformed benchmark JSON"):
            SUMMARY.summarize(text + '\n{"type":"attitude_benchmark",broken}')

    def test_cli_returns_nonzero_on_incomplete_suite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "serial.log"
            path.write_text(log(suite()[:3]))
            result = subprocess.run([sys.executable, str(HERE / "summarize.py"), str(path)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("requires platform", result.stderr)


if __name__ == "__main__":
    unittest.main()
