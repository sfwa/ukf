"""Executable boundary checks: reject invalid requests before writing estimates."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = Path(os.environ.get(
    "EUROC_BENCHMARK_EXECUTABLE",
    ROOT / "build/examples/euroc_benchmark/euroc_benchmark",
)).resolve()
HEADER = "time,gx,gy,gz,ax,ay,az,mx,my,mz\n"
ROW0 = "0,0,0,0,9.80665,0,0,0,0,0\n"
ROW1 = "0.005,0,0,0,9.80665,0,0,0,0,0\n"


class MAVCommandLineTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.source = self.directory / "sensors.csv"
        self.output = self.directory / "attitude.csv"
        self.source.write_text(HEADER+ROW0+ROW1)

    def command(self, extra=(), algorithm="mav-srukf"):
        return [str(EXECUTABLE), "--input", str(self.source), "--output", str(self.output),
                "--filter", algorithm, "--axes", "6", *extra]

    def rejects(self, command):
        original = self.source.read_bytes()
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("euroc_benchmark:", result.stderr)
        self.assertEqual(self.source.read_bytes(), original)
        self.assertFalse(self.output.exists(), "Rejected input must not produce an estimate")

    def test_required_filter_and_six_axis_selection(self):
        command = self.command()
        for omitted in ("--filter", "--axes"):
            with self.subTest(omitted=omitted):
                changed = command.copy()
                index = changed.index(omitted)
                del changed[index:index+2]
                self.rejects(changed)
        changed = command.copy()
        changed[changed.index("--axes")+1] = "9"
        self.rejects(changed)

    def test_unknown_options_filters_and_missing_values_are_rejected(self):
        # The tuning is fixed in mav_estimator.h; it cannot be changed here.
        for extra in (("--drag-rate", "0.3"), ("--unknown", "1"), ("--diagnostics",)):
            with self.subTest(extra=extra):
                self.rejects(self.command(extra))
        self.rejects(self.command(algorithm="srukf"))
        self.rejects([str(EXECUTABLE), "--describe-filter", "srukf"])

    def test_output_and_diagnostic_paths_cannot_overwrite_inputs(self):
        same_output = self.command()
        same_output[same_output.index("--output")+1] = str(self.source)
        self.rejects(same_output)
        for diagnostic in (self.source, self.output):
            with self.subTest(diagnostic=str(diagnostic)):
                self.rejects(self.command(("--diagnostics", str(diagnostic))))

    def test_malformed_csv_and_unobservable_initial_tilt_are_rejected(self):
        cases = {
            "header": "time,gyro\n0,0\n0.005,0\n",
            "one_sample": HEADER+ROW0,
            "duplicate_time": HEADER+ROW0+ROW0,
            "nonfinite_gyro": HEADER+ROW0+ROW1.replace("0.005,0,", "0.005,nan,", 1),
            "zero_initial_acceleration": HEADER+ROW0.replace("9.80665", "0")+ROW1,
            "missing_column": HEADER+ROW0+"0.005,0,0,0,9.80665,0,0,0,0\n",
            "extra_column": HEADER+ROW0+ROW1.rstrip()+",7\n",
        }
        for label, contents in cases.items():
            with self.subTest(case=label):
                self.source.write_text(contents)
                self.rejects(self.command())

    def test_successful_control_is_identified_in_configuration_and_timing(self):
        for algorithm in ("mav-srukf", "mav-srukf-gravity"):
            with self.subTest(algorithm=algorithm):
                description = subprocess.run([str(EXECUTABLE), "--describe-filter", algorithm],
                                             capture_output=True, text=True, check=True, timeout=30)
                self.assertEqual(json.loads(description.stdout)["stage"], algorithm)
                result = subprocess.run(self.command(algorithm=algorithm), capture_output=True,
                                        text=True, check=True, timeout=30)
                self.assertEqual(json.loads(result.stderr)["filter"], algorithm)
                values = np.loadtxt(self.output, delimiter=",", skiprows=1, ndmin=2)
                self.assertEqual(values.shape, (2, 5))
                np.testing.assert_allclose(np.linalg.norm(values[:, 1:], axis=1), 1, atol=1e-10, rtol=0)
                self.output.unlink()


if __name__ == "__main__":
    unittest.main()
