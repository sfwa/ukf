"""Executable-level tests with analytically generated sensor measurements.

Run after building attitude_benchmark. Override its path with
ATTITUDE_BENCHMARK_EXECUTABLE when using a different build directory.
"""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.spatial.transform import Rotation


ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = Path(os.environ.get(
    "ATTITUDE_BENCHMARK_EXECUTABLE",
    ROOT / "build/examples/attitude_benchmark/attitude_benchmark",
)).resolve()
HEADER = "time,gx,gy,gz,ax,ay,az,mx,my,mz"


class EstimatorTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name)

    def sensors(self, time, truth, gyro):
        """Specific force at zero translation and a fixed world magnetic field."""
        acc = truth.inv().apply(np.tile([0, 0, 9.80665], (len(time), 1)))
        mag = truth.inv().apply(np.tile([0.6, 0, -0.8], (len(time), 1)))
        values = np.column_stack((time, gyro, acc, mag))
        np.savetxt(self.path / "sensors.csv", values, delimiter=",", header=HEADER, comments="")

    def run_filter(self, algorithm, axes, expect_success=True):
        output = self.path / f"{algorithm}-{axes}.csv"
        result = subprocess.run([
            str(EXECUTABLE), "--input", str(self.path / "sensors.csv"),
            "--output", str(output), "--filter", algorithm, "--axes", str(axes),
        ], capture_output=True, text=True)
        if not expect_success:
            self.assertNotEqual(result.returncode, 0)
            return
        self.assertEqual(result.returncode, 0, result.stderr)
        timing = json.loads(result.stderr)
        data = np.loadtxt(output, delimiter=",", skiprows=1)
        self.assertEqual(timing["samples"], len(data))
        self.assertEqual(timing["updates"], len(data) - 1)
        self.assertGreater(timing["elapsed_seconds"], 0)
        np.testing.assert_allclose(np.linalg.norm(data[:, 1:], axis=1), 1, atol=1e-12)
        return Rotation.from_quat(data[:, [2, 3, 4, 1]])

    def test_stationary_tilt(self):
        time = np.arange(501) * 0.01
        initial = Rotation.from_euler("xyz", [35, -20, 60], degrees=True)
        truth = Rotation.from_quat(np.tile(initial.as_quat(), (len(time), 1)))
        self.sensors(time, truth, np.zeros((len(time), 3)))
        for axes in (6, 9):
            for algorithm in ("ukf", "srukf"):
                with self.subTest(axes=axes, algorithm=algorithm):
                    estimated = self.run_filter(algorithm, axes)
                    # Only a constant initial world yaw is unobservable in six-axis mode.
                    if axes == 6:
                        estimated = (truth[0] * estimated[0].inv()) * estimated
                    error = np.rad2deg((truth.inv() * estimated).magnitude())
                    self.assertLess(error.max(), 0.01)

    def test_multiaxis_rotation_and_core_equivalence(self):
        # Nonzero roll, pitch and yaw, rotating around a non-coordinate body axis.
        # There is no translational acceleration: accelerometer truth is exact.
        time = np.arange(1001) * 0.01
        omega = np.array([0.3, -0.5, 0.8])
        initial = Rotation.from_euler("xyz", [35, -20, 60], degrees=True)
        truth = initial * Rotation.from_rotvec(time[:, None] * omega)
        self.sensors(time, truth, np.tile(omega, (len(time), 1)))
        for axes in (6, 9):
            estimates = []
            for algorithm in ("ukf", "srukf"):
                with self.subTest(axes=axes, algorithm=algorithm):
                    estimated = self.run_filter(algorithm, axes)
                    estimates.append(estimated)
                    if axes == 6:
                        estimated = (truth[0] * estimated[0].inv()) * estimated
                    error = np.rad2deg((truth.inv() * estimated).magnitude())
                    self.assertLess(error.max(), 0.02)
            difference = np.rad2deg((estimates[0].inv() * estimates[1]).magnitude())
            self.assertLess(difference.max(), 1e-5)

    def test_stationary_gyro_bias_convergence(self):
        time = np.arange(6001) * 0.01
        initial = Rotation.from_euler("xyz", [35, -20, 60], degrees=True)
        truth = Rotation.from_quat(np.tile(initial.as_quat(), (len(time), 1)))
        self.sensors(time, truth, np.tile([0.01, -0.015, 0.02], (len(time), 1)))
        for axes in (6, 9):
            for algorithm in ("ukf", "srukf"):
                with self.subTest(axes=axes, algorithm=algorithm):
                    estimated = self.run_filter(algorithm, axes)
                    gravity_est = estimated.inv().apply(np.tile([0, 0, 1], (len(time), 1)))
                    gravity_true = initial.inv().apply([0, 0, 1])
                    inclination = np.rad2deg(np.arccos(np.clip(gravity_est @ gravity_true, -1, 1)))
                    self.assertLess(inclination[-1000:].max(), 0.1)
                    if axes == 9:
                        error = np.rad2deg((truth.inv() * estimated).magnitude())
                        self.assertLess(error[-1000:].max(), 0.1)
                        # Convergence removes almost all of the injected ~1 deg/s drift.
                        drift = np.rad2deg((estimated[-1001].inv() * estimated[-1]).magnitude()) / 10
                        self.assertLess(drift, 0.002)

    def test_bad_input_fails(self):
        row = "0,0,0,0,0,0,9.81,1,0,0"
        cases = {
            "empty": "",
            "wrong_header": "t,x\n0,0\n1,0\n",
            "single_sample": HEADER + "\n" + row + "\n",
            "duplicate_time": HEADER + "\n" + row + "\n" + row + "\n",
            "nan": HEADER + "\n" + row + "\n1,nan,0,0,0,0,9.81,1,0,0\n",
            # The first accelerometer reading sets the initial tilt.
            "zero_initial_acc": HEADER + "\n0,0,0,0,0,0,0,1,0,0\n1,0,0,0,0,0,9.81,1,0,0\n",
            "extra_column": HEADER + "\n" + row + "\n1,0,0,0,0,0,9.81,1,0,0,2\n",
        }
        for name, contents in cases.items():
            with self.subTest(name=name):
                (self.path / "sensors.csv").write_text(contents)
                self.run_filter("ukf", 9, expect_success=False)


if __name__ == "__main__":
    unittest.main()
