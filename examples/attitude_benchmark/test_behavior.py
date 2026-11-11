"""Synthetic behavioral checks for the UKF attitude estimator.

These recordings have independently known motion, sensor bias and disturbances.
No recorded dataset or ground truth is passed to the executable. Both filter
cores run the same model, so every check applies to `ukf` and `srukf` alike.
"""
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
FILTERS = ("ukf", "srukf")


class EstimatorBehaviorTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.run_number = 0

    @staticmethod
    def recording(time, omega=None, bias=None, translation=None):
        """Constant body angular rate and optional world-frame acceleration."""
        omega = np.zeros(3) if omega is None else np.asarray(omega)
        bias = np.zeros(3) if bias is None else np.asarray(bias)
        elapsed = time - time[0]
        initial = Rotation.from_euler("xyz", [20, -15, 40], degrees=True)
        truth = initial * Rotation.from_rotvec(elapsed[:, None] * omega)
        acceleration = np.tile([0., 0., 9.80665], (len(time), 1))
        if translation is not None:
            acceleration += translation
        acc = truth.inv().apply(acceleration)
        mag = truth.inv().apply(np.tile([0.6, 0., -0.8], (len(time), 1)))
        gyro = np.tile(omega + bias, (len(time), 1))
        return np.column_stack((time, gyro, acc, mag)), truth

    def replay(self, samples, axes=6, algorithm="srukf"):
        self.run_number += 1
        stem = self.directory / str(self.run_number)
        source, output, diagnostics = (stem.with_suffix(ext) for ext in (".in.csv", ".out.csv", ".diag.csv"))
        np.savetxt(source, samples, delimiter=",", header=HEADER, comments="", fmt="%.17g")
        result = subprocess.run([
            str(EXECUTABLE), "--input", str(source), "--output", str(output),
            "--filter", algorithm, "--axes", str(axes), "--diagnostics", str(diagnostics),
        ], capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        values = np.loadtxt(output, delimiter=",", skiprows=1, ndmin=2)
        self.assertEqual(values.shape, (len(samples), 5))
        np.testing.assert_allclose(values[:, 0], samples[:, 0], atol=1e-10, rtol=0)
        self.assertTrue(np.isfinite(values).all())
        np.testing.assert_allclose(np.linalg.norm(values[:, 1:], axis=1), 1, atol=1e-10, rtol=0)
        details = np.genfromtxt(diagnostics, delimiter=",", names=True, ndmin=1)
        self.assertEqual(len(details), len(samples))
        return Rotation.from_quat(values[:, 1:], scalar_first=True), details

    @staticmethod
    def gravity(quaternions):
        return quaternions.inv().apply(np.tile([0., 0., 1.], (len(quaternions), 1)))

    def tilt_error(self, estimate, truth):
        cosine = np.sum(self.gravity(estimate) * self.gravity(truth), axis=1)
        return np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))

    @staticmethod
    def bias(details):
        return np.column_stack([details[name] for name in ("bias_x", "bias_y", "bias_z")])

    def test_rest_observes_all_three_gyro_bias_components(self):
        time = np.arange(2001) / 100
        bias = np.array([0.01, -0.015, 0.02])
        samples, _ = self.recording(time, bias=bias)
        for algorithm in FILTERS:
            with self.subTest(algorithm=algorithm):
                _, details = self.replay(samples, algorithm=algorithm)
                detected = np.flatnonzero(details["rest"])
                self.assertGreater(len(detected), 0)
                self.assertGreaterEqual(time[detected[0]], 1.49)
                self.assertLess(time[detected[0]], 3)
                np.testing.assert_allclose(self.bias(details)[-1], bias, atol=5e-4, rtol=0)
                self.assertTrue((details["rest"][time > 5] == 1).all())

    def test_constant_rotation_is_not_mistaken_for_rest(self):
        time = np.arange(2001) / 100
        # Rotating around the world vertical keeps acceleration constant. A
        # detector based solely on variance would incorrectly declare rest.
        initial = Rotation.from_euler("xyz", [20, -15, 40], degrees=True)
        omega = initial.inv().apply([0., 0., np.deg2rad(20)])
        samples, _ = self.recording(time, omega=omega)
        _, details = self.replay(samples)
        self.assertFalse(details["rest"].any())
        np.testing.assert_allclose(self.bias(details), 0, atol=1e-4, rtol=0)

    def test_translation_does_not_poison_tilt_or_gyro_bias(self):
        time = np.arange(4001) / 100
        bias = np.array([0.01, -0.015, 0.02])
        moving = (time >= 6) & (time < 26)
        translation = np.zeros((len(time), 3))
        translation[moving, 0] = 3 * np.sin(2 * np.pi * 0.7 * (time[moving] - 6))
        translation[moving, 1] = 2 * np.sin(2 * np.pi * 1.1 * (time[moving] - 6))
        samples, truth = self.recording(time, bias=bias, translation=translation)
        for algorithm in FILTERS:
            with self.subTest(algorithm=algorithm):
                estimate, details = self.replay(samples, algorithm=algorithm)
                # Rest before the motion measures the bias; the motion must
                # neither be mistaken for rest nor disturb that estimate.
                np.testing.assert_allclose(self.bias(details)[time >= 6],
                                           np.broadcast_to(bias, (np.count_nonzero(time >= 6), 3)),
                                           atol=8e-4, rtol=0)
                self.assertLess(np.mean(details["rest"][(time > 8) & (time < 24)]), 0.1)
                # The accelerometer points up to 18 degrees away from vertical
                # during this motion. Treating it as gravity would tilt the
                # estimate by a similar amount; integrating it must not.
                apparent = np.rad2deg(np.arctan2(np.linalg.norm(translation[:, :2], axis=1), 9.80665))
                self.assertGreater(apparent.max(), 15)
                tilt = self.tilt_error(estimate, truth)
                self.assertLess(np.sqrt(np.mean(tilt[moving] ** 2)), 1)
                self.assertLess(tilt[time > 38].max(), 0.5)

    def test_prefix_causality(self):
        time = np.arange(4001) / 100
        acceleration = np.zeros((len(time), 3))
        acceleration[time > 20, 0] = 3 * np.sin(2 * np.pi * time[time > 20])
        samples, _ = self.recording(time, bias=[0.004, -0.006, 0.008], translation=acceleration)
        # A future magnetic perturbation must not alter any earlier estimate.
        samples[time > 25, 7:10] += [1., -2., 3.]
        full, _ = self.replay(samples, axes=9)
        prefix, _ = self.replay(samples[:1501], axes=9)
        np.testing.assert_allclose(full[:1501].as_quat(), prefix.as_quat(), atol=1e-12, rtol=0)

    def test_magnetometer_cannot_change_six_axis_solution_or_nine_axis_tilt(self):
        time = np.arange(3001) / 100
        samples, _ = self.recording(time, omega=[0.07, -0.04, 0.16])
        disturbed = samples.copy()
        # Include initialization: six-axis estimates cannot depend on the first
        # magnetic sample either. Later discontinuities exercise updates.
        disturbed[:, 7:10] *= -2
        disturbed[time > 8, 7:10] = [2., -3., 4.]
        clean6, clean_details = self.replay(samples, axes=6)
        disturbed6, disturbed_details = self.replay(disturbed, axes=6)
        disturbed9, heading_details = self.replay(disturbed, axes=9)
        np.testing.assert_allclose(clean6.as_quat(), disturbed6.as_quat(), atol=1e-12, rtol=0)
        np.testing.assert_allclose(self.gravity(clean6), self.gravity(disturbed9), atol=1e-10, rtol=0)
        np.testing.assert_allclose(self.bias(clean_details), self.bias(disturbed_details), atol=1e-12, rtol=0)
        np.testing.assert_allclose(self.bias(clean_details), self.bias(heading_details), atol=1e-12, rtol=0)

    def test_magnetic_disturbance_is_rejected_and_clean_field_recovers(self):
        time = np.arange(6001) / 100
        samples, _ = self.recording(time)
        changed = samples.copy()
        disturbed = (time >= 10) & (time < 15)
        # A large horizontal change also changes field norm/dip, making this a
        # physically detectable disturbance rather than an ambiguous pure yaw.
        changed[disturbed, 7:10] += [2., -3., 4.]
        clean, _ = self.replay(samples, axes=9)
        estimate, details = self.replay(changed, axes=9)
        self.assertGreater(np.mean(details["magnetometer_rejected"][(time > 11) & (time < 14)]), 0.95)
        self.assertLess(np.mean(details["magnetometer_rejected"][time > 20]), 0.01)
        # No magnetic corruption is allowed to leak into tilt, even transiently.
        np.testing.assert_allclose(self.gravity(clean), self.gravity(estimate), atol=1e-10, rtol=0)
        discrepancy = np.rad2deg((clean.inv() * estimate).magnitude())
        self.assertLess(discrepancy[time > 50].max(), 1)

    def test_nine_axis_heading_points_to_magnetic_north(self):
        time = np.arange(3001) / 100
        samples, truth = self.recording(time, omega=[0.07, -0.04, 0.16])
        # The recording's magnetic field lies in the world X-Z plane, so the
        # true attitude already has magnetic north along world +X.
        estimate, _ = self.replay(samples, axes=9)
        error = np.rad2deg((truth.inv() * estimate).magnitude())
        self.assertLess(error.max(), 0.2)

    def test_absolute_timestamp_origin_does_not_change_state(self):
        time = np.arange(1201) / 100
        samples, _ = self.recording(time, bias=[0.004, -0.006, 0.008])
        shifted = samples.copy()
        shifted[:, 0] += 10000
        original, _ = self.replay(samples, axes=9)
        offset, _ = self.replay(shifted, axes=9)
        np.testing.assert_allclose(original.as_quat(), offset.as_quat(), atol=1e-9, rtol=0)

    def test_variable_sampling_intervals_preserve_physical_rotation(self):
        # Large sampling jitter exercises dt-dependent propagation and noise
        # scaling rather than just rounding perturbations.
        intervals = np.resize(np.array([0.004, 0.017, 0.009, 0.013]), 2000)
        time = np.concatenate(([0.], np.cumsum(intervals)))
        samples, truth = self.recording(time, omega=[0.12, -0.08, 0.3])
        estimate, details = self.replay(samples, axes=6)
        # Initial absolute yaw is unobservable; remove its single constant
        # gauge offset, with no alignment of later samples.
        estimate = (truth[0] * estimate[0].inv()) * estimate
        error = np.rad2deg((truth.inv() * estimate).magnitude())
        self.assertLess(error.max(), 0.2)
        self.assertFalse(details["rest"].any())

    def test_unobservable_heading_does_not_destabilize_the_filter(self):
        # Turning steadily about the vertical is never rest, so the gyroscope
        # bias about that axis is never measured and the filter's uncertainty
        # about its heading keeps growing. The estimate must stay sound; left
        # unbounded, that uncertainty breaks the filter after about 100 s.
        time = np.arange(20001) / 100
        initial = Rotation.from_euler("xyz", [20, -15, 40], degrees=True)
        omega = initial.inv().apply([0., 0., np.deg2rad(20)])
        samples, truth = self.recording(time, omega=omega)
        for algorithm in FILTERS:
            with self.subTest(algorithm=algorithm):
                estimate, details = self.replay(samples, algorithm=algorithm)
                self.assertFalse(details["rest"].any())
                estimate = (truth[0] * estimate[0].inv()) * estimate
                error = np.rad2deg((truth.inv() * estimate).magnitude())
                self.assertLess(error.max(), 0.2)

    def test_accelerometer_dropout_coasts_on_the_gyroscope(self):
        time = np.arange(3001) / 100
        omega = np.array([0.12, -0.08, 0.3])
        samples, truth = self.recording(time, omega=omega)
        dropout = (time >= 5) & (time < 7)
        changed = samples.copy()
        # A zero reading is what free fall looks like: it says nothing about
        # which way is down, so the attitude must simply follow the gyroscope.
        changed[dropout, 4:7] = 0
        estimate, details = self.replay(changed)
        estimate = (truth[0] * estimate[0].inv()) * estimate
        error = np.rad2deg((truth.inv() * estimate).magnitude())
        self.assertLess(error[dropout].max(), 0.2)
        self.assertLess(error[time > 20].max(), 0.2)
        self.assertFalse(details["rest"].any())
        # Same tolerance as the rest check above: the dropout must not be
        # mistaken for a gyroscope bias.
        np.testing.assert_allclose(self.bias(details), 0, atol=5e-4, rtol=0)


if __name__ == "__main__":
    unittest.main()
