"""Geometry and causal-adapter regressions; run with unittest discovery."""
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

from metrics import compute_errors, multiply
from references import REFERENCE_NAMES, initial_attitude, run_reference


def quat(rotvec):
    return Rotation.from_rotvec(rotvec).as_quat(scalar_first=True)


class MetricsTests(unittest.TestCase):
    def setUp(self):
        self.time = np.arange(1001)*0.01
        self.truth = np.tile([1., 0., 0., 0.], (len(self.time), 1))
        self.valid = np.ones(len(self.time), bool)

    def test_quaternion_sign_does_not_change_errors(self):
        estimate = self.truth.copy()
        estimate[::2] *= -1
        summary, _ = compute_errors(estimate, self.truth, self.valid, self.time)
        self.assertLess(summary['total_rmse_deg'], 1e-5)
        self.assertLess(summary['heading_rmse_deg'], 1e-5)

    def test_fixed_heading_alignment_preserves_later_drift(self):
        angle = np.full(len(self.time), np.deg2rad(40.))
        angle[self.time >= 5.] += np.deg2rad(20.)
        estimate = quat(np.column_stack((np.zeros(len(angle)), np.zeros(len(angle)), angle)))
        summary, _ = compute_errors(estimate, self.truth, self.valid, self.time)
        self.assertAlmostEqual(summary['heading_alignment_deg'], -40., places=7)
        self.assertAlmostEqual(summary['heading_rmse_deg'], 20., places=7)
        self.assertAlmostEqual(summary['total_rmse_deg'], 20., places=7)
        self.assertAlmostEqual(summary['inclination_rmse_deg'], 0., places=5)
        self.assertAlmostEqual(summary['raw_heading_rmse_deg'], 60., places=7)

    def test_tilt_cannot_be_aligned_away(self):
        estimate = np.tile(quat([np.deg2rad(30.), 0., 0.]), (len(self.time), 1))
        summary, _ = compute_errors(estimate, self.truth, self.valid, self.time)
        self.assertAlmostEqual(summary['inclination_rmse_deg'], 30., places=7)
        self.assertAlmostEqual(summary['total_rmse_deg'], 30., places=7)
        self.assertAlmostEqual(summary['heading_rmse_deg'], 0., places=7)

    def test_ground_truth_gaps_are_excluded(self):
        truth = self.truth.copy()
        valid = self.valid.copy()
        valid[600:700] = False
        truth[~valid] = np.nan
        summary, _ = compute_errors(self.truth, truth, valid, self.time)
        self.assertEqual(summary['samples'], 401)
        self.assertLess(summary['total_rmse_deg'], 1e-5)

    def test_failed_estimate_is_not_silently_masked(self):
        estimate = self.truth.copy()
        estimate[700] = np.nan
        with self.assertRaisesRegex(ValueError, 'estimate'):
            compute_errors(estimate, self.truth, self.valid, self.time)

    def test_missing_reference_fails(self):
        valid = np.zeros_like(self.valid)
        with self.assertRaisesRegex(ValueError, 'No valid ground truth'):
            compute_errors(self.truth, self.truth, valid, self.time)

    def test_leading_optical_gap_uses_first_valid_window(self):
        valid = self.time >= 6.
        truth = self.truth.copy()
        truth[~valid] = np.nan
        angle = np.full(len(self.time), np.deg2rad(40.))
        angle[self.time >= 7.] += np.deg2rad(20.)
        estimate = quat(np.column_stack((angle*0, angle*0, angle)))
        summary, errors = compute_errors(estimate, truth, valid, self.time)
        self.assertEqual(summary['alignment_start_seconds'], 6.)
        self.assertAlmostEqual(summary['heading_alignment_deg'], -40., places=7)
        self.assertAlmostEqual(summary['heading_rmse_deg'], 20., places=7)
        self.assertEqual(summary['samples'], 301)
        self.assertFalse(np.any(errors['score_mask'][self.time < 7.]))

    def test_movement_mask_does_not_remove_initial_alignment(self):
        movement = self.time >= 7.
        estimate = np.tile(quat([0., 0., .5]), (len(self.time), 1))
        summary, _ = compute_errors(estimate, self.truth, self.valid, self.time, score_mask=movement)
        self.assertEqual(summary["samples"], 301)
        self.assertLess(summary["total_rmse_deg"], 1e-5)

    def test_world_yaw_gauge_invariance(self):
        truth = quat(np.column_stack((.2*np.sin(self.time), .3*np.cos(self.time), .1*self.time)))
        estimate = multiply(quat([0., 0., .7]), truth)
        summary, _ = compute_errors(estimate, truth, self.valid, self.time)
        self.assertLess(summary['inclination_rmse_deg'], 1e-5)
        self.assertLess(summary['total_rmse_deg'], 1e-5)


class ReferenceTests(unittest.TestCase):
    def data(self, duration=2.0):
        time = np.arange(int(duration/.01)+1)*.01
        # Steady spin, tilted sensor, physically consistent fields.
        world_yaw = Rotation.from_rotvec(np.column_stack((time*0, time*0, .4*time)))
        rotation = world_yaw * Rotation.from_euler('xy', [20., -30.], degrees=True)
        acc = rotation.inv().apply(np.tile([0., 0., 9.80665], (len(time), 1)))
        mag = rotation.inv().apply(np.tile([.8, 0., -.6], (len(time), 1)))
        gyr = rotation.inv().apply(np.tile([0., 0., .4], (len(time), 1)))
        return time, gyr, acc, mag, rotation.as_quat(scalar_first=True)

    def test_sensor_initialization(self):
        _, _, acc, mag, truth = self.data()
        q = initial_attitude(acc[0], mag[0])
        self.assertAlmostEqual(abs(np.dot(q, truth[0])), 1., places=12)

    def test_nine_axis_reference_frame_and_propagation(self):
        time, gyr, acc, mag, truth = self.data()
        for name in REFERENCE_NAMES:
            with self.subTest(name=name):
                q = run_reference(name, time, gyr, acc, mag, 9)
                error = 2*np.rad2deg(np.arccos(np.clip(abs(np.sum(q*truth, axis=1)), 0, 1)))
                self.assertLess(float(np.max(error)), 1.0)

    def test_six_axis_never_reads_magnetometer_values(self):
        time, gyr, acc, mag, _ = self.data(.3)
        for name in REFERENCE_NAMES:
            with self.subTest(name=name):
                a = run_reference(name, time, gyr, acc, mag, 6)
                b = run_reference(name, time, gyr, acc, mag*7+2, 6)
                np.testing.assert_array_equal(a, b)

    def test_reference_outputs_are_causal(self):
        time, gyr, acc, mag, _ = self.data(.3)
        for name in REFERENCE_NAMES:
            for axes in (6, 9):
                with self.subTest(name=name, axes=axes):
                    full = run_reference(name, time, gyr, acc, mag, axes)
                    prefix = run_reference(name, time[:15], gyr[:15], acc[:15], mag[:15], axes)
                    np.testing.assert_allclose(full[:15], prefix, atol=1e-14)

    def test_vqf_rejects_irregular_timestamps(self):
        time, gyr, acc, mag, _ = self.data(.3)
        time[-1] += .001
        with self.assertRaisesRegex(ValueError, 'uniformly'):
            run_reference('vqf', time, gyr, acc, mag, 6)


if __name__ == '__main__':
    unittest.main()
