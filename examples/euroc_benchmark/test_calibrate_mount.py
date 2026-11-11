"""Check that development mounting calibration cannot read evaluation data."""
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

import calibrate_mount


class MountCalibrationTests(unittest.TestCase):
    def fixture(self):
        names = (*calibrate_mount.DEVELOPMENT, "V2_03_difficult")
        files = [{"path": f"{name}/mav0/{sensor}/data.csv", "sha256": "fixture"}
                 for name in names for sensor in ("imu0", "state_groundtruth_estimate0")]
        sample = np.arange(2000)
        time = sample * .005
        axis = np.array([np.cos(np.deg2rad(20)), 0., -np.sin(np.deg2rad(20))])
        force = (9.80665 + 1.5 * np.sin(2*np.pi*60*time))[:, None] * axis
        imu = np.column_stack((np.zeros_like(force), force))
        timestamps = 1_403_636_579_758_555_392 + sample.astype(np.int64) * 5_000_000
        return files, timestamps, imu

    def test_only_development_imu_files_are_fetched_and_read(self):
        files, timestamps, imu = self.fixture()
        fetched, read = [], []

        def fetch(cache, entry):
            path = Path(entry["path"])
            self.assertIn(path.parts[0], calibrate_mount.DEVELOPMENT)
            self.assertEqual(path.parts[-2], "imu0")
            fetched.append(str(path))
            return cache / path

        def read_csv(path, columns):
            self.assertEqual(path.parts[-2], "imu0")
            self.assertEqual(columns, 7)
            read.append(str(path))
            return timestamps, imu

        with patch.object(calibrate_mount, "manifest", return_value={"files": files}), \
                patch.object(calibrate_mount, "fetch", side_effect=fetch), \
                patch.object(calibrate_mount, "read_csv", side_effect=read_csv):
            result = calibrate_mount.calibrate(Path("cache"))
        self.assertEqual(len(fetched), 2)
        self.assertEqual(len(read), 2)
        self.assertEqual(result["nominal_thrust_tilt_degrees"], 20)
        self.assertAlmostEqual(result["equal_sequence_mean_pitch_degrees"], 20, places=10)
        self.assertEqual({row["sequence"] for row in result["per_sequence"]}, set(calibrate_mount.DEVELOPMENT))

    def test_real_sampling_gap_cannot_be_treated_as_clock_quantization(self):
        files, timestamps, imu = self.fixture()
        timestamps[100:] += 20_000_000
        with patch.object(calibrate_mount, "manifest", return_value={"files": files}), \
                patch.object(calibrate_mount, "fetch", return_value=Path("unused.csv")), \
                patch.object(calibrate_mount, "read_csv", return_value=(timestamps, imu)):
            with self.assertRaisesRegex(ValueError, "deviate from the nominal"):
                calibrate_mount.calibrate(Path("cache"))


if __name__ == "__main__":
    unittest.main()
