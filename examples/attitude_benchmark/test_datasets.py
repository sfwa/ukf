"""Offline fixtures for source formats, reference gaps and download integrity."""
import hashlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np
from scipy.io import savemat

from datasets import SASSARI_SENSORS, fetch, load_broad, load_sassari


class DatasetFormatTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def test_broad_keeps_sensor_timeline_through_optical_gaps(self):
        path = self.root / "fixture.hdf5"
        acc = np.arange(18, dtype=float).reshape(6, 3) + .125
        gyr = acc / 100
        mag = acc * 3
        quat = np.tile([2., 0., 0., 0.], (6, 1))
        quat[:2] = np.nan  # Leading gap is common in the real recordings.
        quat[4] = 0       # Invalid zero reference, never a sensor dropout.
        movement = np.array([False, False, True, True, False, True])
        with h5py.File(path, "w") as data:
            data.attrs["sampling_rate"] = 1 / .0035
            for name, value in (("imu_acc", acc), ("imu_gyr", gyr),
                                ("imu_mag", mag), ("opt_quat", quat),
                                ("movement", movement)):
                data[name] = value
        trial = load_broad(path, "fixture")
        self.assertEqual(trial.dataset, "broad")
        np.testing.assert_allclose(trial.time, np.arange(6) * .0035)
        np.testing.assert_array_equal(trial.acc, acc)
        np.testing.assert_array_equal(trial.gyr, gyr)
        np.testing.assert_array_equal(trial.mag, mag)
        np.testing.assert_array_equal(trial.valid, [False, False, True, True, False, True])
        np.testing.assert_array_equal(trial.movement, movement)
        np.testing.assert_array_equal(trial.quat_gt[trial.valid],
                                      np.tile([1., 0., 0., 0.], (3, 1)))
        self.assertTrue(np.isnan(trial.quat_gt[~trial.valid]).all())

    def sassari_fixture(self):
        # Distinct columns/units expose accidental axis swaps or deg/s conversion.
        n = 8
        data = {"Qs": np.tile([.5, .5, -.5, .5], (n, 1)),
                "indz": np.array([1, 3]), "indx": np.array([4]),
                "indy": np.array([6]), "indarb": np.array([8])}
        for i, sensor in enumerate(SASSARI_SENSORS):
            samples = np.arange(n * 14, dtype=float).reshape(n, 14) / 100 + i
            samples[:, 0] = np.arange(n)[::-1] + 1234  # Native packet counter.
            samples[:, 10:14] = np.nan  # Unused proprietary attitude must not leak.
            data[sensor] = samples
        return data

    def test_sassari_indices_columns_and_common_clock(self):
        data = self.sassari_fixture()
        path = self.root / "fixture.mat"
        savemat(path, data)
        trials = list(load_sassari(path, "slow"))
        self.assertEqual(len(trials), 6)
        for trial, sensor in zip(trials, SASSARI_SENSORS):
            with self.subTest(sensor=sensor):
                self.assertEqual(trial.name, f"slow_{sensor}")
                self.assertEqual(trial.dataset, "sassari")
                np.testing.assert_array_equal(trial.time, np.arange(8) / 100.)
                np.testing.assert_array_equal(trial.acc, data[sensor][:, 1:4])
                np.testing.assert_array_equal(trial.gyr, data[sensor][:, 4:7])
                np.testing.assert_array_equal(trial.mag, data[sensor][:, 7:10])
                np.testing.assert_array_equal(trial.quat_gt, data["Qs"])
                np.testing.assert_array_equal(trial.valid, np.ones(8, bool))
                np.testing.assert_array_equal(trial.movement,
                                              [True, False, True, True, False, True, False, True])

    def test_sassari_invalid_matlab_index_is_rejected(self):
        data = self.sassari_fixture()
        data["indz"] = np.array([0, 1])
        path = self.root / "fixture.mat"
        savemat(path, data)
        with self.assertRaisesRegex(ValueError, "movement indices"):
            list(load_sassari(path, "slow"))

    def test_sensor_failure_is_not_masked_as_missing_ground_truth(self):
        data = self.sassari_fixture()
        data["XS1"][3, 4] = np.nan
        data["Qs"][3] = np.nan
        path = self.root / "fixture.mat"
        savemat(path, data)
        with self.assertRaisesRegex(ValueError, "gyr measurements"):
            list(load_sassari(path, "slow", ("XS1",)))


class DownloadIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.payload = b"expected dataset bytes"
        self.entry = {"path": "broad/fixture.hdf5", "url": "https://invalid.example/fixture",
                      "bytes": len(self.payload),
                      "sha256": hashlib.sha256(self.payload).hexdigest()}

    def test_valid_cache_is_verified_and_used_without_network(self):
        path = self.root / self.entry["path"]
        path.parent.mkdir()
        path.write_bytes(self.payload)
        with patch("datasets.urlopen", side_effect=AssertionError("network accessed")):
            self.assertEqual(fetch(self.root, self.entry), path)

    def test_same_size_corruption_is_replaced_only_by_verified_download(self):
        path = self.root / self.entry["path"]
        path.parent.mkdir()
        path.write_bytes(b"x" * len(self.payload))
        with patch("datasets.urlopen", return_value=io.BytesIO(self.payload)) as request:
            self.assertEqual(fetch(self.root, self.entry), path)
            request.assert_called_once_with(self.entry["url"], timeout=60)
        self.assertEqual(path.read_bytes(), self.payload)
        self.assertEqual(list(path.parent.iterdir()), [path])

    def test_bad_download_fails_and_leaves_no_partial_cache(self):
        with patch("datasets.urlopen", side_effect=lambda *args, **kwargs:
                   io.BytesIO(b"x" * len(self.payload))) as request, patch("datasets.clock.sleep"):
            with self.assertRaisesRegex(RuntimeError, "SHA-256/size mismatch"):
                fetch(self.root, self.entry)
            self.assertEqual(request.call_count, 3)
        path = self.root / self.entry["path"]
        self.assertFalse(path.exists())
        self.assertEqual(list(path.parent.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
