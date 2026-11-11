"""EuRoC clock precision, reference boundaries and download-integrity checks."""
import hashlib
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

SPEC = importlib.util.spec_from_file_location("euroc_test_datasets", Path(__file__).with_name("datasets.py"))
datasets = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = datasets
SPEC.loader.exec_module(datasets)


class InterpolationTests(unittest.TestCase):
    def test_shortest_arc_sign_invariance_and_no_extrapolation(self):
        epoch = 1_403_636_579_758_555_392
        times = epoch + np.array([0, 10_000_000], dtype=np.int64)
        query = epoch + np.array([-1, 0, 5_000_000, 10_000_000, 10_000_001], dtype=np.int64)
        end = np.array([np.sqrt(.5), 0, 0, np.sqrt(.5)])
        quaternions = np.array([[1., 0, 0, 0], -end])
        result, valid = datasets.interpolate_reference(query, times, quaternions)
        np.testing.assert_array_equal(valid, [False, True, True, True, False])
        expected = np.array([np.cos(np.pi / 8), 0, 0, np.sin(np.pi / 8)])
        np.testing.assert_allclose(result[2], expected, atol=1e-14)
        self.assertTrue(np.isnan(result[~valid]).all())

    def test_large_reference_gaps_are_not_bridged(self):
        times = np.array([0, 5_000_000, 100_000_000, 105_000_000], dtype=np.int64)
        query = np.array([2_000_000, 5_000_000, 50_000_000, 100_000_000, 103_000_000], dtype=np.int64)
        _, valid = datasets.interpolate_reference(query, times, np.tile([1., 0, 0, 0], (4, 1)))
        np.testing.assert_array_equal(valid, [True, True, False, True, True])

    def test_invalid_reference_rows_do_not_disappear_into_interpolation(self):
        times = np.array([0, 5_000_000, 10_000_000], dtype=np.int64)
        reference = np.array([[1., 0, 0, 0], [np.nan] * 4, [1., 0, 0, 0]])
        _, valid = datasets.interpolate_reference(np.array([0, 2_000_000, 5_000_000, 8_000_000, 10_000_000]),
                                                  times, reference)
        np.testing.assert_array_equal(valid, [True, False, False, False, True])


class LoaderTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def test_epoch_nanoseconds_are_not_rounded_through_floats(self):
        epoch = 1_403_636_579_758_555_392
        path = self.directory / "tiny.csv"
        path.write_text(f"#timestamp,value\n{epoch},1\n{epoch + 1},2\n{epoch + 7},3\n")
        timestamps, _ = datasets.read_csv(path, 2)
        np.testing.assert_array_equal(timestamps - epoch, [0, 1, 7])

    def fixture(self):
        epoch = 1_403_636_579_758_555_392
        imu = self.directory / "mav0/imu0"
        reference = self.directory / "mav0/state_groundtruth_estimate0"
        imu.mkdir(parents=True)
        reference.mkdir(parents=True)
        calibration = "T_BS:\n  data: [" + ",".join(map(str, np.eye(4).ravel())) + "]\n"
        (imu / "sensor.yaml").write_text(calibration)
        (reference / "sensor.yaml").write_text(calibration)
        (imu / "data.csv").write_text("#timestamp,gx,gy,gz,ax,ay,az\n" + "".join(
            f"{epoch + i * 5_000_000},0.1,0.2,0.3,0,0,9.80665\n" for i in range(3)))
        # Deliberately extreme reference bias estimates must not alter sensors.
        (reference / "data.csv").write_text("#timestamp,reference\n" + "".join(
            f"{epoch + i * 5_000_000},0,0,0,1,0,0,0,0,0,0,1000,1000,1000,1000,1000,1000\n"
            for i in (1, 2)))
        return imu, reference

    def test_reference_support_does_not_crop_imu_or_calibrate_sensor_bias(self):
        self.fixture()
        trial = datasets.load_sequence(self.directory, "V1_02_medium")
        self.assertEqual(len(trial.time), 3)
        np.testing.assert_allclose(trial.time, [0, .005, .01], atol=1e-16)
        np.testing.assert_allclose(trial.gyr, np.tile([.1, .2, .3], (3, 1)), atol=0)
        np.testing.assert_allclose(trial.acc, np.tile([0, 0, 9.80665], (3, 1)), atol=0)
        np.testing.assert_array_equal(trial.valid, [False, True, True])
        self.assertFalse(trial.metadata["sensor_bias_columns_used"])

    def test_unexpected_sensor_frame_is_rejected(self):
        imu, _ = self.fixture()
        (imu / "sensor.yaml").write_text("T_BS:\n data: [0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1]\n")
        with self.assertRaisesRegex(ValueError, "IMU-body reference frame"):
            datasets.load_sequence(self.directory, "V1_02_medium")

    def test_cached_corruption_is_replaced_only_by_hash_verified_content(self):
        source = self.directory / "source.csv"
        source.write_bytes(b"verified bytes\n")
        entry = dict(path="cache.csv", bytes=source.stat().st_size,
                     sha256=hashlib.sha256(source.read_bytes()).hexdigest(), url=source.as_uri())
        cache = self.directory / "cache"
        cache.mkdir()
        target = cache / entry["path"]
        target.write_bytes(b"corrupt")
        self.assertEqual(datasets.fetch(cache, entry).read_bytes(), source.read_bytes())
        source.write_bytes(b"bad replacement")
        target.write_bytes(b"preserve on failed download")
        with self.assertRaisesRegex(RuntimeError, "SHA-256/size mismatch"):
            datasets.fetch(cache, entry)
        self.assertEqual(target.read_bytes(), b"preserve on failed download")

    def test_sequence_selection_never_fetches_unselected_files(self):
        entries = [dict(sequence=name, path=name + "/example") for name in datasets.SEQUENCES]
        with patch.object(datasets, "manifest", return_value={"files": entries}), \
                patch.object(datasets, "fetch") as fetch, \
                patch.object(datasets, "load_sequence", side_effect=lambda path, name: name):
            selected = list(datasets.load_trials(self.directory, datasets.DEVELOPMENT))
        self.assertEqual(selected, list(datasets.DEVELOPMENT))
        self.assertEqual([call.args[1]["sequence"] for call in fetch.call_args_list], list(datasets.DEVELOPMENT))


if __name__ == "__main__":
    unittest.main()
