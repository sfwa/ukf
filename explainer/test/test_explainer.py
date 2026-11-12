#!/usr/bin/env python3
"""Checks that the explainer shows what the benchmark's estimator does.

The page claims to run the benchmark's estimator on the benchmark's inputs.
These tests hold it to that:

  - the packed recordings decode to the very numbers the benchmark reads;
  - the explainer's C interface, built natively, writes the same attitude
    file as euroc_benchmark, byte for byte;
  - the WebAssembly module, driven by the page's own JavaScript, agrees
    with the native program to rounding error;
  - a filter step taken apart for the figures agrees with the real step;
  - the page's error figures agree with the benchmark's metrics.

Run from the repository root, after building (see explainer/README.md):

    python -m unittest discover -s explainer/test -v

A test is skipped if what it needs has not been built. Set
EXPLAINER_REQUIRE_ALL=1 to make that a failure instead (CI does).
"""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SITE = Path(os.environ.get("EXPLAINER_SITE", ROOT / "explainer" / "_site"))
BUILD = Path(os.environ.get("UKF_BUILD", ROOT / "build"))
CACHE = Path(os.environ.get("EUROC_CACHE", ROOT / ".benchmark-cache" / "euroc"))
BENCHMARK = BUILD / "examples" / "euroc_benchmark" / "euroc_benchmark"
NATIVE_CHECK = BUILD / "explainer" / "explainer_check"
FILTERS = ("mav-srukf", "mav-srukf-gravity")


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


datasets = load_module("euroc_datasets", ROOT / "examples" / "euroc_benchmark" / "datasets.py")
metrics = load_module("attitude_metrics", ROOT / "examples" / "attitude_benchmark" / "metrics.py")
pack_euroc = load_module("pack_euroc", ROOT / "explainer" / "tools" / "pack_euroc.py")


def need(test: unittest.TestCase, path: Path, what: str) -> None:
    if path.exists():
        return
    message = f"{what} not found at {path}"
    if os.environ.get("EXPLAINER_REQUIRE_ALL"):
        test.fail(message)
    test.skipTest(message)


def recordings(test: unittest.TestCase) -> list[dict]:
    need(test, SITE / "data" / "recordings.json", "The packed recordings")
    return json.loads((SITE / "data" / "recordings.json").read_text())["recordings"]


def write_sensor_csv(trial, path: Path) -> None:
    """The estimator input, written as the benchmark's own harness writes it."""
    values = np.column_stack((trial.time, trial.gyr, trial.acc, np.zeros_like(trial.acc)))
    np.savetxt(path, values, delimiter=",", comments="", fmt="%.17g", header="time,gx,gy,gz,ax,ay,az,mx,my,mz")


class PackTest(unittest.TestCase):
    def test_packs_hold_the_benchmark_inputs(self):
        for entry in recordings(self):
            with self.subTest(recording=entry["name"]):
                trial = datasets.load_sequence(CACHE / entry["name"], entry["name"])
                header, arrays = pack_euroc.unpack((SITE / "data" / entry["file"]).read_bytes())
                self.assertEqual(header["samples"], len(trial.time))
                # The IMU samples and their times are exact.
                np.testing.assert_array_equal(arrays["time_ns"][:, 0].astype(float) * 1e-9, trial.time)
                np.testing.assert_array_equal(pack_euroc.decode_gyroscope(arrays["gyroscope"]), trial.gyr)
                np.testing.assert_array_equal(pack_euroc.decode_accelerometer(arrays["accelerometer"]), trial.acc)
                # The reference attitude is rounded, by less than one part in a million.
                np.testing.assert_array_equal(arrays["reference_valid"][:, 0].astype(bool), trial.valid)
                attitude = arrays["reference_attitude"] * header["arrays"]["reference_attitude"]["scale"]
                sign = np.sign(np.sum(attitude * np.nan_to_num(trial.quat_gt), axis=1))[:, None]
                error = np.abs(attitude * sign - trial.quat_gt)[trial.valid]
                self.assertLess(error.max(), 3e-7)


class NativeInterfaceTest(unittest.TestCase):
    def test_interface_writes_the_benchmark_output(self):
        need(self, BENCHMARK, "euroc_benchmark")
        need(self, NATIVE_CHECK, "explainer_check")
        name = recordings(self)[0]["name"]
        trial = datasets.load_sequence(CACHE / name, name)
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            write_sensor_csv(trial, work / "sensors.csv")
            for flt in FILTERS:
                with self.subTest(filter=flt):
                    subprocess.run([BENCHMARK, "--input", work / "sensors.csv", "--output", work / "benchmark.csv",
                                    "--filter", flt, "--axes", "6"], check=True, capture_output=True)
                    result = subprocess.run([NATIVE_CHECK, "--input", work / "sensors.csv", "--output", work / "check.csv",
                                             "--filter", flt, "--explain-every", "193"],
                                            check=True, capture_output=True, text=True)
                    self.assertEqual((work / "check.csv").read_bytes(), (work / "benchmark.csv").read_bytes())
                    summary = json.loads(result.stdout)
                    self.assertGreater(summary["explained_steps"], 50)
                    self.assertLess(summary["worst_step_difference"], 1e-9)


class WasmTest(unittest.TestCase):
    def test_page_code_agrees_with_the_native_estimator(self):
        need(self, BENCHMARK, "euroc_benchmark")
        need(self, SITE / "ukf_explainer.wasm", "The WebAssembly module")
        if shutil.which("node") is None:
            need(self, Path("/nonexistent/node"), "Node.js")
        entry = recordings(self)[0]
        trial = datasets.load_sequence(CACHE / entry["name"], entry["name"])
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            subprocess.run(["node", ROOT / "explainer" / "test" / "run_wasm.mjs", SITE / "ukf_explainer.wasm",
                            SITE / "data" / entry["file"], work], check=True, capture_output=True)
            summary = json.loads((work / "summary.json").read_text())

            # The page's decoder gives the estimator the benchmark's inputs, bit for bit.
            imu = np.fromfile(work / "imu.f64").reshape(-1, 7)
            np.testing.assert_array_equal(imu[:, 0], trial.time)
            np.testing.assert_array_equal(imu[:, 1:4], trial.gyr)
            np.testing.assert_array_equal(imu[:, 4:7], trial.acc)

            write_sensor_csv(trial, work / "sensors.csv")
            for flt in FILTERS:
                with self.subTest(filter=flt):
                    subprocess.run([BENCHMARK, "--input", work / "sensors.csv", "--output", work / "benchmark.csv",
                                    "--filter", flt, "--axes", "6"], check=True, capture_output=True)
                    native = np.loadtxt(work / "benchmark.csv", delimiter=",", skiprows=1)[:, 1:]
                    wasm = np.fromfile(work / f"{flt}.f64").reshape(-1, 4)
                    page = summary["filters"][flt]
                    self.assertEqual(page["valid"], len(trial.time))
                    # Same estimator, different compiler and maths library: equal to rounding error.
                    self.assertLess(np.abs(wasm - native).max(), 1e-9)
                    # A step taken apart for the figures is the step the estimator took.
                    self.assertLess(page["worst_step_difference"], 1e-9)
                    # The page's error figures are the benchmark's.
                    expected, _ = metrics.compute_errors(native, trial.quat_gt, trial.valid, trial.time)
                    self.assertEqual(page["scored_samples"], expected["samples"])
                    for key in ("inclination", "heading", "total"):
                        self.assertAlmostEqual(page["rms_deg"][key], expected[f"{key}_rmse_deg"], delta=1e-4)


if __name__ == "__main__":
    unittest.main()
