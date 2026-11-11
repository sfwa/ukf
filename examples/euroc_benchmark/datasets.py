"""Verified EuRoC IMU files and evaluation-only reference interpolation.

EuRoC publishes angular velocity and specific force at nominally 200 Hz, with
no recorded magnetometer or rotor commands. Its state-groundtruth files are
offline batch estimates: Leica positions + IMU in MH, Vicon poses + IMU in V1/V2.
Neither those states nor their estimated sensor biases enter estimator inputs.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import re
import tempfile
from typing import Iterator, Sequence
from urllib.request import urlopen

import numpy as np

SEQUENCES = (
    "MH_01_easy", "MH_02_easy", "MH_03_medium", "MH_04_difficult", "MH_05_difficult",
    "V1_01_easy", "V1_02_medium", "V1_03_difficult",
    "V2_01_easy", "V2_02_medium", "V2_03_difficult",
)
DEVELOPMENT = ("MH_01_easy", "V1_02_medium")
MAX_REFERENCE_GAP_NS = 20_000_000
MANIFEST_PATH = Path(__file__).with_name("datasets.lock.json")


@dataclass
class Trial:
    name: str
    time: np.ndarray
    gyr: np.ndarray
    acc: np.ndarray
    quat_gt: np.ndarray
    valid: np.ndarray
    timestamp_ns: np.ndarray
    metadata: dict = field(default_factory=dict)
    dataset: str = "euroc"

    @property
    def mag(self) -> np.ndarray:
        # Compatibility placeholder only; EuRoC cannot support a nine-axis run.
        return np.zeros_like(self.acc)

    def validate(self) -> "Trial":
        n = len(self.time)
        if n < 2 or self.time.shape != (n,) or not np.isfinite(self.time).all():
            raise ValueError(f"{self.name}: invalid timestamps")
        if self.timestamp_ns.dtype != np.int64 or self.timestamp_ns.shape != (n,):
            raise ValueError(f"{self.name}: timestamps must retain integer nanoseconds")
        if not np.all(np.diff(self.timestamp_ns) > 0) or not np.all(np.diff(self.time) > 0):
            raise ValueError(f"{self.name}: timestamps must increase strictly")
        for key in ("gyr", "acc"):
            data = getattr(self, key)
            if data.shape != (n, 3) or not np.isfinite(data).all():
                raise ValueError(f"{self.name}: invalid {key}")
        if self.quat_gt.shape != (n, 4) or self.valid.shape != (n,) or not self.valid.any():
            raise ValueError(f"{self.name}: missing valid reference orientations")
        if not np.isfinite(self.quat_gt[self.valid]).all():
            raise ValueError(f"{self.name}: non-finite valid reference")
        return self


def manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text())


def _matches(path: Path, entry: dict) -> bool:
    return (path.is_file() and path.stat().st_size == entry["bytes"]
            and hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"])


def fetch(cache: Path, entry: dict) -> Path:
    """Hash-check cached files and atomically replace only verified downloads."""
    relative = Path(entry["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Dataset path must remain inside the cache")
    path = Path(cache) / relative
    if _matches(path, entry):
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for _ in range(3):
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
                temporary = Path(stream.name)
                with urlopen(entry["url"], timeout=60) as response:
                    remaining = entry["bytes"] + 1
                    while remaining:
                        block = response.read(min(1024 * 1024, remaining))
                        if not block:
                            break
                        stream.write(block)
                        remaining -= len(block)
            if not _matches(temporary, entry):
                raise ValueError(f"SHA-256/size mismatch for {entry['url']}")
            temporary.replace(path)
            return path
        except (OSError, ValueError) as error:
            last_error = error
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    raise RuntimeError(f"Could not fetch verified dataset {entry['url']}: {last_error}") from last_error


def read_csv(path: Path, columns: int) -> tuple[np.ndarray, np.ndarray]:
    # Read the first column as decimal integers, never as epoch-sized doubles.
    values = np.loadtxt(path, delimiter=",", skiprows=1, dtype=str, ndmin=2)
    if values.shape[1] != columns or len(values) < 2:
        raise ValueError(f"{path}: expected at least two rows with {columns} columns")
    timestamps = values[:, 0].astype(np.int64)
    if not np.all(np.diff(timestamps) > 0):
        raise ValueError(f"{path}: timestamps must increase strictly")
    return timestamps, values[:, 1:].astype(float)


def interpolate_reference(imu_ns: np.ndarray, reference_ns: np.ndarray,
                          quaternion: np.ndarray,
                          max_gap_ns: int = MAX_REFERENCE_GAP_NS) -> tuple[np.ndarray, np.ndarray]:
    """Shortest-arc SLERP on bounded intervals, without extrapolation or gap filling."""
    imu_ns = np.asarray(imu_ns, dtype=np.int64)
    reference_ns = np.asarray(reference_ns, dtype=np.int64)
    quaternion = np.asarray(quaternion, dtype=float)
    if (len(reference_ns) < 2 or quaternion.shape != (len(reference_ns), 4)
            or not np.all(np.diff(reference_ns) > 0) or max_gap_ns <= 0):
        raise ValueError("Malformed reference timestamps/quaternions")
    norms = np.linalg.norm(quaternion, axis=1)
    good = np.isfinite(quaternion).all(axis=1) & (norms > 1e-8)
    q = quaternion.copy()
    q[good] /= norms[good, None]
    output = np.full((len(imu_ns), 4), np.nan)
    valid = np.zeros(len(imu_ns), dtype=bool)
    lower = np.searchsorted(reference_ns, imu_ns, side="right") - 1
    candidate = np.flatnonzero((lower >= 0) & (lower < len(reference_ns)))
    exact = candidate[imu_ns[candidate] == reference_ns[lower[candidate]]]
    exact = exact[good[lower[exact]]]
    output[exact] = q[lower[exact]]
    valid[exact] = True
    between = np.flatnonzero((lower >= 0) & (lower + 1 < len(reference_ns)) & ~valid)
    left = lower[between]
    gaps = reference_ns[left + 1] - reference_ns[left]
    keep = good[left] & good[left + 1] & (gaps <= max_gap_ns)
    between, left, gaps = between[keep], left[keep], gaps[keep]
    if len(between):
        alpha = ((imu_ns[between] - reference_ns[left]).astype(float) / gaps)[:, None]
        q0, q1 = q[left], q[left + 1].copy()
        dot = np.sum(q0 * q1, axis=1)
        q1[dot < 0] *= -1
        dot = np.clip(np.abs(dot), 0, 1)
        angle = np.arccos(dot)
        sine = np.sin(angle)
        result = (1 - alpha) * q0 + alpha * q1
        curved = sine > 1e-8
        a = alpha[curved]
        theta = angle[curved, None]
        result[curved] = (np.sin((1 - a) * theta) * q0[curved]
                          + np.sin(a * theta) * q1[curved]) / sine[curved, None]
        result /= np.linalg.norm(result, axis=1)[:, None]
        output[between] = result
        valid[between] = True
    return output, valid


def _require_identity_frame(path: Path) -> None:
    # These locked calibration files use a simple ASL matrix literal. Avoid a
    # YAML dependency and fail instead of silently accepting another rig frame.
    match = re.search(r"T_BS:.*?data:\s*\[([^]]+)\]", path.read_text(), re.S)
    if match is None:
        raise ValueError(f"{path}: missing T_BS calibration")
    matrix = np.fromstring(match[1], sep=",")
    if matrix.size != 16 or not np.allclose(matrix.reshape(4, 4), np.eye(4), atol=1e-12, rtol=0):
        raise ValueError(f"{path}: this benchmark requires the published IMU-body reference frame")


def load_sequence(directory: Path, name: str) -> Trial:
    mav = Path(directory) / "mav0"
    _require_identity_frame(mav / "imu0/sensor.yaml")
    _require_identity_frame(mav / "state_groundtruth_estimate0/sensor.yaml")
    imu_ns, imu = read_csv(mav / "imu0/data.csv", 7)
    gt_ns, gt = read_csv(mav / "state_groundtruth_estimate0/data.csv", 17)
    quat, valid = interpolate_reference(imu_ns, gt_ns, gt[:, 3:7])
    is_machine_hall = name.startswith("MH_")
    warning = ("Machine Hall attitude is part of a Leica-position/IMU batch solution; "
               "the source calibration states that orientation is observed only by the IMU."
               if is_machine_hall else
               "Reference is an offline Vicon-pose/IMU batch solution and shares IMU data with the estimators.")
    if name == "V1_01_easy":
        warning += (" OpenVINS reports inaccurate orientation in the canonical V1_01_easy reference; "
                    "this benchmark retains and flags the original reference instead of substituting a revised trajectory.")
    metadata = {
        "frame": "Published IMU/body frame S; scalar-first q_RS rotates S into z-up world R; both T_BS are identity.",
        "ground_truth_source": "state_groundtruth_estimate0/data.csv (canonical original)",
        "ground_truth_warning": warning,
        "ground_truth_independent_attitude": False,
        "external_attitude_observation": not is_machine_hall,
        "ground_truth_quality_flag": "known_orientation_issue" if name == "V1_01_easy" else "fused_reference",
        "reference_interpolation": "Shortest-arc SLERP; integer-nanosecond brackets; no extrapolation.",
        "max_reference_gap_seconds": MAX_REFERENCE_GAP_NS / 1e9,
        "imu_timestamp_origin_ns": int(imu_ns[0]),
        "nominal_sample_rate_hz": 200,
        "sensor_bias_columns_used": False,
        "magnetometer_available": False,
        "rotor_commands_available": False,
    }
    # Subtract before casting: epoch-sized float conversion loses clock precision.
    return Trial(name, (imu_ns - imu_ns[0]).astype(float) * 1e-9,
                 imu[:, :3], imu[:, 3:6], quat, valid, imu_ns, metadata).validate()


def load_trials(cache: Path, sequences: Sequence[str] | None = None) -> Iterator[Trial]:
    selected = tuple(SEQUENCES if sequences is None else sequences)
    if len(set(selected)) != len(selected) or set(selected) - set(SEQUENCES):
        raise ValueError("Select distinct canonical EuRoC sequence names")
    files = manifest()["files"]
    for name in selected:
        entries = [entry for entry in files if entry["sequence"] == name]
        if not entries:
            raise ValueError(f"Manifest has no files for {name}")
        for entry in entries:
            fetch(Path(cache), entry)
        yield load_sequence(Path(cache) / name, name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path(".benchmark-cache/euroc"))
    parser.add_argument("--sequences", choices=SEQUENCES, nargs="+")
    args = parser.parse_args()
    for trial in load_trials(args.cache, args.sequences):
        print(f"{trial.name}: {len(trial.time)} IMU samples, {trial.time[-1]:.3f} s, "
              f"{np.count_nonzero(trial.valid)} reference samples")
