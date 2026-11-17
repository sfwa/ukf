"""Pinned BROAD/Sassari downloads and conversion to common IMU trial arrays.

No reference orientation enters an estimator. Reference frames are z-up but their
horizontal directions differ: BROAD uses ENU and Sassari the optical lab frame.
The evaluator must explicitly specify any *constant* heading alignment.

BROAD is already synchronized and sensor/optical aligned by its authors. Sassari
v5 is synchronized at 100 Hz; its authors nominally aligned all sensor and optical
local axes by mounting them on the same board. We use Qs unchanged and do not fit
a mounting correction using evaluation data. Residual mounting errors therefore
remain part of the benchmark. The proprietary IMU orientation columns are unused.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import tempfile
import time as clock
from typing import Iterator, Sequence
from urllib.request import urlopen

import h5py
import numpy as np
from scipy.io import loadmat

MANIFEST_PATH = Path(__file__).with_name("datasets.lock.json")
SASSARI_SENSORS = ("XS1", "XS2", "AP1", "AP2", "SH1", "SH2")


@dataclass
class Trial:
    name: str
    dataset: str
    time: np.ndarray
    gyr: np.ndarray
    acc: np.ndarray
    mag: np.ndarray
    quat_gt: np.ndarray
    valid: np.ndarray
    movement: np.ndarray | None = None

    def validate(self) -> "Trial":
        """Reject malformed input instead of silently shortening/skipping a trial."""
        n = len(self.time)
        if n < 2 or self.time.shape != (n,) or not np.isfinite(self.time).all():
            raise ValueError(f"{self.name}: invalid timestamps")
        if not np.all(np.diff(self.time) > 0):
            raise ValueError(f"{self.name}: timestamps must increase strictly")
        for name in ("gyr", "acc", "mag"):
            data = getattr(self, name)
            if data.shape != (n, 3) or not np.isfinite(data).all():
                raise ValueError(f"{self.name}: invalid {name} measurements")
        if self.quat_gt.shape != (n, 4) or self.valid.shape != (n,):
            raise ValueError(f"{self.name}: invalid ground-truth shapes")
        if not self.valid.any():
            raise ValueError(f"{self.name}: no valid ground truth")
        if self.movement is not None and self.movement.shape != (n,):
            raise ValueError(f"{self.name}: invalid movement mask")
        return self


def manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def _matches(path: Path, entry: dict) -> bool:
    if not path.is_file() or path.stat().st_size != entry["bytes"]:
        return False
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest() == entry["sha256"]


def fetch(cache: Path, entry: dict) -> Path:
    """Verify cached content every time; downloads are atomic and hash checked."""
    path = Path(cache) / entry["path"]
    if _matches(path, entry):
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    error = None
    for attempt in range(3):
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as out:
                temporary = Path(out.name)
                with urlopen(entry["url"], timeout=60) as response:
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        out.write(block)
            if not _matches(temporary, entry):
                raise ValueError(f"SHA-256/size mismatch for {entry['url']}")
            temporary.replace(path)
            return path
        except (OSError, ValueError) as exc:
            error = exc
            if attempt < 2:
                clock.sleep(1 + attempt)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    raise RuntimeError(f"Could not fetch verified dataset {entry['url']}: {error}") from error


def _reference(quat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    quat = np.asarray(quat, dtype=float).copy()
    norm = np.linalg.norm(quat, axis=1)
    valid = np.isfinite(quat).all(axis=1) & (norm > 1e-8)
    quat[valid] /= norm[valid, None]
    # Keep invalid reference samples explicit; they are excluded only in scoring.
    quat[~valid] = np.nan
    return quat, valid


def load_broad(path: Path, name: str) -> Trial:
    with h5py.File(path, "r") as data:
        acc = np.asarray(data["imu_acc"], dtype=float)
        gyr = np.asarray(data["imu_gyr"], dtype=float)
        mag = np.asarray(data["imu_mag"], dtype=float)
        quat, valid = _reference(np.asarray(data["opt_quat"]))
        movement = np.asarray(data["movement"], dtype=bool)
        sampling_rate = float(data.attrs["sampling_rate"])
    if not np.isfinite(sampling_rate) or sampling_rate <= 0:
        raise ValueError(f"{name}: invalid BROAD sampling rate")
    return Trial(name, "broad", np.arange(len(acc)) / sampling_rate,
                 gyr, acc, mag, quat, valid, movement).validate()


def load_sassari(path: Path, speed: str, sensors: Sequence[str] = SASSARI_SENSORS) -> Iterator[Trial]:
    data = loadmat(path, squeeze_me=True)
    quat, valid = _reference(data["Qs"])
    n = len(quat)
    movement = np.zeros(n, dtype=bool)
    # Source indices are MATLAB one-based sample indices, not start/end pairs.
    for key in ("indz", "indx", "indy", "indarb"):
        indices = np.atleast_1d(data[key]).astype(np.int64) - 1
        if np.any(indices < 0) or np.any(indices >= n):
            raise ValueError(f"{speed}: invalid Sassari movement indices")
        movement[indices] = True
    for sensor in sensors:
        samples = np.asarray(data[sensor], dtype=float)
        if samples.shape != (n, 14):
            raise ValueError(f"{speed}/{sensor}: unexpected Sassari shape {samples.shape}")
        # Column 0 has incompatible native timestamps/packet counters. The
        # published rows are already synchronized/resampled at exactly 100 Hz.
        yield Trial(f"{speed}_{sensor}", "sassari", np.arange(n) / 100.0,
                    samples[:, 4:7], samples[:, 1:4], samples[:, 7:10],
                    quat.copy(), valid.copy(), movement.copy()).validate()


def load_trials(cache: Path, datasets: Sequence[str] = ("broad", "sassari"),
                smoke: bool = False) -> Iterator[Trial]:
    """Load all 39 BROAD + 18 Sassari trials, or one complete trial per dataset.

    Smoke uses BROAD 01 and Sassari slow/XS1, without truncating their recordings.
    Full runs retain every hardware unit and speed; failures are never skipped.
    """
    unknown = set(datasets) - {"broad", "sassari"}
    if unknown:
        raise ValueError(f"Unknown datasets: {sorted(unknown)}")
    seen = set()
    for entry in manifest()["files"]:
        dataset = entry["dataset"]
        if dataset not in datasets or (smoke and dataset in seen):
            continue
        path = fetch(Path(cache), entry)
        seen.add(dataset)
        if dataset == "broad":
            yield load_broad(path, entry["name"])
        else:
            yield from load_sassari(path, entry["name"], ("XS1",) if smoke else SASSARI_SENSORS)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path(".benchmark-cache"))
    parser.add_argument("--datasets", choices=("broad", "sassari"), nargs="+", default=["broad", "sassari"])
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    for trial in load_trials(args.cache, args.datasets, args.smoke):
        print(f"{trial.dataset}/{trial.name}: {len(trial.time)} samples, "
              f"{trial.time[-1]:.2f} s, {np.count_nonzero(trial.valid)} valid references")


if __name__ == "__main__":
    main()
