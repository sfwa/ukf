#!/usr/bin/env python3
"""Pack EuRoC recordings into small binary files for the explainer page.

The page runs the benchmark's estimator in the browser, so it needs the same
IMU samples as the benchmark, plus the reference trajectory to compare with.
This tool reads both through the benchmark's own loader (verified against
examples/euroc_benchmark/datasets.lock.json) and writes one ".pack" file per
recording.

IMU samples
    Stored as the sensor's integer counts. The published CSV values are
    exact multiples of one count, and these rules give back the very same
    doubles, bit for bit:

        gyroscope      rad/s  = count / 25 / 180 * pi
        accelerometer  m/s^2  = count / 1200 * 9.80665
        time           s      = (nanoseconds since the first sample) * 1e-9

    This tool checks that for every value, and stops if it is ever not true.

Reference
    The attitude is the benchmark's own interpolation of the reference at
    each IMU sample time, so the page can score a run as the benchmark does.
    It is rounded to 2^-22 (an angle of about 0.00003 degrees). Position,
    velocity and the reference's sensor biases are only drawn and plotted,
    so they are kept at every fourth sample and rounded more coarsely.

File layout
    Before compression: "UKFX", uint32 version, uint32 header length, a JSON
    header, then the arrays. Integers are little-endian. The header names
    every array and says how to decode it, so the page never hard-codes an
    offset. The whole file is then gzip-compressed.

    Each array holds "count" rows of "components" integers. To make the
    arrays compress well they are stored as follows, and decoded in the
    reverse order:

      1. "delta" times over, each row is replaced by its difference from
         the row before (the first row keeps its value).
      2. The array is written one component after another.
      3. The bytes are written one byte position after another: first the
         lowest byte of every integer, then the next byte, and so on.

    A row of an array with "every": N belongs to IMU sample N * row.

The dataset is not part of this repository and keeps its own terms; see the
header's "source" entry and explainer/README.md.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import struct
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
EUROC = ROOT / "examples" / "euroc_benchmark"
FORMAT_VERSION = 1
MAGIC = b"UKFX"

# The recordings on the page by default. The first is loaded when the page opens.
DEFAULT_SEQUENCES = ("V1_02_medium", "MH_01_easy", "MH_03_medium", "V1_03_difficult", "V2_01_easy")

ATTITUDE_SCALE = 2.0 ** -22
REFERENCE_EVERY = 4
POSITION_SCALE = 1e-4
VELOCITY_SCALE = 1e-4
BIAS_SCALE = 1e-6


def load_datasets_module():
    """Use the benchmark's loader, so that both read the same verified files."""
    spec = importlib.util.spec_from_file_location("euroc_datasets", EUROC / "datasets.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["euroc_datasets"] = module
    spec.loader.exec_module(module)
    return module


def decode_gyroscope(count: np.ndarray) -> np.ndarray:
    return count.astype(float) / 25 / 180 * math.pi


def decode_accelerometer(count: np.ndarray) -> np.ndarray:
    return count.astype(float) / 1200 * 9.80665


def to_counts(values: np.ndarray, decode, name: str) -> np.ndarray:
    """Recover integer counts, and prove that decoding them is exact."""
    one = decode(np.array([1]))[0]
    counts = np.rint(values / one)
    if np.abs(counts).max() > 32767:
        raise ValueError(f"{name}: counts do not fit in 16 bits")
    counts = counts.astype(np.int64)
    if not np.array_equal(decode(counts), values):
        raise ValueError(f"{name}: values are not exact multiples of one count")
    return counts


def interpolate_vectors(imu_ns, reference_ns, values, max_gap_ns):
    """Linear interpolation in the brackets used by the benchmark's attitude
    interpolation: no extrapolation, and no filling of gaps over max_gap_ns."""
    lower = np.searchsorted(reference_ns, imu_ns, side="right") - 1
    output = np.full((len(imu_ns), values.shape[1]), np.nan)
    inside = (lower >= 0) & (lower < len(reference_ns))
    exact = inside & (imu_ns == reference_ns[np.clip(lower, 0, len(reference_ns) - 1)])
    output[exact] = values[lower[exact]]
    between = np.flatnonzero((lower >= 0) & (lower + 1 < len(reference_ns)) & ~exact)
    left = lower[between]
    gaps = reference_ns[left + 1] - reference_ns[left]
    keep = gaps <= max_gap_ns
    between, left, gaps = between[keep], left[keep], gaps[keep]
    alpha = ((imu_ns[between] - reference_ns[left]).astype(float) / gaps)[:, None]
    output[between] = (1 - alpha) * values[left] + alpha * values[left + 1]
    return output


def hold_invalid(values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Give samples without a reference the nearest valid value, so that the
    stored differences stay small. The "valid" array still marks them."""
    index = np.where(valid, np.arange(len(valid)), -1)
    previous = np.maximum.accumulate(index)
    previous[previous < 0] = int(np.argmax(valid))
    return values[previous]


def encode_array(rows: np.ndarray, kind: str, delta: int, name: str) -> bytes:
    """Difference, then write component by component and byte by byte."""
    steps = np.asarray(rows, dtype=np.int64).reshape(len(rows), -1)
    for _ in range(delta):
        steps = np.diff(steps, axis=0, prepend=np.zeros_like(steps[:1]))
    dtype = np.dtype(kind).newbyteorder("<")
    limits = np.iinfo(dtype)
    if steps.min() < limits.min or steps.max() > limits.max:
        raise ValueError(f"{name}: values do not fit in {kind}")
    by_component = np.ascontiguousarray(steps.T.astype(dtype))
    by_byte = np.frombuffer(by_component.tobytes(), dtype=np.uint8).reshape(-1, dtype.itemsize)
    return np.ascontiguousarray(by_byte.T).tobytes()


def decode_array(raw: bytes, description: dict) -> np.ndarray:
    """The reverse of encode_array; returns integers, one row per sample."""
    dtype = np.dtype(description["type"]).newbyteorder("<")
    count, components = description["count"], description["components"]
    size = count * components * dtype.itemsize
    by_byte = np.frombuffer(raw, dtype=np.uint8, count=size, offset=description["offset"])
    by_component = np.ascontiguousarray(by_byte.reshape(dtype.itemsize, -1).T).view(dtype)
    rows = by_component.reshape(components, count).T.astype(np.int64)
    for _ in range(description["delta"]):
        rows = np.cumsum(rows, axis=0)
    return rows


def pack_sequence(datasets, cache: Path, name: str) -> tuple[bytes, dict]:
    manifest = datasets.manifest()
    entries = [entry for entry in manifest["files"] if entry["sequence"] == name]
    if not entries:
        raise ValueError(f"The manifest has no files for {name}")
    for entry in entries:
        datasets.fetch(cache, entry)
    trial = datasets.load_sequence(cache / name, name)
    samples = len(trial.time)

    # IMU: integer counts, checked to decode to the published doubles.
    gyroscope = to_counts(trial.gyr, decode_gyroscope, f"{name} gyroscope")
    accelerometer = to_counts(trial.acc, decode_accelerometer, f"{name} accelerometer")

    # Time: nanoseconds since the first sample.
    elapsed_ns = trial.timestamp_ns - trial.timestamp_ns[0]
    if not np.array_equal(elapsed_ns.astype(float) * 1e-9, trial.time):
        raise ValueError(f"{name}: time does not decode exactly")

    # Reference at the IMU sample times.
    reference_ns, reference = datasets.read_csv(
        cache / name / "mav0" / "state_groundtruth_estimate0" / "data.csv", 17)
    valid = trial.valid
    columns = [0, 1, 2, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    vectors = interpolate_vectors(trial.timestamp_ns, reference_ns, reference[:, columns],
                                  datasets.MAX_REFERENCE_GAP_NS)
    if not np.array_equal(np.isfinite(vectors).all(axis=1), valid):
        raise ValueError(f"{name}: position and attitude references cover different samples")
    vectors = hold_invalid(vectors, valid)[::REFERENCE_EVERY]

    # q and -q are the same rotation. Keep the sign steady, so differences stay small.
    attitude = hold_invalid(trial.quat_gt, valid).copy()
    for i in range(1, samples):
        if np.dot(attitude[i], attitude[i - 1]) < 0:
            attitude[i] *= -1

    def steps(values, scale):
        return np.rint(values / scale).astype(np.int64)

    arrays = [
        ("time_ns", elapsed_ns, dict(type="int32", delta=2, unit="ns since the first sample")),
        ("gyroscope", gyroscope, dict(type="int16", delta=0, decode="count / 25 / 180 * pi",
                                      unit="rad/s, IMU frame")),
        ("accelerometer", accelerometer, dict(type="int16", delta=0, decode="count / 1200 * 9.80665",
                                              unit="m/s^2, IMU frame")),
        ("reference_valid", valid.astype(np.int64), dict(
            type="uint8", delta=0, unit="1 where the reference covers the sample")),
        ("reference_attitude", steps(attitude, ATTITUDE_SCALE), dict(
            type="int32", delta=1, scale=ATTITUDE_SCALE,
            unit="quaternion w x y z, IMU frame to world frame, +Z up")),
        ("reference_position", steps(vectors[:, 0:3], POSITION_SCALE), dict(
            type="int32", delta=2, scale=POSITION_SCALE, every=REFERENCE_EVERY, unit="m, world frame")),
        ("reference_velocity", steps(vectors[:, 3:6], VELOCITY_SCALE), dict(
            type="int32", delta=1, scale=VELOCITY_SCALE, every=REFERENCE_EVERY, unit="m/s, world frame")),
        ("reference_gyro_bias", steps(vectors[:, 6:9], BIAS_SCALE), dict(
            type="int32", delta=1, scale=BIAS_SCALE, every=REFERENCE_EVERY, unit="rad/s, IMU frame")),
        ("reference_accelerometer_bias", steps(vectors[:, 9:12], BIAS_SCALE), dict(
            type="int32", delta=1, scale=BIAS_SCALE, every=REFERENCE_EVERY, unit="m/s^2, IMU frame")),
    ]

    body = bytearray()
    described = {}
    for key, data, description in arrays:
        data = data.reshape(len(data), -1)
        description.update(components=data.shape[1], count=len(data), offset=len(body))
        description.setdefault("every", 1)
        raw = encode_array(data, description["type"], description["delta"], f"{name} {key}")
        body += raw
        # Decode again, so that a mistake here can never reach the page.
        if not np.array_equal(decode_array(bytes(body), description), data):
            raise ValueError(f"{name} {key}: does not decode to what was encoded")
        described[key] = description

    vicon = name.startswith(("V1_", "V2_"))
    header = {
        "format": FORMAT_VERSION,
        "sequence": name,
        "samples": samples,
        "duration_s": float(trial.time[-1]),
        "split": "development" if name in datasets.DEVELOPMENT else "validation",
        "reference": {
            "kind": "Vicon poses and IMU" if vicon else "Leica positions and IMU",
            "note": trial.metadata["ground_truth_warning"],
            "samples": int(valid.sum()),
        },
        "source": {
            "dataset": manifest["dataset"],
            "page": manifest["source"],
            "citation_doi": manifest["doi"],
            "archive_doi": manifest["archive_doi"],
            "terms": "In Copyright - Non-Commercial Use Permitted (ETH Research Collection)",
            "mirror_repository": manifest["mirror_repository"],
            "mirror_revision": manifest["mirror_revision"],
            "files": [{"path": entry["path"], "bytes": entry["bytes"], "sha256": entry["sha256"]}
                      for entry in entries],
        },
        "arrays": described,
    }
    encoded = json.dumps(header, separators=(",", ":")).encode()
    packed = MAGIC + struct.pack("<II", FORMAT_VERSION, len(encoded)) + encoded + bytes(body)
    # mtime=0 keeps the output identical from one build to the next.
    compressed = gzip.compress(packed, compresslevel=9, mtime=0)
    summary = {
        "name": name,
        "file": f"{name}.pack",
        "samples": samples,
        "duration_s": round(float(trial.time[-1]), 3),
        "split": header["split"],
        "reference": header["reference"]["kind"],
        "bytes": len(compressed),
        "sha256": hashlib.sha256(compressed).hexdigest(),
    }
    return compressed, summary


def unpack(compressed: bytes) -> tuple[dict, dict]:
    """Read a pack file back: the header, and each array as integers."""
    packed = gzip.decompress(compressed)
    if packed[:4] != MAGIC:
        raise ValueError("Not a pack file")
    version, length = struct.unpack("<II", packed[4:12])
    if version != FORMAT_VERSION:
        raise ValueError(f"Unknown pack format {version}")
    header = json.loads(packed[12:12 + length])
    body = packed[12 + length:]
    return header, {key: decode_array(body, description) for key, description in header["arrays"].items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--cache", type=Path, default=ROOT / ".benchmark-cache" / "euroc",
                        help="where verified EuRoC files are kept (default: the benchmark's cache)")
    parser.add_argument("--output", type=Path, required=True, help="directory for the .pack files")
    parser.add_argument("--sequences", nargs="+", default=list(DEFAULT_SEQUENCES),
                        help="recordings to pack; the first is the one the page opens with")
    args = parser.parse_args()

    datasets = load_datasets_module()
    unknown = sorted(set(args.sequences) - set(datasets.SEQUENCES))
    if unknown or len(set(args.sequences)) != len(args.sequences):
        parser.error(f"choose distinct EuRoC recording names; unknown: {unknown}")

    args.output.mkdir(parents=True, exist_ok=True)
    recordings = []
    for name in args.sequences:
        compressed, summary = pack_sequence(datasets, args.cache, name)
        (args.output / summary["file"]).write_bytes(compressed)
        recordings.append(summary)
        print(f"{name}: {summary['samples']} samples, {summary['duration_s']:.1f} s, "
              f"{summary['bytes'] / 1024:.0f} KiB")
    index = {"format": FORMAT_VERSION, "recordings": recordings}
    (args.output / "recordings.json").write_text(json.dumps(index, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
