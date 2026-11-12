#!/usr/bin/env python3
"""Reproduce the fixed thrust-axis calibration from development IMU data only.

This is an offline development calibration, not an estimator update or proof
that the selected samples are airborne. The two sequences are fixed here; no
ground-truth file, orientation, velocity, or sensor-bias column is accessed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from datasets import fetch, manifest, read_csv

ROOT = Path(__file__).resolve().parents[2]
DEVELOPMENT = ("MH_01_easy", "V1_02_medium")
SAMPLE_PERIOD_NS = 5_000_000
# Rotor detector constants, identical to MAV_Tuning::rotors in mav_estimator.h.
VIBRATION_LOWPASS_HZ = 20
VIBRATION_ENERGY_TAU_SECONDS = .5
ROTOR_ON_RMS_M_S2 = .8
ROTOR_OFF_RMS_M_S2 = .5


def rotor_detector(time: np.ndarray, force: np.ndarray) -> np.ndarray:
    """Return the runtime rotor state for each sample.

    This is a sample-by-sample port of RotorDetector in rotor_vibration.h, so
    calibration selects exactly the samples where the estimator treats the
    rotors as running: vibration is the specific force minus its exact
    zero-order-hold second-order Butterworth low-pass (20 Hz cutoff, variable
    sample interval), its squared norm feeds a 0.5 s exponential moving
    average, and the RMS switches the state on above 0.8 m/s^2 and off below
    0.5 m/s^2.
    """
    tau = math.sqrt(2.0) / (2 * math.pi * VIBRATION_LOWPASS_HZ)
    value = [0.0, 0.0, 0.0]
    rate = [0.0, 0.0, 0.0]
    energy = 0.0
    rotors = False
    state = np.zeros(len(time), dtype=bool)
    for index in range(len(time)):
        sample = [float(component) for component in force[index]]
        dt = float(time[index] - time[index - 1]) if index else 0.0
        if index == 0:
            # LowPassFilter::update starts from the first sample.
            value = sample
        else:
            h = dt / tau
            decay, cosine, sine = math.exp(-h), math.cos(h), math.sin(h)
            error = [value[axis] - sample[axis] for axis in range(3)]
            value = [sample[axis] + decay * ((cosine + sine) * error[axis] + tau * sine * rate[axis])
                     for axis in range(3)]
            rate = [decay * ((-2 * sine / tau) * error[axis] + (cosine - sine) * rate[axis])
                    for axis in range(3)]
        vibration = [sample[axis] - value[axis] for axis in range(3)]
        gain = 1.0 - math.exp(-dt / VIBRATION_ENERGY_TAU_SECONDS)
        energy += gain * (sum(component * component for component in vibration) - energy)
        rms = math.sqrt(max(0.0, energy))
        if rms > ROTOR_ON_RMS_M_S2:
            rotors = True
        elif rms < ROTOR_OFF_RMS_M_S2:
            rotors = False
        state[index] = rotors
    return state


def calibrate(cache: Path) -> dict:
    records = []
    files = manifest()["files"]
    for sequence in DEVELOPMENT:
        relative = f"{sequence}/mav0/imu0/data.csv"
        entries = [entry for entry in files if entry["path"] == relative]
        if len(entries) != 1:
            raise ValueError(f"Expected exactly one pinned development IMU file: {relative}")
        entry = entries[0]
        timestamps, imu = read_csv(fetch(cache, entry), 7)
        relative_ns = timestamps - timestamps[0]
        clock_error = relative_ns - np.arange(len(timestamps), dtype=np.int64) * SAMPLE_PERIOD_NS
        if np.max(np.abs(clock_error)) > 1000:
            raise ValueError(f"{sequence}: timestamps deviate from the nominal 200 Hz clock by more than 1 us")
        force = imu[:, 3:6]
        if not np.isfinite(force).all():
            raise ValueError(f"{sequence}: non-finite accelerometer measurements")
        # The estimator receives these same relative times in seconds.
        time = relative_ns.astype(float) * 1e-9
        selected = rotor_detector(time, force) & (relative_ns >= 5_000_000_000)
        if np.count_nonzero(selected) < 200:
            raise ValueError(f"{sequence}: insufficient rotor-vibration samples for calibration")
        median = np.median(force[selected], axis=0)
        pitch = float(np.rad2deg(np.arctan2(-median[2], median[0])))
        lateral = float(np.rad2deg(np.arctan2(median[1], np.hypot(median[0], median[2]))))
        records.append({
            "sequence": sequence,
            "input_path": relative,
            "input_sha256": entry["sha256"],
            "samples": len(timestamps),
            "selected_samples": int(np.count_nonzero(selected)),
            "selected_duration_seconds": float(np.count_nonzero(selected) * 0.005),
            "selected_fraction": float(np.mean(selected)),
            "median_specific_force_m_s2": median.tolist(),
            "median_force_axis_imu": (median / np.linalg.norm(median)).tolist(),
            "pitch_degrees": pitch,
            "lateral_degrees": lateral,
            "max_nominal_clock_deviation_ns": int(np.max(np.abs(clock_error))),
        })
    combined_pitch = float(np.mean([record["pitch_degrees"] for record in records]))
    nominal_pitch = float(5 * np.floor(combined_pitch / 5 + .5))
    radians = np.deg2rad(nominal_pitch)
    return {
        "purpose": "Offline sensor-only development calibration; no ground truth or held-out recordings",
        "assumptions": [
            "Rotor thrust dominates the median specific force during selected motor-vibration samples.",
            "The two development recordings receive equal weight regardless of their duration.",
            "The nominal mounting angle is rounded to the nearest five degrees.",
            "The fixed thrust axis has zero IMU-y component by an approximate bilateral-symmetry assumption.",
            "The result is empirical calibration, not a published exact Firefly-to-IMU transform.",
        ],
        "selection": {
            "detector": "runtime RotorDetector of rotor_vibration.h",
            "vibration": "specific force minus exact zero-order-hold second-order Butterworth low-pass",
            "lowpass_cutoff_hz": VIBRATION_LOWPASS_HZ,
            "nominal_sample_rate_hz": 200,
            "energy_ewma_time_constant_seconds": VIBRATION_ENERGY_TAU_SECONDS,
            "rms_on_threshold_m_s2": ROTOR_ON_RMS_M_S2,
            "rms_off_threshold_m_s2": ROTOR_OFF_RMS_M_S2,
            "initial_filter_transient_excluded_seconds": 5,
        },
        "per_sequence": records,
        "equal_sequence_mean_pitch_degrees": combined_pitch,
        "nominal_thrust_tilt_degrees": nominal_pitch,
        "nominal_thrust_axis_imu": [float(np.cos(radians)), 0., float(-np.sin(radians))],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=ROOT / ".benchmark-cache/euroc")
    args = parser.parse_args()
    print(json.dumps(calibrate(args.cache), indent=2, allow_nan=False))
