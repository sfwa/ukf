"""Quaternion attitude errors with one initial heading alignment.

Quaternion convention: scalar-first, body to world, world +Z up. Inclination
and heading use the swing/twist decomposition of q_est * conjugate(q_gt), as
in the BROAD paper. The fixed yaw gauge is estimated in a one-second window
starting at the first optical-valid sample (recordings may have leading optical
gaps). No roll/pitch adjustment or whole-trajectory alignment is performed.
"""
from __future__ import annotations

import numpy as np


def multiply(a, b):
    aw, ax, ay, az = np.moveaxis(a, -1, 0)
    bw, bx, by, bz = np.moveaxis(b, -1, 0)
    return np.stack((aw*bw-ax*bx-ay*by-az*bz,
                     aw*bx+ax*bw+ay*bz-az*by,
                     aw*by-ax*bz+ay*bw+az*bx,
                     aw*bz+ax*by-ay*bx+az*bw), axis=-1)


def _angles(delta):
    twist_norm = np.hypot(delta[:, 0], delta[:, 3])
    inclination = 2 * np.arccos(np.clip(twist_norm, 0, 1))
    total = 2 * np.arccos(np.clip(np.abs(delta[:, 0]), 0, 1))
    heading = (2*np.arctan2(delta[:, 3], delta[:, 0])+np.pi) % (2*np.pi)-np.pi
    heading[twist_norm < 1e-10] = np.nan  # heading undefined for a 180-degree tilt
    return inclination, heading, total


def compute_errors(q_est, q_gt, valid, time, burn_in=5.0, alignment_seconds=1.0, score_mask=None):
    """Return (summary, per_sample_errors), all angular quantities in degrees.

    Reported total/heading errors are aligned by a SINGLE initial yaw offset.
    Raw total/heading errors are also retained, but arbitrary optical-world yaw
    means they are not generally absolute magnetic-heading accuracy measures.
    Ground-truth gaps are excluded; invalid estimates raise rather than vanish
    from the score. burn_in is relative to time[0]; alignment_seconds is a
    window length starting at the first optical-valid sample. Scoring begins
    only after BOTH the burn-in and alignment window have ended.
    Optional score_mask selects motion samples for scoring only: initial yaw
    alignment still uses all optical-valid samples in the initial window.
    """
    q_est, q_gt = np.asarray(q_est, float), np.asarray(q_gt, float)
    valid, time = np.asarray(valid, bool), np.asarray(time, float)
    n = len(time)
    if q_est.shape != (n, 4) or q_gt.shape != (n, 4) or valid.shape != (n,):
        raise ValueError("Inconsistent metric input shapes")
    if n < 2 or not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0):
        raise ValueError("Expected finite increasing timestamps")
    if not 0 < alignment_seconds <= burn_in:
        raise ValueError("Require 0 < alignment_seconds <= burn_in")
    for label, q in (("estimate", q_est), ("ground truth", q_gt[valid])):
        if not np.all(np.isfinite(q)) or not np.allclose(np.linalg.norm(q, axis=1), 1, atol=1e-4):
            raise ValueError(f"Invalid {label} quaternion")
    # Safe placeholders prevent warnings at invalid optical samples.
    reference = q_gt.copy()
    reference[~valid] = [1., 0., 0., 0.]
    reference[:, 1:] *= -1
    delta = multiply(q_est, reference)
    delta /= np.linalg.norm(delta, axis=1)[:, None]
    inclination, raw_heading, raw_total = _angles(delta)
    if not np.any(valid):
        raise ValueError("No valid ground truth")
    alignment_start = float(time[np.flatnonzero(valid)[0]])
    alignment_end = alignment_start + alignment_seconds
    initial = valid & (time < alignment_end) & np.isfinite(raw_heading)
    score = valid & (time >= max(time[0]+burn_in, alignment_end))
    if score_mask is not None:
        score_mask = np.asarray(score_mask, bool)
        if score_mask.shape != (n,):
            raise ValueError("Inconsistent score mask shape")
        score &= score_mask
    if not np.any(initial) or not np.any(score):
        raise ValueError("No valid ground truth in initialization or scored interval")
    mean = np.mean(np.exp(1j*raw_heading[initial]))
    if abs(mean) < 0.1:
        raise ValueError("Initial heading is too dispersed for a reliable fixed alignment")
    offset = float(np.angle(mean))
    correction = np.array([np.cos(offset/2), 0., 0., -np.sin(offset/2)])
    _, heading, total = _angles(multiply(correction, delta))
    errors = {
        "inclination_deg": np.rad2deg(inclination),
        "heading_deg": np.rad2deg(heading),
        "total_deg": np.rad2deg(total),
        "raw_heading_deg": np.rad2deg(raw_heading),
        "raw_total_deg": np.rad2deg(raw_total),
        "score_mask": score,
    }
    summary = {"samples": int(score.sum()), "heading_alignment_deg": float(np.rad2deg(-offset)),
               "alignment_start_seconds": alignment_start-float(time[0])}
    for key in ("inclination", "heading", "total", "raw_heading", "raw_total"):
        values = errors[key+"_deg"][score]
        values = values[np.isfinite(values)]
        if not len(values):
            raise ValueError(f"No defined {key} errors in scored interval")
        summary[key+"_rmse_deg"] = float(np.sqrt(np.mean(values**2)))
        summary[key+"_p95_deg"] = float(np.percentile(np.abs(values), 95))
    summary["heading_samples"] = int(np.count_nonzero(score & np.isfinite(errors["heading_deg"])))
    return summary, errors
