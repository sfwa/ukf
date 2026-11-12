"""Causal, pinned reference implementations; scalar-first body-to-world output.

All outputs use a right-handed world with +Z up and (for nine-axis filters)
horizontal magnetic north along +X. ahrs.Mahony and VQF use north along +Y
internally; a fixed world-frame rotation adapts their nine-axis outputs.
"""
from __future__ import annotations

import numpy as np
from ahrs.filters import Madgwick, Mahony
from scipy.spatial.transform import Rotation
from vqf import VQF

REFERENCE_NAMES = ("madgwick", "mahony", "vqf")
PARAMETERS = {
    "madgwick": {"implementation": "ahrs 0.4.0", "gain": 0.05},
    "mahony": {"implementation": "ahrs 0.4.0", "k_P": 1.0, "k_I": 0.3},
    "vqf": {"implementation": "vqf 2.0.1", "parameters": "online defaults"},
}


def _yaw_rotate(q: np.ndarray, angle: float) -> np.ndarray:
    """Left-multiply one or many scalar-first quaternions by a world Z rotation."""
    c, s = np.cos(angle / 2), np.sin(angle / 2)
    w, x, y, z = np.moveaxis(q, -1, 0)
    return np.stack((c*w-s*z, c*x-s*y, c*y+s*x, c*z+s*w), axis=-1)


def initial_attitude(acc: np.ndarray, mag: np.ndarray | None = None) -> np.ndarray:
    """Shortest rotation mapping measured up to +Z, then magnetic yaw if supplied."""
    norm = np.linalg.norm(acc)
    if norm < 1e-12:
        raise ValueError("Cannot initialize from a zero acceleration vector")
    up = np.asarray(acc, dtype=float) / norm
    # Same shortest-vector rotation used by Eigen::Quaternion::FromTwoVectors.
    if up[2] < -1 + 1e-12:
        q = np.array([0., 1., 0., 0.])
    else:
        q = np.array([1 + up[2], up[1], -up[0], 0.])
        q /= np.linalg.norm(q)
    if mag is not None:
        m = Rotation.from_quat(q, scalar_first=True).apply(mag)
        if np.linalg.norm(m[:2]) < 1e-12:
            raise ValueError("Cannot initialize from a vertical magnetic field")
        q = _yaw_rotate(q, -np.arctan2(m[1], m[0]))
    return q


def run_reference(name, time, gyr, acc, mag, axes):
    """Run an online reference without ground truth or future sensor samples.

    ahrs filters get each measured dt and midpoint gyro samples. VQF requires a
    constant sampling interval (both datasets provide this); variable dt is
    rejected rather than silently resampled. Initialization consumes the first
    acc/mag sample but does not integrate an extra gyro interval at t=0.
    """
    if name not in REFERENCE_NAMES or axes not in (6, 9):
        raise ValueError("Unknown reference or axis count")
    time = np.asarray(time, dtype=float)
    gyr, acc, mag = (np.ascontiguousarray(x, dtype=float) for x in (gyr, acc, mag))
    n = len(time)
    if n < 2 or time.ndim != 1 or any(x.shape != (n, 3) for x in (gyr, acc, mag)):
        raise ValueError("Expected time[N] and sensors[N,3], with N >= 2")
    if not all(np.all(np.isfinite(x)) for x in (time, gyr, acc, mag)):
        raise ValueError("Non-finite reference input")
    dt = np.diff(time)
    if np.any(dt <= 0):
        raise ValueError("Timestamps must increase")
    if np.any(np.linalg.norm(acc, axis=1) < 1e-12):
        raise ValueError("Zero accelerometer vector")
    if axes == 9 and np.any(np.linalg.norm(mag, axis=1) < 1e-12):
        raise ValueError("Zero magnetometer vector")
    q = np.empty((n, 4))
    if name == "vqf":
        if not np.allclose(dt, dt[0], rtol=1e-5, atol=1e-9):
            raise ValueError("VQF requires uniformly sampled input")
        filt = VQF(float(dt[0]))
        filt.updateAcc(acc[0])
        if axes == 9:
            filt.updateMag(mag[0])
        get_q = filt.getQuat9D if axes == 9 else filt.getQuat6D
        q[0] = get_q()
        for i in range(1, n):
            filt.updateGyr((gyr[i-1] + gyr[i]) * 0.5)
            filt.updateAcc(acc[i])
            if axes == 9:
                filt.updateMag(mag[i])
            q[i] = get_q()
    else:
        filt = Madgwick(gain=0.05) if name == "madgwick" else Mahony(k_P=1.0, k_I=0.3)
        q[0] = initial_attitude(acc[0], mag[0] if axes == 9 else None)
        if name == "mahony" and axes == 9:
            q[0] = _yaw_rotate(q[0], np.pi/2)
        for i in range(1, n):
            omega = (gyr[i-1] + gyr[i]) * 0.5
            if axes == 9:
                q[i] = filt.updateMARG(q[i-1], omega, acc[i], mag[i], dt=float(dt[i-1]))
            else:
                q[i] = filt.updateIMU(q[i-1], omega, acc[i], dt=float(dt[i-1]))
    if axes == 9 and name in ("mahony", "vqf"):
        q = _yaw_rotate(q, -np.pi/2)
    if not np.all(np.isfinite(q)) or not np.allclose(np.linalg.norm(q, axis=1), 1, atol=1e-5):
        raise ValueError(f"{name} produced invalid orientation output")
    return q
