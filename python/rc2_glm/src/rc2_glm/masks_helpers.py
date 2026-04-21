"""Velocity filtering and acceleration helpers.

Replicates lib/fcn/general/filter_trace.m: forward-then-backward
order-3 Butterworth low-pass at 50 Hz (zero-phase). MATLAB applies the
filter twice (forward, then on the reversed signal) instead of using
filtfilt — for our purposes scipy.signal.filtfilt is the standard
zero-phase equivalent.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import butter, filtfilt


def filter_velocity(
    velocity: np.ndarray, fs: float, cutoff: float = 50.0, order: int = 3
) -> np.ndarray:
    if velocity.size < 3 * order + 1:
        return velocity.astype(np.float64, copy=True)
    wn = cutoff / (fs / 2.0)
    b, a = butter(order, wn, btype="low")
    return filtfilt(b, a, velocity).astype(np.float64)


def acceleration_from_velocity(velocity: np.ndarray) -> np.ndarray:
    """Replicate MATLAB: acc = 100 * diff(v); val = [acc(1); midpoint(acc); acc(end)].

    Returns the same length as `velocity`.
    """
    if velocity.size < 2:
        return np.zeros_like(velocity, dtype=np.float64)
    acc = 100.0 * np.diff(velocity.astype(np.float64))
    out = np.empty_like(velocity, dtype=np.float64)
    out[0] = acc[0]
    out[1:-1] = 0.5 * (acc[:-1] + acc[1:])
    out[-1] = acc[-1]
    return out
