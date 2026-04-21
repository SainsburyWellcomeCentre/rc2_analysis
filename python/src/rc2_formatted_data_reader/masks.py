"""Motion and stationary mask computation.

The simple `motion_mask` / `stationary_mask` functions match the API in the
reader: pure velocity thresholding.

`treadmill_motion_mask` replicates the more complete MATLAB
`Trial.treadmill_motion_mask` logic (velocity threshold + acceleration
threshold + minimum stationary-bout duration). It does NOT apply
`analysis_window`; the caller composes that separately if needed.
"""

from __future__ import annotations

import numpy as np


def motion_mask(velocity: np.ndarray, threshold: float = 1.0) -> np.ndarray:
    return np.abs(velocity) >= threshold


def stationary_mask(velocity: np.ndarray, threshold: float = 1.0) -> np.ndarray:
    return np.abs(velocity) < threshold


def treadmill_motion_mask(
    velocity: np.ndarray,
    acceleration: np.ndarray,
    fs: float,
    vel_thresh: float = 1.0,
    acc_thresh: float = 0.5,
    min_dur: float = 0.2,
) -> np.ndarray:
    """Replicate MATLAB Trial.treadmill_motion_mask.

    Stationary samples are those with `velocity < vel_thresh` AND
    `|acceleration| < acc_thresh`. Any contiguous stationary run shorter
    than `min_dur` seconds is promoted to motion. Returns the negation
    (motion mask).
    """
    velocity = np.asarray(velocity, dtype=np.float64)
    stationary = (velocity < vel_thresh) & (np.abs(acceleration) < acc_thresh)

    if not stationary.any():
        return ~stationary

    min_samples = int(round(min_dur * fs))
    diff = np.diff(stationary.astype(np.int8), prepend=0, append=0)
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    durations = ends - starts

    for s, e, d in zip(starts, ends, durations):
        if d < min_samples:
            stationary[s:e] = False

    return ~stationary
