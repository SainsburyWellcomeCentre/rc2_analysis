"""Raised-cosine and onset-kernel bases.

Direct ports of `make_raised_cosine_basis` (line 5506) and
`make_onset_kernel_basis` (line 5531) from
scripts/glm_single_cluster_analysis.m. These must match MATLAB output
to ~1e-12.
"""

from __future__ import annotations

import numpy as np


def raised_cosine_basis(
    x: np.ndarray, n_bases: int, x_min: float, x_max: float
) -> np.ndarray:
    """Log-shifted (Weber-law) raised cosine bases over [x_min, x_max].

    Returns array of shape (len(x), n_bases).
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    epsilon = 0.5
    log_x = np.log(x + epsilon)
    log_min = np.log(x_min + epsilon)
    log_max = np.log(x_max + epsilon)

    centers = np.linspace(log_min, log_max, n_bases)
    delta = (log_max - log_min) / (n_bases - 1) if n_bases > 1 else (log_max - log_min)
    width = delta * 1.5

    z = (log_x[:, None] - centers[None, :]) / width * np.pi
    np.clip(z, -np.pi, np.pi, out=z)
    return 0.5 * (1.0 + np.cos(z))


def onset_kernel_basis(
    t_since_onset: np.ndarray, n_bases: int, t_max: float
) -> np.ndarray:
    """Causal raised cosine bases over [0, t_max] (Park et al. 2014 style).

    Bases are zero for t_since_onset < 0 (stationary periods).
    Returns array of shape (len(t_since_onset), n_bases).
    """
    t = np.asarray(t_since_onset, dtype=np.float64).ravel()
    n = t.size
    B = np.zeros((n, n_bases), dtype=np.float64)

    centers = np.linspace(0.0, t_max, n_bases)
    delta = t_max / (n_bases - 1) if n_bases > 1 else t_max
    width = delta * 1.5

    motion_mask = t >= 0
    t_motion = t[motion_mask]
    if t_motion.size:
        z = (t_motion[:, None] - centers[None, :]) / width * np.pi
        np.clip(z, -np.pi, np.pi, out=z)
        B[motion_mask] = 0.5 * (1.0 + np.cos(z))
    return B
