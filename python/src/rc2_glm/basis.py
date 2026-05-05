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


def raised_cosine_basis_linear(
    x: np.ndarray, n_bases: int, x_min: float, x_max: float
) -> np.ndarray:
    """Linear-spaced raised cosine bases over [x_min, x_max].

    Used for value-axis bases on signals that include negative values
    (e.g. z-scored behavioural covariates like face motion energy). For
    positive-only signals with broad dynamic range, prefer
    ``raised_cosine_basis`` (log-shifted Weber spacing).

    Centres are linearly spaced from x_min to x_max; width is 1.5×delta
    so adjacent bases overlap and the partition-of-unity property holds
    in the interior. Values outside [x_min, x_max] are clipped to the
    nearest edge basis (saturating, not zero).

    Returns array of shape (len(x), n_bases).
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    centers = np.linspace(x_min, x_max, n_bases)
    delta = (x_max - x_min) / (n_bases - 1) if n_bases > 1 else (x_max - x_min)
    width = delta * 1.5

    z = (x[:, None] - centers[None, :]) / width * np.pi
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


def history_basis(
    n_bases: int, t_max_s: float, bin_width_s: float,
) -> np.ndarray:
    """Log-spaced raised-cosine bases over post-spike lags [Δt, t_max_s].

    Each row of the returned array is one lag bin (1..n_lag_bins, where
    ``n_lag_bins = round(t_max_s / bin_width_s)``); each column is one
    basis function. The bases are evaluated on ``log(lag + 0.5)`` to give
    Weber-style spacing — fine resolution near lag 0 (refractory), coarse
    at long lags (slow adaptation). Pillow et al. 2008 convention.

    Bases at lag 0 (the *current* bin's spike) are intentionally absent —
    the convolution must use *past* spikes only to remain causal.

    Returns array of shape (n_lag_bins, n_bases). Each column has unit
    peak, like ``raised_cosine_basis``.
    """
    n_lag_bins = max(1, int(round(t_max_s / bin_width_s)))
    # Lag indices 1..n_lag_bins (lag 0 excluded for causality — see
    # convolve_history's Δt offset).
    lags = np.arange(1, n_lag_bins + 1, dtype=np.float64)
    return raised_cosine_basis(lags, n_bases, x_min=lags[0], x_max=lags[-1])


def convolve_history(
    spike_counts: np.ndarray,
    trial_ids: np.ndarray,
    basis: np.ndarray,
) -> np.ndarray:
    """Per-bin causal convolution of spike counts with a history basis.

    For each basis column k, the output at bin t within trial T is::

        h_k(t) = Σ_{τ=1..n_lag} basis[τ-1, k] · spike_counts[t - τ]

    where the sum is restricted to bins within the same trial T (no
    cross-trial leakage). Bins where the lookback would cross the trial
    start are zero-padded — equivalent to "no spikes before the trial"
    so the bins at the very start of each trial have small history
    features.

    Parameters
    ----------
    spike_counts : (n_total_bins,) int or float
        Per-bin spike count for ONE cluster, ordered as the binned table.
    trial_ids : (n_total_bins,) int
        Per-bin trial identifier; convolution is reset at each trial
        boundary (zero-padded lookback).
    basis : (n_lag_bins, n_bases)
        Output of ``history_basis``.

    Returns
    -------
    (n_total_bins, n_bases) float
    """
    spike_counts = np.asarray(spike_counts, dtype=np.float64)
    trial_ids = np.asarray(trial_ids)
    n_total = spike_counts.size
    n_lag, n_bases = basis.shape
    if n_total == 0:
        return np.zeros((0, n_bases), dtype=np.float64)
    if trial_ids.shape[0] != n_total:
        raise ValueError(
            f"spike_counts ({n_total}) and trial_ids ({trial_ids.shape[0]}) "
            "must have the same length"
        )

    out = np.zeros((n_total, n_bases), dtype=np.float64)
    # Identify trial boundaries — first bin index of each trial.
    trial_starts = np.concatenate(
        [[0], np.where(np.diff(trial_ids) != 0)[0] + 1, [n_total]]
    )
    for t_start, t_end in zip(trial_starts[:-1], trial_starts[1:]):
        # Per-trial spike counts.
        sc = spike_counts[t_start:t_end]
        n_t = sc.size
        # For each lag τ in 1..n_lag, the contribution to bin i within
        # this trial is sc[i - τ] if i - τ >= 0, else 0.
        # Vectorise via shifted slices.
        for k in range(n_bases):
            for lag in range(1, n_lag + 1):
                w = basis[lag - 1, k]
                if w == 0.0:
                    continue
                # bin i within trial uses sc[i - lag]; valid for i >= lag.
                if lag >= n_t:
                    continue
                out[t_start + lag : t_end, k] += w * sc[: n_t - lag]
    return out
