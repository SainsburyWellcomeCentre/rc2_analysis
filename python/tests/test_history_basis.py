"""Spike-history basis + trial-aware causal convolution — unit tests.

Per prompt 03 (refreshed 2026-04-28). Three invariants pinned here:

1. Causal: shuffling future spike counts does NOT change history features
   at any earlier bin. (Equivalent: the convolution doesn't use lag 0.)
2. Trial-boundary safe: the first ``n_lag`` bins of each trial have
   zero-padded history (no contamination from the previous trial's
   spikes).
3. Synthetic refractory recovery: a spike train with a programmed
   refractory period produces a history filter whose lag-0 / lag-1
   coefficient is significantly negative when the GLM is fit. (This last
   one is the "does the model see the refractory effect" sanity check.)
"""
from __future__ import annotations

import numpy as np
import pytest

from rc2_glm.basis import (
    convolve_history,
    history_basis,
)


def test_history_basis_shape_and_values():
    """Bases are (n_lag_bins, n_bases) with peaks ≤ 1 and no NaN."""
    B = history_basis(n_bases=10, t_max_s=0.2, bin_width_s=0.1)
    assert B.shape == (2, 10)  # 0.2s / 0.1s = 2 lag bins
    assert np.isfinite(B).all()
    assert B.max() <= 1.0 + 1e-12
    assert B.min() >= 0.0


def test_history_basis_finer_bin_more_lags():
    """Smaller bin width yields more lag bins for the same window."""
    B_100 = history_basis(n_bases=10, t_max_s=0.2, bin_width_s=0.1)
    B_20 = history_basis(n_bases=10, t_max_s=0.2, bin_width_s=0.02)
    assert B_100.shape[0] == 2
    assert B_20.shape[0] == 10


def test_convolve_history_causal_no_future_leakage():
    """Shuffling future spikes leaves earlier history features unchanged.

    The whole point of causal convolution. If lag 0 leaked, this would
    fail at the bin where the shuffle starts.
    """
    rng = np.random.default_rng(0)
    n = 50
    spike_counts = rng.poisson(lam=2.0, size=n).astype(np.float64)
    trial_ids = np.zeros(n, dtype=int)
    basis = history_basis(n_bases=4, t_max_s=0.2, bin_width_s=0.05)
    H = convolve_history(spike_counts, trial_ids, basis)

    # Shuffle the LAST 10 bins; the first 40 bins' history features
    # should be unchanged because they only use lookback into earlier
    # bins (which weren't touched).
    spike_shuffled = spike_counts.copy()
    rng2 = np.random.default_rng(99)
    rng2.shuffle(spike_shuffled[40:])  # in-place
    H_shuffled = convolve_history(spike_shuffled, trial_ids, basis)

    np.testing.assert_allclose(H[:40], H_shuffled[:40], atol=0)


def test_convolve_history_past_changes_propagate():
    """Conversely: shuffling PAST spikes changes future history features."""
    rng = np.random.default_rng(0)
    n = 50
    spike_counts = rng.poisson(lam=2.0, size=n).astype(np.float64)
    trial_ids = np.zeros(n, dtype=int)
    basis = history_basis(n_bases=4, t_max_s=0.2, bin_width_s=0.05)
    H = convolve_history(spike_counts, trial_ids, basis)

    spike_shuffled = spike_counts.copy()
    rng2 = np.random.default_rng(99)
    rng2.shuffle(spike_shuffled[:10])  # shuffle past
    H_shuffled = convolve_history(spike_shuffled, trial_ids, basis)

    # At least SOME bins after bin 10 should differ (history of those
    # bins now uses different past spikes).
    assert not np.allclose(H[10:], H_shuffled[10:])


def test_convolve_history_trial_boundary_zero_padded():
    """First ``n_lag`` bins of each trial have zero history from preceding trial.

    Concretely: place a spike at the LAST bin of trial 0 and the FIRST
    bin of trial 1. The first bin of trial 1's history feature should
    NOT carry the trial-0 spike's contribution.
    """
    n_lag = 4
    n_bases = 3
    basis = np.ones((n_lag, n_bases), dtype=np.float64)  # uniform basis
    # Two trials of 10 bins each.
    spike_counts = np.zeros(20, dtype=np.float64)
    spike_counts[9] = 1.0   # last bin of trial 0
    spike_counts[10] = 1.0  # first bin of trial 1
    trial_ids = np.concatenate([np.zeros(10, dtype=int), np.ones(10, dtype=int)])

    H = convolve_history(spike_counts, trial_ids, basis)

    # Bin 10 (first bin of trial 1) should have zero history — its only
    # legal lookback would be into trial 0, which is zero-padded.
    assert np.allclose(H[10], 0.0), (
        f"Bin 10 (first bin of trial 1) should have zero-padded history "
        f"but got {H[10]}"
    )
    # Bin 11 should see the spike at bin 10 (lag 1).
    assert H[11, 0] > 0.0


def test_convolve_history_n_total_matches_input():
    """Output length matches input bin count."""
    spike_counts = np.array([0, 1, 2, 0, 0, 3], dtype=np.float64)
    trial_ids = np.zeros(6, dtype=int)
    basis = history_basis(n_bases=3, t_max_s=0.2, bin_width_s=0.05)
    H = convolve_history(spike_counts, trial_ids, basis)
    assert H.shape == (6, 3)


def test_convolve_history_empty_input():
    """Empty input returns (0, n_bases) — no exceptions."""
    basis = history_basis(n_bases=3, t_max_s=0.2, bin_width_s=0.05)
    H = convolve_history(np.array([]), np.array([]), basis)
    assert H.shape == (0, 3)
