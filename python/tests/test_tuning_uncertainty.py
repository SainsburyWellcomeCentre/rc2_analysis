"""Per-trial tuning-curve uncertainty band — unit tests.

Prompt 12 (2026-04-27): adds per-trial spread bands to the model-row
tuning-curve panels. The spread is computed by predicting per-bin
rates, collapsing to per-trial means, then summarising (q25/q75 by
default). Tests live here to keep the bands' invariants pinned.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest

from rc2_glm.config import GLMConfig, MIN_TRIALS_FOR_BAND
from rc2_glm.plots import (
    _band_edges,
    _marginal_curve_with_spread,
    _per_trial_summary,
)


# ----------------------------------------------------------------------
# _per_trial_summary — pure-numerics tests, no GLM machinery
# ----------------------------------------------------------------------


def test_per_trial_summary_zero_heterogeneity_collapses_band():
    """When every trial gives identical predictions, IQR width = 0."""
    # 4 trials, 5 bins each, all bins have identical mu — so per-trial
    # means are also identical and the IQR is exactly zero.
    n_trials = 4
    n_bins_per_trial = 5
    trial_ids = np.repeat(np.arange(n_trials), n_bins_per_trial)
    mu = np.full(trial_ids.size, 2.5, dtype=np.float64)
    s = _per_trial_summary(mu, trial_ids)
    assert s["n_trials"] == n_trials
    assert s["q25"] == s["q75"], "IQR should collapse when trials agree"
    assert s["q05"] == s["q95"], "wide quantiles should also collapse"
    assert s["std"] == 0.0


def test_per_trial_summary_iqr_within_q05_q95():
    """q25/q75 must lie inside q05/q95 by construction."""
    rng = np.random.default_rng(42)
    n_trials = 50
    n_bins_per_trial = 3
    trial_ids = np.repeat(np.arange(n_trials), n_bins_per_trial)
    # per-trial means drawn from a non-degenerate distribution
    per_trial_means = rng.normal(loc=5.0, scale=1.5, size=n_trials)
    mu = np.repeat(per_trial_means, n_bins_per_trial)
    s = _per_trial_summary(mu, trial_ids)
    assert s["q05"] <= s["q25"] <= s["q75"] <= s["q95"]


def test_per_trial_summary_drops_nan_trials():
    """A trial whose bins are all NaN doesn't count toward n_trials."""
    trial_ids = np.array([0, 0, 1, 1, 2, 2])
    mu = np.array([1.0, 1.0, np.nan, np.nan, 3.0, 3.0])
    s = _per_trial_summary(mu, trial_ids)
    # Trials 0 and 2 contribute; trial 1 is dropped (mean is NaN).
    assert s["n_trials"] == 2


def test_per_trial_summary_handles_none_mu():
    """When the predictor returned None (rank-deficient design), report NaNs."""
    s = _per_trial_summary(None, trial_ids=np.array([0, 1]))
    assert np.isnan(s["q25"])
    assert np.isnan(s["q75"])
    assert s["n_trials"] == 0


# ----------------------------------------------------------------------
# _band_edges — mode dispatch + sparse-bin guard
# ----------------------------------------------------------------------


def test_band_edges_none_returns_none():
    summary = {
        "mean": np.array([1.0, 1.0]),
        "q05": np.array([0.5, 0.5]),
        "q25": np.array([0.8, 0.8]),
        "q75": np.array([1.2, 1.2]),
        "q95": np.array([1.5, 1.5]),
        "std": np.array([0.2, 0.2]),
        "n_trials": np.array([10, 10]),
    }
    assert _band_edges(summary, "none") is None


def test_band_edges_iqr_picks_q25_q75():
    summary = {
        "mean": np.array([1.0]),
        "q05": np.array([0.0]),
        "q25": np.array([0.7]),
        "q75": np.array([1.3]),
        "q95": np.array([2.0]),
        "std": np.array([0.5]),
        "n_trials": np.array([10]),
    }
    lo, hi = _band_edges(summary, "iqr")
    assert lo[0] == 0.7
    assert hi[0] == 1.3


def test_band_edges_wide_quantile_picks_q05_q95():
    summary = {
        "mean": np.array([1.0]),
        "q05": np.array([0.0]),
        "q25": np.array([0.7]),
        "q75": np.array([1.3]),
        "q95": np.array([2.0]),
        "std": np.array([0.5]),
        "n_trials": np.array([10]),
    }
    lo, hi = _band_edges(summary, "wide-quantile")
    assert lo[0] == 0.0
    assert hi[0] == 2.0


def test_band_edges_std_picks_mean_pm_std():
    summary = {
        "mean": np.array([1.0]),
        "q05": np.array([0.0]),
        "q25": np.array([0.7]),
        "q75": np.array([1.3]),
        "q95": np.array([2.0]),
        "std": np.array([0.5]),
        "n_trials": np.array([10]),
    }
    lo, hi = _band_edges(summary, "std")
    assert lo[0] == 0.5
    assert hi[0] == 1.5


def test_band_edges_unknown_mode_raises():
    summary = {
        "mean": np.array([1.0]),
        "q05": np.array([0.0]), "q25": np.array([0.7]),
        "q75": np.array([1.3]), "q95": np.array([2.0]),
        "std": np.array([0.5]), "n_trials": np.array([10]),
    }
    with pytest.raises(ValueError):
        _band_edges(summary, "made-up-mode")


def test_band_edges_sparse_bin_masked():
    """Grid points with n_trials < MIN_TRIALS_FOR_BAND are NaN-masked."""
    summary = {
        "mean": np.array([1.0, 1.0, 1.0]),
        "q05": np.array([0.0, 0.0, 0.0]),
        "q25": np.array([0.7, 0.7, 0.7]),
        "q75": np.array([1.3, 1.3, 1.3]),
        "q95": np.array([2.0, 2.0, 2.0]),
        "std": np.array([0.5, 0.5, 0.5]),
        "n_trials": np.array([10, 1, MIN_TRIALS_FOR_BAND]),  # middle is sparse
    }
    lo, hi = _band_edges(summary, "iqr")
    assert lo[0] == 0.7  # dense
    assert np.isnan(lo[1])  # sparse middle masked
    assert lo[2] == 0.7  # at threshold (>= MIN_TRIALS_FOR_BAND) — drawn
    assert np.isnan(hi[1])


# ----------------------------------------------------------------------
# _marginal_curve_with_spread — empty-input guard
# ----------------------------------------------------------------------


def test_marginal_curve_with_spread_empty_dataframe_returns_nans():
    """An empty cond_df (no motion bins for that condition) → all NaN."""
    cfg = GLMConfig()
    sf_ref = np.asarray(cfg.sf_levels, dtype=np.float64)
    or_ref = np.asarray(cfg.or_levels, dtype=np.float64)
    grid = np.linspace(0.0, 50.0, 10)
    empty_df = pd.DataFrame(
        columns=["trial_id", "speed", "tf", "sf", "orientation",
                 "time_since_onset"],
    )
    out = _marginal_curve_with_spread(
        beta=np.zeros(1),
        train_names=[],
        selected_vars=[],
        config=cfg,
        sf_ref=sf_ref,
        or_ref=or_ref,
        cond_df=empty_df,
        sweep_var="Speed",
        grid=grid,
    )
    assert out["n_trials"].sum() == 0
    assert np.all(np.isnan(out["mean"]))
    assert np.all(np.isnan(out["q25"]))
