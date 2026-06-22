"""Tests for rc2_glm.tuning_significance.

Cover the contract: a real monotonic/peaked tuning is detected with small p and
a sensible best model; flat noise is not significant (p ~ U(0,1)); the bootstrap
is deterministic under a fixed seed; degenerate inputs return p=NaN without
crashing. These run on synthetic per-trial-per-bin matrices (no GLM/data deps).
"""

import numpy as np
import pandas as pd
import pytest

from rc2_glm import tuning_significance as ts


def _matrix_from_curve(curve, n_trials=30, noise=0.5, seed=0):
    """(n_trials × n_bins) matrix = mean curve + Gaussian per-trial noise."""
    rng = np.random.default_rng(seed)
    return curve[None, :] + rng.normal(0.0, noise, size=(n_trials, curve.size))


# --------------------------------------------------------------------------- #
# Real tuning is detected.
# --------------------------------------------------------------------------- #
def test_linear_tuning_significant():
    x = np.linspace(0.0, 7.0, 20)
    matrix = _matrix_from_curve(2.0 + 1.5 * x, noise=0.5, seed=1)
    r = ts.tuning_significance(matrix, x, value="tf", condition="V",
                               kind="linear", n_reps=500, seed=1)
    assert r.best_model is not None
    assert r.rsq > 0.9
    assert r.p < 0.05


def test_flat_aggregate_lower_rsq_but_detects_linear():
    """MATLAB 'flat' R² (across all per-trial points) is much smaller than the
    mean-curve R², but a real linear trend is still significant."""
    x = np.linspace(0.0, 7.0, 20)
    matrix = _matrix_from_curve(2.0 + 1.5 * x, noise=2.0, seed=11)
    r_mean = ts.tuning_significance(matrix, x, value="tf", condition="V",
                                    kind="linear", aggregate="mean",
                                    linear_families=("linear",), n_reps=500, seed=1)
    r_flat = ts.tuning_significance(matrix, x, value="tf", condition="V",
                                    kind="linear", aggregate="flat",
                                    linear_families=("linear",), n_reps=500, seed=1)
    assert r_flat.aggregate == "flat"
    assert r_flat.rsq < r_mean.rsq          # trial scatter drags flat R² down
    assert r_flat.p < 0.05                  # real trend still beats the null
    # rsq_mean (fit-vs-mean-curve) is reported regardless of the aggregate and is
    # much higher than the across-trials flat R² for a clean trend.
    assert np.isfinite(r_flat.rsq_mean)
    assert r_flat.rsq_mean > r_flat.rsq
    assert np.isclose(r_mean.rsq_mean, r_mean.rsq, atol=1e-9)  # mean agg: same number


def test_select_criterion_rsq_mean_picks_max_mean_r2():
    """rsq_mean selection chooses the family with the highest R² on the mean
    curve (the MATLAB ModelSelectionTuning rule), not the lowest BIC."""
    x = np.linspace(0.0, 0.12, 20)
    curve = 5.0 * np.exp(-((x - 0.06) ** 2) / (2 * 0.02 ** 2)) + 1.0
    matrix = _matrix_from_curve(curve, noise=0.4, seed=21)
    r = ts.tuning_significance(matrix, x, value="sf", condition="V", kind="linear",
                               aggregate="flat", select_criterion="rsq_mean",
                               n_reps=300, seed=1)
    assert r.select_criterion == "rsq_mean"
    # The chosen family must have the best (max) mean-curve R² among all families.
    fits = ts.fit_tuning(x, np.nanmean(matrix, 0), kind="linear")
    best_mean = max(ts.rsq_against_mean(matrix, x, f["name"], f["params"])
                    for f in fits if f["ok"])
    assert np.isclose(r.rsq_mean, best_mean, atol=1e-9)


def test_gaussian_tuning_detected_and_low_p():
    x = np.linspace(0.0, 0.12, 20)
    curve = 5.0 * np.exp(-((x - 0.06) ** 2) / (2 * 0.02 ** 2)) + 1.0
    matrix = _matrix_from_curve(curve, noise=0.3, seed=2)
    r = ts.tuning_significance(matrix, x, value="sf", condition="V",
                               kind="linear", n_reps=500, seed=1)
    assert r.p < 0.05
    # BIC should not prefer a 1-param-cheaper model that misses the peak badly.
    assert r.rsq > 0.8


def test_orientation_vonmises():
    x = np.linspace(0.0, 170.0, 18)
    curve = 4.0 * np.exp(2.0 * (np.cos(np.deg2rad(2 * (x - 90.0))) - 1.0)) + 1.0
    matrix = _matrix_from_curve(curve, noise=0.3, seed=3)
    r = ts.tuning_significance(matrix, x, value="or", condition="V",
                               kind="circular", n_reps=500, seed=1)
    assert r.best_model == "vonmises_180"
    assert r.p < 0.05


# --------------------------------------------------------------------------- #
# Flat noise is not significant.
# --------------------------------------------------------------------------- #
def test_flat_noise_not_significant():
    x = np.linspace(0.0, 7.0, 20)
    matrix = _matrix_from_curve(np.full(20, 3.0), noise=1.0, seed=4)
    r = ts.tuning_significance(matrix, x, value="tf", condition="V",
                               kind="linear", n_reps=1000, seed=1)
    # No real shape → observed R² sits inside the null → not significant.
    assert np.isnan(r.p) or r.p > 0.05


# --------------------------------------------------------------------------- #
# Determinism.
# --------------------------------------------------------------------------- #
def test_seed_deterministic():
    x = np.linspace(0.0, 7.0, 20)
    matrix = _matrix_from_curve(2.0 + 0.8 * x, noise=0.8, seed=5)
    r1 = ts.tuning_significance(matrix, x, value="tf", condition="V",
                                kind="linear", n_reps=400, seed=1)
    r2 = ts.tuning_significance(matrix, x, value="tf", condition="V",
                                kind="linear", n_reps=400, seed=1)
    assert r1.p == r2.p
    assert np.allclose(r1.null_rsq, r2.null_rsq)


# --------------------------------------------------------------------------- #
# Degenerate inputs → p=NaN, no crash.
# --------------------------------------------------------------------------- #
def test_empty_matrix_nan():
    r = ts.tuning_significance(None, np.linspace(0, 1, 20), value="tf",
                               condition="V", kind="linear")
    assert np.isnan(r.p) and r.best_model is None


def test_too_few_trials_nan():
    x = np.linspace(0.0, 7.0, 20)
    matrix = _matrix_from_curve(2.0 + x, noise=0.2, seed=6, n_trials=2)
    r = ts.tuning_significance(matrix, x, value="tf", condition="V",
                               kind="linear")
    assert np.isnan(r.p) and r.best_model is None


def test_flat_mean_curve_tss_zero_nan():
    x = np.linspace(0.0, 7.0, 20)
    matrix = np.full((10, 20), 3.0)  # identical every bin/trial → TSS=0
    r = ts.tuning_significance(matrix, x, value="tf", condition="V",
                               kind="linear", n_reps=100, seed=1)
    assert np.isnan(r.p)


# --------------------------------------------------------------------------- #
# per_trial_bin_matrix.
# --------------------------------------------------------------------------- #
def test_per_trial_bin_matrix_shapes_and_condition():
    rng = np.random.default_rng(7)
    n = 4000
    df = pd.DataFrame({
        "condition": rng.choice(["V", "VT"], size=n),
        "trial_id": rng.integers(0, 20, size=n),
        "tf": rng.uniform(0.0, 7.0, size=n),
        "spike_count": rng.poisson(1.0, size=n).astype(float),
    })
    mat, cen = ts.per_trial_bin_matrix(df, "tf", bw=0.02, condition="V", n_bins=20)
    assert mat is not None
    assert mat.shape[1] == 20 and cen.size == 20
    # VT path independent; missing column → (None, None).
    assert ts.per_trial_bin_matrix(df, "nope", bw=0.02, condition="V")[0] is None


def test_per_trial_bin_matrix_missing_condition():
    df = pd.DataFrame({
        "condition": ["V"] * 100,
        "trial_id": list(range(100)),
        "tf": np.linspace(0, 7, 100),
        "spike_count": np.ones(100),
    })
    assert ts.per_trial_bin_matrix(df, "tf", bw=0.02, condition="VT")[0] is None
