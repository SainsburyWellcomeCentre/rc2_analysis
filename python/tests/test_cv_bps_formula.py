"""Nail down the CV-bits-per-spike formula with a hand-computed Poisson example.

Both Python and MATLAB compute, per test fold::

    ll_test = sum_i ( y_i * log(mu_i) - mu_i - gammaln(y_i + 1) )

then aggregate pooled::

    cv_bps = (sum_folds ll_test) / (sum_folds n_spikes_test) / log(2)

This keeps the Poisson ``-gammaln(y+1)`` term; absolute bps is therefore
bin-width dependent. This test locks the convention.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

from rc2_glm.cross_validation import cross_validate_glm


def _hand_bps(y, mu_per_fold_per_bin, fold_ids):
    total_ll = 0.0
    total_spikes = 0.0
    for f, mu_vec in mu_per_fold_per_bin.items():
        test = fold_ids == f
        y_t = y[test].astype(float)
        mu_t = mu_vec
        total_ll += float(np.sum(y_t * np.log(mu_t) - mu_t - gammaln(y_t + 1.0)))
        total_spikes += float(y_t.sum())
    return (total_ll / total_spikes) / np.log(2)


def test_cv_bps_matches_hand_computed_poisson_3bins_2folds():
    # 3 rows, 2 folds: row 0 in fold 0, rows 1 and 2 in fold 1.
    # Use a single predictor x so fits are deterministic and easy to reason
    # about; we do not rely on the fitter's exact beta — we check that the
    # aggregation + LL formula recovers the hand-computed value.
    X = np.array([[1.0, 0.0],
                  [1.0, 1.0],
                  [1.0, 2.0]])
    y = np.array([1.0, 2.0, 0.0])
    offset = np.zeros(3)
    fold_ids = np.array([0, 1, 1])

    cv = cross_validate_glm(X, y, offset, fold_ids, lambda_ridge=0.0)

    # Reconstruct mu per fold from cv.cv_predicted_count
    mu_by_fold = {
        0: cv.cv_predicted_count[fold_ids == 0],
        1: cv.cv_predicted_count[fold_ids == 1],
    }
    hand = _hand_bps(y, mu_by_fold, fold_ids)

    assert np.isfinite(cv.cv_bits_per_spike)
    assert np.isclose(cv.cv_bits_per_spike, hand, atol=1e-12), (
        cv.cv_bits_per_spike, hand,
    )
    # Also: the pooled LL ratio equals (sum of fold LLs) / sum(y_test).
    assert np.isclose(
        cv.cv_bits_per_spike,
        (cv.fold_log_likelihood.sum() / y.sum()) / np.log(2),
        atol=1e-12,
    )


def test_cv_bps_keeps_gammaln_term():
    # Same inputs as above, but verify the LL uses -gammaln(y+1).
    # If a future change drops the constant term the absolute bps shifts by
    # sum(gammaln(y+1)) / sum(y) / log(2); this test catches that regression.
    X = np.array([[1.0, 0.0],
                  [1.0, 1.0],
                  [1.0, 2.0]])
    y = np.array([1.0, 2.0, 0.0])
    offset = np.zeros(3)
    fold_ids = np.array([0, 1, 1])

    cv = cross_validate_glm(X, y, offset, fold_ids, lambda_ridge=0.0)
    mu = cv.cv_predicted_count

    # Hand LL WITH gammaln
    ll_with = float(np.sum(y * np.log(mu) - mu - gammaln(y + 1.0)))
    # Hand LL WITHOUT gammaln (the wrong formula we're guarding against)
    ll_without = float(np.sum(y * np.log(mu) - mu))

    bps_with = (ll_with / y.sum()) / np.log(2)
    bps_without = (ll_without / y.sum()) / np.log(2)

    # Implementation matches the "with gammaln" form.
    assert np.isclose(cv.cv_bits_per_spike, bps_with, atol=1e-9)
    # And the two forms do differ here (y=2 contributes gammaln(3)=log(2) ≠ 0),
    # so this is a real assertion rather than a tautology. The exact gap is
    # sum(gammaln(y+1)) / sum(y) / log(2) = log(2)/3/log(2) = 1/3 bits/spike.
    expected_gap = float(np.sum(gammaln(y + 1.0))) / y.sum() / np.log(2)
    assert np.isclose(bps_without - bps_with, expected_gap, atol=1e-9), (
        bps_without - bps_with, expected_gap,
    )
    assert abs(bps_with - bps_without) > 0.1


def test_cv_bps_scales_with_bin_width_convention():
    # Two bins aggregated into one bin (same total spikes, doubled rate) gives
    # a different absolute bps. This pins the "bin-width dependent" property
    # and documents that Python↔MATLAB absolute bps offsets can come purely
    # from a different bin-width choice.
    y_fine = np.array([1.0, 1.0, 0.0, 0.0])          # 4 bins, total 2 spikes
    y_coarse = np.array([2.0, 0.0])                   # 2 bins, total 2 spikes
    mu_fine = np.array([0.5, 0.5, 0.5, 0.5])          # rate 0.5 / fine bin
    mu_coarse = np.array([1.0, 1.0])                  # rate 1.0 / coarse bin

    def bps(y, mu):
        ll = float(np.sum(y * np.log(mu) - mu - gammaln(y + 1.0)))
        return (ll / y.sum()) / np.log(2)

    fine = bps(y_fine, mu_fine)
    coarse = bps(y_coarse, mu_coarse)
    # If bps were bin-width invariant these would be equal; they are not.
    assert not np.isclose(fine, coarse, atol=0.1), (fine, coarse)
