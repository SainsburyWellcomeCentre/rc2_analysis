"""Trial-level k-fold cross-validation and bits-per-spike.

Mirrors MATLAB ``cross_validate_glm`` (line 5861). Splits along
fold_ids (entire trials per fold), fits on train, scores on test, and
reports

    bits_per_spike = (sum_test_log_likelihood / total_test_spikes) / log(2)

This matches the MATLAB convention where the null model is the empty
likelihood (i.e. log(2) normalisation only). The "Δ bps vs null"
comparison used by forward selection then computes
``bps(model) - bps(null_model)`` which is the standard
Hardcastle bits-per-spike improvement.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from rc2_glm.fitting import fit_poisson_glm


@dataclass
class CVResult:
    cv_log_likelihood: float
    cv_bits_per_spike: float
    fold_log_likelihood: np.ndarray
    cv_predicted_count: np.ndarray


def make_trial_folds(
    trial_ids_per_bin: np.ndarray, n_folds: int = 5, seed: int = 0
) -> np.ndarray:
    """Assign each row to a fold based on its trial_id.

    Trials are shuffled with `seed` and dealt round-robin across folds,
    so every fold has comparable trial coverage.
    """
    trial_ids_per_bin = np.asarray(trial_ids_per_bin)
    unique_trials = np.unique(trial_ids_per_bin)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(unique_trials)
    trial_to_fold = {int(t): int(i % n_folds) for i, t in enumerate(perm)}
    return np.array([trial_to_fold[int(t)] for t in trial_ids_per_bin], dtype=np.int64)


def cross_validate_glm(
    X: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float,
    fold_ids: np.ndarray,
    lambda_ridge: float = 0.0,
    backend: str = "irls",
    eta_clip: float = 20.0,
    mu_floor: float = 1e-10,
) -> CVResult:
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    fold_ids = np.asarray(fold_ids).ravel()

    if np.isscalar(offset):
        offset_arr = np.full(y.size, float(offset))
    else:
        offset_arr = np.asarray(offset, dtype=np.float64).ravel()
        if offset_arr.size == 1:
            offset_arr = np.full(y.size, float(offset_arr[0]))

    folds = np.unique(fold_ids)
    fold_ll = np.zeros(folds.size)
    cv_predicted = np.full(y.size, np.nan)
    total_ll = 0.0
    total_spikes = 0.0

    for fi, f in enumerate(folds):
        test_idx = fold_ids == f
        train_idx = ~test_idx

        if not test_idx.any() or not train_idx.any():
            continue

        fit = fit_poisson_glm(
            X[train_idx], y[train_idx], offset_arr[train_idx],
            lambda_ridge=lambda_ridge, backend=backend,
        )

        eta = np.clip(X[test_idx] @ fit.beta + offset_arr[test_idx],
                      -eta_clip, eta_clip)
        mu = np.maximum(np.exp(eta), mu_floor)
        cv_predicted[test_idx] = mu

        y_test = y[test_idx]
        ll_test = float(np.sum(y_test * np.log(mu) - mu - gammaln(y_test + 1.0)))
        fold_ll[fi] = ll_test
        total_ll += ll_test
        total_spikes += float(y_test.sum())

    bps = (total_ll / total_spikes) / np.log(2) if total_spikes > 0 else float("nan")
    return CVResult(
        cv_log_likelihood=total_ll,
        cv_bits_per_spike=bps,
        fold_log_likelihood=fold_ll,
        cv_predicted_count=cv_predicted,
    )
