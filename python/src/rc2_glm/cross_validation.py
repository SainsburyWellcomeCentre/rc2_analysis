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
    trial_ids_per_bin: np.ndarray,
    n_folds: int = 5,
    seed: int = 0,
    condition_labels_per_bin: np.ndarray | None = None,
    strategy: str = "condition-stratified",
    profile_ids_per_bin: np.ndarray | None = None,
) -> np.ndarray:
    """Assign each row to a fold.

    ``strategy="condition-stratified"`` (default): k-fold assignment over
    unique ``(trial_id, condition)`` pairs, stratified independently by
    condition — matches MATLAB's ``trial_fold`` loop in
    ``glm_single_cluster_analysis.m:1479-1484``.

    ``strategy="speed-profile"``: **2-fold** assignment where fold id =
    ``profile_id - 1`` regardless of ``n_folds`` — every bin belonging
    to trials of profile 1 goes to fold 0, profile 2 goes to fold 1.
    Train-on-one-profile / test-on-the-other generalisation check,
    mirroring MATLAB's ``sp_fold`` in
    ``glm_single_cluster_analysis.m:2291-2293``. Requires
    ``profile_ids_per_bin``.
    """
    trial_ids_per_bin = np.asarray(trial_ids_per_bin)

    if strategy == "speed-profile":
        if profile_ids_per_bin is None:
            raise ValueError(
                "strategy='speed-profile' requires profile_ids_per_bin; "
                "ensure rc2-glm is run with a StimulusLookup so TrialData.profile_id "
                "is populated and time_binning threads it through."
            )
        profile_ids = np.asarray(profile_ids_per_bin, dtype=np.int64).ravel()
        if profile_ids.shape[0] != trial_ids_per_bin.shape[0]:
            raise ValueError(
                "profile_ids_per_bin must match trial_ids_per_bin length"
            )
        # Per-trial consistency check: every bin of one trial must carry
        # the same profile_id. If profile_id varies within a trial, the
        # upstream plumbing has a bug.
        for tid in np.unique(trial_ids_per_bin):
            uniq = np.unique(profile_ids[trial_ids_per_bin == tid])
            if uniq.size != 1:
                raise ValueError(
                    f"trial_id={int(tid)} has inconsistent profile_ids "
                    f"{uniq.tolist()} — check time_binning output."
                )
        return (profile_ids - 1).astype(np.int64)

    if strategy != "condition-stratified":
        raise ValueError(
            f"unknown strategy {strategy!r}; "
            "valid: 'condition-stratified', 'speed-profile'"
        )

    rng = np.random.default_rng(seed)

    if condition_labels_per_bin is None:
        unique_trials = np.unique(trial_ids_per_bin)
        perm = rng.permutation(unique_trials)
        trial_to_fold = {int(t): int(i % n_folds) for i, t in enumerate(perm)}
        return np.array([trial_to_fold[int(t)] for t in trial_ids_per_bin], dtype=np.int64)

    conditions = np.asarray(condition_labels_per_bin, dtype=object)
    if conditions.shape[0] != trial_ids_per_bin.shape[0]:
        raise ValueError("condition_labels_per_bin must match trial_ids_per_bin length")

    pair_to_fold: dict[tuple[int, str], int] = {}
    for condition in np.unique(conditions):
        cond_mask = conditions == condition
        unique_trials = np.unique(trial_ids_per_bin[cond_mask])
        perm = rng.permutation(unique_trials)
        for i, trial_id in enumerate(perm):
            pair_to_fold[(int(trial_id), str(condition))] = int(i % n_folds)

    return np.array(
        [pair_to_fold[(int(trial_id), str(condition))]
         for trial_id, condition in zip(trial_ids_per_bin, conditions)],
        dtype=np.int64,
    )


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
