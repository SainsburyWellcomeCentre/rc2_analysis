"""Cross-validation and forward-selection tests on synthetic data."""

from __future__ import annotations

import numpy as np

from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.forward_selection import forward_select


def test_make_trial_folds_partitions_trials_round_robin():
    trial_ids = np.repeat(np.arange(10), 5)
    folds = make_trial_folds(trial_ids, n_folds=5, seed=0)
    # Every trial gets exactly one fold assignment
    by_trial = {int(t): set(folds[trial_ids == t]) for t in np.unique(trial_ids)}
    for t, fset in by_trial.items():
        assert len(fset) == 1, f"trial {t} spans folds {fset}"
    # All 5 folds used
    assert set(folds.tolist()) == {0, 1, 2, 3, 4}


def test_make_trial_folds_stratifies_by_condition_when_provided():
    trial_ids = np.array([1, 1, 2, 2, 3, 3, 4, 4])
    conditions = np.array([
        "stationary", "VT", "stationary", "VT",
        "stationary", "VT", "stationary", "VT",
    ], dtype=object)

    folds = make_trial_folds(
        trial_ids,
        n_folds=2,
        seed=0,
        condition_labels_per_bin=conditions,
    )

    by_pair = {
        (int(t), str(c)): set(folds[(trial_ids == t) & (conditions == c)])
        for t, c in zip(trial_ids, conditions)
    }
    for pair, fset in by_pair.items():
        assert len(fset) == 1, f"pair {pair} spans folds {fset}"

    stationary_folds = {next(iter(fset)) for pair, fset in by_pair.items() if pair[1] == "stationary"}
    vt_folds = {next(iter(fset)) for pair, fset in by_pair.items() if pair[1] == "VT"}
    assert stationary_folds == {0, 1}
    assert vt_folds == {0, 1}


def _synth_poisson_with_trials(n_trials=20, n_per_trial=50, seed=0):
    rng = np.random.default_rng(seed)
    n = n_trials * n_per_trial
    X = np.column_stack([np.ones(n), rng.normal(size=n)])
    eta = -1.0 + 0.7 * X[:, 1]
    y = rng.poisson(np.exp(eta))
    trial_ids = np.repeat(np.arange(n_trials), n_per_trial)
    return X, y, trial_ids


def test_cv_bps_above_null_for_good_predictor():
    X, y, trial_ids = _synth_poisson_with_trials(seed=1)
    folds = make_trial_folds(trial_ids, n_folds=5, seed=0)
    cv = cross_validate_glm(X, y, 0.0, folds, lambda_ridge=0.0)
    cv_null = cross_validate_glm(X[:, :1], y, 0.0, folds, lambda_ridge=0.0)
    assert cv.cv_bits_per_spike > cv_null.cv_bits_per_spike


def test_forward_selection_selects_speed_when_only_speed_drives_y():
    """Synthesise spike counts driven by speed-bases only; selection
    should pick Speed and reject pure noise predictors (TF/SF/OR)."""
    rng = np.random.default_rng(7)
    n_trials, n_per = 40, 60
    n = n_trials * n_per
    speed = rng.uniform(0, 50, n)
    tf = rng.uniform(0, 5, n)
    sf = rng.choice([0.003, 0.006, 0.012], n)
    or_v = rng.choice([0.0, np.pi / 4, np.pi / 2], n)
    onset_t = np.tile(np.linspace(0, 2, n_per), n_trials)

    from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis
    B_speed = raised_cosine_basis(speed, 5, 0.0, 50.0)
    B_tf = raised_cosine_basis(tf, 5, 0.0, 7.3)
    B_onset = onset_kernel_basis(onset_t, 6, 2.0)

    # True log-rate depends only on speed (via first speed basis)
    eta = -1.5 + 1.2 * B_speed[:, 0] + 0.8 * B_speed[:, 1]
    y = rng.poisson(np.exp(eta)).astype(np.float64)
    trial_ids = np.repeat(np.arange(n_trials), n_per)
    folds = make_trial_folds(trial_ids, n_folds=5, seed=0)

    res = forward_select(
        B_speed, B_tf, B_onset, sf, or_v, y, 0.0, folds,
    )
    assert "Speed" in res.selected_vars, res.selected_vars
    assert res.final_cv_bps > res.null_cv_bps
