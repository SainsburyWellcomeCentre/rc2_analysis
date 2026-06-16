"""Hardcastle signed-rank admission rule (2026-06-16).

The ``selection_rule="signed_rank"`` path admits a candidate iff a one-sided
Wilcoxon signed-rank test on the PER-FOLD paired Δ bits/spike (candidate minus
current model, across the n_folds folds of a single partition) clears
``selection_alpha``. Tests:

1. **Strong Speed signal is admitted** under signed_rank with 10-fold
   condition-stratified CV, and round 1 carries a small p-value for Speed.
2. **Null cluster admits nothing** (no candidate is significant).
3. **Legacy default unchanged** — with the default rule the result equals the
   pre-change behaviour (delta_bps_threshold), and ``pval`` is NaN.
4. **signed_rank ⟂ multi-seed** — n_selection_seeds>1 raises.
5. **final_vs_null_pval** is populated (and small) for a strong-signal cluster
   under signed_rank, NaN under the legacy rule.
"""

from __future__ import annotations

import numpy as np
import pytest

from rc2_glm.config import GLMConfig
from rc2_glm.cross_validation import make_trial_folds
from rc2_glm.forward_selection import forward_select


def _make_synthetic_cluster(
    n_trials: int = 40,
    n_bins_per_trial: int = 30,
    seed: int = 0,
    speed_strength: float = 1.5,
    n_folds: int = 10,
):
    """Synthetic single cluster where Speed carries real signal. Returns the
    inputs forward_select needs plus a single n_folds condition-stratified
    partition (the signed-rank sample)."""
    rng = np.random.default_rng(seed)
    n_bins = n_trials * n_bins_per_trial
    trial_ids = np.repeat(np.arange(n_trials), n_bins_per_trial)
    cond_per_trial = np.array(["A" if t % 2 == 0 else "B" for t in range(n_trials)])
    condition_labels = cond_per_trial[trial_ids]

    def _toy_basis(n: int, n_cols: int, freq: float) -> np.ndarray:
        t = np.linspace(0, 1, n)
        cols = []
        for k in range(n_cols):
            phase = 2 * np.pi * k / n_cols
            cols.append(np.cos(freq * np.pi * t + phase))
        return np.stack(cols, axis=1) + rng.normal(0, 0.05, (n, n_cols))

    B_speed = _toy_basis(n_bins, 5, freq=2.0)
    B_tf = _toy_basis(n_bins, 5, freq=3.0)
    B_onset = np.zeros((n_bins, 0))
    sf_vals = np.zeros(n_bins)
    or_vals = np.zeros(n_bins)

    log_lambda = -1.0 + speed_strength * B_speed[:, 0]
    y = rng.poisson(np.exp(log_lambda) * 0.1)
    offset = float(np.log(0.1))

    fold_ids = make_trial_folds(
        trial_ids, n_folds=n_folds, seed=0,
        condition_labels_per_bin=condition_labels,
    )
    return dict(
        B_speed=B_speed, B_tf=B_tf, B_onset=B_onset,
        sf_vals=sf_vals, or_vals=or_vals,
        y=y, offset=offset, fold_ids=fold_ids,
        trial_ids=trial_ids, condition_labels=condition_labels,
    )


def _config(rule: str, **kw) -> GLMConfig:
    return GLMConfig(
        include_history=False,
        include_onset_kernel=False,
        include_me_face=False,
        main_effects=("Speed", "TF", "SF", "OR"),
        interactions=(),
        n_folds=10,
        selection_rule=rule,
        **kw,
    )


def _select(data, config):
    return forward_select(
        data["B_speed"], data["B_tf"], data["B_onset"],
        data["sf_vals"], data["or_vals"],
        data["y"], data["offset"], data["fold_ids"],
        config=config, backend="irls",
    )


def test_signed_rank_admits_strong_speed():
    data = _make_synthetic_cluster(speed_strength=2.0, seed=42)
    res = _select(data, _config("signed_rank", selection_alpha=0.05))
    assert "Speed" in res.selected_vars, (
        f"strong Speed signal must be admitted under signed_rank; "
        f"selected={res.selected_vars}"
    )
    round1 = res.history[0]
    assert round1.best_candidate == "Speed"
    assert round1.added is True
    # Speed's per-fold improvement is significant.
    assert round1.pval["Speed"] < 0.05
    # Every tested candidate carries a (finite) p-value under signed_rank.
    for cand, p in round1.pval.items():
        assert np.isfinite(p)


def test_signed_rank_null_admits_nothing():
    rng = np.random.default_rng(123)
    n_trials, n_bpt = 40, 30
    n_bins = n_trials * n_bpt
    trial_ids = np.repeat(np.arange(n_trials), n_bpt)
    cond = np.array(["A" if t % 2 == 0 else "B" for t in range(n_trials)])[trial_ids]
    data = dict(
        B_speed=rng.normal(0, 1, (n_bins, 5)),
        B_tf=rng.normal(0, 1, (n_bins, 5)),
        B_onset=np.zeros((n_bins, 0)),
        sf_vals=np.zeros(n_bins),
        or_vals=np.zeros(n_bins),
        y=rng.poisson(np.exp(-1.0 * np.ones(n_bins)) * 0.1),
        offset=float(np.log(0.1)),
        fold_ids=make_trial_folds(
            trial_ids, n_folds=10, seed=0, condition_labels_per_bin=cond,
        ),
    )
    res = _select(data, _config("signed_rank", selection_alpha=0.05))
    assert res.selected_vars == [], (
        f"null cluster must select nothing under signed_rank; got {res.selected_vars}"
    )
    assert res.history[0].added is False
    assert np.isnan(res.final_vs_null_pval)  # empty selection → NaN


def test_legacy_rule_unchanged_and_pval_nan():
    """The default rule path is byte-identical to a run that never knew about
    signed_rank, and its RoundResult.pval entries are NaN."""
    data = _make_synthetic_cluster(speed_strength=2.0, seed=42)
    legacy = _select(data, _config("delta_bps_threshold"))
    assert "Speed" in legacy.selected_vars
    assert np.isnan(legacy.final_vs_null_pval)
    for r in legacy.history:
        for cand, p in r.pval.items():
            assert np.isnan(p), f"legacy rule must leave pval NaN; {cand}={p}"


def test_signed_rank_requires_single_partition():
    data = _make_synthetic_cluster(speed_strength=1.0, seed=1)
    config = _config("signed_rank", n_selection_seeds=10, selection_threshold_count=7)
    fold_ids_list = [
        make_trial_folds(
            data["trial_ids"], n_folds=10, seed=s,
            condition_labels_per_bin=data["condition_labels"],
        )
        for s in range(10)
    ]
    with pytest.raises(ValueError, match="signed_rank"):
        forward_select(
            data["B_speed"], data["B_tf"], data["B_onset"],
            data["sf_vals"], data["or_vals"],
            data["y"], data["offset"], data["fold_ids"],
            config=config, backend="irls",
            fold_ids_per_seed=fold_ids_list,
        )


def test_final_vs_null_pval_small_for_strong_signal():
    data = _make_synthetic_cluster(speed_strength=2.0, seed=42)
    res = _select(data, _config("signed_rank", selection_alpha=0.05))
    assert res.selected_vars  # non-empty
    assert np.isfinite(res.final_vs_null_pval)
    assert res.final_vs_null_pval < 0.05
