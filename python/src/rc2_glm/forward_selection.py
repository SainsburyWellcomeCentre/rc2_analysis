"""Hardcastle-style hierarchical forward selection.

Two phases (mirrors MATLAB ``forward_select_model`` at line 5903):

- **Phase 1 — Main effects.** Test {Speed, TF, SF, OR} one at a time.
  Add the best if Δ CV bps > threshold. Repeat until no main effect
  passes.
- **Phase 2 — Interactions.** Only test interactions whose BOTH
  parent main effects were selected in Phase 1. Same Δ-bps rule.

Returns the selected variable list, a per-round history, and the
final / null CV bits-per-spike.

Multi-seed admission (added 2026-05-08, prompt 13). When
``fold_ids_per_seed`` is provided with N > 1 partitions, each
candidate's Δ-bps is computed under each of the N partitions and a
candidate is admitted iff at least ``selection_threshold_count`` of N
exceed ``delta_bps_threshold``. The "best" candidate of a passing
round is the admitted one with the highest mean Δ across partitions.
N=1 (the default) reduces exactly to single-seed behaviour.

Signed-rank admission (added 2026-06-16). When
``config.selection_rule == "signed_rank"`` a candidate is admitted iff a
one-sided Wilcoxon signed-rank test on the PER-FOLD paired Δ bits/spike
(candidate minus current model, across the n_folds folds of a single
partition) gives p < ``config.selection_alpha`` — Hardcastle et al. 2017,
Neuron. This replaces the fixed Δ-bps threshold with a per-fold
significance test and is mutually exclusive with multi-seed voting
(the folds ARE the test sample, so N must be 1).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy.stats import wilcoxon

from rc2_glm.config import GLMConfig, INTERACTION_PARENTS
from rc2_glm.cross_validation import CVResult, cross_validate_glm
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.penalty import build_penalty_matrix

logger = logging.getLogger(__name__)


def _signed_rank_greater(per_fold_delta: np.ndarray) -> float:
    """One-sided Wilcoxon signed-rank p-value for H1: median per-fold Δ > 0.

    The Hardcastle et al. 2017 admission statistic. NaN folds (zero held-out
    spikes) are dropped. Returns 1.0 when the test is undefined (all
    differences ~0, or too few non-zero paired differences for scipy)."""
    d = np.asarray(per_fold_delta, dtype=np.float64)
    d = d[np.isfinite(d)]
    if d.size < 1 or np.allclose(d, 0.0):
        return 1.0
    try:
        _, p = wilcoxon(d, alternative="greater", zero_method="wilcox")
    except ValueError:
        # scipy raises when every difference is zero or n is too small.
        return 1.0
    return float(p)


def _penalty_for(
    names: list[str],
    config: GLMConfig,
    history_basis_mat: np.ndarray | None,
) -> np.ndarray | None:
    """Build the smoothness penalty for a model's columns, or None.

    Returns None (→ caller falls back to scalar ridge) unless
    ``config.history_smooth_lambda`` is set and a history basis is available.
    """
    smooth = getattr(config, "history_smooth_lambda", None)
    if smooth is None or history_basis_mat is None:
        return None
    return build_penalty_matrix(
        names, config.lambda_ridge,
        lambda_min=config.lambda_ridge_min,
        history_smooth_lambda=smooth,
        history_basis_mat=history_basis_mat,
    )


@dataclass
class RoundResult:
    round: int
    phase: int
    tested: dict[str, float]              # candidate → mean cv_bps across seeds
    delta_bps: dict[str, float]           # candidate → mean Δ across seeds
    best_candidate: str | None
    best_delta_bps: float                 # mean Δ of best candidate
    added: bool
    cv_bps_after: float                   # mean across seeds of new current cv_bps
    # Multi-seed evidence (added 2026-05-08). For single-seed runs each
    # candidate's lists / counts are length 1 / value ∈ {0, 1}.
    delta_bps_per_seed: dict[str, list[float]] = field(default_factory=dict)
    admitted_count: dict[str, int] = field(default_factory=dict)
    n_seeds: int = 1
    # Per-candidate one-sided signed-rank p-value (added 2026-06-16). NaN for
    # the legacy delta_bps_threshold rule; populated only under signed_rank.
    pval: dict[str, float] = field(default_factory=dict)


@dataclass
class SelectionResult:
    selected_vars: list[str]
    history: list[RoundResult]
    null_cv_bps: float
    final_cv_bps: float
    null_cv: CVResult           # canonical (seed-0) null CV result
    final_cv: CVResult | None   # canonical (seed-0) final CV result
    # Signed-rank p-value of the assembled final model vs the null, under the
    # single fold partition (added 2026-06-16). NaN for the legacy rule or an
    # empty selection.
    final_vs_null_pval: float = float("nan")
    # The admission rule that produced this selection (added 2026-06-16), so
    # downstream CSV rows self-describe the gate without re-reading config.
    selection_rule: str = "delta_bps_threshold"


def forward_select(
    B_speed: np.ndarray,
    B_tf: np.ndarray,
    B_onset: np.ndarray,
    sf_vals: np.ndarray,
    or_vals: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float,
    fold_ids: np.ndarray,
    config: GLMConfig | None = None,
    backend: str = "irls",
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
    *,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
    B_accel: np.ndarray | None = None,
    B_sf: np.ndarray | None = None,
    B_or: np.ndarray | None = None,
    fold_ids_per_seed: list[np.ndarray] | None = None,
) -> SelectionResult:
    """Hardcastle-style hierarchical forward selection.

    Keyword-only basis arguments:

    - ``B_history``: per-cluster spike-history feature matrix. When
      provided AND ``config.include_history`` is True, ``"History"`` is
      added as a Phase-1 candidate (NOT through ``config.main_effects``
      — History bypasses that tuple by historical convention).
    - ``B_me_face``: per-bin ME_face raised-cosine basis. When ``None``
      (camera absent or ``--no-me-face``), ``"ME_face"`` is dropped from
      the Phase-1 candidate list even if it appears in
      ``config.main_effects``. ``"ME_face_x_Speed"`` Phase-2 eligibility
      auto-resolves via ``INTERACTION_PARENTS`` (it requires both parents
      to have been selected in Phase 1).
    - ``fold_ids_per_seed``: optional list of N independent fold-id
      arrays for multi-seed admission. When None or length 1, falls
      back to single-seed (the standard Hardcastle behaviour) using
      ``fold_ids``. When length > 1, each candidate is evaluated under
      every partition and admitted iff
      ``config.selection_threshold_count`` of N partitions clear
      ``config.delta_bps_threshold``. The canonical ``null_cv`` /
      ``final_cv`` returned in ``SelectionResult`` are the seed-0
      results (downstream coefficient extraction is single-seed).

    The onset-kernel inclusion is gated by ``config.include_onset_kernel``.
    """
    config = config or GLMConfig()
    include_onset = getattr(config, "include_onset_kernel", True)
    include_history = (
        getattr(config, "include_history", False) and B_history is not None
    )

    # Rebuild the (n_lag, n_bases) history basis for the smoothness penalty.
    # Same args as pipeline's convolve_history input, so columns align.
    history_basis_mat = None
    if getattr(config, "history_smooth_lambda", None) is not None and include_history:
        from rc2_glm.basis import history_basis
        history_basis_mat = history_basis(
            n_bases=config.n_history_bases,
            t_max_s=config.history_window_s,
            bin_width_s=config.time_bin_width,
        )

    # Resolve seed list. Single-seed (back-compat): wrap fold_ids in a
    # one-element list. Multi-seed: use the provided list as-is.
    if fold_ids_per_seed is None or len(fold_ids_per_seed) <= 1:
        seeds = [fold_ids if fold_ids_per_seed is None else fold_ids_per_seed[0]]
        threshold_count = 1
    else:
        seeds = list(fold_ids_per_seed)
        threshold_count = int(getattr(config, "selection_threshold_count", 1))
        n_seeds_cfg = int(getattr(config, "n_selection_seeds", 1))
        if len(seeds) != n_seeds_cfg:
            # Trust the explicit fold_ids_per_seed length over the config
            # value; pipeline.py is the canonical builder.
            n_seeds_cfg = len(seeds)
        if not (1 <= threshold_count <= len(seeds)):
            raise ValueError(
                f"selection_threshold_count={threshold_count} must lie in "
                f"[1, n_selection_seeds={len(seeds)}]"
            )
    n_seeds = len(seeds)

    selection_rule = getattr(config, "selection_rule", "delta_bps_threshold")
    if selection_rule == "signed_rank" and n_seeds > 1:
        raise ValueError(
            "selection_rule='signed_rank' requires a single fold partition "
            f"(n_selection_seeds=1); got n_seeds={n_seeds}. The folds are the "
            "signed-rank test sample, so multi-seed voting is redundant."
        )

    common_assembler_kwargs = dict(
        B_history=B_history if include_history else None,
        B_me_face=B_me_face,
        B_accel=B_accel,
        B_sf=B_sf,
        B_or=B_or,
        include_onset_kernel=include_onset,
    )

    # Null model: intercept (+ onset kernel if include_onset_kernel)
    X_null, null_names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, [],
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        **common_assembler_kwargs,
    )
    null_penalty = _penalty_for(null_names, config, history_basis_mat)
    null_cv_per_seed = [
        cross_validate_glm(
            X_null, y, offset, seed_folds,
            lambda_ridge=config.lambda_ridge, backend=backend,
            penalty_matrix=null_penalty,
        )
        for seed_folds in seeds
    ]
    null_cv = null_cv_per_seed[0]
    null_bps_per_seed = np.array(
        [cv.cv_bits_per_spike for cv in null_cv_per_seed]
    )

    selected: list[str] = []
    current_cv_per_seed = null_cv_per_seed
    current_bps_per_seed = null_bps_per_seed
    history: list[RoundResult] = []
    round_num = 0

    # ----- Phase 1: main effects -----
    remaining = list(config.main_effects)
    if "Acceleration" in remaining and B_accel is None:
        # include_acceleration off, or no acceleration column — drop it from
        # the candidate list (mirrors the ME_face guard below).
        remaining.remove("Acceleration")
    if "ME_face" in remaining and B_me_face is None:
        # Camera data absent for this probe / cluster, or --no-me-face
        # was set at the pipeline level. Drop ME_face from the candidate
        # list rather than testing-then-discarding (cheaper, and the
        # selection_history reflects only candidates we actually had data
        # for).
        remaining.remove("ME_face")
    if include_history:
        remaining.append("History")
    while remaining:
        round_num += 1
        round_result = _try_candidates(
            remaining, selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
            y, offset, seeds, current_bps_per_seed,
            phase=1, round_num=round_num, backend=backend, config=config,
            threshold_count=threshold_count, current_cv_list=current_cv_per_seed,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
            history_basis_mat=history_basis_mat,
            **common_assembler_kwargs,
        )
        history.append(round_result)
        if round_result.added and round_result.best_candidate is not None:
            selected.append(round_result.best_candidate)
            remaining.remove(round_result.best_candidate)
            # Refit & store CV per seed for the new model
            current_cv_per_seed = _cv_for_selected_per_seed(
                selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
                y, offset, seeds, backend, config=config,
                sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
                history_basis_mat=history_basis_mat,
                **common_assembler_kwargs,
            )
            current_bps_per_seed = np.array(
                [cv.cv_bits_per_spike for cv in current_cv_per_seed]
            )
        else:
            break

    # ----- Phase 2: interactions (only those with both parents selected) -----
    # History interactions are excluded by default
    # (config.allow_history_interactions=False); interaction list only
    # contains stimulus-variable interactions.
    eligible: list[str] = [
        name for name in config.interactions
        if all(parent in selected for parent in INTERACTION_PARENTS[name])
    ]
    while eligible:
        round_num += 1
        round_result = _try_candidates(
            eligible, selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
            y, offset, seeds, current_bps_per_seed,
            phase=2, round_num=round_num, backend=backend, config=config,
            threshold_count=threshold_count, current_cv_list=current_cv_per_seed,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
            history_basis_mat=history_basis_mat,
            **common_assembler_kwargs,
        )
        history.append(round_result)
        if round_result.added and round_result.best_candidate is not None:
            selected.append(round_result.best_candidate)
            eligible.remove(round_result.best_candidate)
            current_cv_per_seed = _cv_for_selected_per_seed(
                selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
                y, offset, seeds, backend, config=config,
                sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
                history_basis_mat=history_basis_mat,
                **common_assembler_kwargs,
            )
            current_bps_per_seed = np.array(
                [cv.cv_bits_per_spike for cv in current_cv_per_seed]
            )
        else:
            break

    final_cv = current_cv_per_seed[0] if current_cv_per_seed is not None else None
    final_bps = (
        float(np.mean(current_bps_per_seed))
        if n_seeds > 1
        else float(current_bps_per_seed[0])
    )

    # Hardcastle final-vs-null check: the assembled model must clear the
    # signed-rank test against the null. Forward steps each passed vs the
    # running model, but signed-rank is not transitive, so we record (and
    # warn on) the final-vs-null p-value without auto-clearing the selection.
    final_vs_null_pval = float("nan")
    if selection_rule == "signed_rank" and selected and final_cv is not None:
        final_vs_null_pval = _signed_rank_greater(
            final_cv.fold_bits_per_spike - null_cv.fold_bits_per_spike
        )
        if final_vs_null_pval >= config.selection_alpha:
            logger.warning(
                "signed_rank: final model %s does NOT clear signed-rank vs "
                "null (p=%.3g ≥ alpha=%.3g) despite per-step admission",
                selected, final_vs_null_pval, config.selection_alpha,
            )

    return SelectionResult(
        selected_vars=selected,
        history=history,
        null_cv_bps=float(np.mean(null_bps_per_seed))
        if n_seeds > 1
        else float(null_bps_per_seed[0]),
        final_cv_bps=final_bps,
        null_cv=null_cv,
        final_cv=final_cv,
        final_vs_null_pval=final_vs_null_pval,
        selection_rule=selection_rule,
    )


def _try_candidates(
    candidates: Sequence[str],
    already_selected: Sequence[str],
    B_speed: np.ndarray,
    B_tf: np.ndarray,
    B_onset: np.ndarray,
    sf_vals: np.ndarray,
    or_vals: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float,
    fold_ids_list: Sequence[np.ndarray],
    current_bps_per_seed: np.ndarray,
    *,
    phase: int,
    round_num: int,
    backend: str,
    config: GLMConfig,
    threshold_count: int,
    current_cv_list: Sequence[CVResult],
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
    history_basis_mat: np.ndarray | None = None,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
    B_accel: np.ndarray | None = None,
    B_sf: np.ndarray | None = None,
    B_or: np.ndarray | None = None,
    include_onset_kernel: bool = True,
) -> RoundResult:
    """Evaluate each candidate under the cv-fold partition(s) and admit by
    ``config.selection_rule``:

    - ``"delta_bps_threshold"`` (legacy): admit candidates whose Δ cv_bps
      clears ``delta_bps_threshold`` in ≥ ``threshold_count`` of N partitions;
      "best" = highest mean Δ. Reduces to single-seed when N == 1.
    - ``"signed_rank"`` (Hardcastle): a single partition; admit candidates
      whose per-fold paired Δ bits/spike (vs ``current_cv_list[0]``) is
      significant by a one-sided Wilcoxon signed-rank test at
      ``config.selection_alpha``; "best" = highest **median** per-fold Δ.
    """
    threshold = config.delta_bps_threshold
    selection_rule = getattr(config, "selection_rule", "delta_bps_threshold")
    n_seeds = len(fold_ids_list)
    current_fold_bps = (
        current_cv_list[0].fold_bits_per_spike
        if selection_rule == "signed_rank" and len(current_cv_list) > 0
        else None
    )

    tested: dict[str, float] = {}
    deltas_mean: dict[str, float] = {}
    delta_bps_per_seed: dict[str, list[float]] = {}
    admitted_count: dict[str, int] = {}
    pvals: dict[str, float] = {}
    rank_key: dict[str, float] = {}     # the value the "best" pick maximises

    for cand in candidates:
        test_vars = list(already_selected) + [cand]
        X_test, test_names = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals, test_vars,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
            B_history=B_history,
            B_me_face=B_me_face,
            B_accel=B_accel,
            B_sf=B_sf,
            B_or=B_or,
            include_onset_kernel=include_onset_kernel,
        )
        if X_test.shape[1] >= y.size:
            tested[cand] = -np.inf
            deltas_mean[cand] = -np.inf
            delta_bps_per_seed[cand] = [-np.inf] * n_seeds
            admitted_count[cand] = 0
            pvals[cand] = 1.0
            rank_key[cand] = -np.inf
            continue

        penalty = _penalty_for(test_names, config, history_basis_mat)
        cv_bps_per_seed: list[float] = []
        deltas_for_cand: list[float] = []
        cand_cv_seed0: CVResult | None = None
        for seed_idx, seed_folds in enumerate(fold_ids_list):
            cv = cross_validate_glm(
                X_test, y, offset, seed_folds,
                lambda_ridge=config.lambda_ridge, backend=backend,
                penalty_matrix=penalty,
            )
            if seed_idx == 0:
                cand_cv_seed0 = cv
            cv_bps_per_seed.append(cv.cv_bits_per_spike)
            deltas_for_cand.append(
                cv.cv_bits_per_spike - float(current_bps_per_seed[seed_idx])
            )

        tested[cand] = float(np.mean(cv_bps_per_seed))
        deltas_mean[cand] = float(np.mean(deltas_for_cand))
        delta_bps_per_seed[cand] = deltas_for_cand
        admitted_count[cand] = int(sum(d > threshold for d in deltas_for_cand))

        if selection_rule == "signed_rank" and current_fold_bps is not None:
            per_fold_delta = cand_cv_seed0.fold_bits_per_spike - current_fold_bps
            pvals[cand] = _signed_rank_greater(per_fold_delta)
            valid = per_fold_delta[np.isfinite(per_fold_delta)]
            rank_key[cand] = float(np.median(valid)) if valid.size else -np.inf
        else:
            pvals[cand] = float("nan")
            rank_key[cand] = deltas_mean[cand]

    # Admission depends on the rule. Among passing candidates, pick the one
    # that maximises rank_key (mean Δ for the legacy rule; median per-fold Δ
    # for signed_rank).
    if selection_rule == "signed_rank":
        alpha = config.selection_alpha
        passing = [c for c in candidates if pvals.get(c, 1.0) < alpha]
    else:
        passing = [
            c for c in candidates
            if admitted_count.get(c, 0) >= threshold_count
        ]
    if passing:
        best_candidate: str | None = max(passing, key=lambda c: rank_key[c])
        best_delta = deltas_mean[best_candidate]
        added = True
        cv_after = tested[best_candidate]
    else:
        best_candidate = None
        best_delta = -np.inf
        added = False
        cv_after = float(np.mean(current_bps_per_seed))

    return RoundResult(
        round=round_num,
        phase=phase,
        tested=tested,
        delta_bps=deltas_mean,
        best_candidate=best_candidate,
        best_delta_bps=best_delta,
        added=added,
        cv_bps_after=cv_after,
        delta_bps_per_seed=delta_bps_per_seed,
        admitted_count=admitted_count,
        n_seeds=n_seeds,
        pval=pvals,
    )


def _cv_for_selected_per_seed(
    selected: Sequence[str],
    B_speed: np.ndarray,
    B_tf: np.ndarray,
    B_onset: np.ndarray,
    sf_vals: np.ndarray,
    or_vals: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float,
    fold_ids_list: Sequence[np.ndarray],
    backend: str,
    *,
    config: GLMConfig,
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
    history_basis_mat: np.ndarray | None = None,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
    B_accel: np.ndarray | None = None,
    B_sf: np.ndarray | None = None,
    B_or: np.ndarray | None = None,
    include_onset_kernel: bool = True,
) -> list[CVResult]:
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, list(selected),
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        B_history=B_history,
        B_me_face=B_me_face,
        B_accel=B_accel,
        B_sf=B_sf,
        B_or=B_or,
        include_onset_kernel=include_onset_kernel,
    )
    penalty = _penalty_for(names, config, history_basis_mat)
    return [
        cross_validate_glm(
            X, y, offset, seed_folds,
            lambda_ridge=config.lambda_ridge, backend=backend,
            penalty_matrix=penalty,
        )
        for seed_folds in fold_ids_list
    ]
