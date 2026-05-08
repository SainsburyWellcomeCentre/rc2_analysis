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
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from rc2_glm.config import GLMConfig, INTERACTION_PARENTS
from rc2_glm.cross_validation import CVResult, cross_validate_glm
from rc2_glm.design_matrix import assemble_design_matrix_selected


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


@dataclass
class SelectionResult:
    selected_vars: list[str]
    history: list[RoundResult]
    null_cv_bps: float
    final_cv_bps: float
    null_cv: CVResult           # canonical (seed-0) null CV result
    final_cv: CVResult | None   # canonical (seed-0) final CV result


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

    common_assembler_kwargs = dict(
        B_history=B_history if include_history else None,
        B_me_face=B_me_face,
        include_onset_kernel=include_onset,
    )

    # Null model: intercept (+ onset kernel if include_onset_kernel)
    X_null, _ = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, [],
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        **common_assembler_kwargs,
    )
    null_cv_per_seed = [
        cross_validate_glm(
            X_null, y, offset, seed_folds,
            lambda_ridge=config.lambda_ridge, backend=backend,
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
            threshold_count=threshold_count,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
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
            threshold_count=threshold_count,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
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

    return SelectionResult(
        selected_vars=selected,
        history=history,
        null_cv_bps=float(np.mean(null_bps_per_seed))
        if n_seeds > 1
        else float(null_bps_per_seed[0]),
        final_cv_bps=final_bps,
        null_cv=null_cv,
        final_cv=final_cv,
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
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
    include_onset_kernel: bool = True,
) -> RoundResult:
    """Evaluate each candidate under N cv-fold partitions, admit
    candidates that pass threshold_count of N, pick the passing
    candidate with the highest mean Δ.

    Reduces to single-seed when ``len(fold_ids_list) == 1`` and
    ``threshold_count == 1``.
    """
    threshold = config.delta_bps_threshold
    n_seeds = len(fold_ids_list)

    tested: dict[str, float] = {}
    deltas_mean: dict[str, float] = {}
    delta_bps_per_seed: dict[str, list[float]] = {}
    admitted_count: dict[str, int] = {}

    for cand in candidates:
        test_vars = list(already_selected) + [cand]
        X_test, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals, test_vars,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
            B_history=B_history,
            B_me_face=B_me_face,
            include_onset_kernel=include_onset_kernel,
        )
        if X_test.shape[1] >= y.size:
            tested[cand] = -np.inf
            deltas_mean[cand] = -np.inf
            delta_bps_per_seed[cand] = [-np.inf] * n_seeds
            admitted_count[cand] = 0
            continue

        cv_bps_per_seed: list[float] = []
        deltas_for_cand: list[float] = []
        for seed_idx, seed_folds in enumerate(fold_ids_list):
            cv = cross_validate_glm(
                X_test, y, offset, seed_folds,
                lambda_ridge=config.lambda_ridge, backend=backend,
            )
            cv_bps_per_seed.append(cv.cv_bits_per_spike)
            deltas_for_cand.append(
                cv.cv_bits_per_spike - float(current_bps_per_seed[seed_idx])
            )

        tested[cand] = float(np.mean(cv_bps_per_seed))
        deltas_mean[cand] = float(np.mean(deltas_for_cand))
        delta_bps_per_seed[cand] = deltas_for_cand
        admitted_count[cand] = int(sum(d > threshold for d in deltas_for_cand))

    # Admit only candidates that clear threshold in ≥k of N seeds.
    # Among those, pick the one with the highest mean Δ.
    passing = [
        c for c in candidates
        if admitted_count.get(c, 0) >= threshold_count
    ]
    if passing:
        best_candidate: str | None = max(passing, key=lambda c: deltas_mean[c])
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
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
    include_onset_kernel: bool = True,
) -> list[CVResult]:
    X, _ = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, list(selected),
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        B_history=B_history,
        B_me_face=B_me_face,
        include_onset_kernel=include_onset_kernel,
    )
    return [
        cross_validate_glm(
            X, y, offset, seed_folds,
            lambda_ridge=config.lambda_ridge, backend=backend,
        )
        for seed_folds in fold_ids_list
    ]
