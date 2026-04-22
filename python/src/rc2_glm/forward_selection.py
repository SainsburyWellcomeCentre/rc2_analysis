"""Hardcastle-style hierarchical forward selection.

Two phases (mirrors MATLAB ``forward_select_model`` at line 5903):

- **Phase 1 — Main effects.** Test {Speed, TF, SF, OR} one at a time.
  Add the best if Δ CV bps > threshold. Repeat until no main effect
  passes.
- **Phase 2 — Interactions.** Only test interactions whose BOTH
  parent main effects were selected in Phase 1. Same Δ-bps rule.

Returns the selected variable list, a per-round history, and the
final / null CV bits-per-spike.
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
    tested: dict[str, float]              # candidate → cv_bps
    delta_bps: dict[str, float]
    best_candidate: str | None
    best_delta_bps: float
    added: bool
    cv_bps_after: float


@dataclass
class SelectionResult:
    selected_vars: list[str]
    history: list[RoundResult]
    null_cv_bps: float
    final_cv_bps: float
    null_cv: CVResult
    final_cv: CVResult | None


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
) -> SelectionResult:
    config = config or GLMConfig()

    # Null model: intercept + onset kernel
    X_null, _ = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, [],
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )
    null_cv = cross_validate_glm(
        X_null, y, offset, fold_ids,
        lambda_ridge=config.lambda_ridge, backend=backend,
    )

    selected: list[str] = []
    current_cv = null_cv
    history: list[RoundResult] = []
    round_num = 0

    # ----- Phase 1: main effects -----
    remaining = list(config.main_effects)
    while remaining:
        round_num += 1
        round_result = _try_candidates(
            remaining, selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
            y, offset, fold_ids, current_cv.cv_bits_per_spike,
            phase=1, round_num=round_num, backend=backend, config=config,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        )
        history.append(round_result)
        if round_result.added and round_result.best_candidate is not None:
            selected.append(round_result.best_candidate)
            remaining.remove(round_result.best_candidate)
            # Refit & store CV for the new model so subsequent rounds keep going
            current_cv = _cv_for_selected(
                selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
                y, offset, fold_ids, backend, config=config,
                sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
            )
        else:
            break

    # ----- Phase 2: interactions (only those with both parents selected) -----
    eligible: list[str] = [
        name for name in config.interactions
        if all(parent in selected for parent in INTERACTION_PARENTS[name])
    ]
    while eligible:
        round_num += 1
        round_result = _try_candidates(
            eligible, selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
            y, offset, fold_ids, current_cv.cv_bits_per_spike,
            phase=2, round_num=round_num, backend=backend, config=config,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        )
        history.append(round_result)
        if round_result.added and round_result.best_candidate is not None:
            selected.append(round_result.best_candidate)
            eligible.remove(round_result.best_candidate)
            current_cv = _cv_for_selected(
                selected, B_speed, B_tf, B_onset, sf_vals, or_vals,
                y, offset, fold_ids, backend, config=config,
                sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
            )
        else:
            break

    return SelectionResult(
        selected_vars=selected,
        history=history,
        null_cv_bps=null_cv.cv_bits_per_spike,
        final_cv_bps=current_cv.cv_bits_per_spike,
        null_cv=null_cv,
        final_cv=current_cv,
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
    fold_ids: np.ndarray,
    current_bps: float,
    *,
    phase: int,
    round_num: int,
    backend: str,
    config: GLMConfig,
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
) -> RoundResult:
    threshold = config.delta_bps_threshold

    tested: dict[str, float] = {}
    deltas: dict[str, float] = {}
    best_candidate: str | None = None
    best_delta = -np.inf

    for cand in candidates:
        test_vars = list(already_selected) + [cand]
        X_test, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals, test_vars,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        )
        if X_test.shape[1] >= y.size:
            tested[cand] = -np.inf
            deltas[cand] = -np.inf
            continue
        cv = cross_validate_glm(
            X_test, y, offset, fold_ids,
            lambda_ridge=config.lambda_ridge, backend=backend,
        )
        tested[cand] = cv.cv_bits_per_spike
        deltas[cand] = cv.cv_bits_per_spike - current_bps
        if deltas[cand] > best_delta:
            best_delta = deltas[cand]
            best_candidate = cand

    added = best_candidate is not None and best_delta > threshold
    cv_after = (
        tested[best_candidate] if (added and best_candidate is not None) else current_bps
    )
    return RoundResult(
        round=round_num,
        phase=phase,
        tested=tested,
        delta_bps=deltas,
        best_candidate=best_candidate,
        best_delta_bps=best_delta,
        added=added,
        cv_bps_after=cv_after,
    )


def _cv_for_selected(
    selected: Sequence[str],
    B_speed: np.ndarray,
    B_tf: np.ndarray,
    B_onset: np.ndarray,
    sf_vals: np.ndarray,
    or_vals: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float,
    fold_ids: np.ndarray,
    backend: str,
    *,
    config: GLMConfig,
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
) -> CVResult:
    X, _ = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, list(selected),
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )
    return cross_validate_glm(
        X, y, offset, fold_ids,
        lambda_ridge=config.lambda_ridge, backend=backend,
    )
