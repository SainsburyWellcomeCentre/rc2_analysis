"""End-to-end GLM pipeline: probe .mat → forward selection → CSVs.

Mirrors the orchestration in scripts/glm_single_cluster_analysis.m
(Section 4 onwards). For each VISp cluster:

1. (Optional) Prefilter via stationary-vs-motion Wilcoxon test.
2. Time-bin the trials into a long-format table.
3. Build raised-cosine bases for speed and TF, plus a causal onset
   kernel basis. Reference-code SF and OR.
4. Run trial-level k-fold cross-validation with Hardcastle-style
   forward selection.
5. Emit three CSVs:
   - ``glm_model_comparison.csv``  one row per (cluster, glm_type)
   - ``glm_selection_history.csv`` one row per forward-selection round
   - ``glm_coefficients.csv``      one row per fitted coefficient

The pipeline only runs the ``time`` glm_type (the equivalent of
``glm_type='time'`` in MATLAB — uses time-binned spike counts with the
``log(bin_width)`` offset).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis
from rc2_glm.config import GLMConfig
from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix, assemble_design_matrix_selected
from rc2_glm.fitting import fit_poisson_glm
from rc2_glm.forward_selection import SelectionResult, forward_select
from rc2_glm.io import ProbeData, load_probe_data
from rc2_glm.prefilter import prefilter_probe
from rc2_glm.time_binning import bin_cluster


GLM_TYPE = "time"   # only one glm_type in the Python port


@dataclass
class ClusterFit:
    cluster_id: int
    selection: SelectionResult
    additive_cv_bps: float
    full_interaction_cv_bps: float
    n_bins: int
    n_spikes: int
    coefficients: pd.DataFrame   # columns: model, coefficient, estimate, se


@dataclass
class PipelineResult:
    probe_id: str
    prefilter: pd.DataFrame
    model_comparison: pd.DataFrame
    selection_history: pd.DataFrame
    coefficients: pd.DataFrame


def run_pipeline(
    mat_path: str | Path,
    config: GLMConfig | None = None,
    output_dir: str | Path | None = None,
    stimulus_lookup: StimulusLookup | None = None,
    backend: str = "irls",
    visp_only: bool = True,
) -> PipelineResult:
    """Load a probe .mat, fit GLMs cluster-by-cluster, write CSVs."""
    config = config or GLMConfig()

    probe = load_probe_data(
        mat_path,
        config=config,
        stimulus_lookup=stimulus_lookup,
        visp_only=visp_only,
    )

    if config.apply_prefilter:
        prefilter_df = prefilter_probe(probe, config=config)
        keep_ids = set(prefilter_df.loc[prefilter_df["should_run_glm"], "cluster_id"])
        clusters = [c for c in probe.clusters if c.cluster_id in keep_ids]
    else:
        prefilter_df = pd.DataFrame()
        clusters = probe.clusters

    comparison_rows: list[dict] = []
    history_rows: list[dict] = []
    coefficient_rows: list[dict] = []

    for cluster in clusters:
        df = bin_cluster(probe, cluster)
        if df.empty:
            continue
        fit = _fit_one_cluster(probe, df, cluster.cluster_id, config, backend)
        if fit is None:
            continue
        comparison_rows.append(_comparison_row(probe.probe_id, df, fit))
        history_rows.extend(_history_rows(probe.probe_id, fit))
        coefficient_rows.extend(_coefficient_rows(probe.probe_id, fit))

    comparison_df = pd.DataFrame(comparison_rows)
    history_df = pd.DataFrame(history_rows)
    coef_df = pd.DataFrame(coefficient_rows)

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        comparison_df.to_csv(out / "glm_model_comparison.csv", index=False)
        history_df.to_csv(out / "glm_selection_history.csv", index=False)
        coef_df.to_csv(out / "glm_coefficients.csv", index=False)
        if not prefilter_df.empty:
            prefilter_df.to_csv(out / "prefilter_decision_tree.csv", index=False)

    return PipelineResult(
        probe_id=probe.probe_id,
        prefilter=prefilter_df,
        model_comparison=comparison_df,
        selection_history=history_df,
        coefficients=coef_df,
    )


def _fit_one_cluster(
    probe: ProbeData,
    df: pd.DataFrame,
    cluster_id: int,
    config: GLMConfig,
    backend: str,
) -> ClusterFit | None:
    motion = df[df["condition"] != "stationary"].copy()
    if motion.empty or motion["spike_count"].sum() == 0:
        return None

    speed = motion["speed"].to_numpy(dtype=np.float64)
    tf = motion["tf"].to_numpy(dtype=np.float64)
    onset = motion["time_since_onset"].to_numpy(dtype=np.float64)
    sf_vals = motion["sf"].to_numpy(dtype=np.float64)
    or_vals = motion["orientation"].to_numpy(dtype=np.float64)
    y = motion["spike_count"].to_numpy(dtype=np.float64)
    trial_ids = motion["trial_id"].to_numpy(dtype=np.int64)

    B_speed = raised_cosine_basis(
        speed, config.n_speed_bases, *config.speed_range
    )
    B_tf = raised_cosine_basis(
        tf, config.n_tf_bases, *config.tf_range
    )
    B_onset = onset_kernel_basis(
        onset, config.n_onset_bases, config.onset_range[1]
    )

    offset = float(np.log(config.time_bin_width))
    fold_ids = make_trial_folds(trial_ids, config.n_folds, config.cv_seed)

    selection = forward_select(
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids,
        config=config, backend=backend,
    )

    additive_cv = _cv_for_label(
        "Additive", B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids, backend,
    )
    full_int_cv = _cv_for_label(
        "FullInteraction", B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids, backend,
    )

    coef_df = _coef_table_for_models(
        selection.selected_vars,
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, config, backend,
    )

    return ClusterFit(
        cluster_id=cluster_id,
        selection=selection,
        additive_cv_bps=additive_cv,
        full_interaction_cv_bps=full_int_cv,
        n_bins=int(motion.shape[0]),
        n_spikes=int(y.sum()),
        coefficients=coef_df,
    )


def _cv_for_label(
    label: str,
    B_speed: np.ndarray, B_tf: np.ndarray, B_onset: np.ndarray,
    sf_vals: np.ndarray, or_vals: np.ndarray,
    y: np.ndarray, offset: float,
    fold_ids: np.ndarray, backend: str,
) -> float:
    X, _ = assemble_design_matrix(
        B_speed, B_tf, B_onset, sf_vals, or_vals, label,
    )
    if X.shape[1] >= y.size:
        return float("nan")
    cv = cross_validate_glm(X, y, offset, fold_ids, lambda_ridge=0.0, backend=backend)
    return float(cv.cv_bits_per_spike)


def _coef_table_for_models(
    selected_vars: list[str],
    B_speed: np.ndarray, B_tf: np.ndarray, B_onset: np.ndarray,
    sf_vals: np.ndarray, or_vals: np.ndarray,
    y: np.ndarray, offset: float,
    config: GLMConfig, backend: str,
) -> pd.DataFrame:
    rows: list[dict] = []
    for label, vars_for_label in (
        ("Selected", selected_vars),
        ("Additive", ["Speed", "TF", "SF", "OR"]),
    ):
        X, names = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals, vars_for_label,
        )
        if X.shape[1] == 0 or X.shape[1] >= y.size:
            continue
        fit = fit_poisson_glm(
            X, y, offset, lambda_ridge=config.lambda_ridge, backend=backend,
        )
        for name, beta, se in zip(names, fit.beta, fit.se):
            rows.append({
                "model": label, "coefficient": name,
                "estimate": float(beta), "se": float(se),
            })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# CSV row builders
# --------------------------------------------------------------------------- #


def _comparison_row(probe_id: str, df: pd.DataFrame, fit: ClusterFit) -> dict:
    sel = fit.selection
    selected = sel.selected_vars
    selected_set = set(selected)
    delta = sel.final_cv_bps - sel.null_cv_bps
    has_int = any(v for v in selected if "_x_" in v)

    return {
        "probe_id": probe_id,
        "cluster_id": fit.cluster_id,
        "n_trials": int(df["trial_id"].nunique()),
        f"{GLM_TYPE}_n_bins": fit.n_bins,
        f"{GLM_TYPE}_n_spikes": fit.n_spikes,
        f"{GLM_TYPE}_selected_vars": "+".join(selected) if selected else "Null",
        f"{GLM_TYPE}_n_selected_vars": len(selected),
        f"{GLM_TYPE}_selection_rounds": len(sel.history),
        f"{GLM_TYPE}_Null_cv_bps": sel.null_cv_bps,
        f"{GLM_TYPE}_Selected_cv_bps": sel.final_cv_bps,
        f"{GLM_TYPE}_Additive_cv_bps": fit.additive_cv_bps,
        f"{GLM_TYPE}_FullInteraction_cv_bps": fit.full_interaction_cv_bps,
        f"{GLM_TYPE}_delta_selected_vs_null": delta,
        f"{GLM_TYPE}_delta_additive_vs_null": fit.additive_cv_bps - sel.null_cv_bps,
        f"{GLM_TYPE}_delta_interaction": fit.full_interaction_cv_bps - fit.additive_cv_bps,
        f"{GLM_TYPE}_delta_selected_vs_additive": sel.final_cv_bps - fit.additive_cv_bps,
        f"{GLM_TYPE}_is_speed_tuned": "Speed" in selected_set,
        f"{GLM_TYPE}_is_tf_tuned": "TF" in selected_set,
        f"{GLM_TYPE}_is_sf_tuned": "SF" in selected_set,
        f"{GLM_TYPE}_is_or_tuned": "OR" in selected_set,
        f"{GLM_TYPE}_has_interaction": has_int,
        f"{GLM_TYPE}_has_speed_x_tf": "Speed_x_TF" in selected_set,
        f"{GLM_TYPE}_has_speed_x_sf": "Speed_x_SF" in selected_set,
        f"{GLM_TYPE}_has_speed_x_or": "Speed_x_OR" in selected_set,
        f"{GLM_TYPE}_has_tf_x_sf": "TF_x_SF" in selected_set,
        f"{GLM_TYPE}_has_tf_x_or": "TF_x_OR" in selected_set,
        f"{GLM_TYPE}_has_sf_x_or": "SF_x_OR" in selected_set,
    }


def _history_rows(probe_id: str, fit: ClusterFit) -> list[dict]:
    return [
        {
            "probe_id": probe_id,
            "cluster_id": fit.cluster_id,
            "round": r.round,
            "phase": r.phase,
            "best_candidate": r.best_candidate or "",
            "delta_bps": r.best_delta_bps,
            "added": r.added,
            "cv_bps_after": r.cv_bps_after,
        }
        for r in fit.selection.history
    ]


def _coefficient_rows(probe_id: str, fit: ClusterFit) -> list[dict]:
    rows: list[dict] = []
    for _, row in fit.coefficients.iterrows():
        rows.append({
            "probe_id": probe_id,
            "cluster_id": fit.cluster_id,
            "glm_type": GLM_TYPE,
            "model": row["model"],
            "coefficient": row["coefficient"],
            "estimate": row["estimate"],
            "se": row["se"],
        })
    return rows
