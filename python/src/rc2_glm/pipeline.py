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
6. (Optional) Write matplotlib figures under ``<out_dir>/figs/`` —
   basis functions, forward-selection summary, per-cluster selection
   waterfall, per-cluster coefficients, and the MATLAB Fig 5
   reconstructed tuning curves. PDF by default; PNG also supported.

The pipeline only runs the ``time`` glm_type (the equivalent of
``glm_type='time'`` in MATLAB — uses time-binned spike counts with the
``log(bin_width)`` offset).
"""

from __future__ import annotations

import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from rc2_formatted_data_reader import FormattedDataReader, StimulusLookup
from rc2_formatted_data_reader.reader import PROTOCOL_VELOCITY_CHANNEL
from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis
from rc2_glm.config import GLMConfig
from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix, assemble_design_matrix_selected
from rc2_glm.fitting import fit_poisson_glm
from rc2_glm.forward_selection import SelectionResult, forward_select
from rc2_glm.io import ProbeData, load_probe_data
from rc2_glm.prefilter import per_trial_firing_rate, prefilter_probe
from rc2_glm.time_binning import bin_cluster


GLM_TYPE = "time"   # only one glm_type in the Python port

logger = logging.getLogger("rc2_glm")


@dataclass
class ClusterFit:
    cluster_id: int
    selection: SelectionResult
    additive_cv_bps: float
    full_interaction_cv_bps: float
    n_bins: int
    n_spikes: int
    coefficients: pd.DataFrame   # columns: model, coefficient, estimate, se
    cluster_df: pd.DataFrame     # full cluster table (motion + stationary rows)
    model_betas: dict[str, np.ndarray]           # 'Null'|'Selected'|'Additive'|'FullInteraction' -> β
    model_col_names: dict[str, list[str]]        # same keys -> design column names
    model_predictions: dict[str, np.ndarray]     # same keys -> per-motion-row predicted FR (Hz)


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
    make_plots: bool = True,
    plot_format: str = "pdf",
    plot_clusters: int | None = None,
) -> PipelineResult:
    """Load a probe .mat, fit GLMs cluster-by-cluster, write CSVs + plots."""
    config = config or GLMConfig()

    out_path = Path(output_dir) if output_dir is not None else None
    log_handlers: list[logging.Handler] = []
    if out_path is not None:
        out_path.mkdir(parents=True, exist_ok=True)
        log_handlers = _setup_logging(out_path / "run.log")

    try:
        return _run_pipeline_inner(
            mat_path=mat_path,
            config=config,
            output_dir=out_path,
            stimulus_lookup=stimulus_lookup,
            backend=backend,
            visp_only=visp_only,
            make_plots=make_plots,
            plot_format=plot_format,
            plot_clusters=plot_clusters,
        )
    finally:
        for h in log_handlers:
            logger.removeHandler(h)
            h.close()


def _run_pipeline_inner(
    mat_path: str | Path,
    config: GLMConfig,
    output_dir: Path | None,
    stimulus_lookup: StimulusLookup | None,
    backend: str,
    visp_only: bool,
    make_plots: bool,
    plot_format: str,
    plot_clusters: int | None,
) -> PipelineResult:
    _banner("Loading data")
    logger.info("mat file: %s", mat_path)
    logger.info("backend: %s | visp_only: %s | prefilter: %s",
                backend, visp_only, config.apply_prefilter)

    probe = load_probe_data(
        mat_path,
        config=config,
        stimulus_lookup=stimulus_lookup,
        visp_only=visp_only,
    )
    logger.info("probe_id=%s | n_trials=%d | n_clusters=%d",
                probe.probe_id, len(probe.trials), len(probe.clusters))

    _log_trial_channel_summary(probe)
    if output_dir is not None:
        _emit_trial_channel_report(mat_path, probe, output_dir / "diagnostics")
        _emit_stationary_motion_fr(probe, output_dir / "diagnostics")

    if config.apply_prefilter:
        _banner("Prefilter")
        prefilter_df = prefilter_probe(probe, config=config)
        keep_ids = set(prefilter_df.loc[prefilter_df["should_run_glm"], "cluster_id"])
        _log_prefilter_summary(prefilter_df, keep_ids)
        clusters = [c for c in probe.clusters if c.cluster_id in keep_ids]
    else:
        prefilter_df = pd.DataFrame()
        clusters = probe.clusters
        logger.info("prefilter skipped (apply_prefilter=False)")

    _banner("Forward selection")
    logger.info("fitting %d clusters", len(clusters))
    comparison_rows: list[dict] = []
    history_rows: list[dict] = []
    coefficient_rows: list[dict] = []
    cluster_fits: list[tuple[ClusterFit, pd.DataFrame]] = []

    t0 = time.perf_counter()
    for idx, cluster in enumerate(clusters, start=1):
        df = bin_cluster(probe, cluster)
        if df.empty:
            logger.info("[%d/%d] cluster %d: no bins, skipped",
                        idx, len(clusters), cluster.cluster_id)
            continue
        if len(df) < 10:
            logger.info("[%d/%d] cluster %d: fewer than 10 bins, skipped",
                        idx, len(clusters), cluster.cluster_id)
            continue
        t_cluster = time.perf_counter()
        fit = _fit_one_cluster(probe, df, cluster.cluster_id, config, backend)
        if fit is None:
            logger.info("[%d/%d] cluster %d: no usable spikes, skipped",
                        idx, len(clusters), cluster.cluster_id)
            continue
        comparison_rows.append(_comparison_row(probe.probe_id, df, fit))
        history_rows.extend(_history_rows(probe.probe_id, fit))
        coefficient_rows.extend(_coefficient_rows(probe.probe_id, fit))
        cluster_fits.append((fit, df))
        logger.info(
            "[%d/%d] cluster %d: selected=%s | null %.3f → sel %.3f bps (%.1fs)",
            idx, len(clusters), cluster.cluster_id,
            "+".join(fit.selection.selected_vars) or "Null",
            fit.selection.null_cv_bps, fit.selection.final_cv_bps,
            time.perf_counter() - t_cluster,
        )
    logger.info("forward selection done in %.1fs", time.perf_counter() - t0)

    comparison_df = pd.DataFrame(comparison_rows)
    history_df = pd.DataFrame(history_rows)
    coef_df = pd.DataFrame(coefficient_rows)

    if output_dir is not None:
        _banner("CSV export")
        comparison_df.to_csv(output_dir / "glm_model_comparison.csv", index=False)
        logger.info("wrote glm_model_comparison.csv (%d rows)", len(comparison_df))
        history_df.to_csv(output_dir / "glm_selection_history.csv", index=False)
        logger.info("wrote glm_selection_history.csv (%d rows)", len(history_df))
        coef_df.to_csv(output_dir / "glm_coefficients.csv", index=False)
        logger.info("wrote glm_coefficients.csv (%d rows)", len(coef_df))
        if cluster_fits:
            _emit_tuning_curves(
                probe.probe_id, cluster_fits, config,
                output_dir / "diagnostics",
            )
        if not prefilter_df.empty:
            prefilter_df.to_csv(output_dir / "prefilter_decision_tree.csv", index=False)
            logger.info("wrote prefilter_decision_tree.csv (%d rows)", len(prefilter_df))

        if make_plots and cluster_fits:
            _banner("Plots")
            _write_all_plots(
                probe.probe_id, comparison_df, history_df, coef_df,
                cluster_fits, config, output_dir / "figs",
                plot_format=plot_format, plot_clusters=plot_clusters,
                backend=backend,
            )

    _banner("Done")
    logger.info("probe=%s | retained clusters=%d", probe.probe_id, len(comparison_rows))

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
    motion_mask = df["condition"] != "stationary"
    motion = df.loc[motion_mask].copy()
    if motion.empty or df["spike_count"].sum() == 0:
        return None

    speed = df["speed"].to_numpy(dtype=np.float64)
    tf = df["tf"].to_numpy(dtype=np.float64)
    onset = df["time_since_onset"].to_numpy(dtype=np.float64)
    sf_vals = df["sf"].to_numpy(dtype=np.float64)
    or_vals = df["orientation"].to_numpy(dtype=np.float64)
    y = df["spike_count"].to_numpy(dtype=np.float64)
    trial_ids = df["trial_id"].to_numpy(dtype=np.int64)
    condition_labels = df["condition"].to_numpy(dtype=object)

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
    fold_ids = make_trial_folds(
        trial_ids,
        config.n_folds,
        config.cv_seed,
        condition_labels_per_bin=condition_labels,
    )

    sf_valid = sf_vals[(sf_vals != 0.0) & ~np.isnan(sf_vals)]
    sf_ref_levels = np.sort(np.unique(sf_valid)).tolist() if sf_valid.size > 0 else []

    or_valid = or_vals[(or_vals != 0.0) & ~np.isnan(or_vals)]
    or_ref_levels = np.sort(np.unique(or_valid)).tolist() if or_valid.size > 0 else []

    selection = forward_select(
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids,
        config=config, backend=backend,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )

    additive_cv = _cv_for_label(
        "Additive", B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids, backend, config,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )
    full_int_cv = _cv_for_label(
        "FullInteraction", B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids, backend, config,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )

    coef_df, betas, col_names_by_model, preds = _fit_plot_models(
        selection.selected_vars,
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, config, backend,
        motion_row_mask=motion_mask.to_numpy(dtype=bool),
        cluster_id=cluster_id,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )

    return ClusterFit(
        cluster_id=cluster_id,
        selection=selection,
        additive_cv_bps=additive_cv,
        full_interaction_cv_bps=full_int_cv,
        n_bins=int(df.shape[0]),
        n_spikes=int(y.sum()),
        coefficients=coef_df,
        cluster_df=df,
        model_betas=betas,
        model_col_names=col_names_by_model,
        model_predictions=preds,
    )


def _cv_for_label(
    label: str,
    B_speed: np.ndarray, B_tf: np.ndarray, B_onset: np.ndarray,
    sf_vals: np.ndarray, or_vals: np.ndarray,
    y: np.ndarray, offset: float,
    fold_ids: np.ndarray, backend: str, config: GLMConfig,
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
) -> float:
    X, _ = assemble_design_matrix(
        B_speed, B_tf, B_onset, sf_vals, or_vals, label,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )
    if X.shape[1] >= y.size:
        return float("nan")
    lambda_ridge = config.full_interaction_lambda if label == "FullInteraction" else 0.0
    cv = cross_validate_glm(X, y, offset, fold_ids, lambda_ridge=lambda_ridge, backend=backend)
    return float(cv.cv_bits_per_spike)


_ADDITIVE_VARS = ["Speed", "TF", "SF", "OR"]
_FULL_INTERACTION_VARS = [
    "Speed", "TF", "SF", "OR",
    "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
    "TF_x_SF", "TF_x_OR", "SF_x_OR",
]


def _fit_plot_models(
    selected_vars: list[str],
    B_speed: np.ndarray, B_tf: np.ndarray, B_onset: np.ndarray,
    sf_vals: np.ndarray, or_vals: np.ndarray,
    y: np.ndarray, offset: float,
    config: GLMConfig, backend: str,
    motion_row_mask: np.ndarray,
    cluster_id: int,
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
) -> tuple[
    pd.DataFrame,
    dict[str, np.ndarray],
    dict[str, list[str]],
    dict[str, np.ndarray],
]:
    """Fit Null / Selected / Additive / FullInteraction in-sample.

    Returns the CSV coefficient table (Selected + Additive only, preserving the
    existing export surface) plus the per-model β vectors, design column
    names, and per-motion-row predicted firing rates (Hz) that the
    per-cluster plots need. Each fit decision (fitted, skipped because
    ``X.shape[1] >= y.size``, or empty design) is logged so a run-log
    reader can explain blank rows in ``cluster_<id>_overview.pdf``.
    """
    csv_rows: list[dict] = []
    betas: dict[str, np.ndarray] = {}
    col_names_by_model: dict[str, list[str]] = {}
    preds: dict[str, np.ndarray] = {}

    model_defs: list[tuple[str, list[str]]] = [
        ("Null", []),
        ("Selected", list(selected_vars)),
        ("Additive", _ADDITIVE_VARS),
        ("FullInteraction", _FULL_INTERACTION_VARS),
    ]

    for label, vars_for_label in model_defs:
        X, names = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals, vars_for_label,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        )
        if X.shape[1] == 0:
            logger.warning(
                "cluster %d: model %s has 0 design columns, skipping (vars=%s)",
                cluster_id, label, vars_for_label,
            )
            continue
        if X.shape[1] >= y.size:
            logger.warning(
                "cluster %d: model %s skipped — design has %d cols but only "
                "%d bins (X.shape[1] >= y.size). Row will render as N/A.",
                cluster_id, label, X.shape[1], y.size,
            )
            continue
        lambda_ridge = (
            config.full_interaction_lambda
            if label == "FullInteraction"
            else config.lambda_ridge
        )
        fit = fit_poisson_glm(
            X, y, offset, lambda_ridge=lambda_ridge, backend=backend,
        )
        betas[label] = np.asarray(fit.beta, dtype=np.float64).copy()
        col_names_by_model[label] = list(names)
        preds_full = np.exp(np.clip(X @ fit.beta, -20.0, 20.0))
        preds[label] = preds_full[motion_row_mask]
        logger.info(
            "cluster %d: fit %s | %d coefs | pred FR range [%.2f, %.2f] Hz",
            cluster_id, label, X.shape[1],
            float(preds[label].min()), float(preds[label].max()),
        )
        for name, beta_v, se_v in zip(names, fit.beta, fit.se):
            csv_rows.append({
                "model": label, "coefficient": name,
                "estimate": float(beta_v), "se": float(se_v),
            })
    return pd.DataFrame(csv_rows), betas, col_names_by_model, preds


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


# --------------------------------------------------------------------------- #
# Tuning-curve export (prediction-space parity)
# --------------------------------------------------------------------------- #


# Number of grid points used for Speed / TF tuning curves. 100 is plenty
# for Pearson r to stabilise but small enough to keep the CSV modest even
# across all probes / clusters / models.
_TUNING_CURVE_GRID_N: int = 100

# Variables we export tuning curves for. Only the two that live on a
# continuous raised-cosine basis (Speed, TF) — SF / OR are categorical
# dummies, so their "tuning curve" is just the per-level β and is already
# captured by the coefficient-sign check.
_TUNING_CURVE_VARIABLES: tuple[tuple[str, str, str], ...] = (
    # (variable_label, coef_prefix, config_range_attr)
    ("Speed", "Speed_", "speed_range"),
    ("TF", "TF_", "tf_range"),
)


def _tuning_curve_grid(
    variable: str, config: GLMConfig
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(x_grid, B_grid)`` used to evaluate a tuning curve.

    The grid covers the configured ``speed_range`` / ``tf_range`` with
    ``_TUNING_CURVE_GRID_N`` linearly-spaced points (the raised-cosine
    basis already encodes log-spacing internally via ``log(x + 0.5)``,
    so a linear grid in x is the natural axis for plotting / parity).
    """
    if variable == "Speed":
        lo, hi = config.speed_range
        n_b = config.n_speed_bases
    elif variable == "TF":
        lo, hi = config.tf_range
        n_b = config.n_tf_bases
    else:
        raise ValueError(f"unsupported tuning-curve variable: {variable!r}")

    x_grid = np.linspace(float(lo), float(hi), _TUNING_CURVE_GRID_N)
    B_grid = raised_cosine_basis(x_grid, n_b, float(lo), float(hi))
    return x_grid, B_grid


def _extract_variable_betas(
    betas: np.ndarray, col_names: list[str], coef_prefix: str,
) -> np.ndarray | None:
    """Return the subvector of ``betas`` whose column names start with ``coef_prefix``.

    Order is taken from ``col_names`` (Speed_1, Speed_2, …). Returns None
    if the variable isn't present in this model's design — the caller
    should skip the row.
    """
    idx = [
        i for i, name in enumerate(col_names)
        if name.startswith(coef_prefix) and "_x_" not in name
    ]
    if not idx:
        return None
    return np.asarray(betas, dtype=np.float64)[idx]


def _tuning_curves_for_fit(
    probe_id: str, fit: ClusterFit, config: GLMConfig,
) -> list[dict]:
    """Emit tuning-curve rows for one cluster, one row per (model, variable, grid point).

    Only the main effect of each variable is exported — interactions are
    left out so the curve is the marginal ``exp(B_grid @ β_var)`` gain
    factor. This is a prediction-space quantity that is invariant to
    basis-rotation of β, which is what we want for parity.
    """
    rows: list[dict] = []
    for variable, coef_prefix, _range_attr in _TUNING_CURVE_VARIABLES:
        x_grid, B_grid = _tuning_curve_grid(variable, config)
        for model_label, beta_vec in fit.model_betas.items():
            col_names = fit.model_col_names.get(model_label, [])
            beta_var = _extract_variable_betas(beta_vec, col_names, coef_prefix)
            if beta_var is None:
                continue
            if beta_var.size != B_grid.shape[1]:
                # Defensive: a downstream change that renames or drops a
                # basis column would land here. Skip rather than crash.
                logger.warning(
                    "cluster %d / %s / %s: have %d β for %d bases — skipping tuning row",
                    fit.cluster_id, model_label, variable,
                    int(beta_var.size), int(B_grid.shape[1]),
                )
                continue
            log_gain = B_grid @ beta_var
            gain = np.exp(np.clip(log_gain, -20.0, 20.0))
            for x, g, lg in zip(x_grid, gain, log_gain):
                rows.append({
                    "probe_id": probe_id,
                    "cluster_id": fit.cluster_id,
                    "glm_type": GLM_TYPE,
                    "model": model_label,
                    "variable": variable,
                    "grid_x": float(x),
                    "log_gain": float(lg),
                    "gain": float(g),
                })
    return rows


def _emit_tuning_curves(
    probe_id: str,
    cluster_fits: list[tuple[ClusterFit, pd.DataFrame]],
    config: GLMConfig,
    out_dir: Path,
) -> None:
    """Write ``<out_dir>/tuning_curves.csv`` with marginal Speed / TF curves."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for fit, _df in cluster_fits:
        rows.extend(_tuning_curves_for_fit(probe_id, fit, config))
    if not rows:
        logger.info("no tuning-curve rows to write (no Speed/TF βs?)")
        return
    df = pd.DataFrame(rows)
    path = out_dir / "tuning_curves.csv"
    df.to_csv(path, index=False)
    logger.info("wrote %s (%d rows)", path.name, len(df))


# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #


def _setup_logging(log_path: Path) -> list[logging.Handler]:
    """Attach a StreamHandler(stdout) + FileHandler(log_path) to the logger.

    Returns the handlers so callers can remove them when the run ends.
    """
    logger.setLevel(logging.INFO)
    logger.propagate = False

    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(fmt)
    stream.setLevel(logging.INFO)
    logger.addHandler(stream)

    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(fmt)
    file_handler.setLevel(logging.INFO)
    logger.addHandler(file_handler)

    return [stream, file_handler]


def _banner(title: str) -> None:
    line = "=" * 72
    logger.info(line)
    logger.info(title)
    logger.info(line)


def _log_prefilter_summary(prefilter_df: pd.DataFrame, keep_ids: set[int]) -> None:
    counts = prefilter_df["category"].value_counts()
    for cat, n in counts.items():
        logger.info("  %s: %d", cat, int(n))
    logger.info("total clusters: %d | keep for GLM: %d",
                int(len(prefilter_df)), len(keep_ids))


# --------------------------------------------------------------------------- #
# Trial-channel diagnostics (protocol velocity selection)
# --------------------------------------------------------------------------- #


# All channels we audit per trial. The "chosen" channel is derived from
# protocol via PROTOCOL_VELOCITY_CHANNEL; the rest are audited so that a
# wrong choice or a dead channel is visible in the report.
_AUDIT_CHANNELS: tuple[str, ...] = (
    "filtered_teensy",
    "filtered_teensy_2",
    "stage",
    "multiplexer_output",
)

# A motion fraction below this on the chosen channel is flagged as a
# likely-bad trial (matches the prompt's 1% heuristic).
_LOW_MOTION_WARN_FRACTION: float = 0.01


def _log_trial_channel_summary(probe: ProbeData) -> None:
    """Log protocol counts and which velocity channel each protocol resolved to."""
    proto_counts: dict[str, int] = {}
    channel_counts: dict[str, int] = {}
    for t in probe.trials:
        proto_counts[t.protocol] = proto_counts.get(t.protocol, 0) + 1
        channel_counts[t.velocity_channel] = channel_counts.get(t.velocity_channel, 0) + 1
    logger.info("protocol distribution: %s",
                ", ".join(f"{k}: {v}" for k, v in sorted(proto_counts.items())))
    logger.info("velocity channel distribution: %s",
                ", ".join(f"{k}: {v}" for k, v in sorted(channel_counts.items())))


def _emit_trial_channel_report(
    mat_path: str | Path,
    probe: ProbeData,
    out_dir: Path,
) -> None:
    """Write per-trial motion fractions across all candidate channels.

    One row per trial, with motion fraction (|v| >= 1 cm/s) computed on each
    of ``filtered_teensy``, ``filtered_teensy_2``, ``stage``, and
    ``multiplexer_output``, plus the chosen channel's motion fraction. Trials
    whose chosen-channel motion fraction is < 1% get a WARNING line — on
    passive paradigms this is the canary for "channel selection is still
    wrong somewhere."
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    vel_thresh = probe.config.velocity_threshold
    rows: list[dict] = []
    with FormattedDataReader(mat_path) as reader:
        for t in probe.trials:
            s, e = reader.trial_bounds(t.trial_idx)
            duration_s = (e - s) / reader.fs if reader.fs else float("nan")
            row: dict = {
                "trial_idx": t.trial_idx,
                "trial_id": t.trial_id,
                "protocol": t.protocol,
                "condition": t.condition,
                "duration_s": duration_s,
                "chosen_channel": t.velocity_channel,
            }
            for ch in _AUDIT_CHANNELS:
                trace = reader.session_channel(ch, s, e)
                if trace is None:
                    row[f"motion_frac_{ch}"] = float("nan")
                else:
                    row[f"motion_frac_{ch}"] = float(
                        (np.abs(trace) >= vel_thresh).mean()
                    )
            row["motion_frac_chosen"] = row.get(f"motion_frac_{t.velocity_channel}", float("nan"))
            rows.append(row)

    df = pd.DataFrame(rows)
    report_path = out_dir / "trial_channel_report.csv"
    df.to_csv(report_path, index=False)
    logger.info("wrote %s (%d rows)",
                report_path.relative_to(out_dir.parent), len(df))

    low = df[df["motion_frac_chosen"] < _LOW_MOTION_WARN_FRACTION]
    for _, row in low.iterrows():
        logger.warning(
            "trial_idx=%d trial_id=%d %s: chosen channel %s motion fraction "
            "%.3f%% < %.1f%% threshold (likely bad trial)",
            int(row["trial_idx"]), int(row["trial_id"]), row["protocol"],
            row["chosen_channel"], 100.0 * row["motion_frac_chosen"],
            100.0 * _LOW_MOTION_WARN_FRACTION,
        )


def _emit_stationary_motion_fr(probe: ProbeData, out_dir: Path) -> None:
    """Write Python-native per-(cluster,trial) stationary/motion FR.

    Parallel to the MATLAB-generated ``csvs/stationary_vs_motion_fr/<probe>.csv``.
    Required by ``rc2-glm-compare`` check 7 (stationary-vs-motion FR parity)
    — the velocity fix must bring the Python and MATLAB numbers within 1e-3
    median relative error.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for cluster in probe.clusters:
        for trial in probe.trials:
            stat_fr, mot_fr = per_trial_firing_rate(cluster, trial, probe.fs)
            stat_time = float(trial.stationary_mask.sum()) / probe.fs
            mot_time = float(trial.motion_mask.sum()) / probe.fs
            rows.append({
                "probe_id": probe.probe_id,
                "trial_id": trial.trial_id,
                "trial_group_label": trial.condition,
                "cluster_id": cluster.cluster_id,
                "stationary_fr": stat_fr,
                "stationary_time": stat_time,
                "motion_fr": mot_fr,
                "motion_time": mot_time,
            })
    df = pd.DataFrame(rows)
    path = out_dir / "stationary_vs_motion_fr_python.csv"
    df.to_csv(path, index=False)
    logger.info("wrote %s (%d rows)",
                path.relative_to(out_dir.parent), len(df))


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #


def _write_all_plots(
    probe_id: str,
    comparison_df: pd.DataFrame,
    history_df: pd.DataFrame,
    coef_df: pd.DataFrame,
    cluster_fits: list[tuple["ClusterFit", pd.DataFrame]],
    config: GLMConfig,
    figs_dir: Path,
    *,
    plot_format: str,
    plot_clusters: int | None,
    backend: str,
) -> None:
    from rc2_glm import plots

    figs_dir.mkdir(parents=True, exist_ok=True)

    for path in plots.save_figure(
        plots.plot_basis_functions(config),
        figs_dir / "basis_functions",
        fmt=plot_format,
    ):
        logger.info("wrote %s", path.relative_to(figs_dir.parent))

    for path in plots.save_figure(
        plots.plot_forward_selection_summary(comparison_df, history_df),
        figs_dir / "forward_selection_summary",
        fmt=plot_format,
    ):
        logger.info("wrote %s", path.relative_to(figs_dir.parent))

    fits_to_plot = (
        cluster_fits if plot_clusters is None else cluster_fits[:plot_clusters]
    )
    for fit, _ in fits_to_plot:
        cid = fit.cluster_id
        for fn, name in (
            (plots.plot_cluster_model_overview(
                probe_id=probe_id,
                cluster_id=cid,
                cluster_df=fit.cluster_df,
                selected_vars=fit.selection.selected_vars,
                model_betas=fit.model_betas,
                model_col_names=fit.model_col_names,
                model_predictions=fit.model_predictions,
                config=config,
             ), f"cluster_{cid}_overview"),
            (plots.plot_tuning_curves(
                probe_id=probe_id,
                cluster_id=cid,
                cluster_df=fit.cluster_df,
                model_betas=fit.model_betas,
                model_col_names=fit.model_col_names,
                config=config,
             ), f"cluster_{cid}_tuning"),
        ):
            for path in plots.save_figure(fn, figs_dir / name, fmt=plot_format):
                logger.info("wrote %s", path.relative_to(figs_dir.parent))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _resolve_mat_path(arg_mat: str | None, env_path: str | None, env_dir: str | None) -> Path:
    """Resolve the formatted .mat file from CLI arg + env vars.

    Precedence:
      1. CLI arg that contains a path separator or exists on disk.
      2. CLI arg as a bare filename inside ``RC2_FORMATTED_DATA_DIR``.
      3. ``RC2_FORMATTED_DATA_PATH`` env var.
      4. First ``*.mat`` inside ``RC2_FORMATTED_DATA_DIR``.
    """
    if arg_mat:
        p = Path(arg_mat).expanduser()
        if p.exists():
            return p
        if env_dir and "/" not in arg_mat and "\\" not in arg_mat:
            candidate = Path(env_dir).expanduser() / arg_mat
            if candidate.exists():
                return candidate
        raise SystemExit(f"rc2-glm: no such file: {arg_mat}")
    if env_path:
        p = Path(env_path).expanduser()
        if p.exists():
            return p
    if env_dir:
        d = Path(env_dir).expanduser()
        if d.is_dir():
            mats = sorted(d.glob("*.mat"))
            if mats:
                return mats[0]
    raise SystemExit(
        "rc2-glm: no .mat specified. Pass one as a positional arg or set "
        "RC2_FORMATTED_DATA_PATH / RC2_FORMATTED_DATA_DIR in python/.env."
    )


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: ``rc2-glm [mat] [out_dir] [--backend irls|nemos]``.

    Reads defaults from ``python/.env`` via ``python-dotenv``:

    - ``RC2_FORMATTED_DATA_PATH`` — single .mat file
    - ``RC2_FORMATTED_DATA_DIR``  — directory (first *.mat chosen)
    - ``RC2_MC_SEQUENCE_PATH``    — motion-cloud sequence .mat
    - ``RC2_MC_FOLDERS_PATH``     — image_folders.mat
    - ``RC2_GLM_OUTPUT_DIR``      — default output directory
    """
    import argparse
    import os
    from dataclasses import replace

    from dotenv import load_dotenv

    load_dotenv(_find_env_file())

    env_path = os.environ.get("RC2_FORMATTED_DATA_PATH")
    env_dir = os.environ.get("RC2_FORMATTED_DATA_DIR")
    env_out = os.environ.get("RC2_GLM_OUTPUT_DIR")
    env_mc_seq = os.environ.get("RC2_MC_SEQUENCE_PATH")
    env_mc_fol = os.environ.get("RC2_MC_FOLDERS_PATH")

    parser = argparse.ArgumentParser(
        prog="rc2-glm",
        description=(
            "Run the rc2 Poisson GLM pipeline. Paths default to values in "
            "python/.env (RC2_FORMATTED_DATA_DIR, RC2_MC_*_PATH, "
            "RC2_GLM_OUTPUT_DIR); CLI flags override."
        ),
    )
    parser.add_argument(
        "mat", nargs="?", default=None,
        help="Formatted v7.3 .mat file (path, or bare filename inside "
             "RC2_FORMATTED_DATA_DIR). Optional if .env provides one.",
    )
    parser.add_argument(
        "out_dir", nargs="?", default=env_out,
        help="Directory for CSV outputs (default: $RC2_GLM_OUTPUT_DIR)",
    )
    parser.add_argument(
        "--backend", choices=("irls", "nemos"), default="irls",
        help="Fitting backend (default: irls)",
    )
    parser.add_argument(
        "--all-clusters", action="store_true",
        help="Include non-VISp clusters (default: VISp only)",
    )
    parser.add_argument(
        "--no-prefilter", dest="prefilter", action="store_false",
        help="Skip the stationary/motion Wilcoxon prefilter before GLM",
    )
    parser.add_argument(
        "--mc-sequence", default=env_mc_seq,
        help="Motion-cloud sequence .mat (default: $RC2_MC_SEQUENCE_PATH)",
    )
    parser.add_argument(
        "--mc-folders", default=env_mc_fol,
        help="image_folders.mat (default: $RC2_MC_FOLDERS_PATH)",
    )
    parser.add_argument(
        "--no-plots", dest="make_plots", action="store_false",
        help="Skip all matplotlib figures (CSVs + log still written)",
    )
    parser.add_argument(
        "--plot-format", choices=("pdf", "png", "both"), default="pdf",
        help="Format for plot output (default: pdf).",
    )
    parser.add_argument(
        "--plot-clusters", type=int, default=None,
        help="Cap per-cluster plots to the first N retained clusters (debug).",
    )
    parser.set_defaults(make_plots=True, prefilter=True)
    args = parser.parse_args(argv)

    mat_path = _resolve_mat_path(args.mat, env_path, env_dir)
    if not args.out_dir:
        raise SystemExit(
            "rc2-glm: no output directory. Pass as second positional arg "
            "or set RC2_GLM_OUTPUT_DIR in python/.env."
        )

    config = replace(GLMConfig(), apply_prefilter=args.prefilter)

    lookup = None
    if args.mc_sequence and args.mc_folders:
        from rc2_formatted_data_reader import StimulusLookup
        lookup = StimulusLookup(args.mc_sequence, args.mc_folders)

    result = run_pipeline(
        mat_path=mat_path,
        config=config,
        output_dir=args.out_dir,
        stimulus_lookup=lookup,
        backend=args.backend,
        visp_only=not args.all_clusters,
        make_plots=args.make_plots,
        plot_format=args.plot_format,
        plot_clusters=args.plot_clusters,
    )
    # run_pipeline has already torn down its handlers; a bare print is fine.
    print(f"rc2-glm: wrote {len(result.model_comparison)} cluster rows to {args.out_dir}")
    return 0


def _find_env_file() -> str | None:
    """Locate python/.env by walking up from the CWD, then from this file."""
    from dotenv import find_dotenv

    found = find_dotenv(usecwd=True)
    if found:
        return found
    here = Path(__file__).resolve().parent
    for parent in (here, *here.parents):
        candidate = parent / ".env"
        if candidate.is_file():
            return str(candidate)
    return None


if __name__ == "__main__":
    raise SystemExit(main())
