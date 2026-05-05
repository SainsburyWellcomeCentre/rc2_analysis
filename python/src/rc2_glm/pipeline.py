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


_FIT_STATUS_FITTED = "fitted"
_FIT_STATUS_P_GE_N = "skipped_p_ge_n"
_FIT_STATUS_EMPTY = "skipped_empty_design"


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
    # Status trail for each model: one of the _FIT_STATUS_* constants above.
    # ``cv_status`` tracks the cross-validated Additive / FullInteraction
    # columns in ``glm_model_comparison.csv`` (NaN values become
    # explicable); ``refit_status`` tracks the full-data refits behind
    # ``model_betas`` / ``model_predictions`` (blank rows in per-cluster
    # plots become explicable). The two can disagree on the same cluster
    # because the CV folds have ~80% of the bins that the refit has, so
    # ``p>=n`` can fire for the CV but not the refit.
    cv_status: dict[str, str]        # keys: "Additive", "FullInteraction"
    refit_status: dict[str, str]     # keys: "Null", "Selected", "Additive", "FullInteraction"
    # Post-hoc speed-profile CV (MATLAB glm_single_cluster_analysis.m:2246-2367).
    # Populated only when config.profile_cv_diagnostic is True; otherwise
    # every entry is NaN. Keys: "Null", "Selected", "NoSpeed". "NoSpeed"
    # is the Selected model re-fit with Speed dropped, for the speed-unique
    # contribution diagnostic — NaN if Speed wasn't in the Selected set.
    profile_cv_bps: dict[str, float]


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
    n_jobs: int = -1,
    cluster_filter: set[int] | None = None,
) -> PipelineResult:
    """Load a probe .mat, fit GLMs cluster-by-cluster, write CSVs + plots.

    ``n_jobs`` controls the per-cluster fit parallelism (joblib loky /
    process backend). -1 uses all cores; 1 forces serial execution (also
    useful for debugging — lets per-cluster worker log lines reach
    run.log). Output is deterministic up to BLAS-level FP rounding
    (~1e-11 absolute in CV-bps); cluster order and variable selection
    are preserved exactly.
    """
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
            cluster_filter=cluster_filter,
            n_jobs=n_jobs,
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
    n_jobs: int,
    cluster_filter: set[int] | None = None,
) -> PipelineResult:
    _banner("Loading data")
    logger.info("mat file: %s", mat_path)
    logger.info("backend: %s | visp_only: %s | prefilter: %s",
                backend, visp_only, config.apply_prefilter)
    if backend == "nemos":
        from rc2_glm.fitting import configure_jax_device
        resolved = configure_jax_device(config.device)
        logger.info("nemos → JAX backend=%s (requested device=%s)",
                    resolved, config.device)

    probe = load_probe_data(
        mat_path,
        config=config,
        stimulus_lookup=stimulus_lookup,
        visp_only=visp_only,
    )
    if cluster_filter is not None:
        # Upstream cluster allow-list. Restrict the probe's cluster set
        # to ids that appear in ``cluster_filter`` (typically loaded from
        # MATLAB's prefilter_decision_tree.csv so Python fits the same
        # ~100 clusters MATLAB does instead of all VISp). Clusters absent
        # from the filter are dropped here — the prefilter + GLM loop
        # never see them.
        before = len(probe.clusters)
        kept = [c for c in probe.clusters if c.cluster_id in cluster_filter]
        probe.clusters = kept
        logger.info(
            "cluster_filter applied: %d → %d clusters "
            "(dropped %d absent from the filter list)",
            before, len(kept), before - len(kept),
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
    logger.info("fitting %d clusters (n_jobs=%d)", len(clusters), n_jobs)
    comparison_rows: list[dict] = []
    history_rows: list[dict] = []
    coefficient_rows: list[dict] = []
    cluster_fits: list[tuple[ClusterFit, pd.DataFrame]] = []

    t0 = time.perf_counter()
    # Phase 1: bin every cluster on the main thread (cheap, touches `probe`).
    # Drop empties / too-short runs up-front so joblib dispatches real work only.
    fit_tasks: list[tuple[int, int, pd.DataFrame]] = []
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
        fit_tasks.append((idx, cluster.cluster_id, df))

    # Phase 2: fit (serial or process-parallel via joblib loky).
    # Processes avoid BLAS thread oversubscription on Accelerate / OpenBLAS
    # (threading backend was break-even here because each IRLS matmul
    # already saturates all cores). Trade-off: per-cluster info lines
    # emitted inside workers stay in the worker's stderr and do NOT
    # land in run.log; summary lines below are the user-visible trace.
    if n_jobs == 1 or len(fit_tasks) <= 1:
        fit_results = [
            _fit_one_cluster(df, cid, config, backend)
            for (_, cid, df) in fit_tasks
        ]
    else:
        from joblib import Parallel, delayed
        fit_results = Parallel(n_jobs=n_jobs)(
            delayed(_fit_one_cluster)(df, cid, config, backend)
            for (_, cid, df) in fit_tasks
        )

    # Phase 3: collect in deterministic cluster order.
    for (idx, cid, df), fit in zip(fit_tasks, fit_results):
        if fit is None:
            logger.info("[%d/%d] cluster %d: no usable spikes, skipped",
                        idx, len(clusters), cid)
            continue
        comparison_rows.append(_comparison_row(probe.probe_id, df, fit))
        history_rows.extend(_history_rows(probe.probe_id, fit))
        coefficient_rows.extend(_coefficient_rows(probe.probe_id, fit))
        cluster_fits.append((fit, df))
        logger.info(
            "[%d/%d] cluster %d: selected=%s | null %.3f → sel %.3f bps",
            idx, len(clusters), cid,
            "+".join(fit.selection.selected_vars) or "Null",
            fit.selection.null_cv_bps, fit.selection.final_cv_bps,
        )
    logger.info("forward selection done in %.1fs", time.perf_counter() - t0)
    _log_fit_status_summary(cluster_fits)

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
            _emit_trial_level_metrics(
                probe.probe_id, cluster_fits, config,
                output_dir / "diagnostics",
            )
            from rc2_glm.precomputed_bins import load_precomputed_bin_edges
            try:
                pre_bins = load_precomputed_bin_edges(probe.mat_path)
            except Exception as e:  # noqa: BLE001 — best-effort cache, fallback OK
                logger.info(
                    "precomputed bin edges unavailable (%s); per-trial "
                    "uncertainty CSV will use fallback grids", e,
                )
                pre_bins = None
            _emit_tuning_curves_uncertainty(
                probe.probe_id, cluster_fits, config,
                output_dir / "diagnostics",
                precomputed_bins=pre_bins,
            )
        if not prefilter_df.empty:
            prefilter_df.to_csv(output_dir / "prefilter_decision_tree.csv", index=False)
            logger.info("wrote prefilter_decision_tree.csv (%d rows)", len(prefilter_df))

        if make_plots and cluster_fits:
            _banner("Plots")
            _write_all_plots(
                probe.probe_id, comparison_df, history_df, coef_df,
                cluster_fits, config, output_dir / "figs",
                formatted_mat_path=probe.mat_path,
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


def _log_fit_status_summary(
    cluster_fits: list[tuple[ClusterFit, pd.DataFrame]],
) -> None:
    """Log per-model skip counts so the run log explains NaN CV-bps rows."""
    if not cluster_fits:
        return
    cv_counts: dict[str, dict[str, int]] = {}
    refit_counts: dict[str, dict[str, int]] = {}
    for fit, _ in cluster_fits:
        for label, status in fit.cv_status.items():
            cv_counts.setdefault(label, {}).setdefault(status, 0)
            cv_counts[label][status] += 1
        for label, status in fit.refit_status.items():
            refit_counts.setdefault(label, {}).setdefault(status, 0)
            refit_counts[label][status] += 1
    for label in sorted(cv_counts):
        parts = ", ".join(
            f"{st}={cv_counts[label][st]}" for st in sorted(cv_counts[label])
        )
        logger.info("cv-status %s: %s", label, parts)
    for label in sorted(refit_counts):
        parts = ", ".join(
            f"{st}={refit_counts[label][st]}" for st in sorted(refit_counts[label])
        )
        logger.info("refit-status %s: %s", label, parts)


def _fit_one_cluster(
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
    profile_ids = (
        df["profile_id"].to_numpy(dtype=np.int64)
        if "profile_id" in df.columns else None
    )

    B_speed = raised_cosine_basis(
        speed, config.n_speed_bases, *config.speed_range
    )
    B_tf = raised_cosine_basis(
        tf, config.n_tf_bases, *config.tf_range
    )
    B_onset = onset_kernel_basis(
        onset, config.n_onset_bases, config.onset_range[1]
    )
    # Per-cluster spike-history features (causal, trial-aware). Built
    # only when config.include_history is True (prompt 03 opt-in);
    # otherwise None — every downstream caller treats None as "no history".
    if getattr(config, "include_history", False):
        from rc2_glm.basis import history_basis, convolve_history
        h_basis = history_basis(
            n_bases=config.n_history_bases,
            t_max_s=config.history_window_s,
            bin_width_s=config.time_bin_width,
        )
        B_history = convolve_history(y, trial_ids, h_basis)
    else:
        B_history = None

    # Face motion energy basis (Speed-style raised cosines over the
    # z-scored ME range, prompt 06, 2026-04-30). Z-scoring stats are
    # computed on motion-condition rows of THIS cluster's df — same as
    # any other probe-cluster of this session, since me_face_raw is a
    # session-level signal aligned to the same per-bin grid.
    # B_me_face is None when the camera trace is absent (all-NaN
    # me_face_raw column from time_binning) or when fewer than 10 finite
    # motion-row samples exist for z-score statistics — downstream
    # forward_select / assemble_design_matrix branches treat None as
    # "no ME term, drop ME_face from Phase-1 candidates".
    if "me_face_raw" in df.columns and getattr(config, "include_me_face", True):
        me_raw = df["me_face_raw"].to_numpy(dtype=np.float64)
        motion_rows = (df["condition"] != "stationary").to_numpy()
        me_motion_finite = np.isfinite(me_raw) & motion_rows
        if int(me_motion_finite.sum()) >= 10:
            me_mean = float(me_raw[me_motion_finite].mean())
            me_std = float(me_raw[me_motion_finite].std(ddof=0))
            if me_std <= 0.0:
                me_std = 1.0  # degenerate session; basis = constant column
            me_z = np.where(np.isfinite(me_raw), (me_raw - me_mean) / me_std, 0.0)
            from rc2_glm.basis import raised_cosine_basis_linear
            B_me_face = raised_cosine_basis_linear(
                me_z,
                config.n_me_face_bases,
                config.me_face_range[0],
                config.me_face_range[1],
            )
        else:
            B_me_face = None
    else:
        B_me_face = None

    offset = float(np.log(config.time_bin_width))
    fold_ids = make_trial_folds(
        trial_ids,
        config.n_folds,
        config.cv_seed,
        condition_labels_per_bin=condition_labels,
        strategy=config.cv_strategy,
        profile_ids_per_bin=profile_ids,
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
        B_history=B_history,
        B_me_face=B_me_face,
    )

    additive_cv, additive_status = _cv_for_label(
        "Additive", B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids, backend, config,
        cluster_id=cluster_id,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        B_history=B_history,
        B_me_face=B_me_face,
    )
    full_int_cv, full_int_status = _cv_for_label(
        "FullInteraction", B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, fold_ids, backend, config,
        cluster_id=cluster_id,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        B_history=B_history,
        B_me_face=B_me_face,
    )

    coef_df, betas, col_names_by_model, preds, refit_status = _fit_plot_models(
        selection.selected_vars,
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        y, offset, config, backend,
        motion_row_mask=motion_mask.to_numpy(dtype=bool),
        cluster_id=cluster_id,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        B_history=B_history,
        B_me_face=B_me_face,
    )

    profile_cv_bps = _profile_cv_diagnostic(
        config=config,
        selection=selection,
        trial_ids=trial_ids,
        profile_ids=profile_ids,
        standard_fold_ids=fold_ids,
        B_speed=B_speed, B_tf=B_tf, B_onset=B_onset,
        sf_vals=sf_vals, or_vals=or_vals,
        y=y, offset=offset, backend=backend,
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
        cv_status={"Additive": additive_status, "FullInteraction": full_int_status},
        refit_status=refit_status,
        profile_cv_bps=profile_cv_bps,
    )


def _profile_cv_diagnostic(
    *,
    config: GLMConfig,
    selection: SelectionResult,
    trial_ids: np.ndarray,
    profile_ids: np.ndarray | None,
    standard_fold_ids: np.ndarray,
    B_speed: np.ndarray, B_tf: np.ndarray, B_onset: np.ndarray,
    sf_vals: np.ndarray, or_vals: np.ndarray,
    y: np.ndarray, offset: float, backend: str,
    cluster_id: int,
    sf_ref_levels: list[float] | None,
    or_ref_levels: list[float] | None,
) -> dict[str, float]:
    """Post-hoc speed-profile CV on Null / Selected / Selected-without-Speed.

    When ``config.profile_cv_diagnostic`` is False (or no profile_ids
    are available), returns all-NaN so downstream serialisation is
    uniform. Otherwise mirrors MATLAB glm_single_cluster_analysis.m:2289-2362:
    - Null = no variables (intercept + onset kernel only).
    - Selected = the vars the condition-stratified forward selection picked.
    - NoSpeed = Selected minus Speed (if Speed was selected).

    MATLAB runs each model under BOTH fold schemes so the paired Δ-bps
    plot in panel 2 (speed unique contribution under standard vs profile
    CV) has paired numbers. We mirror that: keys with suffix ``_standard``
    are under the condition-stratified fold (the primary fit's fold_ids),
    keys without suffix are under the profile-fold scheme. CV is skipped
    (NaN) when either profile group has < 5 bins or the design is
    rank-deficient — matches MATLAB's viability gate.
    """
    nan_keys = (
        "Null", "Selected", "NoSpeed",
        "Null_standard", "Selected_standard", "NoSpeed_standard",
    )
    nan_result = {k: float("nan") for k in nan_keys}
    if not config.profile_cv_diagnostic:
        return nan_result
    if profile_ids is None:
        logger.info(
            "cluster %d: profile-CV diagnostic skipped — no profile_ids "
            "(run rc2-glm with --mc-sequence + --mc-folders so "
            "StimulusLookup populates TrialData.profile_id)",
            cluster_id,
        )
        return nan_result

    # Profile folds (deterministic). Fails loud if inconsistent.
    try:
        profile_fold_ids = make_trial_folds(
            trial_ids,
            strategy="speed-profile",
            profile_ids_per_bin=profile_ids,
        )
    except ValueError as exc:
        logger.warning("cluster %d: profile-CV diagnostic skipped — %s",
                       cluster_id, exc)
        return nan_result

    # Viability: at least 5 bins in each profile.
    for pid in (1, 2):
        if int((profile_ids == pid).sum()) < 5:
            return nan_result

    out = dict(nan_result)

    def _cv_with_vars(vars_: list[str], fold_ids: np.ndarray) -> float:
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals, vars_,
            sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        )
        if X.shape[1] == 0 or X.shape[1] >= y.size:
            return float("nan")
        cv = cross_validate_glm(
            X, y, offset, fold_ids,
            lambda_ridge=config.lambda_ridge, backend=backend,
        )
        return float(cv.cv_bits_per_spike)

    selected_list = list(selection.selected_vars)
    no_speed = [v for v in selected_list if "Speed" not in v]
    has_speed = "Speed" in selected_list

    out["Null"] = _cv_with_vars([], profile_fold_ids)
    out["Selected"] = _cv_with_vars(selected_list, profile_fold_ids)
    out["Null_standard"] = _cv_with_vars([], standard_fold_ids)
    out["Selected_standard"] = _cv_with_vars(selected_list, standard_fold_ids)
    if has_speed:
        out["NoSpeed"] = _cv_with_vars(no_speed, profile_fold_ids)
        out["NoSpeed_standard"] = _cv_with_vars(no_speed, standard_fold_ids)
    return out


def _cv_for_label(
    label: str,
    B_speed: np.ndarray, B_tf: np.ndarray, B_onset: np.ndarray,
    sf_vals: np.ndarray, or_vals: np.ndarray,
    y: np.ndarray, offset: float,
    fold_ids: np.ndarray, backend: str, config: GLMConfig,
    cluster_id: int,
    sf_ref_levels: list[float] | None = None,
    or_ref_levels: list[float] | None = None,
    *,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
) -> tuple[float, str]:
    include_onset = getattr(config, "include_onset_kernel", True)
    X, _ = assemble_design_matrix(
        B_speed, B_tf, B_onset, sf_vals, or_vals, label,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        B_history=B_history,
        B_me_face=B_me_face,
        include_onset_kernel=include_onset,
    )
    if X.shape[1] == 0:
        logger.warning(
            "cluster %d: CV(%s) skipped — 0 design columns",
            cluster_id, label,
        )
        return float("nan"), _FIT_STATUS_EMPTY
    if X.shape[1] >= y.size:
        logger.warning(
            "cluster %d: CV(%s) skipped — design has %d cols vs %d bins "
            "(X.shape[1] >= y.size); comparison CSV row will be NaN.",
            cluster_id, label, X.shape[1], y.size,
        )
        return float("nan"), _FIT_STATUS_P_GE_N
    lambda_ridge = (
        config.full_interaction_lambda
        if label == "FullInteraction"
        else config.lambda_ridge
    )
    cv = cross_validate_glm(X, y, offset, fold_ids, lambda_ridge=lambda_ridge, backend=backend)
    return float(cv.cv_bits_per_spike), _FIT_STATUS_FITTED


_ADDITIVE_VARS = ["Speed", "TF", "SF", "OR"]
_FULL_INTERACTION_VARS = [
    "Speed", "TF", "SF", "OR", "ME_face",
    "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
    "TF_x_SF", "TF_x_OR", "SF_x_OR", "ME_face_x_Speed",
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
    *,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
) -> tuple[
    pd.DataFrame,
    dict[str, np.ndarray],
    dict[str, list[str]],
    dict[str, np.ndarray],
    dict[str, str],
]:
    """Fit Null / Selected / Additive / FullInteraction in-sample.

    Returns the CSV coefficient table (Selected + Additive only, preserving the
    existing export surface) plus the per-model β vectors, design column
    names, per-motion-row predicted firing rates (Hz) for the cluster
    plots, and a ``refit_status`` dict keyed by model label with one of
    ``_FIT_STATUS_FITTED`` / ``_FIT_STATUS_P_GE_N`` /
    ``_FIT_STATUS_EMPTY`` so blank rows in ``cluster_<id>_overview.pdf``
    can be explained from the run log or the comparison CSV without
    re-reading the warnings stream.
    """
    csv_rows: list[dict] = []
    betas: dict[str, np.ndarray] = {}
    col_names_by_model: dict[str, list[str]] = {}
    preds: dict[str, np.ndarray] = {}
    refit_status: dict[str, str] = {}

    include_onset = getattr(config, "include_onset_kernel", True)
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
            B_history=B_history,
            B_me_face=B_me_face,
            include_onset_kernel=include_onset,
        )
        if X.shape[1] == 0:
            logger.warning(
                "cluster %d: refit %s has 0 design columns, skipping (vars=%s)",
                cluster_id, label, vars_for_label,
            )
            refit_status[label] = _FIT_STATUS_EMPTY
            continue
        if X.shape[1] >= y.size:
            logger.warning(
                "cluster %d: refit %s skipped — design has %d cols but only "
                "%d bins (X.shape[1] >= y.size). Row will render as N/A.",
                cluster_id, label, X.shape[1], y.size,
            )
            refit_status[label] = _FIT_STATUS_P_GE_N
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
        # For the Selected model, also compute the "no-history" variant —
        # same fit, but with History_* coefficients zeroed before
        # multiplication. This lets the post-hoc α-sweep
        # (scripts/run_history_alpha_sweep.py) interpolate
        # log-predictions between α=0 (no history) and α=1 (full history),
        # quantifying the trial-level vs tuning-curve trade-off the
        # spike-history filter induces.
        if label == "Selected":
            history_mask = np.array(
                [n.startswith("History_") for n in names], dtype=bool,
            )
            if history_mask.any():
                beta_no_h = fit.beta.copy()
                beta_no_h[history_mask] = 0.0
                preds_no_h_full = np.exp(np.clip(X @ beta_no_h, -20.0, 20.0))
                preds["Selected_no_history"] = preds_no_h_full[motion_row_mask]
                col_names_by_model["Selected_no_history"] = list(names)
                betas["Selected_no_history"] = beta_no_h
            # Selected_no_me_face: zero out ME_face main-effect and
            # ME_face × Speed interaction columns at prediction time.
            # Mirrors the Selected_no_history pattern; downstream
            # variance partitioning (unique_ME_face_bps) uses this.
            # ME_face main columns: "ME_face_1..5"; interaction:
            # "MEf{i}_x_Spd{j}".
            me_mask = np.array(
                [n.startswith("ME_face_") or n.startswith("MEf") for n in names],
                dtype=bool,
            )
            if me_mask.any():
                beta_no_me = fit.beta.copy()
                beta_no_me[me_mask] = 0.0
                preds_no_me_full = np.exp(np.clip(X @ beta_no_me, -20.0, 20.0))
                preds["Selected_no_me_face"] = preds_no_me_full[motion_row_mask]
                col_names_by_model["Selected_no_me_face"] = list(names)
                betas["Selected_no_me_face"] = beta_no_me
        refit_status[label] = _FIT_STATUS_FITTED
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
    return pd.DataFrame(csv_rows), betas, col_names_by_model, preds, refit_status


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
        f"{GLM_TYPE}_Additive_cv_status": fit.cv_status.get("Additive", _FIT_STATUS_FITTED),
        f"{GLM_TYPE}_FullInteraction_cv_bps": fit.full_interaction_cv_bps,
        f"{GLM_TYPE}_FullInteraction_cv_status": fit.cv_status.get(
            "FullInteraction", _FIT_STATUS_FITTED,
        ),
        f"{GLM_TYPE}_FullInteraction_refit_status": fit.refit_status.get(
            "FullInteraction", _FIT_STATUS_FITTED,
        ),
        f"{GLM_TYPE}_delta_selected_vs_null": delta,
        f"{GLM_TYPE}_delta_additive_vs_null": fit.additive_cv_bps - sel.null_cv_bps,
        f"{GLM_TYPE}_delta_interaction": fit.full_interaction_cv_bps - fit.additive_cv_bps,
        f"{GLM_TYPE}_delta_selected_vs_additive": sel.final_cv_bps - fit.additive_cv_bps,
        f"{GLM_TYPE}_is_speed_tuned": "Speed" in selected_set,
        f"{GLM_TYPE}_is_tf_tuned": "TF" in selected_set,
        f"{GLM_TYPE}_is_sf_tuned": "SF" in selected_set,
        f"{GLM_TYPE}_is_or_tuned": "OR" in selected_set,
        f"{GLM_TYPE}_is_history_tuned": "History" in selected_set,
        f"{GLM_TYPE}_is_me_face_tuned": "ME_face" in selected_set,
        f"{GLM_TYPE}_has_interaction": has_int,
        f"{GLM_TYPE}_has_speed_x_tf": "Speed_x_TF" in selected_set,
        f"{GLM_TYPE}_has_speed_x_sf": "Speed_x_SF" in selected_set,
        f"{GLM_TYPE}_has_speed_x_or": "Speed_x_OR" in selected_set,
        f"{GLM_TYPE}_has_tf_x_sf": "TF_x_SF" in selected_set,
        f"{GLM_TYPE}_has_tf_x_or": "TF_x_OR" in selected_set,
        f"{GLM_TYPE}_has_sf_x_or": "SF_x_OR" in selected_set,
        f"{GLM_TYPE}_has_me_face_x_speed": "ME_face_x_Speed" in selected_set,
        # Post-hoc speed-profile CV (MATLAB glm_single_cluster_analysis.m:2246-2367).
        # Both fold schemes so panel 2 of the comparison plot can compute
        # Selected − NoSpeed under each. NaN when the diagnostic is
        # disabled; NoSpeed* additionally NaN when Speed wasn't selected.
        f"{GLM_TYPE}_profile_cv_bps_Null": fit.profile_cv_bps.get("Null", float("nan")),
        f"{GLM_TYPE}_profile_cv_bps_Selected": fit.profile_cv_bps.get("Selected", float("nan")),
        f"{GLM_TYPE}_profile_cv_bps_NoSpeed": fit.profile_cv_bps.get("NoSpeed", float("nan")),
        f"{GLM_TYPE}_profile_cv_bps_Null_standard": fit.profile_cv_bps.get("Null_standard", float("nan")),
        f"{GLM_TYPE}_profile_cv_bps_Selected_standard": fit.profile_cv_bps.get("Selected_standard", float("nan")),
        f"{GLM_TYPE}_profile_cv_bps_NoSpeed_standard": fit.profile_cv_bps.get("NoSpeed_standard", float("nan")),
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


def _emit_trial_level_metrics(
    probe_id: str,
    cluster_fits: list[tuple[ClusterFit, pd.DataFrame]],
    config: GLMConfig,
    out_dir: Path,
) -> None:
    """Per-(cluster, model) trial-level + per-trial Pearson r.

    Trial-level Pearson is correlation between the cluster's full
    sequence of motion-bin observed counts and the model's predicted
    rate (× bin width). Per-trial Pearson is the same correlation
    computed *within each trial* and then averaged (mean / median).
    Together these say "does the model track within-trial firing-rate
    fluctuations" and "across-trial-and-time, does the model's
    prediction line up with what we observed?".
    """
    from scipy.stats import pearsonr
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for fit, df in cluster_fits:
        motion_mask = (df["condition"] != "stationary").to_numpy()
        if motion_mask.sum() < 5:
            continue
        df_motion = df.loc[motion_mask].reset_index(drop=True)
        observed = df_motion["spike_count"].to_numpy(dtype=np.float64)
        trial_ids = df_motion["trial_id"].to_numpy()
        for model, pred_fr in fit.model_predictions.items():
            if pred_fr is None:
                continue
            # pred_fr is rates (Hz) per motion row; convert to expected
            # counts via Δt.
            expected = np.asarray(pred_fr, dtype=np.float64) * config.time_bin_width
            if expected.shape != observed.shape:
                continue
            try:
                r_total, _ = pearsonr(observed, expected)
            except (ValueError, RuntimeWarning):
                r_total = float("nan")
            # Per-trial Pearson — drop trials with < 3 bins or zero
            # variance in either side.
            per_trial_r: list[float] = []
            for tid in np.unique(trial_ids):
                m = trial_ids == tid
                if m.sum() < 3:
                    continue
                o = observed[m]; e = expected[m]
                if np.std(o) == 0 or np.std(e) == 0:
                    continue
                try:
                    r, _ = pearsonr(o, e)
                    if np.isfinite(r):
                        per_trial_r.append(float(r))
                except (ValueError, RuntimeWarning):
                    continue
            rows.append({
                "probe_id": probe_id,
                "cluster_id": fit.cluster_id,
                "model": model,
                "n_motion_bins": int(observed.size),
                "n_trials": int(np.unique(trial_ids).size),
                "pearson_r_overall": float(r_total) if np.isfinite(r_total) else float("nan"),
                "pearson_r_per_trial_mean": float(np.mean(per_trial_r)) if per_trial_r else float("nan"),
                "pearson_r_per_trial_median": float(np.median(per_trial_r)) if per_trial_r else float("nan"),
                "n_per_trial_used": len(per_trial_r),
            })
    if not rows:
        logger.info("no trial-level metric rows to write")
        return
    out = pd.DataFrame(rows)
    path = out_dir / "trial_level_metrics.csv"
    out.to_csv(path, index=False)
    logger.info("wrote %s (%d rows)", path.name, len(out))

    # History α-sweep: for clusters with both "Selected" and
    # "Selected_no_history" predictions, interpolate in log-prediction
    # space (η_α = η_no_history + α·(η_full − η_no_history); λ_α = exp η_α)
    # and compute trial-level Pearson r at each α. Quantifies the
    # potentiate/depotentiate question: at what α does the
    # autoregressive gain crossover the stimulus-tuning loss?
    alpha_grid = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
    alpha_rows: list[dict] = []
    for fit, df in cluster_fits:
        if "Selected" not in fit.model_predictions:
            continue
        if "Selected_no_history" not in fit.model_predictions:
            continue
        motion_mask = (df["condition"] != "stationary").to_numpy()
        if motion_mask.sum() < 5:
            continue
        df_motion = df.loc[motion_mask].reset_index(drop=True)
        observed = df_motion["spike_count"].to_numpy(dtype=np.float64)
        pred_full = np.asarray(fit.model_predictions["Selected"], dtype=np.float64)
        pred_no_h = np.asarray(fit.model_predictions["Selected_no_history"], dtype=np.float64)
        if pred_full.shape != observed.shape or pred_no_h.shape != observed.shape:
            continue
        # log-rate space
        log_full = np.log(np.maximum(pred_full, 1e-12))
        log_no_h = np.log(np.maximum(pred_no_h, 1e-12))
        for alpha in alpha_grid:
            log_alpha = log_no_h + alpha * (log_full - log_no_h)
            pred_alpha = np.exp(np.clip(log_alpha, -20.0, 20.0))
            expected_alpha = pred_alpha * config.time_bin_width
            try:
                r, _ = pearsonr(observed, expected_alpha)
            except (ValueError, RuntimeWarning):
                r = float("nan")
            alpha_rows.append({
                "probe_id": probe_id,
                "cluster_id": fit.cluster_id,
                "alpha": float(alpha),
                "pearson_r_overall": float(r) if np.isfinite(r) else float("nan"),
            })
    if alpha_rows:
        alpha_path = out_dir / "history_alpha_sweep.csv"
        pd.DataFrame(alpha_rows).to_csv(alpha_path, index=False)
        logger.info("wrote %s (%d rows)", alpha_path.name, len(alpha_rows))


# Variable → list of (condition, sweep_kwargs) pairs. Mirrors the
# layout in ``plots._predict_model_row_trial_averaged`` so the CSV is a
# faithful textual mirror of the plotted curves.
_PER_TRIAL_PANELS: dict[str, list[tuple[str, str, dict]]] = {
    "Speed": [
        ("T_Vstatic", "Speed", {"fix_tf": 0.0}),
        ("VT", "Speed", {}),
    ],
    "TF": [
        ("V", "TF", {"fix_speed": 0.0}),
        ("VT", "TF", {}),
    ],
    "SF": [
        ("V", "SF", {"fix_speed": 0.0, "fix_tf": None}),
        ("VT", "SF", {"fix_speed": None, "fix_tf": None}),
    ],
    "OR": [
        ("V", "OR", {"fix_speed": 0.0, "fix_tf": None}),
        ("VT", "OR", {"fix_speed": None, "fix_tf": None}),
    ],
}


def _emit_tuning_curves_uncertainty(
    probe_id: str,
    cluster_fits: list[tuple[ClusterFit, pd.DataFrame]],
    config: GLMConfig,
    out_dir: Path,
    *,
    precomputed_bins: "PrecomputedBinEdges | None" = None,
) -> None:
    """Write ``<out_dir>/tuning_curves_uncertainty.csv``.

    One row per (probe, cluster, model, variable, condition, grid_x).
    Columns: gain_mean (== plot line value), gain_q05/q25/q75/q95,
    gain_std, n_trials_in_bin. Sibling to ``tuning_curves.csv`` (which
    carries the theoretical, condition-agnostic marginal curve used for
    parity); this file is the per-condition trial-averaged curve with
    spread that the plot bands are read off.

    Only emitted when ``config.tuning_curve_mode == "trial-averaged"``;
    steady-state mode evaluates one prediction per condition with no
    per-trial axis to summarise.
    """
    if config.tuning_curve_mode != "trial-averaged":
        return
    # Lazy import to avoid a circular at module load.
    from rc2_glm.plots import (
        _marginal_curve_with_spread,
        _marginal_levels_with_spread,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    sf_ref = np.asarray(config.sf_levels, dtype=np.float64)
    or_ref = np.asarray(config.or_levels, dtype=np.float64)
    sf_plot = sf_ref.copy()
    or_plot = or_ref.copy()

    fallback_spd = np.linspace(*config.speed_range, 11)
    fallback_tf = np.linspace(*config.tf_range, 11)

    def _spd_centres(cond: str) -> np.ndarray:
        edges = (
            precomputed_bins.speed_edges(cond)
            if precomputed_bins is not None else None
        )
        if edges is None:
            edges = fallback_spd
        return 0.5 * (edges[:-1] + edges[1:])

    def _tf_centres(cond: str) -> np.ndarray:
        edges = (
            precomputed_bins.tf_edges(cond)
            if precomputed_bins is not None else None
        )
        if edges is None:
            edges = fallback_tf
        return 0.5 * (edges[:-1] + edges[1:])

    rows: list[dict] = []
    for fit, df in cluster_fits:
        motion = df[df["condition"] != "stationary"]
        if motion.empty:
            continue
        cond_bins: dict[str, pd.DataFrame] = {}
        for cond in ("T_Vstatic", "V", "VT"):
            sub = motion[motion["condition"] == cond]
            cond_bins[cond] = sub if not sub.empty else None

        for model_label, beta in fit.model_betas.items():
            train_names = fit.model_col_names.get(model_label)
            if beta is None or train_names is None:
                continue
            sel_vars = _vars_from_names_for_model(model_label, train_names)
            for variable, panels in _PER_TRIAL_PANELS.items():
                for cond, _kind, kwargs in panels:
                    cd = cond_bins.get(cond)
                    if cd is None:
                        continue
                    if variable in ("Speed", "TF"):
                        grid = (
                            _spd_centres(cond) if variable == "Speed"
                            else _tf_centres(cond)
                        )
                        summary = _marginal_curve_with_spread(
                            beta, train_names, sel_vars, config, sf_ref, or_ref,
                            cd, sweep_var=variable, grid=grid, **kwargs,
                        )
                    else:
                        levels = sf_plot if variable == "SF" else or_plot
                        # SF/OR panels need fix_or / fix_sf for the held-
                        # constant categorical. Use the same convention as
                        # _predict_model_row_trial_averaged.
                        held = (
                            {"fix_or": float(or_plot[0])}
                            if variable == "SF"
                            else {"fix_sf": float(sf_plot[0])}
                        )
                        summary = _marginal_levels_with_spread(
                            beta, train_names, sel_vars, config, sf_ref, or_ref,
                            cd, level_kind=variable, levels=levels,
                            **kwargs, **held,
                        )
                        grid = np.arange(levels.size, dtype=np.float64)
                    n_grid = grid.size
                    for j in range(n_grid):
                        rows.append({
                            "probe_id": probe_id,
                            "cluster_id": fit.cluster_id,
                            "glm_type": GLM_TYPE,
                            "model": model_label,
                            "variable": variable,
                            "condition": cond,
                            "grid_x": float(grid[j]),
                            "gain_mean": float(summary["mean"][j]),
                            "gain_q05": float(summary["q05"][j]),
                            "gain_q25": float(summary["q25"][j]),
                            "gain_q75": float(summary["q75"][j]),
                            "gain_q95": float(summary["q95"][j]),
                            "gain_std": float(summary["std"][j]),
                            "n_trials_in_bin": int(summary["n_trials"][j]),
                        })
    if not rows:
        logger.info("no per-trial uncertainty rows to write")
        return
    out_df = pd.DataFrame(rows)
    path = out_dir / "tuning_curves_uncertainty.csv"
    out_df.to_csv(path, index=False)
    sparse = int((out_df["n_trials_in_bin"] < 3).sum())
    logger.info(
        "wrote %s (%d rows; %d grid-points have <3 trials and will plot no band)",
        path.name, len(out_df), sparse,
    )


def _vars_from_names_for_model(
    model_label: str, train_names: list[str]
) -> list[str]:
    """Variables present on the design columns of ``train_names``.

    Mirrors ``plots._vars_from_names`` and ``plots`` model_vars dict but
    avoids importing those private helpers here.
    """
    vars_present: list[str] = []
    for prefix in ("Speed", "TF", "SF", "OR",
                   "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
                   "TF_x_SF", "TF_x_OR", "SF_x_OR"):
        if any(n.startswith(prefix + "_") or n == prefix for n in train_names):
            vars_present.append(prefix)
    return vars_present


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
    formatted_mat_path: Path,
    plot_format: str,
    plot_clusters: int | None,
    backend: str,
) -> None:
    from rc2_glm import plots
    from rc2_glm.precomputed_bins import load_precomputed_bin_edges

    figs_dir.mkdir(parents=True, exist_ok=True)

    precomputed_bins = load_precomputed_bin_edges(formatted_mat_path)

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

    for path in plots.save_figure(
        plots.plot_cv_bps_deltas(comparison_df),
        figs_dir / "cv_bps_deltas",
        fmt=plot_format,
    ):
        logger.info("wrote %s", path.relative_to(figs_dir.parent))

    # Speed-profile CV comparison: only populated when the pipeline ran
    # with --profile-cv-diagnostic. The plot function renders an
    # explanatory stub when columns are all-NaN, so this call is safe
    # to make unconditionally.
    if config.profile_cv_diagnostic:
        for path in plots.save_figure(
            plots.plot_speed_profile_cv_comparison(comparison_df),
            figs_dir / "speed_profile_cv_comparison",
            fmt=plot_format,
        ):
            logger.info("wrote %s", path.relative_to(figs_dir.parent))

    fits_to_plot = (
        cluster_fits if plot_clusters is None else cluster_fits[:plot_clusters]
    )
    for fit, _ in fits_to_plot:
        cid = fit.cluster_id
        # Build a small per-cluster coefficient DataFrame matching the
        # `glm_coefficients.csv` schema so plot_cluster_kernels can read
        # one cluster's worth without hitting disk.
        kernel_rows: list[dict] = []
        for model_name, betas in fit.model_betas.items():
            if betas is None:
                continue
            for col_name, beta in zip(fit.model_col_names[model_name], betas):
                kernel_rows.append({
                    "cluster_id": cid, "model": model_name,
                    "coefficient": col_name, "estimate": float(beta),
                })
        cluster_coef_df = pd.DataFrame(kernel_rows)

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
                precomputed_bins=precomputed_bins,
                model_predictions=fit.model_predictions,
             ), f"cluster_{cid}_tuning"),
            (plots.plot_cluster_kernels(
                probe_id=probe_id,
                cluster_id=cid,
                coef_df_cluster=cluster_coef_df,
                config=config,
             ), f"cluster_{cid}_kernels"),
        ):
            for path in plots.save_figure(fn, figs_dir / name, fmt=plot_format):
                logger.info("wrote %s", path.relative_to(figs_dir.parent))

    # Trial-level predictions: MATLAB Section 8b picks the top 4 clusters with
    # ≥2 selected vars by Selected cv_bps. Fewer than 4 is fine if the session
    # doesn't produce enough rich models — we just render what we have.
    trial_level_candidates = [
        fit for fit, _ in cluster_fits
        if len(fit.selection.selected_vars) >= 2
        and np.isfinite(fit.selection.final_cv_bps)
    ]
    trial_level_candidates.sort(
        key=lambda f: f.selection.final_cv_bps, reverse=True,
    )
    for fit in trial_level_candidates[:4]:
        cid = fit.cluster_id
        for path in plots.save_figure(
            plots.plot_trial_level_predictions(
                probe_id=probe_id,
                cluster_id=cid,
                cluster_df=fit.cluster_df,
                model_predictions=fit.model_predictions,
                config=config,
                cv_bps=fit.selection.final_cv_bps,
                selected_vars=list(fit.selection.selected_vars),
                seed=0,
            ),
            figs_dir / f"trial_level_cluster_{cid}",
            fmt=plot_format,
        ):
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
             "RC2_FORMATTED_DATA_DIR). Optional: if omitted, runs every "
             "*.mat in RC2_FORMATTED_DATA_DIR (per-probe subdirs + "
             "concatenated top-level CSVs).",
    )
    parser.add_argument(
        "out_dir", nargs="?", default=env_out,
        help="Directory for CSV outputs (default: $RC2_GLM_OUTPUT_DIR)",
    )
    parser.add_argument(
        "--single-probe", action="store_true",
        help="In multi-probe mode (no positional mat + RC2_FORMATTED_DATA_DIR), "
             "this flag restricts the run to a single probe — resolved via "
             "the same rules as before (RC2_FORMATTED_DATA_PATH, first *.mat). "
             "Useful for debugging without changing the env.",
    )
    parser.add_argument(
        "--backend", choices=("irls", "nemos"), default="irls",
        help="Fitting backend (default: irls)",
    )
    parser.add_argument(
        "--device", choices=("auto", "cpu", "gpu"), default="auto",
        help="JAX device for the nemos backend. 'auto' leaves JAX to "
             "pick (uses GPU if visible, else CPU). Ignored when "
             "--backend=irls. Applied before any JAX import — don't "
             "import jax in a caller before invoking rc2-glm.",
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
        help=(
            "Cap per-cluster plots to the first N retained clusters "
            "(default: no cap, plot every retained cluster). The retained "
            "set is controlled by the prefilter + --cluster-filter-csv; "
            "those are the right knobs to limit figure volume. Use "
            "--plot-clusters 0 to skip per-cluster plots entirely "
            "(aggregate figures still render)."
        ),
    )
    parser.add_argument(
        "--n-jobs", type=int, default=-1,
        help="Parallel fit workers for the per-cluster loop "
             "(joblib loky / process backend). -1 = all cores (default), "
             "1 = serial (use this to see per-cluster worker logs).",
    )
    parser.add_argument(
        "--tuning-curve-mode",
        choices=("trial-averaged", "steady-state"),
        default="trial-averaged",
        help=(
            "Tuning-curve rendering. 'trial-averaged' (default) marginalises "
            "the onset kernel and unfixed covariates over each cluster's "
            "observed per-condition motion bins. 'steady-state' pins the "
            "onset kernel at t=1.5s (MATLAB convention), kept for back-compat."
        ),
    )
    parser.add_argument(
        "--no-history", dest="include_history", action="store_false",
        help=(
            "Disable the spike-history term (default ON since 2026-04-29). "
            "Without history, the model is intercept + onset (if enabled) "
            "+ forward-selected stimulus terms — i.e. the legacy MATLAB-"
            "equivalent model. Use for one-off comparisons or legacy reruns."
        ),
    )
    # Deprecated alias — kept for back-compat with scripts/notebooks that
    # still pass --include-history. No-op since the new default is True.
    parser.add_argument(
        "--include-history", dest="include_history", action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--with-onset-kernel", dest="include_onset_kernel", action="store_true",
        help=(
            "Re-enable the onset kernel (default OFF since 2026-04-29). "
            "The prompt-03 ablation showed onset adds ~0 CV-bps once "
            "spike history is in; default flipped to False. Use this flag "
            "for legacy MATLAB-parity reruns or ablation experiments."
        ),
    )
    # Deprecated alias — kept for back-compat with scripts that pass
    # --no-onset-kernel. Default flipped True/False/True on
    # 2026-04-23 / 2026-04-29 / 2026-04-30; both flags remain valid.
    parser.add_argument(
        "--no-onset-kernel", dest="include_onset_kernel", action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--no-me-face", dest="include_me_face", action="store_false",
        help=(
            "Disable the face motion energy term ME_face (default ON since "
            "2026-04-30, prompt 06). Without ME_face, Phase-1 candidates "
            "are {Speed, TF, SF, OR}; ME_face_x_Speed is automatically "
            "ineligible in Phase 2 because its parent ME_face was never "
            "selected. Use this flag to compare 'with ME' vs 'without ME' "
            "runs on the same cluster set."
        ),
    )
    # Deprecated alias — kept for back-compat / scripts that pass
    # --include-me-face explicitly. No-op since the default is True.
    parser.add_argument(
        "--include-me-face", dest="include_me_face", action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--motion-energy-camera", choices=("camera0", "camera1"),
        default=None,
        help=(
            "Which camera trace to use for ME_face / ME_body. 'camera0' "
            "(default, face cam) is the prompt-06 v1 choice. 'camera1' "
            "(body cam) is the alternative for ablation runs. Overrides "
            "config.motion_energy_camera."
        ),
    )
    # Defaults updated 2026-04-30 (prompt 06) to match the new GLMConfig
    # defaults: history OFF, onset ON, ME_face ON. Ensures argparse-level
    # defaults stay in sync with config-level defaults so a bare
    # ``rc2-glm <probe> <out>`` call uses the new world without surprises.
    parser.set_defaults(
        include_history=False,
        include_onset_kernel=True,
        include_me_face=True,
    )
    parser.add_argument(
        "--tuning-curve-uncertainty",
        choices=("none", "simulated", "covariate-spread", "wide-quantile",
                 "std", "iqr"),
        default="simulated",
        help=(
            "Per-trial uncertainty band on the model rows of the "
            "tuning-curve panels (cluster_<id>_tuning.pdf). Modes: "
            "'simulated' (default, since 2026-04-28) = parametric "
            "bootstrap — predict λ at the 100 ms training granularity, "
            "draw y_sim ~ Poisson(λ·Δt), collapse simulated rates to the "
            "cache's display bins, band = IQR across trials per display "
            "bin (mean across MC iterations). Directly comparable to the "
            "Observed row's whiskers. "
            "'covariate-spread' = previous default — IQR across trials "
            "of the predicted rate at fixed sweep-x given each trial's "
            "actual non-target covariates (model-structural sensitivity). "
            "'wide-quantile' / 'std' / 'none' as before. 'iqr' is a "
            "deprecated alias for 'covariate-spread'. Only applies in "
            "trial-averaged mode."
        ),
    )
    parser.add_argument(
        "--cv-strategy",
        choices=("condition-stratified", "speed-profile"),
        default="condition-stratified",
        help=(
            "CV fold-assignment strategy. 'condition-stratified' (default) "
            "= k-fold over (trial_id, condition) pairs, MATLAB parity. "
            "'speed-profile' = 2-fold train-on-profile-1 / test-on-profile-2 "
            "(generalisation across the two reproduced velocity trajectories, "
            "MATLAB glm_single_cluster_analysis.m:2291-2293)."
        ),
    )
    parser.add_argument(
        "--profile-cv-diagnostic", dest="profile_cv_diagnostic",
        action="store_true",
        help=(
            "Also run a post-hoc 2-fold speed-profile CV on each cluster's "
            "Null, Selected, and Selected-without-Speed model. Keeps the "
            "primary condition-stratified forward selection unchanged; "
            "adds 3 extra profile_cv_bps_* columns to glm_model_comparison.csv "
            "and renders figs/speed_profile_cv_comparison.pdf. Mirrors "
            "MATLAB glm_single_cluster_analysis.m:2246-2367."
        ),
    )
    parser.add_argument(
        "--bin-width", type=float, default=None,
        help=(
            "Time-bin width in seconds (overrides config.time_bin_width, "
            "default 0.1 = 100 ms). Exploration knob for the Phase-B bin-"
            "width sweep. Absolute cv_bps values scale with log(bin_width); "
            "Spearman / classification / tuning-shape comparisons are "
            "invariant. 0.02 matches the MATLAB reference."
        ),
    )
    parser.add_argument(
        "--cv-seed", type=int, default=None,
        help=(
            "Seed for the condition-stratified fold assignment (overrides "
            "config.cv_seed, default 0). Exploration knob for the Phase-C "
            "pseudo-random-stability sweep: re-running with different seeds "
            "quantifies how much the fits depend on a specific trial-to-fold "
            "permutation vs. being stable under re-shuffling."
        ),
    )
    parser.add_argument(
        "--lambda-ridge", type=float, default=None,
        help=(
            "Ridge L2 penalty on non-intercept columns (overrides "
            "config.lambda_ridge, default 1e-3). Exploration knob for "
            "regularisation sweeps. λ=10⁻³ was tuned for 100 ms; at finer "
            "bin widths (e.g. 20 ms = 5× more bins per cluster) a stronger "
            "λ may be needed to suppress noise in the History-filter "
            "coefficients. Try λ ∈ {1e-3, 1e-2, 1e-1, 1.0}."
        ),
    )
    parser.add_argument(
        "--n-history-bases", type=int, default=None,
        help=(
            "Number of log-spaced raised-cosine history bases (overrides "
            "config.n_history_bases, default 10). At 100 ms only 2 lag "
            "bins are distinct so 10 bases over-parametrise the filter; "
            "5 may be better. At 20 ms 10 bases were data-limited per the "
            "2026-04-29 ridge sweep — fewer may identify cleaner shapes."
        ),
    )
    parser.add_argument(
        "--cluster-filter-csv", type=str, default=None,
        help=(
            "Restrict the pipeline to the (probe_id, cluster_id) pairs in "
            "this CSV (must have those two columns). Typical use: point at "
            "MATLAB's upstream-filtered prefilter_decision_tree.csv so "
            "Python processes the same ~100 clusters MATLAB did instead of "
            "all VISp. Clusters absent from the CSV are silently skipped "
            "at load time; prefilter + GLM never see them. See "
            "results/motion-clouds/findings-so-far.md for context."
        ),
    )
    parser.set_defaults(
        make_plots=True, prefilter=True, profile_cv_diagnostic=False,
    )
    args = parser.parse_args(argv)

    if not args.out_dir:
        raise SystemExit(
            "rc2-glm: no output directory. Pass as second positional arg "
            "or set RC2_GLM_OUTPUT_DIR in python/.env."
        )

    # Build config via replace() so None overrides fall through to the
    # GLMConfig defaults unchanged. Both --bin-width and --cv-seed are
    # exploration knobs; passing them None is a no-op.
    config_overrides: dict[str, object] = {
        "apply_prefilter": args.prefilter,
        "device": args.device,
        "tuning_curve_mode": args.tuning_curve_mode,
        "tuning_curve_uncertainty": args.tuning_curve_uncertainty,
        "include_history": args.include_history,
        "include_onset_kernel": args.include_onset_kernel,
        "include_me_face": args.include_me_face,
        "cv_strategy": args.cv_strategy,
        "profile_cv_diagnostic": args.profile_cv_diagnostic,
    }
    if args.motion_energy_camera is not None:
        config_overrides["motion_energy_camera"] = args.motion_energy_camera
    if args.bin_width is not None:
        config_overrides["time_bin_width"] = float(args.bin_width)
    if args.cv_seed is not None:
        config_overrides["cv_seed"] = int(args.cv_seed)
    if args.lambda_ridge is not None:
        config_overrides["lambda_ridge"] = float(args.lambda_ridge)
    if args.n_history_bases is not None:
        config_overrides["n_history_bases"] = int(args.n_history_bases)
    config = replace(GLMConfig(), **config_overrides)

    lookup = None
    if args.mc_sequence and args.mc_folders:
        from rc2_formatted_data_reader import StimulusLookup
        lookup = StimulusLookup(args.mc_sequence, args.mc_folders)

    # Parse the cluster-filter CSV once; look up per-probe cluster sets
    # inside the multi-probe / single-probe paths below. CSV must have
    # ``probe_id`` and ``cluster_id`` columns; other columns are ignored.
    cluster_filter_by_probe: dict[str, set[int]] | None = None
    if args.cluster_filter_csv:
        cf_df = pd.read_csv(args.cluster_filter_csv)
        if "cluster_id" not in cf_df.columns or "probe_id" not in cf_df.columns:
            raise SystemExit(
                f"rc2-glm: --cluster-filter-csv {args.cluster_filter_csv!r} "
                "must have 'probe_id' and 'cluster_id' columns"
            )
        cluster_filter_by_probe = {
            str(pid): {int(c) for c in grp["cluster_id"].dropna()}
            for pid, grp in cf_df.groupby("probe_id")
        }
        print(f"rc2-glm: cluster_filter active — "
              f"{sum(len(v) for v in cluster_filter_by_probe.values())} clusters "
              f"across {len(cluster_filter_by_probe)} probes from "
              f"{args.cluster_filter_csv}")

    # Multi-probe default: with no positional mat argument and
    # RC2_FORMATTED_DATA_DIR set, run every *.mat in that directory into
    # a per-probe subdir of out_dir, then concatenate into top-level CSVs
    # so downstream (rc2-glm-compare, the notebook) sees one consolidated
    # run. Laura 2026-04-23: "all probes should be always the default."
    # Empty string is also treated as "no mat argument" so users can
    # override out_dir positionally (``rc2-glm "" /path/to/out``).
    no_mat = args.mat is None or args.mat == ""
    multi_probe = (
        no_mat
        and not args.single_probe
        and env_dir
        and Path(env_dir).expanduser().is_dir()
        and len(sorted(Path(env_dir).expanduser().glob("*.mat"))) > 1
    )
    if multi_probe:
        mats = sorted(Path(env_dir).expanduser().glob("*.mat"))
        out_root = Path(args.out_dir)
        out_root.mkdir(parents=True, exist_ok=True)
        runs_root = out_root / "_runs"
        runs_root.mkdir(exist_ok=True)
        total_rows = 0
        for i, mat_path in enumerate(mats, start=1):
            probe_dir = runs_root / mat_path.stem
            probe_dir.mkdir(exist_ok=True)
            print(f"rc2-glm: [{i}/{len(mats)}] running {mat_path.name} → {probe_dir}")
            per_probe_filter = (
                cluster_filter_by_probe.get(mat_path.stem)
                if cluster_filter_by_probe is not None else None
            )
            result = run_pipeline(
                mat_path=mat_path,
                config=config,
                output_dir=probe_dir,
                stimulus_lookup=lookup,
                backend=args.backend,
                visp_only=not args.all_clusters,
                make_plots=args.make_plots,
                plot_format=args.plot_format,
                plot_clusters=args.plot_clusters,
                n_jobs=args.n_jobs,
                cluster_filter=per_probe_filter,
            )
            total_rows += len(result.model_comparison)
        _aggregate_probe_runs(runs_root, out_root)
        print(f"rc2-glm: wrote {total_rows} cluster rows across "
              f"{len(mats)} probes to {out_root}")
        return 0

    mat_path = _resolve_mat_path(args.mat, env_path, env_dir)
    per_probe_filter = (
        cluster_filter_by_probe.get(mat_path.stem)
        if cluster_filter_by_probe is not None else None
    )
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
        n_jobs=args.n_jobs,
        cluster_filter=per_probe_filter,
    )
    # run_pipeline has already torn down its handlers; a bare print is fine.
    print(f"rc2-glm: wrote {len(result.model_comparison)} cluster rows to {args.out_dir}")
    return 0


def _aggregate_probe_runs(runs_root: Path, out_root: Path) -> None:
    """Concatenate per-probe CSVs into top-level aggregate CSVs.

    The per-probe subdirs under ``runs_root`` each carry a full set of
    pipeline outputs (glm_model_comparison.csv, glm_coefficients.csv,
    glm_selection_history.csv, prefilter_decision_tree.csv, plus the
    diagnostics/ subdir). We concatenate the top-level CSVs and the
    ``diagnostics/stationary_vs_motion_fr_python.csv`` trial-level CSV
    so ``rc2-glm-compare`` sees an aggregated run without needing a
    separate aggregate step.
    """
    top_level = (
        "glm_model_comparison.csv",
        "glm_coefficients.csv",
        "glm_selection_history.csv",
        "prefilter_decision_tree.csv",
    )
    for basename in top_level:
        dfs = [
            pd.read_csv(probe / basename)
            for probe in sorted(runs_root.iterdir())
            if (probe / basename).is_file()
        ]
        if not dfs:
            continue
        pd.concat(dfs, ignore_index=True).to_csv(out_root / basename, index=False)

    diag = out_root / "diagnostics"
    diag.mkdir(exist_ok=True)
    svm_dfs = [
        pd.read_csv(probe / "diagnostics" / "stationary_vs_motion_fr_python.csv")
        for probe in sorted(runs_root.iterdir())
        if (probe / "diagnostics" / "stationary_vs_motion_fr_python.csv").is_file()
    ]
    if svm_dfs:
        pd.concat(svm_dfs, ignore_index=True).to_csv(
            diag / "stationary_vs_motion_fr_python.csv", index=False,
        )
    tl_dfs = [
        pd.read_csv(probe / "diagnostics" / "trial_level_metrics.csv")
        for probe in sorted(runs_root.iterdir())
        if (probe / "diagnostics" / "trial_level_metrics.csv").is_file()
    ]
    if tl_dfs:
        pd.concat(tl_dfs, ignore_index=True).to_csv(
            diag / "trial_level_metrics.csv", index=False,
        )
    alpha_dfs = [
        pd.read_csv(probe / "diagnostics" / "history_alpha_sweep.csv")
        for probe in sorted(runs_root.iterdir())
        if (probe / "diagnostics" / "history_alpha_sweep.csv").is_file()
    ]
    if alpha_dfs:
        pd.concat(alpha_dfs, ignore_index=True).to_csv(
            diag / "history_alpha_sweep.csv", index=False,
        )

    # Aggregate forward-selection summary across all probes. Used to live
    # under validation/ via rc2-glm-compare; that tool was retired with
    # MATLAB parity (2026-04-29). Now generated here so the production
    # output dir always carries one canonical "what fraction of clusters
    # selected each variable" plot at the 4-probe level.
    agg_cmp_path = out_root / "glm_model_comparison.csv"
    if agg_cmp_path.is_file():
        try:
            from rc2_glm.plots import plot_forward_selection_summary
            agg_df = pd.read_csv(agg_cmp_path)
            figs_dir = out_root / "figs"
            figs_dir.mkdir(exist_ok=True)
            fig = plot_forward_selection_summary(agg_df)
            fig.suptitle(
                f"Forward-selection summary — {len(sorted(runs_root.iterdir()))} probes, "
                f"{len(agg_df)} clusters",
                fontsize=11,
            )
            fig.savefig(figs_dir / "forward_selection_summary.pdf")
            fig.savefig(figs_dir / "forward_selection_summary.png", dpi=150)
            import matplotlib.pyplot as _plt
            _plt.close(fig)
            logger.info("wrote aggregate %s",
                        (figs_dir / "forward_selection_summary.pdf").relative_to(out_root))
        except Exception as exc:
            logger.warning("aggregate forward-selection summary failed: %s", exc)


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
