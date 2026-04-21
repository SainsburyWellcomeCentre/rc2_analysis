"""Smoke tests for the matplotlib figure builders.

These don't verify visual output — they just confirm the builders run
against tiny inputs, return a :class:`matplotlib.figure.Figure`, and
that :func:`save_figure` writes to disk in the requested format.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from rc2_glm.config import GLMConfig
from rc2_glm.plots import (
    plot_basis_functions,
    plot_cluster_model_overview,
    plot_forward_selection_summary,
    plot_tuning_curves,
    save_figure,
)


def test_plot_basis_functions_returns_figure():
    fig = plot_basis_functions(GLMConfig())
    assert isinstance(fig, Figure)
    n_axes = sum(1 for ax in fig.axes if ax.has_data() or ax.lines)
    assert n_axes >= 1
    plt.close(fig)


def test_plot_forward_selection_summary_minimal_inputs():
    comparison = pd.DataFrame([
        {
            "probe_id": "p0", "cluster_id": 1, "n_trials": 10,
            "time_Null_cv_bps": 0.0, "time_Selected_cv_bps": 0.12,
            "time_n_selected_vars": 2,
        },
        {
            "probe_id": "p0", "cluster_id": 2, "n_trials": 10,
            "time_Null_cv_bps": 0.0, "time_Selected_cv_bps": 0.05,
            "time_n_selected_vars": 1,
        },
    ])
    history = pd.DataFrame([
        {"probe_id": "p0", "cluster_id": 1, "round": 1, "phase": 1,
         "best_candidate": "Speed", "delta_bps": 0.1, "added": True,
         "cv_bps_after": 0.1},
        {"probe_id": "p0", "cluster_id": 2, "round": 1, "phase": 1,
         "best_candidate": "TF", "delta_bps": 0.05, "added": True,
         "cv_bps_after": 0.05},
    ])
    fig = plot_forward_selection_summary(comparison, history)
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 3
    plt.close(fig)


def _synthetic_cluster_df(n_trials: int = 8, n_bins: int = 20) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    cfg = GLMConfig()
    rows = []
    conditions = ["T_Vstatic", "V", "VT"]
    for tid in range(n_trials):
        cond = conditions[tid % 3]
        # stationary prelude: 3 bins before motion onset
        for b in range(3):
            rows.append({
                "probe_id": "p0", "cluster_id": 1, "trial_id": tid,
                "condition": "stationary", "speed": 0.0, "tf": 0.0,
                "sf": np.nan, "orientation": np.nan, "batch_gain": 1.0,
                "spike_count": int(rng.poisson(0.5)),
                "time_in_trial": b * cfg.time_bin_width,
                "time_since_onset": -((3 - b) * cfg.time_bin_width),
            })
        for b in range(n_bins):
            speed = rng.uniform(0.0, cfg.speed_range[1])
            tf_v = 0.0 if cond == "T_Vstatic" else rng.uniform(0.0, cfg.tf_range[1])
            sf = np.nan if cond == "T_Vstatic" else cfg.sf_levels[tid % len(cfg.sf_levels)]
            ori = np.nan if cond == "T_Vstatic" else cfg.or_levels[tid % len(cfg.or_levels)]
            rows.append({
                "probe_id": "p0", "cluster_id": 1, "trial_id": tid,
                "condition": cond, "speed": speed, "tf": tf_v,
                "sf": sf, "orientation": ori, "batch_gain": 1.0,
                "spike_count": int(rng.poisson(1.0)),
                "time_in_trial": b * cfg.time_bin_width,
                "time_since_onset": b * cfg.time_bin_width,
            })
    return pd.DataFrame(rows)


def _synthetic_model_dicts(cluster_df: pd.DataFrame) -> tuple[
    dict[str, np.ndarray], dict[str, list[str]], dict[str, np.ndarray],
]:
    """Fake trained GLM output aligned with the motion rows of ``cluster_df``.

    The column layout matches what ``assemble_design_matrix_selected`` would
    produce for the given ``selected_vars`` — tiny betas let the caller verify
    the plots wire up β, col names, and per-row predictions without running a
    real fit. ``n_motion`` is the count of rows where
    ``condition != 'stationary'``.
    """
    n_motion = int((cluster_df["condition"] != "stationary").sum())
    rng = np.random.default_rng(42)
    cols_null = ["Intercept", "Onset_1", "Onset_2"]
    cols_selected = cols_null + ["Speed_1", "Speed_2", "TF_1", "TF_2"]
    cols_additive = cols_selected + ["SF_0.0060", "SF_0.0120",
                                     "OR_0.000", "OR_0.785", "OR_1.571"]
    cols_full = cols_additive + ["Spd1_x_TF1", "Spd1_x_TF2"]
    betas = {
        "Null": rng.normal(0.0, 0.2, size=len(cols_null)),
        "Selected": rng.normal(0.0, 0.2, size=len(cols_selected)),
        "Additive": rng.normal(0.0, 0.2, size=len(cols_additive)),
        "FullInteraction": rng.normal(0.0, 0.2, size=len(cols_full)),
    }
    col_names = {
        "Null": cols_null,
        "Selected": cols_selected,
        "Additive": cols_additive,
        "FullInteraction": cols_full,
    }
    preds = {
        label: rng.uniform(0.5, 5.0, size=n_motion)
        for label in ("Null", "Selected", "Additive", "FullInteraction")
    }
    return betas, col_names, preds


def test_plot_tuning_curves_shape_and_finiteness():
    cluster_df = _synthetic_cluster_df()
    betas, col_names, _ = _synthetic_model_dicts(cluster_df)
    fig = plot_tuning_curves(
        probe_id="p0", cluster_id=1, cluster_df=cluster_df,
        model_betas=betas, model_col_names=col_names,
        config=GLMConfig(), n_sweep=10,
    )
    assert isinstance(fig, Figure)
    assert len(fig.axes) >= 20
    for ax in fig.axes:
        for line in ax.lines:
            ys = np.asarray(line.get_ydata(), dtype=float)
            ys = ys[~np.isnan(ys)]
            if ys.size:
                assert np.all(np.isfinite(ys))
    plt.close(fig)


def test_plot_cluster_model_overview_shape_and_finiteness():
    cluster_df = _synthetic_cluster_df()
    betas, col_names, preds = _synthetic_model_dicts(cluster_df)
    fig = plot_cluster_model_overview(
        probe_id="p0", cluster_id=1,
        cluster_df=cluster_df,
        selected_vars=["Speed", "TF"],
        model_betas=betas,
        model_col_names=col_names,
        model_predictions=preds,
        config=GLMConfig(),
    )
    assert isinstance(fig, Figure)
    # 4 rows × 6 cols = 24 main axes (colorbars on cols 2-3 add a few more).
    assert len(fig.axes) >= 24
    for ax in fig.axes:
        for line in ax.lines:
            ys = np.asarray(line.get_ydata(), dtype=float)
            ys = ys[~np.isnan(ys)]
            if ys.size:
                assert np.all(np.isfinite(ys))
    plt.close(fig)


def test_plot_cluster_model_overview_handles_missing_model():
    cluster_df = _synthetic_cluster_df()
    betas, col_names, preds = _synthetic_model_dicts(cluster_df)
    for d in (betas, col_names, preds):
        d.pop("FullInteraction")
    fig = plot_cluster_model_overview(
        probe_id="p0", cluster_id=1,
        cluster_df=cluster_df,
        selected_vars=["Speed", "TF"],
        model_betas=betas,
        model_col_names=col_names,
        model_predictions=preds,
        config=GLMConfig(),
    )
    assert isinstance(fig, Figure)
    plt.close(fig)


def test_save_figure_pdf(tmp_path):
    fig = plot_basis_functions(GLMConfig())
    written = save_figure(fig, tmp_path / "basis", fmt="pdf")
    assert len(written) == 1
    assert written[0].suffix == ".pdf"
    assert written[0].exists()


def test_save_figure_both_formats(tmp_path):
    fig = plot_basis_functions(GLMConfig())
    written = save_figure(fig, tmp_path / "basis", fmt="both")
    assert len(written) == 2
    suffixes = sorted(p.suffix for p in written)
    assert suffixes == [".pdf", ".png"]
    for p in written:
        assert p.exists()
