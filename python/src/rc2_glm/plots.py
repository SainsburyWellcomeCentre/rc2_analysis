"""Matplotlib figures mirroring (a subset of) the MATLAB GLM script.

Covers summary + cheap per-cluster panels and the MATLAB Fig 5
reconstructed tuning curves (lines 3298-3980 of
``scripts/glm_single_cluster_analysis.m``). Figs 3 / 6 / 7 / 8 are not
ported — the Python pipeline does not yet export the per-fold /
per-trial traces those plots need.

All figure helpers return a :class:`matplotlib.figure.Figure`. Use
:func:`save_figure` to write PDF and optional PNG to disk.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless backend for CLI use
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis
from rc2_glm.config import GLMConfig
from rc2_glm.design_matrix import assemble_design_matrix_selected

logger = logging.getLogger("rc2_glm")

MODEL_LABELS = ("Null", "Selected", "Additive", "FullInteraction")

# Matches MATLAB cond_colors_rc in glm_single_cluster_analysis.m line 3335.
COND_COLORS: dict[str, tuple[float, float, float]] = {
    "T_Vstatic": (0.2, 0.7, 0.2),
    "V": (0.9, 0.6, 0.1),
    "VT": (0.0, 0.4, 0.8),
    "stationary": (0.0, 0.0, 0.0),
}

MAX_PREDICTED_FR = 400.0   # Hz; matches MATLAB max_predicted_fr soft ceiling
OBS_FR_SMOOTH_WIDTH = 0.1  # seconds; matches MATLAB obs_fr_smooth_width (line 129)
MIN_TRIALS_PER_CELL = 1    # divergence from MATLAB (≥2): show data even with one trial


# --------------------------------------------------------------------------- #
# Save helper
# --------------------------------------------------------------------------- #


def save_figure(fig: Figure, path_base: Path, fmt: str = "pdf") -> list[Path]:
    """Write ``fig`` to ``path_base`` in the requested format(s).

    ``fmt`` is one of ``"pdf"``, ``"png"``, ``"both"``.
    """
    path_base = Path(path_base)
    path_base.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    if fmt in ("pdf", "both"):
        pdf_path = path_base.with_suffix(".pdf")
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
        written.append(pdf_path)

    if fmt in ("png", "both"):
        png_path = path_base.with_suffix(".png")
        fig.savefig(png_path, format="png", dpi=150, bbox_inches="tight")
        written.append(png_path)

    plt.close(fig)
    return written


# --------------------------------------------------------------------------- #
# Fig 1 — basis functions
# --------------------------------------------------------------------------- #


def plot_basis_functions(config: GLMConfig) -> Figure:
    """Overview of the four basis families used by the GLM."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)

    speed_grid = np.linspace(config.speed_range[0], config.speed_range[1], 400)
    B_speed = raised_cosine_basis(speed_grid, config.n_speed_bases, *config.speed_range)
    for i in range(B_speed.shape[1]):
        axes[0, 0].plot(speed_grid, B_speed[:, i], label=f"basis {i + 1}")
    axes[0, 0].set_title(f"Speed raised cosine ({config.n_speed_bases} bases)")
    axes[0, 0].set_xlabel("speed (cm/s)")
    axes[0, 0].set_ylabel("activation")

    tf_grid = np.linspace(config.tf_range[0], config.tf_range[1], 400)
    B_tf = raised_cosine_basis(tf_grid, config.n_tf_bases, *config.tf_range)
    for i in range(B_tf.shape[1]):
        axes[0, 1].plot(tf_grid, B_tf[:, i], label=f"basis {i + 1}")
    axes[0, 1].set_title(f"TF raised cosine ({config.n_tf_bases} bases)")
    axes[0, 1].set_xlabel("TF (Hz)")
    axes[0, 1].set_ylabel("activation")

    sf_levels = np.asarray(config.sf_levels, dtype=np.float64)
    or_levels = np.asarray(config.or_levels, dtype=np.float64)
    axes[1, 0].vlines(sf_levels, 0, 1, colors="#1f77b4", label="SF levels")
    axes[1, 0].vlines(or_levels, 0, 0.5, colors="#d62728", label="OR levels (rad)")
    axes[1, 0].scatter(sf_levels, np.ones_like(sf_levels), marker="D", color="#1f77b4")
    axes[1, 0].scatter(or_levels, np.full_like(or_levels, 0.5),
                       marker="s", color="#d62728")
    axes[1, 0].set_title("SF / OR indicator levels")
    axes[1, 0].set_xlabel("value")
    axes[1, 0].set_yticks([])
    axes[1, 0].legend(loc="upper right", fontsize=8)

    onset_grid = np.linspace(0.0, config.onset_range[1], 400)
    B_onset = onset_kernel_basis(onset_grid, config.n_onset_bases, config.onset_range[1])
    for i in range(B_onset.shape[1]):
        axes[1, 1].plot(onset_grid, B_onset[:, i], label=f"basis {i + 1}")
    axes[1, 1].set_title(f"Onset kernel ({config.n_onset_bases} bases)")
    axes[1, 1].set_xlabel("t since onset (s)")
    axes[1, 1].set_ylabel("activation")

    fig.suptitle("GLM basis functions")
    return fig


# --------------------------------------------------------------------------- #
# Fig 2 — forward-selection summary across clusters
# --------------------------------------------------------------------------- #


def plot_forward_selection_summary(
    comparison_df: pd.DataFrame, history_df: pd.DataFrame,
) -> Figure:
    """Three panels: null-vs-selected scatter, #selected histogram, var counts."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)

    null_col = "time_Null_cv_bps"
    sel_col = "time_Selected_cv_bps"
    n_sel_col = "time_n_selected_vars"

    ax = axes[0]
    if not comparison_df.empty and null_col in comparison_df and sel_col in comparison_df:
        ax.scatter(
            comparison_df[null_col], comparison_df[sel_col],
            s=25, color="#1f77b4", alpha=0.7, edgecolors="none",
        )
        lo = float(min(comparison_df[null_col].min(), comparison_df[sel_col].min()))
        hi = float(max(comparison_df[null_col].max(), comparison_df[sel_col].max()))
        ax.plot([lo, hi], [lo, hi], linestyle=":", color="#888")
    ax.set_title("Null vs Selected CV bps")
    ax.set_xlabel("null CV bps")
    ax.set_ylabel("selected CV bps")

    ax = axes[1]
    if not comparison_df.empty and n_sel_col in comparison_df:
        counts = comparison_df[n_sel_col].value_counts().sort_index()
        ax.bar(counts.index.astype(int), counts.values, color="#2ca02c")
    ax.set_title("Selected variable count per cluster")
    ax.set_xlabel("# selected vars")
    ax.set_ylabel("# clusters")

    ax = axes[2]
    if not history_df.empty and "added" in history_df.columns:
        picked = history_df[history_df["added"]]["best_candidate"]
        var_counts = picked.value_counts()
        ax.bar(range(len(var_counts)), var_counts.values, color="#d62728")
        ax.set_xticks(range(len(var_counts)))
        ax.set_xticklabels(var_counts.index, rotation=45, ha="right")
    ax.set_title("How often each variable was picked")
    ax.set_xlabel("variable")
    ax.set_ylabel("# times picked")

    fig.suptitle("Forward-selection summary across clusters")
    return fig


# --------------------------------------------------------------------------- #
# Fig 5 — reconstructed tuning curves (MATLAB Section 8b)
# --------------------------------------------------------------------------- #


def plot_tuning_curves(
    probe_id: str,
    cluster_id: int,
    cluster_df: pd.DataFrame,
    model_betas: dict[str, np.ndarray],
    model_col_names: dict[str, list[str]],
    config: GLMConfig,
    n_sweep: int = 20,
) -> Figure:
    """Reconstructed tuning curves matching MATLAB Fig 5 (5×4 grid).

    Rows: Observed / Null / Selected / Additive / FullInteraction.
    Columns: Speed / TF / SF / OR.

    Uses the trained β vectors (from ``ClusterFit.model_betas``) — no
    refit. For each sweep we build a prediction design matrix with
    ``is_prediction=True`` (no zero-variance drop) and align its columns
    by name with the training design stored in ``model_col_names``:
    missing columns are filled with zeros and extras dropped. This
    guarantees ``X_aligned @ β`` is well-defined regardless of which
    columns the training fit happened to drop.

    Each panel shows per-condition traces (T_Vstatic green, V orange,
    VT blue; black stationary dot at speed=0), following the MATLAB
    layout in ``scripts/glm_single_cluster_analysis.m`` lines
    3298-3980. Held-constant covariates follow the MATLAB convention:

    - T_Vstatic: TF=0, SF=reference, OR=reference (visual off).
    - V: Speed=0, SF=mode, OR=mode (visual only).
    - VT: Speed=mean, TF=mean, SF=mode, OR=mode (multi-sensory).

    The onset kernel is evaluated at ``steady_state_time = 1.5s`` for
    all motion predictions (MATLAB line 3397) and at all-zeros for the
    stationary dot (MATLAB line 3743).

    Observed row: per-trial aggregation matching MATLAB (lines
    2873, 2958, 3031, 3108, 3210). For each bin/level, rows within that
    condition are first grouped by trial; the mean rate per trial is
    computed; the bin's displayed value is the mean ± std across trials
    (≥2 trials) or a single trial's value (no error bar) if only one
    trial contributed.
    """
    fig, axes = plt.subplots(
        5, 4, figsize=(16, 14), sharex="col", constrained_layout=True,
    )
    row_labels = ("Observed", "Null", "Selected", "Additive", "FullInteraction")
    col_labels = ("Speed", "TF", "SF", "OR")
    for r, lbl in enumerate(row_labels):
        axes[r, 0].set_ylabel(f"{lbl}\nrate (Hz)", fontsize=9)
    for c, lbl in enumerate(col_labels):
        axes[0, c].set_title(lbl, fontsize=10)
    for c, xlbl in enumerate(
        ("Speed (cm/s)", "TF (Hz)", "SF (cpd)", "Orientation (rad)")
    ):
        axes[4, c].set_xlabel(xlbl)

    motion = cluster_df[cluster_df["condition"] != "stationary"]
    if motion.empty:
        fig.suptitle(
            f"No motion bins for {probe_id} cluster {cluster_id}",
            fontsize=12,
        )
        return fig

    speed_train = motion["speed"].to_numpy(dtype=np.float64)
    tf_train = motion["tf"].to_numpy(dtype=np.float64)
    sf_train = motion["sf"].to_numpy(dtype=np.float64)
    or_train = motion["orientation"].to_numpy(dtype=np.float64)

    sf_ref = np.asarray(config.sf_levels, dtype=np.float64)
    or_ref = np.asarray(config.or_levels, dtype=np.float64)

    mean_speed = float(np.nanmean(speed_train)) if speed_train.size else 0.0
    mean_tf_motion = float(np.nanmean(tf_train[tf_train > 0])) \
        if np.any(tf_train > 0) else 0.0
    nonnan_sf = sf_train[~np.isnan(sf_train)]
    mode_sf = float(_mode(nonnan_sf)) if nonnan_sf.size else float(sf_ref[0])
    nonnan_or = or_train[~np.isnan(or_train)]
    mode_or = float(_mode(nonnan_or)) if nonnan_or.size else float(or_ref[0])

    spd_centres = np.linspace(config.speed_range[0], config.speed_range[1], n_sweep)
    tf_centres = np.linspace(config.tf_range[0], config.tf_range[1], n_sweep)

    sf_plot = sf_ref.copy()
    or_plot = or_ref.copy()

    B_onset_steady = onset_kernel_basis(
        np.array([1.5]), config.n_onset_bases, config.onset_range[1],
    )  # shape (1, n_onset)
    B_onset_stat = np.zeros((1, config.n_onset_bases))

    _plot_observed_row(
        axes[0, :], cluster_df, config,
        sf_plot=sf_plot, or_plot=or_plot,
    )

    model_vars = {
        "Null": [],
        "Selected": _vars_from_names(model_col_names.get("Selected", [])),
        "Additive": ["Speed", "TF", "SF", "OR"],
        "FullInteraction": ["Speed", "TF", "SF", "OR",
                            "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
                            "TF_x_SF", "TF_x_OR", "SF_x_OR"],
    }

    for row_idx, label in enumerate(MODEL_LABELS, start=1):
        row_axes = axes[row_idx, :]
        beta = model_betas.get(label)
        train_names = model_col_names.get(label)
        if beta is None or train_names is None:
            for ax in row_axes:
                ax.text(0.5, 0.5, "not fitted", ha="center", va="center",
                        transform=ax.transAxes, color="#888", fontsize=10)
            continue
        _predict_model_row(
            row_axes, label, beta, train_names, model_vars[label], config,
            B_onset_steady=B_onset_steady, B_onset_stat=B_onset_stat,
            sf_ref=sf_ref, or_ref=or_ref,
            sf_plot=sf_plot, or_plot=or_plot,
            mean_speed=mean_speed, mean_tf=mean_tf_motion,
            mode_sf=mode_sf, mode_or=mode_or,
            spd_centres=spd_centres, tf_centres=tf_centres,
        )

    for col in range(4):
        hi = 0.0
        for row in range(5):
            yl = axes[row, col].get_ylim()
            if np.isfinite(yl[1]):
                hi = max(hi, yl[1])
        if hi > 0:
            for row in range(5):
                axes[row, col].set_ylim(0.0, hi * 1.05)

    handles = [
        plt.Line2D([0], [0], marker="o", color=COND_COLORS["T_Vstatic"],
                   linestyle="-", label="T_Vstatic"),
        plt.Line2D([0], [0], marker="o", color=COND_COLORS["V"],
                   linestyle="-", label="V"),
        plt.Line2D([0], [0], marker="o", color=COND_COLORS["VT"],
                   linestyle="-", label="VT"),
        plt.Line2D([0], [0], marker="o", color=COND_COLORS["stationary"],
                   linestyle="", label="stationary"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))

    sel_names = model_col_names.get("Selected", [])
    sel_vars = _vars_from_names(sel_names)
    fig.suptitle(
        f"Reconstructed tuning curves — {probe_id} cluster {cluster_id}"
        f"    selected: {'+'.join(sel_vars) if sel_vars else 'Null'}",
        fontsize=12,
    )
    return fig


# --------------------------------------------------------------------------- #
# Fig 4 — per-cluster model overview with scatter panels
# (MATLAB lines 2727-3290 of glm_single_cluster_analysis.m)
# --------------------------------------------------------------------------- #


# Beta-swarm group layout — matches MATLAB compute_grouped_params (line 6328).
# Keys are MATLAB group tags; values are the palette row and tick label.
_BETA_GROUP_ORDER: tuple[str, ...] = (
    "Intercept", "Speed", "TF", "SF", "OR", "Time",
    "Spd x TF", "Spd x SF", "Spd x OR", "Other",
)
_BETA_GROUP_COLORS: dict[str, tuple[float, float, float]] = {
    "Intercept":  (0.5, 0.5, 0.5),
    "Speed":      (0.17, 0.63, 0.17),
    "TF":         (1.0, 0.50, 0.05),
    "SF":         (0.95, 0.85, 0.10),
    "OR":         (0.84, 0.15, 0.16),
    "Time":       (0.3, 0.8, 0.8),
    "Spd x TF":   (0.6, 0.35, 0.05),
    "Spd x SF":   (0.55, 0.50, 0.05),
    "Spd x OR":   (0.50, 0.10, 0.10),
    "Other":      (0.7, 0.7, 0.7),
}

# Discrete SF / OR palettes — copied verbatim from MATLAB lines 2661-2671.
_SF_PALETTE_LEVELS: tuple[float, ...] = (0.003, 0.006, 0.012)
_SF_PALETTE_COLORS: tuple[tuple[float, float, float], ...] = (
    (0.2, 0.6, 0.2),   # 0.003 cpd = green
    (0.2, 0.4, 0.9),   # 0.006 cpd = blue
    (0.9, 0.2, 0.2),   # 0.012 cpd = red
)
_OR_PALETTE_LEVELS: tuple[float, ...] = (
    -np.pi / 4, 0.0, np.pi / 4, np.pi / 2,
)
_OR_PALETTE_COLORS: tuple[tuple[float, float, float], ...] = (
    (0.9, 0.1, 0.1),    # -45° bright red
    (0.1, 0.8, 0.9),    #   0° cyan
    (0.95, 0.55, 0.0),  #  45° orange
    (0.6, 0.1, 0.85),   #  90° purple
)


def plot_cluster_model_overview(
    probe_id: str,
    cluster_id: int,
    cluster_df: pd.DataFrame,
    selected_vars: list[str],
    model_betas: dict[str, np.ndarray],
    model_col_names: dict[str, list[str]],
    model_predictions: dict[str, np.ndarray],
    config: GLMConfig,
) -> Figure:
    """MATLAB Fig 4 — 4×6 per-cluster model overview with scatter panels.

    Rows: Null / Selected / Additive / FullInteraction.
    Cols: β swarm | cond R² | speed R² | TF R² | SF R² | OR R².

    Observed firing rate is per-trial boxcar-smoothed (width
    ``OBS_FR_SMOOTH_WIDTH``) before binning; per-bin/level cells require
    ≥2 contributing trials (trial-first aggregation) and display the
    mean ± std across trials. The Selected-model row is highlighted in
    red; the column-0 title gets a red ``★`` prefix. Missing models
    (e.g. FullInteraction dropped because ``X.shape[1] >= y.size``)
    render "N/A" in every panel of that row. Mirrors MATLAB
    ``scripts/glm_single_cluster_analysis.m`` lines 2727-3290.
    """
    fig, axes = plt.subplots(4, 6, figsize=(18, 11), constrained_layout=True)

    motion_mask = cluster_df["condition"] != "stationary"
    motion_df = cluster_df.loc[motion_mask].reset_index(drop=True)
    if motion_df.empty:
        fig.suptitle(
            f"No motion bins for {probe_id} cluster {cluster_id}", fontsize=12,
        )
        for ax in axes.ravel():
            ax.set_axis_off()
        return fig

    obs_fr = _smoothed_obs_per_trial(cluster_df, config.time_bin_width,
                                     OBS_FR_SMOOTH_WIDTH)
    obs_fr_motion = obs_fr[motion_mask.to_numpy()]

    selected_label = "Selected"

    vt_sub = motion_df[motion_df["condition"] == "VT"]
    logger.info(
        "overview: cluster %d | motion_rows=%d (T_Vstatic=%d V=%d VT=%d) | "
        "VT speed=[%.2f,%.2f] tf=[%.3f,%.3f] | obs_fr=[%.1f,%.1f] Hz",
        cluster_id, len(motion_df),
        int((motion_df["condition"] == "T_Vstatic").sum()),
        int((motion_df["condition"] == "V").sum()),
        int((motion_df["condition"] == "VT").sum()),
        float(vt_sub["speed"].min()) if len(vt_sub) else float("nan"),
        float(vt_sub["speed"].max()) if len(vt_sub) else float("nan"),
        float(vt_sub["tf"].min()) if len(vt_sub) else float("nan"),
        float(vt_sub["tf"].max()) if len(vt_sub) else float("nan"),
        float(np.nanmin(obs_fr_motion)) if obs_fr_motion.size else float("nan"),
        float(np.nanmax(obs_fr_motion)) if obs_fr_motion.size else float("nan"),
    )

    for r, label in enumerate(MODEL_LABELS):
        row_axes = axes[r, :]
        beta = model_betas.get(label)
        col_names = model_col_names.get(label)
        preds = model_predictions.get(label)
        is_winner = (label == selected_label)

        if beta is None or col_names is None or preds is None:
            for ax in row_axes:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=10, color="#888")
                ax.set_xticks([])
                ax.set_yticks([])
            _apply_row_style(row_axes, label, is_winner, n_coefs=None)
            continue

        _plot_beta_swarm(row_axes[0], np.asarray(beta), list(col_names))
        _plot_condition_scatter(row_axes[1], motion_df, obs_fr_motion, preds)
        _plot_vt_speed_scatter(row_axes[2], motion_df, obs_fr_motion, preds,
                               config.speed_range)
        _plot_vt_tf_scatter(row_axes[3], motion_df, obs_fr_motion, preds,
                            config.tf_range)
        _plot_vt_level_scatter(
            row_axes[4], motion_df, obs_fr_motion, preds,
            value_col="sf", levels=np.asarray(_SF_PALETTE_LEVELS),
            palette=[np.asarray(c) for c in _SF_PALETTE_COLORS],
            reference_levels=np.asarray(_SF_PALETTE_LEVELS),
            axis_title="SF", legend_fmt=lambda v: f"{v:.3f}",
            tol=1e-4,
        )
        _plot_vt_level_scatter(
            row_axes[5], motion_df, obs_fr_motion, preds,
            value_col="orientation", levels=np.asarray(_OR_PALETTE_LEVELS),
            palette=[np.asarray(c) for c in _OR_PALETTE_COLORS],
            reference_levels=np.asarray(_OR_PALETTE_LEVELS),
            axis_title="OR", legend_fmt=lambda v: f"{int(round(np.degrees(v)))}°",
            tol=1e-3,
        )
        _apply_row_style(row_axes, label, is_winner, n_coefs=len(col_names))

    for ax in axes[:, 1:].ravel():
        ax.set_xlabel("Pred", fontsize=7)
        ax.set_ylabel("Obs", fontsize=7)
    for ax in axes.ravel():
        ax.tick_params(labelsize=6)

    fig.suptitle(
        f"Probe: {probe_id}  Cluster: {cluster_id}  "
        f"|  Selected: {'+'.join(selected_vars) if selected_vars else 'Null'}",
        fontsize=10, fontweight="bold",
    )
    return fig


def _apply_row_style(
    row_axes, label: str, is_winner: bool, n_coefs: int | None,
) -> None:
    """Apply red spines + title color for the Selected row, plain for others."""
    swarm_ax = row_axes[0]
    if n_coefs is None:
        title_str = label
    else:
        title_str = f"{label} ({n_coefs} coefs)"
    if is_winner:
        title_str = "★ " + title_str
        swarm_ax.set_title(title_str, fontsize=8, color="red", fontweight="bold")
        for ax in row_axes:
            for spine in ax.spines.values():
                spine.set_edgecolor("red")
                spine.set_linewidth(1.0)
            ax.tick_params(colors="red", labelsize=6)
    else:
        swarm_ax.set_title(title_str, fontsize=8)


def _plot_beta_swarm(ax, beta: np.ndarray, col_names: list[str]) -> None:
    """Grouped β swarm: basis-number text labels + group-mean lines.

    Groups match MATLAB ``compute_grouped_params`` (line 6328). Basis
    number extraction matches ``extract_basis_label`` (line 6391): main
    effect → last digits of the name; interaction → ``"{i}-{j}"``.
    """
    groups: dict[str, list[tuple[float, str]]] = {g: [] for g in _BETA_GROUP_ORDER}
    for b, cn in zip(beta, col_names):
        tag = _beta_group_tag(cn)
        groups[tag].append((float(b), _extract_basis_label(cn)))

    nonempty = [g for g in _BETA_GROUP_ORDER if groups[g]]
    if not nonempty:
        ax.text(0.5, 0.5, "no coefficients", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="#888")
        return

    for xi, tag in enumerate(nonempty, start=1):
        members = groups[tag]
        color = _BETA_GROUP_COLORS[tag]
        bs = np.array([m[0] for m in members], dtype=np.float64)
        labels = [m[1] or "-" for m in members]
        if bs.size > 1:
            jitter = np.linspace(-0.3, 0.3, bs.size)
        else:
            jitter = np.array([0.0])
        # Draw a visible dot per coefficient so every basis shows up
        # regardless of how cramped the text labels end up.
        ax.scatter(
            xi + jitter, bs, s=14, c=[color], edgecolors="black",
            linewidths=0.3, zorder=3,
        )
        for jx, b, lbl in zip(jitter, bs, labels):
            ax.text(
                xi + jx + 0.03, b, lbl, ha="left", va="center",
                fontsize=6, color=color, fontweight="bold", zorder=4,
            )
        ax.plot([xi - 0.35, xi + 0.35], [bs.mean(), bs.mean()],
                "-", color=color, linewidth=1.4, zorder=2)

    ax.axhline(0.0, color="black", linewidth=0.3)
    ax.set_xlim(0.3, len(nonempty) + 0.8)
    ax.set_xticks(np.arange(1, len(nonempty) + 1))
    ax.set_xticklabels(nonempty, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("β", fontsize=7)
    all_betas = np.concatenate(
        [np.array([b for b, _ in groups[g]], dtype=np.float64) for g in nonempty]
    )
    if all_betas.size:
        pad = 0.1 * (np.abs(all_betas).max() + 1e-6)
        ax.set_ylim(float(all_betas.min() - pad), float(all_betas.max() + pad))


def _beta_group_tag(col_name: str) -> str:
    """Map a design column name onto one of ``_BETA_GROUP_ORDER``."""
    if col_name == "Intercept":
        return "Intercept"
    if "_x_TF" in col_name:
        return "Spd x TF"
    if "_x_SF" in col_name:
        return "Spd x SF"
    if "_x_OR" in col_name:
        return "Spd x OR"
    if col_name.startswith("Speed_"):
        return "Speed"
    if col_name.startswith("TF_"):
        return "TF"
    if col_name.startswith("SF_"):
        return "SF"
    if col_name.startswith("OR_"):
        return "OR"
    if col_name.startswith(("Onset_", "Time_")):
        return "Time"
    return "Other"


def _extract_basis_label(col_name: str) -> str:
    """Return the basis index string (MATLAB ``extract_basis_label``)."""
    if col_name == "Intercept":
        return ""
    if "_x_" in col_name:
        a, b = col_name.split("_x_", 1)
        na = re.findall(r"\d+", a)
        nb = re.findall(r"\d+", b)
        if na and nb:
            return f"{na[-1]}-{nb[-1]}"
        return ""
    nums = re.findall(r"\d+", col_name)
    return nums[-1] if nums else ""


def _smoothed_obs_per_trial(
    cluster_df: pd.DataFrame, bin_width: float, smooth_width: float,
) -> np.ndarray:
    """Per-row boxcar-smoothed observed firing rate (Hz).

    Matches MATLAB lines 2738-2762: for each trial, sort by
    ``time_in_trial``, convolve ``spike_count / bin_width`` with a
    uniform boxcar of length ``round(smooth_width / bin_width)`` (``mode='same'``),
    then unsort back to original row order.
    """
    n_smooth = max(1, int(round(smooth_width / bin_width)))
    kernel = np.ones(n_smooth, dtype=np.float64) / n_smooth
    out = np.full(cluster_df.shape[0], np.nan, dtype=np.float64)
    rate = cluster_df["spike_count"].to_numpy(dtype=np.float64) / bin_width
    time_in_trial = cluster_df["time_in_trial"].to_numpy(dtype=np.float64)
    trial_ids = cluster_df["trial_id"].to_numpy()
    row_index = np.arange(cluster_df.shape[0])
    for tid in np.unique(trial_ids):
        sel = np.where(trial_ids == tid)[0]
        if sel.size == 0:
            continue
        order = sel[np.argsort(time_in_trial[sel])]
        smoothed = np.convolve(rate[order], kernel, mode="same")
        out[order] = smoothed
    return out


def _r_squared(obs: np.ndarray, pred: np.ndarray) -> float:
    """Coefficient of determination. Returns NaN if obs variance is too small
    for the statistic to be meaningful (arbitrary prediction then dominates
    ss_tot and R² blows to −1e3 or more)."""
    ss_res = float(np.sum((obs - pred) ** 2))
    ss_tot = float(np.sum((obs - obs.mean()) ** 2))
    if ss_tot <= 1e-8 or obs.size < 2:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def _r2_title(prefix: str, obs: np.ndarray, pred: np.ndarray) -> str:
    """Format R² for a panel title, falling back to a readable string when
    the obs mean-variation is too small to anchor a meaningful R²."""
    r2 = _r_squared(obs, pred)
    if np.isnan(r2):
        return f"{prefix} obs≈const"
    if r2 < -99.0:
        return f"{prefix} R²<<0"
    return f"{prefix} R²={r2:.2f}"


def _plot_condition_scatter(
    ax, motion_df: pd.DataFrame, obs_fr: np.ndarray, preds: np.ndarray,
) -> None:
    """Predicted vs observed mean FR per condition (T_Vstatic / V / VT).

    Shows any condition with ≥``MIN_TRIALS_PER_CELL`` trials. Single-trial
    conditions render with ``std=0`` (no error bar). R² computed across
    the contributing condition means; the title falls back to a plain
    label when all condition means are nearly equal.
    """
    order = ("T_Vstatic", "V", "VT")
    mean_obs: list[float] = []
    mean_pred: list[float] = []
    std_obs: list[float] = []
    std_pred: list[float] = []
    colors: list[tuple[float, float, float]] = []
    labels: list[str] = []
    cond = motion_df["condition"].to_numpy()
    trial_ids = motion_df["trial_id"].to_numpy()
    for c in order:
        mask = cond == c
        if not mask.any():
            continue
        tr_obs, tr_pred = _trial_means_within_mask(
            mask, trial_ids, obs_fr, preds,
        )
        if tr_obs.size < MIN_TRIALS_PER_CELL:
            continue
        mean_obs.append(float(tr_obs.mean()))
        mean_pred.append(float(tr_pred.mean()))
        std_obs.append(float(tr_obs.std(ddof=0)) if tr_obs.size >= 2 else 0.0)
        std_pred.append(float(tr_pred.std(ddof=0)) if tr_pred.size >= 2 else 0.0)
        colors.append(COND_COLORS[c])
        labels.append(c)

    if not mean_obs:
        ax.text(0.5, 0.5, "no trials", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#888")
        return

    for xp, yo, xerr, yerr, clr, lbl in zip(
        mean_pred, mean_obs, std_pred, std_obs, colors, labels,
    ):
        ax.errorbar(xp, yo, yerr=yerr, xerr=xerr, fmt="o",
                    color=clr, markerfacecolor=clr, markersize=5,
                    linewidth=0.8, capsize=3, label=lbl)

    mx = max(max(mean_obs) + max(std_obs), max(mean_pred) + max(std_pred))
    if mx > 0:
        ax.plot([0.0, mx], [0.0, mx], "k--", linewidth=0.5)

    ax.set_title(
        _r2_title("Cond", np.asarray(mean_obs), np.asarray(mean_pred)),
        fontsize=7,
    )
    ax.legend(loc="best", fontsize=4, frameon=True)


def _adaptive_bin_edges(
    x: np.ndarray, config_range: tuple[float, float], n_bins: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """Pick speed/TF bin edges that actually span the data.

    The basis functions use ``config_range`` (e.g. (0, 50) cm/s for speed)
    but real recordings often occupy a small subset of that range (mice
    rarely hit 50 cm/s). Falling back to fixed config-range bins leaves
    9 of 10 bins empty and R² nan. Use the data range instead when it is
    meaningfully tighter than the config range.
    """
    good = x[~np.isnan(x)]
    if good.size == 0:
        edges = np.linspace(config_range[0], config_range[1], n_bins + 1)
        return edges, 0.5 * (edges[:-1] + edges[1:])
    data_lo = float(good.min())
    data_hi = float(good.max())
    cfg_span = config_range[1] - config_range[0]
    data_span = data_hi - data_lo
    if data_span <= 0:
        edges = np.linspace(config_range[0], config_range[1], n_bins + 1)
    elif data_span < 0.3 * cfg_span:
        pad = 0.05 * data_span
        edges = np.linspace(max(config_range[0], data_lo - pad),
                            min(config_range[1], data_hi + pad),
                            n_bins + 1)
    else:
        edges = np.linspace(config_range[0], config_range[1], n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return edges, centres


def _plot_vt_speed_scatter(
    ax,
    motion_df: pd.DataFrame,
    obs_fr: np.ndarray,
    preds: np.ndarray,
    config_range: tuple[float, float],
) -> None:
    """VT-only predicted vs observed scatter across speed bins."""
    vt_mask = motion_df["condition"].to_numpy() == "VT"
    if not vt_mask.any():
        ax.text(0.5, 0.5, "no VT bins", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#888")
        return
    x = motion_df.loc[vt_mask, "speed"].to_numpy(dtype=np.float64)
    bin_edges, bin_centres = _adaptive_bin_edges(x, config_range)
    _plot_binned_scatter(
        ax, x, obs_fr[vt_mask], preds[vt_mask],
        motion_df.loc[vt_mask, "trial_id"].to_numpy(),
        bin_edges, bin_centres,
        cmap=plt.get_cmap("viridis"),
        title_prefix="Spd", colorbar=True,
    )


def _plot_vt_tf_scatter(
    ax,
    motion_df: pd.DataFrame,
    obs_fr: np.ndarray,
    preds: np.ndarray,
    config_range: tuple[float, float],
) -> None:
    """VT-only predicted vs observed scatter across TF bins."""
    vt_mask = motion_df["condition"].to_numpy() == "VT"
    if not vt_mask.any():
        ax.text(0.5, 0.5, "no VT bins", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#888")
        return
    x = motion_df.loc[vt_mask, "tf"].to_numpy(dtype=np.float64)
    bin_edges, bin_centres = _adaptive_bin_edges(x, config_range)
    _plot_binned_scatter(
        ax, x, obs_fr[vt_mask], preds[vt_mask],
        motion_df.loc[vt_mask, "trial_id"].to_numpy(),
        bin_edges, bin_centres,
        cmap=plt.get_cmap("hot"),
        title_prefix="TF", colorbar=True,
    )


def _plot_binned_scatter(
    ax,
    x: np.ndarray,
    obs: np.ndarray,
    pred: np.ndarray,
    trial_ids: np.ndarray,
    bin_edges: np.ndarray,
    bin_centres: np.ndarray,
    *,
    cmap,
    title_prefix: str,
    colorbar: bool,
) -> None:
    obs_mean, obs_std, pred_mean, pred_std, valid = _trial_stats_by_bin(
        x, obs, pred, trial_ids, bin_edges,
    )
    if not valid.any():
        good = ~np.isnan(x)
        msg = (
            f"empty: x∈[{x[good].min():.3f},{x[good].max():.3f}]"
            f" bins∈[{bin_edges[0]:.2f},{bin_edges[-1]:.2f}]"
            if good.any() else "no x values"
        )
        logger.info("%s scatter empty — %s", title_prefix, msg)
        ax.text(0.5, 0.5, msg, ha="center", va="center",
                transform=ax.transAxes, fontsize=6, color="#888")
        return
    lo, hi = float(bin_centres.min()), float(bin_centres.max())
    norm = plt.Normalize(vmin=lo, vmax=hi if hi > lo else lo + 1.0)
    for i in np.flatnonzero(valid):
        color = cmap(norm(bin_centres[i]))
        ax.errorbar(pred_mean[i], obs_mean[i],
                    xerr=pred_std[i], yerr=obs_std[i],
                    fmt="o", color=color, markerfacecolor=color,
                    markersize=4, linewidth=0.8, capsize=2)
    mx = float(max(
        (obs_mean[valid] + obs_std[valid]).max(),
        (pred_mean[valid] + pred_std[valid]).max(),
    ))
    if mx > 0:
        ax.plot([0.0, mx], [0.0, mx], "k--", linewidth=0.5)
    if colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = plt.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
        cb.ax.tick_params(labelsize=5)
    ax.set_title(
        _r2_title(title_prefix, obs_mean[valid], pred_mean[valid]),
        fontsize=7,
    )


def _plot_vt_level_scatter(
    ax,
    motion_df: pd.DataFrame,
    obs_fr: np.ndarray,
    preds: np.ndarray,
    *,
    value_col: str,
    levels: np.ndarray,
    palette: list[np.ndarray],
    reference_levels: np.ndarray,
    axis_title: str,
    legend_fmt,
    tol: float,
) -> None:
    """VT-only scatter across discrete SF / OR levels."""
    vt_mask = motion_df["condition"].to_numpy() == "VT"
    if not vt_mask.any():
        ax.text(0.5, 0.5, "no VT bins", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#888")
        return
    x = motion_df.loc[vt_mask, value_col].to_numpy(dtype=np.float64)
    obs_vt = obs_fr[vt_mask]
    pred_vt = preds[vt_mask]
    trial_ids_vt = motion_df.loc[vt_mask, "trial_id"].to_numpy()

    obs_mean, obs_std, pred_mean, pred_std, valid = _trial_stats_by_level(
        x, obs_vt, pred_vt, trial_ids_vt, levels, tol,
    )
    if not valid.any():
        ax.text(0.5, 0.5, "insufficient trials", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#888")
        return

    for i in np.flatnonzero(valid):
        # Map the observed level to the palette index via closest
        # reference-level match (MATLAB lines 3129, 3231).
        ci = int(np.argmin(np.abs(reference_levels - levels[i])))
        color = palette[ci]
        ax.errorbar(pred_mean[i], obs_mean[i],
                    xerr=pred_std[i], yerr=obs_std[i],
                    fmt="o", color=color, markerfacecolor=color,
                    markersize=5, linewidth=0.8, capsize=3,
                    label=legend_fmt(float(levels[i])))
    mx = float(max(
        (obs_mean[valid] + obs_std[valid]).max(),
        (pred_mean[valid] + pred_std[valid]).max(),
    ))
    if mx > 0:
        ax.plot([0.0, mx], [0.0, mx], "k--", linewidth=0.5)
    if valid.sum() >= 2:
        ax.set_title(
            _r2_title(axis_title, obs_mean[valid], pred_mean[valid]),
            fontsize=7,
        )
    else:
        ax.set_title(axis_title, fontsize=7)
    ax.legend(loc="best", fontsize=4, frameon=True)


def _trial_means_within_mask(
    mask: np.ndarray,
    trial_ids: np.ndarray,
    obs: np.ndarray,
    pred: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Mean obs / pred FR per trial inside ``mask``; drops trials with all-NaN."""
    tr_obs_list: list[float] = []
    tr_pred_list: list[float] = []
    for tid in np.unique(trial_ids[mask]):
        sub = mask & (trial_ids == tid)
        if not sub.any():
            continue
        o = obs[sub]
        p = pred[sub]
        o_m = np.nanmean(o) if np.any(~np.isnan(o)) else np.nan
        p_m = np.nanmean(p) if np.any(~np.isnan(p)) else np.nan
        if np.isnan(o_m) or np.isnan(p_m):
            continue
        tr_obs_list.append(float(o_m))
        tr_pred_list.append(float(p_m))
    return np.asarray(tr_obs_list), np.asarray(tr_pred_list)


def _trial_stats_by_bin(
    x: np.ndarray,
    obs: np.ndarray,
    pred: np.ndarray,
    trial_ids: np.ndarray,
    bin_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Trial-first aggregation per ``bin_edges`` bucket.

    Uses the ``MIN_TRIALS_PER_CELL`` threshold (default 1) so every cell
    with any contributing trial is plotted — a divergence from MATLAB's
    ≥2-trial gate. Single-trial cells render with ``std=0`` (no error bar).
    Returns ``(obs_mean, obs_std, pred_mean, pred_std, valid)`` each of
    length ``len(bin_edges) - 1``.
    """
    n_bins = bin_edges.size - 1
    obs_mean = np.full(n_bins, np.nan)
    obs_std = np.full(n_bins, np.nan)
    pred_mean = np.full(n_bins, np.nan)
    pred_std = np.full(n_bins, np.nan)
    valid = np.zeros(n_bins, dtype=bool)
    good = ~np.isnan(x) & ~np.isnan(obs) & ~np.isnan(pred)
    x_g = x[good]
    obs_g = obs[good]
    pred_g = pred[good]
    tid_g = trial_ids[good]
    if x_g.size == 0:
        return obs_mean, obs_std, pred_mean, pred_std, valid
    inside = (x_g >= bin_edges[0]) & (x_g < bin_edges[-1])
    bin_idx = np.where(inside,
                       np.clip(np.digitize(x_g, bin_edges) - 1, 0, n_bins - 1),
                       -1)
    for b in range(n_bins):
        sel = bin_idx == b
        if not sel.any():
            continue
        tr_obs, tr_pred = _trial_means_within_mask(
            sel, tid_g, obs_g, pred_g,
        )
        if tr_obs.size >= MIN_TRIALS_PER_CELL:
            obs_mean[b] = tr_obs.mean()
            obs_std[b] = tr_obs.std(ddof=0) if tr_obs.size >= 2 else 0.0
            pred_mean[b] = tr_pred.mean()
            pred_std[b] = tr_pred.std(ddof=0) if tr_pred.size >= 2 else 0.0
            valid[b] = True
    return obs_mean, obs_std, pred_mean, pred_std, valid


def _trial_stats_by_level(
    x: np.ndarray,
    obs: np.ndarray,
    pred: np.ndarray,
    trial_ids: np.ndarray,
    levels: np.ndarray,
    tol: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Trial-first aggregation at each discrete level (SF / OR)."""
    n = levels.size
    obs_mean = np.full(n, np.nan)
    obs_std = np.full(n, np.nan)
    pred_mean = np.full(n, np.nan)
    pred_std = np.full(n, np.nan)
    valid = np.zeros(n, dtype=bool)
    good = ~np.isnan(x) & ~np.isnan(obs) & ~np.isnan(pred)
    x_g = x[good]
    obs_g = obs[good]
    pred_g = pred[good]
    tid_g = trial_ids[good]
    if x_g.size == 0:
        return obs_mean, obs_std, pred_mean, pred_std, valid
    for i, level in enumerate(levels):
        sel = np.abs(x_g - level) < tol
        if not sel.any():
            continue
        tr_obs, tr_pred = _trial_means_within_mask(
            sel, tid_g, obs_g, pred_g,
        )
        if tr_obs.size >= MIN_TRIALS_PER_CELL:
            obs_mean[i] = tr_obs.mean()
            obs_std[i] = tr_obs.std(ddof=0) if tr_obs.size >= 2 else 0.0
            pred_mean[i] = tr_pred.mean()
            pred_std[i] = tr_pred.std(ddof=0) if tr_pred.size >= 2 else 0.0
            valid[i] = True
    return obs_mean, obs_std, pred_mean, pred_std, valid


# --------------------------------------------------------------------------- #
# Observed + prediction helpers (Fig 5)
# --------------------------------------------------------------------------- #


def _plot_observed_row(
    row_axes,
    cluster_df: pd.DataFrame,
    config: GLMConfig,
    *,
    sf_plot: np.ndarray,
    or_plot: np.ndarray,
) -> None:
    """Observed rates per condition, aggregated trial-first (matches MATLAB).

    Rates for each bin/level are computed by first grouping within-condition
    rows by ``trial_id`` and taking the mean ``spike_count / bin_width`` per
    trial within that bin, then averaging across trials with std for the
    error bar. Single-trial bins render as a dot with no error bar; empty
    bins are skipped.
    """
    bw = config.time_bin_width
    ax_spd, ax_tf, ax_sf, ax_or = row_axes

    spd_edges = np.linspace(config.speed_range[0], config.speed_range[1], 11)
    spd_centres = 0.5 * (spd_edges[:-1] + spd_edges[1:])
    tf_edges = np.linspace(config.tf_range[0], config.tf_range[1], 11)
    tf_centres = 0.5 * (tf_edges[:-1] + tf_edges[1:])

    # Speed: T_Vstatic + VT (MATLAB skips V for speed since its speed=0)
    for cond in ("T_Vstatic", "VT"):
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        mean, std = _trial_mean_std_by_bin(
            sub, value_col="speed", bin_edges=spd_edges, bin_width=bw,
        )
        good = ~np.isnan(mean)
        if good.any():
            ax_spd.errorbar(spd_centres[good], mean[good], yerr=std[good],
                            fmt="o-", color=COND_COLORS[cond],
                            markersize=4, linewidth=1, capsize=3)

    # Stationary dot at speed=0 — per-trial mean then across trials
    stat = cluster_df[cluster_df["condition"] == "stationary"]
    if not stat.empty:
        per_trial = (
            stat.groupby("trial_id")["spike_count"]
            .mean() / bw
        ).to_numpy(dtype=np.float64)
        per_trial = per_trial[~np.isnan(per_trial)]
        if per_trial.size >= 1:
            mean_stat = float(per_trial.mean())
            std_stat = float(per_trial.std(ddof=0)) if per_trial.size >= 2 else 0.0
            ax_spd.errorbar(0.0, mean_stat, yerr=std_stat,
                            fmt="o", color=COND_COLORS["stationary"],
                            markersize=6, capsize=4)

    # TF: V + VT (skip T_Vstatic since TF=0)
    for cond in ("V", "VT"):
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        mean, std = _trial_mean_std_by_bin(
            sub, value_col="tf", bin_edges=tf_edges, bin_width=bw,
        )
        good = ~np.isnan(mean)
        if good.any():
            ax_tf.errorbar(tf_centres[good], mean[good], yerr=std[good],
                           fmt="o-", color=COND_COLORS[cond],
                           markersize=4, linewidth=1, capsize=3)

    # SF: V + VT
    for cond in ("V", "VT"):
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        mean, std = _trial_mean_std_by_level(
            sub, value_col="sf", levels=sf_plot, tol=0.001, bin_width=bw,
        )
        good = ~np.isnan(mean)
        if good.any():
            ax_sf.errorbar(np.arange(sf_plot.size)[good], mean[good],
                           yerr=std[good],
                           fmt="o-", color=COND_COLORS[cond],
                           markersize=5, linewidth=1, capsize=4)
    ax_sf.set_xticks(np.arange(sf_plot.size))
    ax_sf.set_xticklabels([f"{v:.3f}" for v in sf_plot], fontsize=7)

    # OR: V + VT
    for cond in ("V", "VT"):
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        mean, std = _trial_mean_std_by_level(
            sub, value_col="orientation", levels=or_plot, tol=0.01, bin_width=bw,
        )
        good = ~np.isnan(mean)
        if good.any():
            ax_or.errorbar(np.arange(or_plot.size)[good], mean[good],
                           yerr=std[good],
                           fmt="o-", color=COND_COLORS[cond],
                           markersize=5, linewidth=1, capsize=4)
    ax_or.set_xticks(np.arange(or_plot.size))
    ax_or.set_xticklabels([_format_radians(v) for v in or_plot], fontsize=7)


def _predict_model_row(
    row_axes,
    label: str,
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    *,
    B_onset_steady: np.ndarray,
    B_onset_stat: np.ndarray,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    sf_plot: np.ndarray,
    or_plot: np.ndarray,
    mean_speed: float,
    mean_tf: float,
    mode_sf: float,
    mode_or: float,
    spd_centres: np.ndarray,
    tf_centres: np.ndarray,
) -> None:
    ax_spd, ax_tf, ax_sf, ax_or = row_axes
    n_spd = spd_centres.size
    n_tf = tf_centres.size

    B_onset_spd = np.broadcast_to(B_onset_steady, (n_spd, B_onset_steady.shape[1])).copy()
    B_onset_tf = np.broadcast_to(B_onset_steady, (n_tf, B_onset_steady.shape[1])).copy()

    # --- Speed sweep ---
    _plot_sweep(
        ax_spd, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec=spd_centres,
        tf_vec=np.zeros(n_spd),
        onset_mat=B_onset_spd,
        sf_vec=np.full(n_spd, np.nan),
        or_vec=np.full(n_spd, np.nan),
        color=COND_COLORS["T_Vstatic"],
        xs=spd_centres,
    )
    _plot_sweep(
        ax_spd, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec=spd_centres,
        tf_vec=np.full(n_spd, mean_tf),
        onset_mat=B_onset_spd,
        sf_vec=np.full(n_spd, mode_sf),
        or_vec=np.full(n_spd, mode_or),
        color=COND_COLORS["VT"],
        xs=spd_centres,
    )
    _plot_point(
        ax_spd, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed=0.0, tf=0.0, onset_row=B_onset_stat[0],
        sf=np.nan, orientation=np.nan,
        color=COND_COLORS["stationary"],
        x=0.0,
    )

    # --- TF sweep ---
    _plot_sweep(
        ax_tf, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec=np.zeros(n_tf),
        tf_vec=tf_centres,
        onset_mat=B_onset_tf,
        sf_vec=np.full(n_tf, mode_sf),
        or_vec=np.full(n_tf, mode_or),
        color=COND_COLORS["V"],
        xs=tf_centres,
    )
    _plot_sweep(
        ax_tf, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec=np.full(n_tf, mean_speed),
        tf_vec=tf_centres,
        onset_mat=B_onset_tf,
        sf_vec=np.full(n_tf, mode_sf),
        or_vec=np.full(n_tf, mode_or),
        color=COND_COLORS["VT"],
        xs=tf_centres,
    )

    # --- SF levels ---
    _plot_levels(
        ax_sf, beta, train_names, selected_vars, config, sf_ref, or_ref,
        level_kind="SF", levels=sf_plot,
        speed_const=0.0, tf_const=mean_tf,
        sf_fill=None, or_fill=mode_or,
        onset_row=B_onset_steady[0],
        color=COND_COLORS["V"],
    )
    _plot_levels(
        ax_sf, beta, train_names, selected_vars, config, sf_ref, or_ref,
        level_kind="SF", levels=sf_plot,
        speed_const=mean_speed, tf_const=mean_tf,
        sf_fill=None, or_fill=mode_or,
        onset_row=B_onset_steady[0],
        color=COND_COLORS["VT"],
    )
    ax_sf.set_xticks(np.arange(sf_plot.size))
    ax_sf.set_xticklabels([f"{v:.3f}" for v in sf_plot], fontsize=7)

    # --- OR levels ---
    _plot_levels(
        ax_or, beta, train_names, selected_vars, config, sf_ref, or_ref,
        level_kind="OR", levels=or_plot,
        speed_const=0.0, tf_const=mean_tf,
        sf_fill=mode_sf, or_fill=None,
        onset_row=B_onset_steady[0],
        color=COND_COLORS["V"],
    )
    _plot_levels(
        ax_or, beta, train_names, selected_vars, config, sf_ref, or_ref,
        level_kind="OR", levels=or_plot,
        speed_const=mean_speed, tf_const=mean_tf,
        sf_fill=mode_sf, or_fill=None,
        onset_row=B_onset_steady[0],
        color=COND_COLORS["VT"],
    )
    ax_or.set_xticks(np.arange(or_plot.size))
    ax_or.set_xticklabels([_format_radians(v) for v in or_plot], fontsize=7)


def _plot_sweep(
    ax,
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    *,
    speed_vec: np.ndarray,
    tf_vec: np.ndarray,
    onset_mat: np.ndarray,
    sf_vec: np.ndarray,
    or_vec: np.ndarray,
    color: tuple[float, float, float],
    xs: np.ndarray,
) -> None:
    rates = _predict_rate(
        beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec, tf_vec, onset_mat, sf_vec, or_vec,
    )
    if rates is None:
        return
    ax.plot(xs, rates, "o-", color=color,
            markersize=4, linewidth=1.2, markerfacecolor=color)


def _plot_point(
    ax,
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    *,
    speed: float,
    tf: float,
    onset_row: np.ndarray,
    sf: float,
    orientation: float,
    color: tuple[float, float, float],
    x: float,
) -> None:
    rates = _predict_rate(
        beta, train_names, selected_vars, config, sf_ref, or_ref,
        np.array([speed]), np.array([tf]), onset_row[None, :],
        np.array([sf]), np.array([orientation]),
    )
    if rates is None:
        return
    ax.plot(x, rates[0], "o", color=color, markersize=6,
            markerfacecolor=color)


def _plot_levels(
    ax,
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    *,
    level_kind: str,
    levels: np.ndarray,
    speed_const: float,
    tf_const: float,
    sf_fill: float | None,
    or_fill: float | None,
    onset_row: np.ndarray,
    color: tuple[float, float, float],
) -> None:
    n = levels.size
    if level_kind == "SF":
        sf_vec = levels.copy()
        or_vec = np.full(n, or_fill)
    else:
        sf_vec = np.full(n, sf_fill)
        or_vec = levels.copy()
    rates = _predict_rate(
        beta, train_names, selected_vars, config, sf_ref, or_ref,
        np.full(n, speed_const), np.full(n, tf_const),
        np.broadcast_to(onset_row, (n, onset_row.size)).copy(),
        sf_vec, or_vec,
    )
    if rates is None:
        return
    ax.plot(np.arange(n), rates, "o-", color=color,
            markersize=5, linewidth=1.2, markerfacecolor=color)


def _predict_rate(
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    speed_vec: np.ndarray,
    tf_vec: np.ndarray,
    onset_mat: np.ndarray,
    sf_vec: np.ndarray,
    or_vec: np.ndarray,
) -> np.ndarray | None:
    """Build a sweep design and return predicted firing rate in Hz.

    Training used ``offset = log(bin_width)`` so ``exp(X @ beta + offset)``
    is expected count per bin. Dividing by ``bin_width`` converts to Hz,
    which equals ``exp(X @ beta)`` by algebra.

    The prediction design is built with ``is_prediction=True`` (no
    zero-variance drop) and then aligned by column name with the training
    design: missing columns are filled with zero, extras are dropped. This
    makes the matrix-vector product safe regardless of which columns the
    training fit happened to drop.
    """
    B_speed = raised_cosine_basis(speed_vec, config.n_speed_bases, *config.speed_range)
    B_tf = raised_cosine_basis(tf_vec, config.n_tf_bases, *config.tf_range)
    X_pred, pred_names = assemble_design_matrix_selected(
        B_speed, B_tf, onset_mat, sf_vec, or_vec, selected_vars,
        sf_ref_levels=sf_ref, or_ref_levels=or_ref,
    )
    X_aligned = _align_prediction_columns(X_pred, pred_names, train_names)
    if X_aligned.shape[1] != beta.size:
        logger.debug(
            "tuning-curve predict: aligned %d cols vs %d β — skipping",
            X_aligned.shape[1], beta.size,
        )
        return None
    eta = np.clip(X_aligned @ beta, -20.0, 20.0)
    return np.minimum(np.exp(eta), MAX_PREDICTED_FR)


def _align_prediction_columns(
    X_pred: np.ndarray, pred_names: list[str], train_names: list[str],
) -> np.ndarray:
    """Reshape ``X_pred`` to have the same columns (in order) as training.

    Columns in ``train_names`` that are absent from ``pred_names`` are
    filled with zeros (reference-level behaviour). Columns in
    ``pred_names`` that are absent from ``train_names`` are dropped.
    """
    n_rows = X_pred.shape[0]
    name_to_pred_idx = {n: i for i, n in enumerate(pred_names)}
    out = np.zeros((n_rows, len(train_names)), dtype=np.float64)
    for j, name in enumerate(train_names):
        i = name_to_pred_idx.get(name)
        if i is not None:
            out[:, j] = X_pred[:, i]
    return out


_MAIN_PREFIX: dict[str, str] = {
    "Speed": "Speed_", "TF": "TF_", "SF": "SF_", "OR": "OR_",
}
_INTERACTION_PREFIX: dict[str, str] = {
    "Speed": "Spd", "TF": "TF", "SF": "SF", "OR": "OR",
}
_INTERACTION_VARS: tuple[str, ...] = (
    "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
    "TF_x_SF", "TF_x_OR", "SF_x_OR",
)


def _vars_from_names(col_names: list[str]) -> list[str]:
    """Infer which main-effect / interaction groups a column list contains.

    ``assemble_design_matrix_selected`` names main-effect columns with an
    underscore (``Speed_1``, ``TF_2``, ``SF_0.0060``, ``OR_0.000``) and
    interaction columns with ``_x_`` in them (``Spd1_x_TF2``,
    ``SF0.0060_x_OR0.000`` etc.). Ordering follows the canonical MATLAB
    order so suptitle labels read like ``Speed+TF+SF``.
    """
    present: list[str] = []
    for var, prefix in _MAIN_PREFIX.items():
        if any(n.startswith(prefix) and "_x_" not in n for n in col_names):
            present.append(var)
    for var in _INTERACTION_VARS:
        if any(_interaction_matches(n, var) for n in col_names):
            present.append(var)
    return present


def _interaction_matches(col_name: str, var: str) -> bool:
    """Return True if ``col_name`` belongs to interaction group ``var``."""
    parts = var.split("_x_")
    if len(parts) != 2:
        return False
    pa = _INTERACTION_PREFIX.get(parts[0], parts[0])
    pb = _INTERACTION_PREFIX.get(parts[1], parts[1])
    return col_name.startswith(pa) and f"_x_{pb}" in col_name


def _trial_mean_std_by_bin(
    sub_df: pd.DataFrame,
    *,
    value_col: str,
    bin_edges: np.ndarray,
    bin_width: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-trial mean rate inside each ``bin_edges`` bucket of ``value_col``.

    Divergence from MATLAB: we plot every cell with ≥``MIN_TRIALS_PER_CELL``
    contributing trials (default 1) so the observed data is always visible.
    MATLAB required ≥2 trials and silently dropped single-trial cells.
    Single-trial cells render with ``std=0`` (no error bar).
    """
    x = sub_df[value_col].to_numpy(dtype=np.float64)
    rate = sub_df["spike_count"].to_numpy(dtype=np.float64) / bin_width
    trial_ids = sub_df["trial_id"].to_numpy(dtype=np.int64)
    n_bins = bin_edges.size - 1
    mean = np.full(n_bins, np.nan)
    std = np.full(n_bins, np.nan)
    good = ~np.isnan(x) & ~np.isnan(rate)
    x, rate, trial_ids = x[good], rate[good], trial_ids[good]
    if x.size == 0:
        return mean, std
    bin_idx = np.clip(np.digitize(x, bin_edges) - 1, 0, n_bins - 1)
    inside = (x >= bin_edges[0]) & (x < bin_edges[-1])
    bin_idx = np.where(inside, bin_idx, -1)
    for b in range(n_bins):
        sel = bin_idx == b
        if not sel.any():
            continue
        trial_rates = [
            float(rate[sel & (trial_ids == tid)].mean())
            for tid in np.unique(trial_ids[sel])
            if (sel & (trial_ids == tid)).any()
        ]
        if len(trial_rates) >= MIN_TRIALS_PER_CELL:
            arr = np.asarray(trial_rates, dtype=np.float64)
            mean[b] = float(arr.mean())
            std[b] = float(arr.std(ddof=0)) if arr.size >= 2 else 0.0
    return mean, std


def _trial_mean_std_by_level(
    sub_df: pd.DataFrame,
    *,
    value_col: str,
    levels: np.ndarray,
    tol: float,
    bin_width: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-trial mean rate at each discrete level (SF / OR)."""
    x = sub_df[value_col].to_numpy(dtype=np.float64)
    rate = sub_df["spike_count"].to_numpy(dtype=np.float64) / bin_width
    trial_ids = sub_df["trial_id"].to_numpy(dtype=np.int64)
    mean = np.full(levels.size, np.nan)
    std = np.full(levels.size, np.nan)
    good = ~np.isnan(x) & ~np.isnan(rate)
    x, rate, trial_ids = x[good], rate[good], trial_ids[good]
    if x.size == 0:
        return mean, std
    for i, level in enumerate(levels):
        sel = np.abs(x - level) < tol
        if not sel.any():
            continue
        trial_rates = [
            float(rate[sel & (trial_ids == tid)].mean())
            for tid in np.unique(trial_ids[sel])
            if (sel & (trial_ids == tid)).any()
        ]
        if len(trial_rates) >= MIN_TRIALS_PER_CELL:
            arr = np.asarray(trial_rates, dtype=np.float64)
            mean[i] = float(arr.mean())
            std[i] = float(arr.std(ddof=0)) if arr.size >= 2 else 0.0
    return mean, std


# --------------------------------------------------------------------------- #
# Small numeric helpers
# --------------------------------------------------------------------------- #


def _mode(values: np.ndarray) -> float:
    """Return the most frequent value; ties broken by first appearance."""
    uniq, counts = np.unique(values, return_counts=True)
    return float(uniq[int(np.argmax(counts))])


def _format_radians(value: float) -> str:
    """Human-readable π-fraction for orientation tick labels."""
    pi = np.pi
    candidates = [
        (-pi / 2, "-π/2"), (-pi / 4, "-π/4"), (0.0, "0"),
        (pi / 4, "π/4"), (pi / 2, "π/2"), (3 * pi / 4, "3π/4"), (pi, "π"),
    ]
    for v, label in candidates:
        if abs(value - v) < 1e-3:
            return label
    return f"{value:.2f}"
