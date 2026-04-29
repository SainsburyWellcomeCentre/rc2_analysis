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

from rc2_glm.basis import history_basis, onset_kernel_basis, raised_cosine_basis
from rc2_glm.config import GLMConfig, MIN_TRIALS_FOR_BAND
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.precomputed_bins import PrecomputedBinEdges

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
    """Overview of the continuous basis families used by the GLM.

    Always shows Speed and TF (the two stimulus bases). The third panel
    shows whichever temporal basis the active config selects:
      - History (default since 2026-04-29 — `include_history=True`)
      - Onset kernel (when `include_onset_kernel=True`, opt-in)
      - Both, side-by-side, when both are active for an ablation run.
    """
    show_history = config.include_history
    show_onset = config.include_onset_kernel
    n_panels = 2 + (1 if show_history else 0) + (1 if show_onset else 0)
    if n_panels == 2:
        # Edge case: someone disabled both. Still render an explanatory panel.
        n_panels = 3
    fig, axes = plt.subplots(1, n_panels, figsize=(4.5 * n_panels, 4),
                             constrained_layout=True)

    speed_grid = np.linspace(config.speed_range[0], config.speed_range[1], 400)
    B_speed = raised_cosine_basis(speed_grid, config.n_speed_bases, *config.speed_range)
    for i in range(B_speed.shape[1]):
        axes[0].plot(speed_grid, B_speed[:, i], label=f"basis {i + 1}")
    axes[0].set_title(f"Speed raised cosine ({config.n_speed_bases} bases)")
    axes[0].set_xlabel("speed (cm/s)")
    axes[0].set_ylabel("activation")

    tf_grid = np.linspace(config.tf_range[0], config.tf_range[1], 400)
    B_tf = raised_cosine_basis(tf_grid, config.n_tf_bases, *config.tf_range)
    for i in range(B_tf.shape[1]):
        axes[1].plot(tf_grid, B_tf[:, i], label=f"basis {i + 1}")
    axes[1].set_title(f"TF raised cosine ({config.n_tf_bases} bases)")
    axes[1].set_xlabel("TF (Hz)")
    axes[1].set_ylabel("activation")

    panel_idx = 2
    if show_history:
        n_lag_bins = max(1, int(round(config.history_window_s / config.time_bin_width)))
        h_basis = history_basis(
            config.n_history_bases, config.history_window_s, config.time_bin_width,
        )
        # Lag in milliseconds (lag bin 1..N → time post-spike)
        lag_ms = np.arange(1, n_lag_bins + 1) * config.time_bin_width * 1000.0
        for i in range(h_basis.shape[1]):
            axes[panel_idx].plot(lag_ms, h_basis[:, i], "o-",
                                 markersize=3, label=f"basis {i + 1}")
        axes[panel_idx].set_title(
            f"History basis ({config.n_history_bases} log-spaced over "
            f"{config.history_window_s * 1000:.0f} ms — {n_lag_bins} lag bins resolved)"
        )
        axes[panel_idx].set_xlabel("lag (ms post-spike)")
        axes[panel_idx].set_ylabel("activation")
        axes[panel_idx].grid(True, alpha=0.3)
        panel_idx += 1
    if show_onset:
        onset_grid = np.linspace(0.0, config.onset_range[1], 400)
        B_onset = onset_kernel_basis(onset_grid, config.n_onset_bases, config.onset_range[1])
        for i in range(B_onset.shape[1]):
            axes[panel_idx].plot(onset_grid, B_onset[:, i], label=f"basis {i + 1}")
        axes[panel_idx].set_title(f"Onset kernel ({config.n_onset_bases} bases)")
        axes[panel_idx].set_xlabel("t since onset (s)")
        axes[panel_idx].set_ylabel("activation")
        panel_idx += 1
    if not show_history and not show_onset:
        axes[2].text(0.5, 0.5,
                     "no temporal basis active\n(include_history=False,\ninclude_onset_kernel=False)",
                     ha="center", va="center", transform=axes[2].transAxes)
        axes[2].set_axis_off()

    fig.suptitle("GLM basis functions")
    return fig


# --------------------------------------------------------------------------- #
# Per-cluster kernel-shape plot — basis × β contribution to log λ.
# Sits between the per-cluster β swarm (raw coefficients) and the
# per-cluster tuning curves (predicted firing rate marginalised over
# covariates). Shows the SHAPE of the GLM's learned kernel for each
# variable. Especially informative for History at the various bin widths.
# --------------------------------------------------------------------------- #


def _kernel_for_var(
    coef_rows: pd.DataFrame,
    var: str,
    config: GLMConfig,
) -> tuple[np.ndarray, np.ndarray, str] | None:
    """Return (x_grid, kernel_at_grid, x_label) for one variable, or None
    if no coefficients for that variable are present in ``coef_rows``.

    ``coef_rows`` is the slice of ``glm_coefficients.csv`` for one
    (cluster, model) pair (one row per coefficient).
    """
    if var == "Speed":
        sub = coef_rows[coef_rows["coefficient"].str.match(r"^Speed_\d+$")]
        if sub.empty:
            return None
        sub = sub.assign(
            lag=sub["coefficient"].str.extract(r"(\d+)$").astype(int).iloc[:, 0]
        ).sort_values("lag")
        betas = sub["estimate"].to_numpy()
        x_grid = np.linspace(*config.speed_range, 200)
        B = raised_cosine_basis(x_grid, config.n_speed_bases, *config.speed_range)
        return x_grid, B @ betas, "speed (cm/s)"

    if var == "TF":
        sub = coef_rows[coef_rows["coefficient"].str.match(r"^TF_\d+$")]
        if sub.empty:
            return None
        sub = sub.assign(
            lag=sub["coefficient"].str.extract(r"(\d+)$").astype(int).iloc[:, 0]
        ).sort_values("lag")
        betas = sub["estimate"].to_numpy()
        x_grid = np.linspace(*config.tf_range, 200)
        B = raised_cosine_basis(x_grid, config.n_tf_bases, *config.tf_range)
        return x_grid, B @ betas, "TF (Hz)"

    if var == "History":
        sub = coef_rows[coef_rows["coefficient"].str.match(r"^History_\d+$")]
        if sub.empty:
            return None
        sub = sub.assign(
            lag=sub["coefficient"].str.extract(r"(\d+)$").astype(int).iloc[:, 0]
        ).sort_values("lag")
        betas = sub["estimate"].to_numpy()
        # Reconstruct the lag-bin grid the basis was evaluated on
        n_lag_bins = max(int(round(config.history_window_s / config.time_bin_width)), 1)
        n_bases = sub.shape[0]
        h_basis = history_basis(
            n_bases, config.history_window_s, config.time_bin_width,
        )
        # Lag in milliseconds (lag bin 1..N → time post-spike)
        x_grid = np.arange(1, n_lag_bins + 1) * config.time_bin_width * 1000.0
        return x_grid, h_basis @ betas, "lag (ms post-spike)"

    if var == "Onset":
        sub = coef_rows[coef_rows["coefficient"].str.match(r"^Onset_\d+$")]
        if sub.empty:
            return None
        sub = sub.assign(
            lag=sub["coefficient"].str.extract(r"(\d+)$").astype(int).iloc[:, 0]
        ).sort_values("lag")
        betas = sub["estimate"].to_numpy()
        x_grid = np.linspace(0.0, config.onset_range[1], 200)
        B = onset_kernel_basis(x_grid, config.n_onset_bases, config.onset_range[1])
        return x_grid, B @ betas, "t since onset (s)"

    if var == "SF":
        sub = coef_rows[coef_rows["coefficient"].str.startswith("SF_")]
        if sub.empty:
            return None
        # Categorical levels — coefficient names look like "SF_0.0060".
        levels = sub["coefficient"].str.replace("SF_", "").astype(float).to_numpy()
        return levels, sub["estimate"].to_numpy(), "SF (cpd)"

    if var == "OR":
        sub = coef_rows[coef_rows["coefficient"].str.startswith("OR_")]
        if sub.empty:
            return None
        levels = sub["coefficient"].str.replace("OR_", "").astype(float).to_numpy()
        # Render in degrees for readability.
        return np.degrees(levels), sub["estimate"].to_numpy(), "orientation (deg)"

    return None


def _kernel_for_interaction(
    coef_rows: pd.DataFrame, config: GLMConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Speed × TF interaction kernel as a 2D surface.

    Returns ``(speed_grid, tf_grid, kernel)`` with shape
    ``(len(tf_grid), len(speed_grid))``, or None if no Spd_x_TF
    coefficients are present.
    """
    sub = coef_rows[coef_rows["coefficient"].str.contains("_x_TF")]
    sub = sub[sub["coefficient"].str.startswith("Spd")]
    if sub.empty:
        return None
    coef_map: dict[tuple[int, int], float] = {}
    for _, row in sub.iterrows():
        # "Spd5_x_TF3" → (5, 3)
        m = re.match(r"^Spd(\d+)_x_TF(\d+)$", row["coefficient"])
        if not m:
            continue
        coef_map[(int(m.group(1)), int(m.group(2)))] = float(row["estimate"])
    if not coef_map:
        return None
    n_speed = config.n_speed_bases
    n_tf = config.n_tf_bases
    speed_grid = np.linspace(*config.speed_range, 60)
    tf_grid = np.linspace(*config.tf_range, 60)
    Bs = raised_cosine_basis(speed_grid, n_speed, *config.speed_range)
    Bt = raised_cosine_basis(tf_grid, n_tf, *config.tf_range)
    kernel = np.zeros((tf_grid.size, speed_grid.size), dtype=np.float64)
    for (i, j), beta in coef_map.items():
        if 1 <= i <= n_speed and 1 <= j <= n_tf:
            kernel += beta * np.outer(Bt[:, j - 1], Bs[:, i - 1])
    return speed_grid, tf_grid, kernel


def plot_cluster_kernels(
    probe_id: str,
    cluster_id: int,
    coef_df_cluster: pd.DataFrame,
    config: GLMConfig,
) -> Figure:
    """Per-cluster kernel shapes — what the GLM actually learned.

    For each of the four models (Null / Selected / Additive /
    FullInteraction) and each of six variables (Speed / TF / SF / OR /
    History / Spd × TF), plots the kernel shape ``Σ_k β_k · basis_k(x)``
    that contributes to the log firing rate. Categorical variables (SF /
    OR) render as bar charts; the Speed × TF interaction renders as a
    2D heatmap. Empty panels mean the variable wasn't part of that
    model's fit.

    ``coef_df_cluster`` is the slice of ``glm_coefficients.csv`` for
    this cluster, ``glm_type == "time"`` (one row per (model,
    coefficient) pair).
    """
    models = ("Null", "Selected", "Additive", "FullInteraction")
    var_cols = ("Speed", "TF", "SF", "OR", "History", "Spd × TF")

    fig, axes = plt.subplots(
        len(models), len(var_cols),
        figsize=(20, 11), constrained_layout=True,
    )

    # Header banner: which model is selected for this cluster
    sel_row = coef_df_cluster[coef_df_cluster["model"] == "Selected"]
    selected_vars = "—"
    if not sel_row.empty:
        # Recover the selected_vars string from the coefficient list
        coef_names = sel_row["coefficient"].tolist()
        groups = sorted({_beta_group_tag(c) for c in coef_names if c != "Intercept"})
        selected_vars = "+".join(g for g in groups if g != "Other")

    fig.suptitle(
        f"Per-cluster kernels — {probe_id} · cluster {cluster_id} · "
        f"Selected: {selected_vars}",
        fontsize=12,
    )

    for row_idx, model in enumerate(models):
        sub = coef_df_cluster[coef_df_cluster["model"] == model]
        for col_idx, var in enumerate(var_cols):
            ax = axes[row_idx, col_idx]
            if col_idx == 0:
                ax.set_ylabel(f"{model}\n\nlog λ contribution",
                              fontsize=9, rotation=90)
            if row_idx == 0:
                ax.set_title(var, fontsize=11, fontweight="bold")
            if sub.empty:
                ax.text(0.5, 0.5, "no fit",
                        ha="center", va="center",
                        transform=ax.transAxes,
                        fontsize=10, color=(0.5, 0.5, 0.5))
                ax.set_xticks([])
                ax.set_yticks([])
                continue

            if var == "Spd × TF":
                inter = _kernel_for_interaction(sub, config)
                if inter is None:
                    ax.text(0.5, 0.5, "not selected",
                            ha="center", va="center",
                            transform=ax.transAxes,
                            fontsize=9, color=(0.5, 0.5, 0.5))
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue
                speed_grid, tf_grid, kernel = inter
                vmax = np.abs(kernel).max() or 1.0
                im = ax.pcolormesh(
                    speed_grid, tf_grid, kernel,
                    cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto",
                )
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                ax.set_xlabel("speed (cm/s)", fontsize=8)
                if col_idx == 0:
                    ax.set_ylabel("TF (Hz)", fontsize=8)
                continue

            kdata = _kernel_for_var(sub, var, config)
            if kdata is None:
                ax.text(0.5, 0.5, "not selected",
                        ha="center", va="center",
                        transform=ax.transAxes,
                        fontsize=9, color=(0.5, 0.5, 0.5))
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            x, kernel, xlabel = kdata
            colour = _BETA_GROUP_COLORS.get(var, (0.2, 0.2, 0.2))
            if var in ("SF", "OR"):
                ax.bar(x, kernel, color=colour, edgecolor="black",
                       linewidth=0.5, width=(x.max() - x.min()) / max(len(x), 1) * 0.4)
                ax.axhline(0, color="grey", linewidth=0.6, alpha=0.7)
            else:
                ax.plot(x, kernel, color=colour, linewidth=2)
                ax.fill_between(x, 0, kernel, color=colour, alpha=0.15)
                ax.axhline(0, color="grey", linewidth=0.6, alpha=0.7)
            ax.set_xlabel(xlabel, fontsize=8)
            ax.tick_params(axis="both", labelsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)

    return fig


# --------------------------------------------------------------------------- #
# Fig 2 — forward-selection summary across clusters
# --------------------------------------------------------------------------- #


# MATLAB colour scheme (glm_single_cluster_analysis.m lines 1862-1996).
_COLOR_SPEED = (0.17, 0.63, 0.17)
_COLOR_TF = (1.00, 0.50, 0.05)
_COLOR_SF = (0.95, 0.85, 0.10)
_COLOR_OR = (0.84, 0.15, 0.16)
_COLOR_NULL = (0.75, 0.75, 0.75)
_COLOR_4PLUS = (0.35, 0.35, 0.35)
_COLOR_INT_BLUE = (0.20, 0.40, 0.80)

_MAIN_EFFECT_NAMES = ("Speed", "TF", "SF", "OR")
_MAIN_EFFECT_COLORS: dict[str, tuple[float, float, float]] = {
    "Speed": _COLOR_SPEED, "TF": _COLOR_TF,
    "SF": _COLOR_SF, "OR": _COLOR_OR,
}
_MAIN_EFFECT_LABELS = {"Speed": "Spd", "TF": "TF", "SF": "SF", "OR": "OR"}
_INTERACTION_LABELS = {
    "Speed_x_TF": "SxT", "Speed_x_SF": "SxSF", "Speed_x_OR": "SxO",
    "TF_x_SF": "TxSF", "TF_x_OR": "TxO", "SF_x_OR": "SFxO",
}
# Keys are sorted "+"-joined main-effect lists (matches MATLAB main_key).
_COMBO_COLORS: dict[str, tuple[float, float, float]] = {
    "Speed": _COLOR_SPEED, "TF": _COLOR_TF,
    "SF": _COLOR_SF, "OR": _COLOR_OR,
    "Speed+TF": (0.12, 0.47, 0.71),
    "SF+Speed": (0.26, 0.58, 0.85),
    "OR+Speed": (0.05, 0.30, 0.55),
    "SF+TF": (0.40, 0.60, 0.80),
    "OR+TF": (0.15, 0.40, 0.65),
    "OR+SF": (0.30, 0.50, 0.75),
    "SF+Speed+TF": (0.50, 0.70, 0.90),
    "OR+Speed+TF": (0.35, 0.55, 0.80),
    "OR+SF+Speed": (0.45, 0.65, 0.85),
    "OR+SF+TF": (0.55, 0.75, 0.92),
    "OR+SF+Speed+TF": _COLOR_4PLUS,
}


def _darken(c: tuple[float, float, float], factor: float = 0.65) -> tuple[float, float, float]:
    return tuple(max(x * factor, 0.0) for x in c)  # type: ignore[return-value]


def _parse_selected_vars(s: str) -> list[str]:
    """Split a selected-vars token string on either separator.

    Python writes ``"Speed+TF"``; MATLAB writes ``"Speed, TF"``. Accepting
    both means this plot renders truthfully from either side's CSV —
    without this, a MATLAB-side aggregate collapses every multi-variable
    model into a single-variable bucket (the full token ``"Speed, TF"``
    parses as one unit).
    """
    if not isinstance(s, str) or not s or s == "Null":
        return []
    norm = s.replace(", ", "+").replace(",", "+")
    return [t.strip() for t in norm.split("+") if t.strip()]


def plot_forward_selection_summary(
    comparison_df: pd.DataFrame, history_df: pd.DataFrame | None = None,
) -> Figure:
    """MATLAB Fig 2 port — forward-selection classification summary.

    Three panels mirroring ``glm_single_cluster_analysis.m`` lines
    1828-2209:

    1. Stacked bar by model complexity (Null / 1 var / 2 / 3 / 4+).
       Each stack segment is a distinct model type; colour encodes
       composition (pure colour for single-variable models, blue
       variations for multi-variable combos, darker if the model
       contains interactions). Segments with ≥2 main effects are
       overlaid with horizontal bands cycling through the constituent
       main-effect colours — this replaces MATLAB's diagonal stripe
       pattern with the same horizontal-band scheme MATLAB actually
       draws (`patch` bands, lines 2100-2117).
    2. Variable-inclusion rates (horizontal bar: Speed/TF/SF/OR/Any
       Interact.), same colour scheme as Panel 1.
    3. Pairwise interaction breakdown (horizontal bar for the six
       interactions).

    ``history_df`` is accepted for call-site compatibility and unused —
    Panel 1's complexity view already covers what the old
    null-vs-selected-scatter / var-pick-counts plots showed.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    if (
        comparison_df.empty
        or "time_selected_vars" not in comparison_df.columns
    ):
        fig.suptitle("Forward Selection Summary — no data")
        return fig

    n_total = len(comparison_df)

    # --- Build catalog: sorted-key -> metadata for each distinct model
    catalog: dict[str, dict] = {}
    for raw in comparison_df["time_selected_vars"].astype(str):
        sv = _parse_selected_vars(raw)
        key = "+".join(sorted(sv)) if sv else "Null"
        if key not in catalog:
            main_effects = sorted(v for v in sv if v in _MAIN_EFFECT_NAMES)
            interactions = sorted(v for v in sv if v not in _MAIN_EFFECT_NAMES)
            catalog[key] = {
                "count": 0,
                "complexity": len(sv),
                "main_effects": main_effects,
                "interactions": interactions,
            }
        catalog[key]["count"] += 1

    for key, info in catalog.items():
        if key == "Null":
            info["color"] = _COLOR_NULL
            info["label"] = "Null"
            continue
        main_key = "+".join(info["main_effects"])
        base = _COMBO_COLORS.get(main_key, (0.5, 0.5, 0.5))
        info["color"] = _darken(base) if info["interactions"] else base
        parts = [_MAIN_EFFECT_LABELS[m] for m in info["main_effects"]]
        parts += [_INTERACTION_LABELS.get(i, i) for i in info["interactions"]]
        info["label"] = "+".join(parts)

    # Complexity → bar position: 0 vars → 1, 1 → 2, …, ≥4 → 5.
    by_level: dict[int, list[dict]] = {i: [] for i in range(1, 6)}
    for info in catalog.values():
        level = min(info["complexity"] + 1, 5)
        by_level[level].append(info)
    for level in by_level:
        by_level[level].sort(key=lambda d: -d["count"])

    # --- Panel 1: Model complexity stacked bar ---
    ax = axes[0]
    bar_width = 0.8
    legend_entries: list[dict] = []
    totals = [sum(d["count"] for d in by_level[lv]) for lv in range(1, 6)]
    max_total = max(totals) if totals else 0

    for level, segments in by_level.items():
        y_bottom = 0.0
        for info in segments:
            height = info["count"]
            if height <= 0:
                continue
            ax.bar(
                level, height, width=bar_width, bottom=y_bottom,
                color=info["color"], edgecolor=(0.3, 0.3, 0.3), linewidth=0.5,
            )
            me = info["main_effects"]
            if len(me) >= 2:
                n_eff = len(me)
                band_h = height / (n_eff * 3)
                band_colors = [_MAIN_EFFECT_COLORS[m] for m in me]
                y_pos = y_bottom
                band_idx = 0
                while y_pos < y_bottom + height - 1e-9:
                    band_top = min(y_pos + band_h, y_bottom + height)
                    c = band_colors[band_idx % n_eff]
                    ax.add_patch(Rectangle(
                        (level - bar_width / 2, y_pos),
                        bar_width, band_top - y_pos,
                        facecolor=c, alpha=0.85, edgecolor="none",
                    ))
                    y_pos = band_top
                    band_idx += 1
            y_bottom += height
            legend_entries.append(info)

    for level, total in enumerate(totals, start=1):
        if total <= 0:
            continue
        ax.text(
            level, total + max_total * 0.03,
            f"{total}\n({100 * total / n_total:.0f}%)",
            ha="center", va="bottom", fontsize=7,
        )
    ax.set_xticks(range(1, 6))
    ax.set_xticklabels(["Null", "1 var", "2 vars", "3 vars", "4+ vars"], rotation=30)
    ax.set_ylabel("# Neurons")
    ax.set_title("Model Complexity (Forward Selection)", fontsize=10)
    ax.set_ylim(0, max(max_total * 1.25, 1))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    legend_handles: list[Line2D] = []
    legend_labels: list[str] = []
    for info in legend_entries:
        me = info["main_effects"]
        if len(me) <= 1:
            handle = Line2D(
                [0], [0], marker="s", linestyle="none",
                markerfacecolor=info["color"],
                markeredgecolor=(0.3, 0.3, 0.3),
                markersize=8,
            )
        else:
            handle = Line2D(
                [0], [0], marker="s", linestyle="none",
                markerfacecolor=_MAIN_EFFECT_COLORS[me[0]],
                markeredgecolor=_MAIN_EFFECT_COLORS[me[1]],
                markeredgewidth=2, markersize=8,
            )
        legend_handles.append(handle)
        legend_labels.append(f"{info['label']} ({info['count']})")
    if legend_handles:
        ax.legend(
            legend_handles, legend_labels,
            loc="upper right", fontsize=5, frameon=True,
        )

    # --- Panel 2: Variable-inclusion rates ---
    ax = axes[1]
    # History selection isn't a column in comparison_df — derive it from
    # the selected_vars string (which lists "History" when the spike-history
    # term made it through forward selection).
    history_count = int(
        comparison_df["time_selected_vars"].astype(str)
        .apply(lambda s: "History" in s.split("+")).sum()
    )
    me_counts = [
        int(comparison_df["time_is_speed_tuned"].sum()),
        int(comparison_df["time_is_tf_tuned"].sum()),
        int(comparison_df["time_is_sf_tuned"].sum()),
        int(comparison_df["time_is_or_tuned"].sum()),
        history_count,
        int(comparison_df["time_has_interaction"].sum()),
    ]
    me_labels = ["Speed", "TF", "SF", "OR", "History", "Any Interact."]
    me_colors = [
        _COLOR_SPEED, _COLOR_TF, _COLOR_SF, _COLOR_OR,
        _BETA_GROUP_COLORS["History"], _COLOR_INT_BLUE,
    ]
    y_pos = np.arange(1, len(me_counts) + 1)
    ax.barh(y_pos, me_counts, color=me_colors, edgecolor="none")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(me_labels)
    ax.set_xlabel("# Neurons")
    ax.set_title(f"Variable Inclusion (n={n_total})", fontsize=10)
    max_me = max(me_counts) if any(me_counts) else 1
    for yi, c in zip(y_pos, me_counts):
        ax.text(
            c + max_me * 0.02, yi,
            f"{c} ({100 * c / n_total:.0f}%)",
            va="center", fontsize=8,
        )
    ax.set_xlim(0, max_me * 1.35)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # --- Panel 3: Interaction breakdown ---
    ax = axes[2]
    int_counts = [
        int(comparison_df["time_has_speed_x_tf"].sum()),
        int(comparison_df["time_has_speed_x_sf"].sum()),
        int(comparison_df["time_has_speed_x_or"].sum()),
        int(comparison_df["time_has_tf_x_sf"].sum()),
        int(comparison_df["time_has_tf_x_or"].sum()),
        int(comparison_df["time_has_sf_x_or"].sum()),
    ]
    int_labels = ["Spd x TF", "Spd x SF", "Spd x OR", "TF x SF", "TF x OR", "SF x OR"]
    int_colors = [
        (0.9, 0.2, 0.2), (0.8, 0.5, 0.2), (0.6, 0.2, 0.6),
        (0.5, 0.8, 0.2), (0.2, 0.5, 0.8), (0.9, 0.6, 0.6),
    ]
    y_pos = np.arange(1, len(int_counts) + 1)
    ax.barh(y_pos, int_counts, color=int_colors, edgecolor="none")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(int_labels)
    ax.set_xlabel("# Neurons")
    ax.set_title("Interaction Breakdown", fontsize=10)
    max_int = max(max(int_counts), 1)
    for yi, c in zip(y_pos, int_counts):
        ax.text(
            c + max_int * 0.02, yi,
            f"{c} ({100 * c / n_total:.0f}%)",
            va="center", fontsize=8,
        )
    ax.set_xlim(0, max_int * 1.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    thr = GLMConfig().delta_bps_threshold
    fig.suptitle(
        f"Forward Selection Results (threshold: Δ bps > {thr:.3f})",
        fontsize=12, fontweight="bold",
    )
    return fig


# --------------------------------------------------------------------------- #
# Fig 2b — CV bps deltas (replaces the old null-vs-selected scatter)
# --------------------------------------------------------------------------- #


def plot_speed_profile_cv_comparison(comparison_df: pd.DataFrame) -> Figure:
    """Standard vs speed-profile CV, per-cluster paired comparison.

    Reproduces MATLAB ``speed_profile_cross_validation.png``
    (``glm_single_cluster_analysis.m:2467-2627``). Three panels:

    1. Δ-bps = CV(Selected) − CV(Null), standard vs profile CV, per cluster.
    2. Speed unique contribution = CV(Selected) − CV(Selected-without-Speed),
       standard vs profile CV, per cluster, speed-tuned only.
    3. Stacked bar: median per-cluster ratio retained (profile / standard)
       vs profile-specific (1 − retained), for total Δ-bps and for speed
       contribution.

    The figure is only meaningful when the pipeline ran with
    ``--profile-cv-diagnostic`` — otherwise the profile_cv_bps_* columns
    are NaN and the function returns a short explanatory stub.
    """
    from scipy.stats import wilcoxon

    # Discover the two sets of cv-bps columns regardless of the ``time_``
    # prefix convention used by the caller's CSV load.
    def _col(*candidates: str) -> str | None:
        for c in candidates:
            if c in comparison_df.columns:
                return c
        return None

    col_null = _col("time_Null_cv_bps", "Null_cv_bps")
    col_sel = _col("time_Selected_cv_bps", "Selected_cv_bps")
    col_pnull = _col("time_profile_cv_bps_Null", "profile_cv_bps_Null")
    col_psel = _col("time_profile_cv_bps_Selected", "profile_cv_bps_Selected")
    col_pnospd = _col("time_profile_cv_bps_NoSpeed", "profile_cv_bps_NoSpeed")

    fig, axes = plt.subplots(
        1, 3, figsize=(14.0, 5.5), constrained_layout=True,
    )

    missing_required = any(c is None for c in
                           (col_null, col_sel, col_pnull, col_psel))
    if missing_required or comparison_df[col_pnull].isna().all():
        fig.suptitle(
            "Speed-profile CV comparison — no data "
            "(run with --profile-cv-diagnostic to populate)",
            fontsize=11,
        )
        for ax in axes:
            ax.set_axis_off()
        return fig

    # ---- Per-cluster paired deltas ----
    d_std = (comparison_df[col_sel] - comparison_df[col_null]).to_numpy(dtype=np.float64)
    d_spr = (comparison_df[col_psel] - comparison_df[col_pnull]).to_numpy(dtype=np.float64)
    has_pair = np.isfinite(d_std) & np.isfinite(d_spr)
    d_std_p = d_std[has_pair]
    d_spr_p = d_spr[has_pair]
    n_pair = int(has_pair.sum())
    ylim = (-0.5, 0.5)
    n_clipped = int(
        ((d_std_p < ylim[0]) | (d_std_p > ylim[1])
         | (d_spr_p < ylim[0]) | (d_spr_p > ylim[1])).sum()
    )

    col_std = (0.55, 0.55, 0.55)
    col_pro = (0.25, 0.45, 0.75)

    ax1 = axes[0]
    for i in range(n_pair):
        ax1.plot([1, 2], [d_std_p[i], d_spr_p[i]],
                 color=(0.75, 0.75, 0.75, 0.35), linewidth=0.7)
    ax1.scatter(np.ones(n_pair), d_std_p, s=28, color=col_std, alpha=0.65)
    ax1.scatter(np.full(n_pair, 2), d_spr_p, s=28, color=col_pro, alpha=0.65)
    med_std = float(np.median(d_std_p)) if n_pair else float("nan")
    med_spr = float(np.median(d_spr_p)) if n_pair else float("nan")
    ax1.plot([1, 2], [med_std, med_spr], "k-", linewidth=2.5)
    ax1.plot(1, med_std, "ko", markerfacecolor=col_std, markersize=9)
    ax1.plot(2, med_spr, "ko", markerfacecolor=col_pro, markersize=9)
    ax1.axhline(0.0, color=(0.4, 0.4, 0.4), linestyle="--", linewidth=0.8)
    ax1.set_xlim(0.5, 2.5)
    ax1.set_ylim(*ylim)
    ax1.set_xticks([1, 2])
    ax1.set_xticklabels(["Standard CV", "Speed-Profile CV"])
    ax1.set_ylabel("Δ bps = CV(Selected) − CV(Null)\nbits / spike", fontsize=9)
    ax1.set_title(
        f"Total model information gain  (n={n_pair}, {n_clipped} clipped)",
        fontsize=10,
    )
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    # Wilcoxon one-tailed: H1 = standard CV > profile CV.
    if n_pair >= 1:
        try:
            p_t = float(wilcoxon(d_std_p, d_spr_p, alternative="greater").pvalue)
        except ValueError:
            p_t = float("nan")
    else:
        p_t = float("nan")
    ax1.text(1.5, ylim[1] * 0.88,
             f"Wilcoxon: p {'< 0.001' if p_t < 0.001 else f'= {p_t:.3f}'}",
             ha="center", fontsize=9, color=(0.2, 0.2, 0.2))

    # ---- Panel 2: speed unique contribution ----
    # MATLAB computes this under BOTH fold schemes: the standard-CV
    # Selected vs Selected-without-Speed numbers come from
    # profile_cv_bps_*_standard columns (populated by
    # _profile_cv_diagnostic when --profile-cv-diagnostic is set).
    col_psel_std = _col("time_profile_cv_bps_Selected_standard",
                        "profile_cv_bps_Selected_standard")
    col_pnospd_std = _col("time_profile_cv_bps_NoSpeed_standard",
                          "profile_cv_bps_NoSpeed_standard")
    ax2 = axes[1]
    if col_pnospd is not None and col_psel_std is not None and col_pnospd_std is not None:
        s_std = (comparison_df[col_psel_std] - comparison_df[col_pnospd_std]).to_numpy(
            dtype=np.float64,
        )
        s_spr = (comparison_df[col_psel] - comparison_df[col_pnospd]).to_numpy(dtype=np.float64)
        pair_s = np.isfinite(s_std) & np.isfinite(s_spr)
        s_std_p = s_std[pair_s]
        s_spr_p = s_spr[pair_s]
        n_pair_s = int(pair_s.sum())
    else:
        n_pair_s = 0
    if n_pair_s > 0:
        col_pro_s = (0.15, 0.58, 0.15)
        for i in range(n_pair_s):
            ax2.plot([1, 2], [s_std_p[i], s_spr_p[i]],
                     color=(0.75, 0.75, 0.75, 0.35), linewidth=0.7)
        ax2.scatter(np.ones(n_pair_s), s_std_p, s=28, color=col_std, alpha=0.65)
        ax2.scatter(np.full(n_pair_s, 2), s_spr_p, s=28, color=col_pro_s, alpha=0.65)
        med_s_std = float(np.median(s_std_p))
        med_s_spr = float(np.median(s_spr_p))
        ax2.plot([1, 2], [med_s_std, med_s_spr], "k-", linewidth=2.5)
        ax2.plot(1, med_s_std, "ko", markerfacecolor=col_std, markersize=9)
        ax2.plot(2, med_s_spr, "ko", markerfacecolor=col_pro_s, markersize=9)
        ax2.axhline(0.0, color=(0.4, 0.4, 0.4), linestyle="--", linewidth=0.8)
        ax2.set_xlim(0.5, 2.5)
        ax2.set_ylim(*ylim)
        ax2.set_xticks([1, 2])
        ax2.set_xticklabels(["Standard CV", "Speed-Profile CV"])
        ax2.set_ylabel("Δ bps = CV(Sel) − CV(Sel − Speed)\nbits / spike", fontsize=9)
        ax2.set_title(
            f"Speed unique contribution  (n={n_pair_s} speed-tuned)",
            fontsize=10,
        )
        for sp in ("top", "right"):
            ax2.spines[sp].set_visible(False)
        try:
            p_s = float(wilcoxon(s_std_p, s_spr_p, alternative="greater").pvalue)
        except ValueError:
            p_s = float("nan")
        ax2.text(1.5, ylim[1] * 0.88,
                 f"Wilcoxon: p {'< 0.001' if p_s < 0.001 else f'= {p_s:.3f}'}",
                 ha="center", fontsize=9, color=(0.2, 0.2, 0.2))
    else:
        ax2.text(0.5, 0.5, "no speed-tuned clusters\nwith paired profile/standard Δ",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=9, color="#888")
        ax2.set_axis_off()

    # ---- Panel 3: retained vs profile-specific stacked bar ----
    ax3 = axes[2]
    pos_t = d_std_p > 0
    ret_t = (
        min(100 * float(np.median(d_spr_p[pos_t] / d_std_p[pos_t])), 100)
        if pos_t.any() else float("nan")
    )
    if n_pair_s > 0:
        pos_s = s_std_p > 0
        ret_s = (
            min(100 * float(np.median(s_spr_p[pos_s] / s_std_p[pos_s])), 100)
            if pos_s.any() else float("nan")
        )
    else:
        ret_s = float("nan")
    bar_x = np.array([1.0, 2.0])
    retained = np.array([ret_t, ret_s])
    not_ret = 100 - retained
    retained_safe = np.where(np.isfinite(retained), retained, 0)
    not_ret_safe = np.where(np.isfinite(not_ret), not_ret, 0)
    ax3.bar(bar_x, retained_safe, width=0.55, color=[col_pro, (0.15, 0.58, 0.15)],
            edgecolor="none", label="Generalises")
    ax3.bar(bar_x, not_ret_safe, width=0.55, bottom=retained_safe,
            color=(0.88, 0.88, 0.88), edgecolor="none", label="Profile-specific")
    for bi in range(2):
        if np.isfinite(retained[bi]):
            ax3.text(bar_x[bi], retained[bi] / 2, f"{retained[bi]:.0f}%",
                     ha="center", va="center", fontsize=11, fontweight="bold",
                     color="white")
            ax3.text(bar_x[bi], retained[bi] + not_ret[bi] / 2,
                     f"{not_ret[bi]:.0f}%",
                     ha="center", va="center", fontsize=11, fontweight="bold",
                     color=(0.4, 0.4, 0.4))
    ax3.set_xticks(bar_x)
    ax3.set_xticklabels(["Total model", "Speed"])
    ax3.set_ylim(0, 100)
    ax3.set_ylabel("Median per-cluster ratio: profile / standard (% retained)",
                   fontsize=9)
    ax3.set_title("Generalisation across speed profiles", fontsize=10)
    ax3.legend(loc="upper right", fontsize=8, frameon=False)
    for sp in ("top", "right"):
        ax3.spines[sp].set_visible(False)

    fig.suptitle(
        "Speed-Profile Cross-Validation: generalisation across speed profiles",
        fontsize=11, fontweight="bold",
    )
    return fig


def plot_cv_bps_deltas(comparison_df: pd.DataFrame) -> Figure:
    """Per-cluster Δ CV-bps distributions across model comparisons.

    Raw CV bps is dominated by a cluster-level baseline (its null log-
    likelihood), so plotting absolute CV bps makes it nearly impossible
    to see model-comparison effects — every point sits at the cluster's
    own scale. The informative quantity is the Δ between two nested
    models on the same cluster, which isolates the information gain
    attributable to the added / removed structure.

    Four panels (strip plots, one dot per cluster, jittered):

    - ΔSelected − Null       : info gain from forward selection vs null.
    - ΔAdditive − Null       : info gain of the saturated main-effects
                               model vs null (upper envelope of what
                               main effects alone can contribute).
    - ΔSelected − Additive   : does forward selection beat the saturated
                               additive model? Should typically be ≥ 0
                               for selection to earn its keep, but can
                               be slightly negative when selection drops
                               marginal predictors that the CV actually
                               wanted.
    - ΔFullInteraction − Additive : do all pairwise interactions help?
                               Mostly ≤ 0 on this dataset (p≫n region),
                               which is the motivation for keeping the
                               per-interaction forward-selection loop.

    Each panel shows per-cluster dots, a horizontal median marker, the
    y=0 reference line, and the selection threshold (±config default)
    as a shaded band. An annotation in the corner reports
    ``N>0 | N<0 | median``.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    deltas = [
        ("time_delta_selected_vs_null",
         "Δ Selected − Null",
         "info gain: forward selection vs null",
         (0.20, 0.60, 0.20)),
        ("time_delta_additive_vs_null",
         "Δ Additive − Null",
         "info gain: saturated main effects vs null",
         (0.25, 0.45, 0.75)),
        ("time_delta_selected_vs_additive",
         "Δ Selected − Additive",
         "does selection beat saturated additive?",
         (0.70, 0.40, 0.10)),
        ("time_delta_interaction",
         "Δ FullInteraction − Additive",
         "do pairwise interactions help?",
         (0.60, 0.20, 0.50)),
    ]

    threshold = GLMConfig().delta_bps_threshold
    rng = np.random.default_rng(0)  # deterministic jitter

    for ax, (col, title, subtitle, color) in zip(axes.flat, deltas):
        if col not in comparison_df.columns or comparison_df.empty:
            ax.set_title(f"{title}\n(column missing)", fontsize=10)
            continue
        vals = comparison_df[col].to_numpy(dtype=float)
        vals = vals[~np.isnan(vals)]
        if vals.size == 0:
            ax.set_title(f"{title}\n(no data)", fontsize=10)
            continue
        n_pos = int((vals > 0).sum())
        n_neg = int((vals < 0).sum())
        median = float(np.median(vals))

        jitter = rng.uniform(-0.2, 0.2, size=vals.size)
        ax.scatter(
            jitter, vals, s=25, color=color, alpha=0.65, edgecolors="none",
        )
        ax.axhline(0.0, color=(0.4, 0.4, 0.4), linestyle="--", linewidth=0.8)
        ax.axhspan(-threshold, threshold, color=(0.85, 0.85, 0.85), alpha=0.4,
                   label=f"±{threshold:g} threshold")
        ax.hlines(
            median, -0.35, 0.35,
            colors="black", linewidth=2.5, zorder=3,
        )
        ax.scatter([0], [median], s=60, color="black", zorder=4,
                   label=f"median = {median:+.4f}")

        ax.set_xlim(-0.6, 0.6)
        ax.set_xticks([])
        ax.set_ylabel("Δ bits / spike")
        ax.set_title(title, fontsize=11)
        ax.text(
            0.5, -0.18, subtitle,
            transform=ax.transAxes,
            ha="center", va="top", fontsize=9, color=(0.3, 0.3, 0.3),
        )
        ax.text(
            0.02, 0.98,
            f"N>0: {n_pos}\nN<0: {n_neg}\nmedian: {median:+.4f}",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor=(0.7, 0.7, 0.7),
                      alpha=0.9),
        )
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    n_total = int(len(comparison_df))
    fig.suptitle(
        f"CV bps model-comparison Δ per cluster  (n={n_total}, "
        f"threshold ±{threshold:g})",
        fontsize=12, fontweight="bold",
    )
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
    precomputed_bins: PrecomputedBinEdges | None = None,
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

    The onset-kernel treatment depends on ``config.tuning_curve_mode``:

    - ``"trial-averaged"`` (default, Laura 2026-04-23): at each sweep
      grid point, average predicted rate across the cluster's observed
      ``time_since_onset`` values *within the condition's motion bins*.
      This marginalises out the onset basis in a way that reflects the
      trial-time distribution the neuron actually experienced, rather
      than pinning onset to an arbitrary steady-state time.
    - ``"steady-state"``: evaluate the onset kernel at ``t=1.5s``
      (MATLAB line 3397) — the pre-2026-04-23 behaviour, kept for
      back-compat and side-by-side comparisons.

    The stationary dot is computed from observed stationary bins in
    both modes (SF/OR at observed values, onset=0).

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

    # Trial-averaged mode needs the condition-specific observed onset
    # distributions. Collect them now so each model row reuses the same
    # bin sets. ``None`` when the condition has no motion bins — the
    # row-plotting code falls back to the mean/steady-state scalar.
    cond_bins: dict[str, pd.DataFrame] = {}
    if config.tuning_curve_mode == "trial-averaged":
        motion_cluster = cluster_df[cluster_df["condition"] != "stationary"]
        for cond in ("T_Vstatic", "V", "VT"):
            sub = motion_cluster[motion_cluster["condition"] == cond]
            cond_bins[cond] = sub if not sub.empty else None

    # Per-condition Speed/TF grids, matching the bin_centres used in the
    # Observed row. Using the SAME grid on both rows is essential: a
    # model line plotted on a uniform linspace beside an observed line
    # plotted at 5%-quantile centres gives a figure where x-positions
    # are not comparable, which was Laura's 2026-04-23 finding. Falls
    # back to a linspace if the MATLAB precomputed-edges cache is
    # unavailable — same behaviour as the observed row. (The fallback
    # has 10 bins to match _plot_observed_row's fallback_spd_edges.)
    fallback_spd = np.linspace(*config.speed_range, 11)
    fallback_tf = np.linspace(*config.tf_range, 11)
    spd_centres_by_cond: dict[str, np.ndarray] = {}
    tf_centres_by_cond: dict[str, np.ndarray] = {}
    for cond in ("T_Vstatic", "V", "VT"):
        spd_edges = (
            precomputed_bins.speed_edges(cond)
            if precomputed_bins is not None else None
        )
        if spd_edges is None:
            spd_edges = fallback_spd
        spd_centres_by_cond[cond] = 0.5 * (spd_edges[:-1] + spd_edges[1:])
        tf_edges = (
            precomputed_bins.tf_edges(cond)
            if precomputed_bins is not None else None
        )
        if tf_edges is None:
            tf_edges = fallback_tf
        tf_centres_by_cond[cond] = 0.5 * (tf_edges[:-1] + tf_edges[1:])

    _plot_observed_row(
        axes[0, :], cluster_df, config,
        sf_plot=sf_plot, or_plot=or_plot,
        precomputed_bins=precomputed_bins,
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
            cond_bins=cond_bins,
            spd_centres_by_cond=spd_centres_by_cond,
            tf_centres_by_cond=tf_centres_by_cond,
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
        f"    selected: {'+'.join(sel_vars) if sel_vars else 'Null'}"
        f"    [{config.tuning_curve_mode}]",
        fontsize=12,
    )
    return fig


# --------------------------------------------------------------------------- #
# Fig 8b — per-trial Speed/TF + observed-vs-predicted FR over time
# (MATLAB lines 4941-5198 of glm_single_cluster_analysis.m)
# --------------------------------------------------------------------------- #


_TRIAL_LEVEL_COND_COLORS: dict[str, tuple[float, float, float]] = {
    "T_Vstatic": (0.2, 0.6, 0.2),  # green
    "V":         (0.8, 0.2, 0.2),  # red
    "VT":        (0.2, 0.2, 0.8),  # blue
}
_TRIAL_LEVEL_COND_CAPTION: dict[str, str] = {
    "T_Vstatic": "T (speed only)",
    "V":         "V (visual only)",
    "VT":        "VT (multi-sensory)",
}


def plot_trial_level_predictions(
    probe_id: str,
    cluster_id: int,
    cluster_df: pd.DataFrame,
    model_predictions: dict[str, np.ndarray],
    config: GLMConfig,
    *,
    selected_label: str = "Selected",
    cv_bps: float | None = None,
    selected_vars: list[str] | None = None,
    seed: int = 0,
) -> Figure:
    """Per-trial Speed / TF / observed-vs-predicted FR traces.

    Ports MATLAB ``Section 8b`` (``glm_single_cluster_analysis.m`` lines
    4941-5198): selects 8 representative trials per cluster (2 T_Vstatic,
    2 V, 4 VT, drawn with a deterministic RNG so re-runs produce the
    same figures) and renders an N×3 grid:

    - Col 0: Speed(t) coloured by condition.
    - Col 1: TF(t) coloured by condition; y-axis shared across trials
      for the cluster so amplitude differences are visible.
    - Col 2: Boxcar-smoothed observed FR (grey) overlaid with the
      Selected-model predicted FR (red).

    Deviation from MATLAB: the Python pipeline currently stores the
    full-data refit predictions (``ClusterFit.model_predictions``), not
    concatenated held-out CV predictions. The trace is therefore
    in-sample — fine for "does the model capture the time-course?"
    inspection, but will mask out-of-sample shortfall compared to
    MATLAB's CV version. The pipeline issue tracking a switch to CV
    predictions is Commit 10 (`FullInteraction skip audit trail + stats
    storage review`).
    """
    motion_mask = (cluster_df["condition"] != "stationary").to_numpy()
    motion_df = cluster_df.loc[motion_mask].copy()
    preds = model_predictions.get(selected_label)

    title_parts = [f"Trial-Level Predictions | {probe_id} cluster {cluster_id}"]
    if cv_bps is not None and np.isfinite(cv_bps):
        title_parts.append(f"cv_bps={cv_bps:.3f}")
    if selected_vars:
        title_parts.append(f"vars: {'+'.join(selected_vars) or 'Null'}")
    title = "  |  ".join(title_parts)

    if motion_df.empty or preds is None or preds.size != len(motion_df):
        fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
        ax.text(
            0.5, 0.5,
            "no motion rows or predictions unavailable",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=11, color=(0.5, 0.5, 0.5),
        )
        ax.set_axis_off()
        fig.suptitle(title, fontsize=11, fontweight="bold")
        return fig

    motion_df = motion_df.reset_index(drop=True)
    motion_df["_pred_fr"] = preds
    # Full-cluster obs FR, then filter to motion (keeps boxcar smoothing
    # inside each trial consistent with Fig 4).
    obs_fr_full = _smoothed_obs_per_trial(
        cluster_df, config.time_bin_width, OBS_FR_SMOOTH_WIDTH,
    )
    motion_df["_obs_fr"] = obs_fr_full[motion_mask]

    rng = np.random.default_rng(seed)

    # Unique trials per condition (one row per trial_id × sf × orientation).
    trial_meta = (
        motion_df[["trial_id", "condition", "sf", "orientation"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    picks: list[tuple[int, str, str]] = []

    def _pick_random(cond: str, k: int) -> list[int]:
        pool = trial_meta.loc[trial_meta["condition"] == cond, "trial_id"].to_numpy()
        if pool.size == 0 or k <= 0:
            return []
        take = min(k, pool.size)
        return list(rng.choice(pool, size=take, replace=False))

    for tid in _pick_random("T_Vstatic", 2):
        picks.append((int(tid), "T_Vstatic", "T (speed only)"))
    for tid in _pick_random("V", 2):
        row = trial_meta.loc[trial_meta["trial_id"] == tid].iloc[0]
        picks.append((int(tid), "V",
                      f"V: SF={row['sf']:.3f}, OR={row['orientation']:.2f}"))

    # VT: try to sample 4 distinct (sf, orientation) combos first, then fill.
    vt_meta = trial_meta.loc[trial_meta["condition"] == "VT"].copy()
    if not vt_meta.empty:
        combos = (
            vt_meta.drop_duplicates(subset=["sf", "orientation"])
            .reset_index(drop=True)
        )
        k_combo = min(4, len(combos))
        if k_combo:
            combo_idx = rng.choice(len(combos), size=k_combo, replace=False)
            for ci in combo_idx:
                combo = combos.iloc[int(ci)]
                matching = vt_meta.loc[
                    (vt_meta["sf"] == combo["sf"])
                    & (vt_meta["orientation"] == combo["orientation"]),
                    "trial_id",
                ].to_numpy()
                if matching.size == 0:
                    continue
                tid = int(rng.choice(matching))
                picks.append((
                    tid, "VT",
                    f"VT: SF={combo['sf']:.3f}, OR={combo['orientation']:.2f}",
                ))
        if len(picks) < 8 and len(vt_meta) > k_combo:
            already = {tid for tid, _, _ in picks}
            extra_pool = vt_meta.loc[~vt_meta["trial_id"].isin(already), "trial_id"].to_numpy()
            if extra_pool.size:
                take = min(8 - len(picks), extra_pool.size)
                for tid in rng.choice(extra_pool, size=take, replace=False):
                    row = vt_meta.loc[vt_meta["trial_id"] == tid].iloc[0]
                    picks.append((
                        int(tid), "VT",
                        f"VT: SF={row['sf']:.3f}, OR={row['orientation']:.2f}",
                    ))

    if not picks:
        fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
        ax.text(
            0.5, 0.5, "no motion trials to plot",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=11, color=(0.5, 0.5, 0.5),
        )
        ax.set_axis_off()
        fig.suptitle(title, fontsize=11, fontweight="bold")
        return fig

    # Pre-compute per-trial slices + cluster-wide TF ceiling.
    trial_slices: list[dict] = []
    tf_max_global = 0.0
    for tid, cond, label in picks:
        rows = motion_df.loc[motion_df["trial_id"] == tid].sort_values("time_in_trial")
        if rows.empty:
            continue
        tf_max_global = max(tf_max_global, float(rows["tf"].max(skipna=True) or 0.0))
        trial_slices.append(dict(
            trial_id=tid, cond=cond, label=label,
            time=rows["time_in_trial"].to_numpy(dtype=float),
            speed=rows["speed"].to_numpy(dtype=float),
            tf=rows["tf"].to_numpy(dtype=float),
            obs_fr=rows["_obs_fr"].to_numpy(dtype=float),
            pred_fr=rows["_pred_fr"].to_numpy(dtype=float),
        ))

    n_rows = len(trial_slices)
    fig, axes = plt.subplots(
        n_rows, 3,
        figsize=(10, max(2.0, 1.5 * n_rows)),
        squeeze=False, constrained_layout=True,
    )
    tf_ylim = (0.0, max(tf_max_global * 1.1, 0.1))

    for r, slc in enumerate(trial_slices):
        cond_clr = _TRIAL_LEVEL_COND_COLORS.get(slc["cond"], (0.3, 0.3, 0.3))
        t = slc["time"]
        xlim = (float(t.min() - 0.05), float(t.max() + 0.05))

        ax_sp = axes[r, 0]
        ax_sp.plot(t, slc["speed"], "-", color=cond_clr, linewidth=1.2)
        ax_sp.set_ylabel("Speed (cm/s)", fontsize=8)
        ax_sp.set_xlim(*xlim)
        ax_sp.text(
            -0.20, 0.5, f"Trial {slc['trial_id']}\n{slc['label']}",
            transform=ax_sp.transAxes, ha="right", va="center",
            fontsize=7, color=(0.2, 0.2, 0.2),
        )
        if r == 0:
            ax_sp.set_title("Speed", fontweight="bold", fontsize=10)

        ax_tf = axes[r, 1]
        ax_tf.plot(t, slc["tf"], "-", color=cond_clr, linewidth=1.2)
        ax_tf.set_ylabel("TF (Hz)", fontsize=8)
        ax_tf.set_xlim(*xlim)
        ax_tf.set_ylim(*tf_ylim)
        if r == 0:
            ax_tf.set_title("Temporal Freq", fontweight="bold", fontsize=10)

        ax_fr = axes[r, 2]
        ax_fr.plot(t, slc["obs_fr"], "-", color=(0.5, 0.5, 0.5),
                   linewidth=1.2, label="Observed")
        ax_fr.plot(t, slc["pred_fr"], "-", color=(0.8, 0.2, 0.2),
                   linewidth=1.5, label="Predicted")
        ax_fr.set_ylabel("FR (Hz)", fontsize=8)
        ax_fr.set_xlim(*xlim)
        fr_max = max(
            float(np.nanmax(slc["obs_fr"]) if slc["obs_fr"].size else 0.0),
            float(np.nanmax(slc["pred_fr"]) if slc["pred_fr"].size else 0.0),
            1.0,
        )
        ax_fr.set_ylim(0.0, fr_max * 1.1)
        if r == 0:
            ax_fr.set_title("Firing Rate", fontweight="bold", fontsize=10)
            ax_fr.legend(loc="upper right", fontsize=6, frameon=False)

        for ax in (ax_sp, ax_tf, ax_fr):
            ax.tick_params(labelsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if r == n_rows - 1:
                ax.set_xlabel("Time (s)", fontsize=8)

    fig.suptitle(title, fontsize=11, fontweight="bold")
    return fig


# --------------------------------------------------------------------------- #
# Fig 4 — per-cluster model overview with scatter panels
# (MATLAB lines 2727-3290 of glm_single_cluster_analysis.m)
# --------------------------------------------------------------------------- #


# Beta-swarm group layout — matches MATLAB compute_grouped_params (line 6328).
# Keys are MATLAB group tags; values are the palette row and tick label.
# "History" added 2026-04-29 with the spike-history term — was previously
# falling through to "Other" in the beta-swarm panel.
_BETA_GROUP_ORDER: tuple[str, ...] = (
    "Intercept", "Speed", "TF", "SF", "OR", "Time", "History",
    "Spd x TF", "Spd x SF", "Spd x OR",
    "TF x SF", "TF x OR", "SF x OR",
    "Other",
)
_BETA_GROUP_COLORS: dict[str, tuple[float, float, float]] = {
    "Intercept":  (0.5, 0.5, 0.5),
    "Speed":      (0.17, 0.63, 0.17),
    "TF":         (1.0, 0.50, 0.05),
    "SF":         (0.95, 0.85, 0.10),
    "OR":         (0.84, 0.15, 0.16),
    "Time":       (0.3, 0.8, 0.8),
    "History":    (0.45, 0.20, 0.65),  # purple — distinct from stimuli + Time
    "Spd x TF":   (0.6, 0.35, 0.05),
    "Spd x SF":   (0.55, 0.50, 0.05),
    "Spd x OR":   (0.50, 0.10, 0.10),
    "TF x SF":    (0.85, 0.60, 0.05),
    "TF x OR":    (0.75, 0.20, 0.15),
    "SF x OR":    (0.55, 0.25, 0.35),
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
    """Map a design column name onto one of ``_BETA_GROUP_ORDER``.

    Interaction columns are ``<A>_x_<B>`` where both sides carry a variable
    prefix (``Spd``, ``TF``, ``SF``, ``OR``). Earlier versions used
    ``"_x_SF" in col_name`` which collapsed ``TF1_x_SF0.003`` into the Spd×SF
    bucket — split on ``_x_`` and read the left prefix so all six pairings
    land in their own group.
    """
    if col_name == "Intercept":
        return "Intercept"
    if "_x_" in col_name:
        left, right = col_name.split("_x_", 1)
        if left.startswith("Spd") and right.startswith("TF"):
            return "Spd x TF"
        if left.startswith("Spd") and right.startswith("SF"):
            return "Spd x SF"
        if left.startswith("Spd") and right.startswith("OR"):
            return "Spd x OR"
        if left.startswith("TF") and right.startswith("SF"):
            return "TF x SF"
        if left.startswith("TF") and right.startswith("OR"):
            return "TF x OR"
        if left.startswith("SF") and right.startswith("OR"):
            return "SF x OR"
        return "Other"
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
    if col_name.startswith("History_"):
        return "History"
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


def _per_bin_median_iqr(tuning: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-bin median + Q1 + Q3 across trials, NaN-tolerant.

    ``tuning`` shape (n_trials, n_bins). Returns three (n_bins,) arrays.
    Bins where all trials are NaN return (NaN, NaN, NaN).
    """
    if tuning.size == 0:
        return (np.full(0, np.nan),) * 3
    median = np.nanmedian(tuning, axis=0)
    q1 = np.nanquantile(tuning, 0.25, axis=0)
    q3 = np.nanquantile(tuning, 0.75, axis=0)
    return median, q1, q3


def _plot_observed_row(
    row_axes,
    cluster_df: pd.DataFrame,
    config: GLMConfig,
    *,
    sf_plot: np.ndarray,
    or_plot: np.ndarray,
    precomputed_bins: PrecomputedBinEdges | None = None,
) -> None:
    """Observed rates per condition, aggregated trial-first (matches MATLAB).

    For Speed and TF panels: when ``precomputed_bins`` is provided AND it
    carries the per-cluster firing-rate cache (Laura's call 2026-04-27), the
    plot reads the MATLAB-cached per-trial-per-bin firing rates directly —
    same numbers MATLAB plots from. This guarantees pixel-level parity with
    the ceph TF tuning PDFs and avoids the Python-side recomputation that
    historically diverged from MATLAB by small amounts (different motion-
    window definition, different per-bin spike-counting convention).

    For SF and OR panels: no MATLAB cache exists for these (they're
    categorical, with handfuls of levels), so per-trial means are computed
    from ``cluster_df`` directly — same as before.

    Fallback: when ``precomputed_bins`` is None (e.g. the cache is missing
    or the cluster isn't in it), Speed and TF fall back to the original
    recompute path. Logged once per cluster.
    """
    bw = config.time_bin_width
    ax_spd, ax_tf, ax_sf, ax_or = row_axes
    cluster_id = (
        int(cluster_df["cluster_id"].iloc[0])
        if not cluster_df.empty else None
    )

    fallback_spd_edges = np.linspace(
        config.speed_range[0], config.speed_range[1], 11,
    )
    fallback_tf_edges = np.linspace(
        config.tf_range[0], config.tf_range[1], 11,
    )

    # ---- Speed panel: T_Vstatic + VT (MATLAB skips V for speed since its speed=0) ----
    for cond in ("T_Vstatic", "VT"):
        cached_tuning = (
            None if precomputed_bins is None or cluster_id is None
            else precomputed_bins.speed_tuning(cond, cluster_id)
        )
        if cached_tuning is not None:
            centres = precomputed_bins.speed_centres(cond)
            median, q1, q3 = _per_bin_median_iqr(cached_tuning)
            good = ~np.isnan(median)
            if good.any():
                yerr = np.vstack([median[good] - q1[good], q3[good] - median[good]])
                ax_spd.errorbar(centres[good], median[good], yerr=yerr,
                                fmt="o-", color=COND_COLORS[cond],
                                markersize=4, linewidth=1, capsize=3)
            continue
        # Fallback: recompute from cluster_df
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        edges = (
            precomputed_bins.speed_edges(cond)
            if precomputed_bins is not None else None
        )
        if edges is None:
            edges = fallback_spd_edges
        centres = 0.5 * (edges[:-1] + edges[1:])
        median, q1, q3 = _trial_quantiles_by_bin(
            sub, value_col="speed", bin_edges=edges, bin_width=bw,
        )
        good = ~np.isnan(median)
        if good.any():
            yerr = np.vstack([median[good] - q1[good], q3[good] - median[good]])
            ax_spd.errorbar(centres[good], median[good], yerr=yerr,
                            fmt="o-", color=COND_COLORS[cond],
                            markersize=4, linewidth=1, capsize=3)

    # ---- Stationary dot at speed=0: prefer cache stationary_fr if available ----
    cached_stat_fr = None
    if precomputed_bins is not None and cluster_id is not None:
        # Cache stationary FR is recorded per (variable cache, condition);
        # any of them works since stationary is condition-independent for
        # the cluster — pick speed/VT or first available.
        for store_var, cond in (("speed", "VT"), ("speed", "T_Vstatic"), ("tf", "V")):
            arr = precomputed_bins.stationary_fr(store_var, cond, cluster_id)
            if arr is not None and arr.size > 0:
                cached_stat_fr = arr
                break
    if cached_stat_fr is not None:
        per_trial = cached_stat_fr[~np.isnan(cached_stat_fr)]
    else:
        stat = cluster_df[cluster_df["condition"] == "stationary"]
        if stat.empty:
            per_trial = np.array([], dtype=np.float64)
        else:
            per_trial = (
                stat.groupby("trial_id")["spike_count"].mean() / bw
            ).to_numpy(dtype=np.float64)
            per_trial = per_trial[~np.isnan(per_trial)]
    if per_trial.size >= 1:
        med_stat = float(np.median(per_trial))
        if per_trial.size >= 2:
            q1_stat = float(np.quantile(per_trial, 0.25))
            q3_stat = float(np.quantile(per_trial, 0.75))
        else:
            q1_stat = q3_stat = med_stat
        yerr_stat = np.array([[med_stat - q1_stat], [q3_stat - med_stat]])
        ax_spd.errorbar(0.0, med_stat, yerr=yerr_stat,
                        fmt="o", color=COND_COLORS["stationary"],
                        markersize=6, capsize=4)

    # ---- TF panel: V + VT (skip T_Vstatic since TF=0) ----
    for cond in ("V", "VT"):
        cached_tuning = (
            None if precomputed_bins is None or cluster_id is None
            else precomputed_bins.tf_tuning(cond, cluster_id)
        )
        if cached_tuning is not None:
            centres = precomputed_bins.tf_centres(cond)
            median, q1, q3 = _per_bin_median_iqr(cached_tuning)
            good = ~np.isnan(median)
            if good.any():
                yerr = np.vstack([median[good] - q1[good], q3[good] - median[good]])
                ax_tf.errorbar(centres[good], median[good], yerr=yerr,
                               fmt="o-", color=COND_COLORS[cond],
                               markersize=4, linewidth=1, capsize=3)
            continue
        # Fallback: recompute
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        edges = (
            precomputed_bins.tf_edges(cond)
            if precomputed_bins is not None else None
        )
        if edges is None:
            edges = fallback_tf_edges
        centres = 0.5 * (edges[:-1] + edges[1:])
        median, q1, q3 = _trial_quantiles_by_bin(
            sub, value_col="tf", bin_edges=edges, bin_width=bw,
        )
        good = ~np.isnan(median)
        if good.any():
            yerr = np.vstack([median[good] - q1[good], q3[good] - median[good]])
            ax_tf.errorbar(centres[good], median[good], yerr=yerr,
                           fmt="o-", color=COND_COLORS[cond],
                           markersize=4, linewidth=1, capsize=3)

    # SF: V + VT
    for cond in ("V", "VT"):
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        median, q1, q3 = _trial_quantiles_by_level(
            sub, value_col="sf", levels=sf_plot, tol=0.001, bin_width=bw,
        )
        good = ~np.isnan(median)
        if good.any():
            yerr = np.vstack([median[good] - q1[good], q3[good] - median[good]])
            ax_sf.errorbar(np.arange(sf_plot.size)[good], median[good],
                           yerr=yerr,
                           fmt="o-", color=COND_COLORS[cond],
                           markersize=5, linewidth=1, capsize=4)
    ax_sf.set_xticks(np.arange(sf_plot.size))
    ax_sf.set_xticklabels([f"{v:.3f}" for v in sf_plot], fontsize=7)

    # OR: V + VT
    for cond in ("V", "VT"):
        sub = cluster_df[cluster_df["condition"] == cond]
        if sub.empty:
            continue
        median, q1, q3 = _trial_quantiles_by_level(
            sub, value_col="orientation", levels=or_plot, tol=0.01, bin_width=bw,
        )
        good = ~np.isnan(median)
        if good.any():
            yerr = np.vstack([median[good] - q1[good], q3[good] - median[good]])
            ax_or.errorbar(np.arange(or_plot.size)[good], median[good],
                           yerr=yerr,
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
    cond_bins: dict | None = None,
    spd_centres_by_cond: dict[str, np.ndarray] | None = None,
    tf_centres_by_cond: dict[str, np.ndarray] | None = None,
) -> None:
    ax_spd, ax_tf, ax_sf, ax_or = row_axes

    # Trial-averaged path: run a dedicated helper per panel and return.
    # Keeps the code simpler than threading mode branches through every
    # low-level helper, at the cost of one duplicated per-panel block.
    if config.tuning_curve_mode == "trial-averaged" and cond_bins:
        _predict_model_row_trial_averaged(
            row_axes, beta, train_names, selected_vars, config,
            B_onset_stat=B_onset_stat,
            sf_ref=sf_ref, or_ref=or_ref,
            sf_plot=sf_plot, or_plot=or_plot,
            mean_speed=mean_speed, mean_tf=mean_tf,
            mode_sf=mode_sf, mode_or=mode_or,
            spd_centres_by_cond=spd_centres_by_cond,
            tf_centres_by_cond=tf_centres_by_cond,
            cond_bins=cond_bins,
        )
        return

    # Steady-state path: evaluate each condition's line at that
    # condition's Speed/TF bin centres (same as the Observed row above)
    # rather than a uniform linspace, so the x-axis aligns between
    # observed and predicted. Default centres for conditions without
    # their own cached edges: the old uniform spd_centres / tf_centres.
    spd_by_cond = spd_centres_by_cond or {}
    tf_by_cond = tf_centres_by_cond or {}
    spd_t = spd_by_cond.get("T_Vstatic", spd_centres)
    spd_vt = spd_by_cond.get("VT", spd_centres)
    tf_v = tf_by_cond.get("V", tf_centres)
    tf_vt = tf_by_cond.get("VT", tf_centres)
    n_spd_t = spd_t.size
    n_spd_vt = spd_vt.size
    n_tf_v = tf_v.size
    n_tf_vt = tf_vt.size

    B_onset_spd_t = np.broadcast_to(
        B_onset_steady, (n_spd_t, B_onset_steady.shape[1]),
    ).copy()
    B_onset_spd_vt = np.broadcast_to(
        B_onset_steady, (n_spd_vt, B_onset_steady.shape[1]),
    ).copy()
    B_onset_tf_v = np.broadcast_to(
        B_onset_steady, (n_tf_v, B_onset_steady.shape[1]),
    ).copy()
    B_onset_tf_vt = np.broadcast_to(
        B_onset_steady, (n_tf_vt, B_onset_steady.shape[1]),
    ).copy()

    # --- Speed sweep ---
    _plot_sweep(
        ax_spd, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec=spd_t,
        tf_vec=np.zeros(n_spd_t),
        onset_mat=B_onset_spd_t,
        sf_vec=np.full(n_spd_t, np.nan),
        or_vec=np.full(n_spd_t, np.nan),
        color=COND_COLORS["T_Vstatic"],
        xs=spd_t,
    )
    _plot_sweep(
        ax_spd, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec=spd_vt,
        tf_vec=np.full(n_spd_vt, mean_tf),
        onset_mat=B_onset_spd_vt,
        sf_vec=np.full(n_spd_vt, mode_sf),
        or_vec=np.full(n_spd_vt, mode_or),
        color=COND_COLORS["VT"],
        xs=spd_vt,
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
        speed_vec=np.zeros(n_tf_v),
        tf_vec=tf_v,
        onset_mat=B_onset_tf_v,
        sf_vec=np.full(n_tf_v, mode_sf),
        or_vec=np.full(n_tf_v, mode_or),
        color=COND_COLORS["V"],
        xs=tf_v,
    )
    _plot_sweep(
        ax_tf, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed_vec=np.full(n_tf_vt, mean_speed),
        tf_vec=tf_vt,
        onset_mat=B_onset_tf_vt,
        sf_vec=np.full(n_tf_vt, mode_sf),
        or_vec=np.full(n_tf_vt, mode_or),
        color=COND_COLORS["VT"],
        xs=tf_vt,
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


def _per_trial_summary(
    mu: np.ndarray, trial_ids: np.ndarray,
) -> dict[str, float | int]:
    """Collapse per-bin predicted rates ``mu`` to per-trial means and summarise.

    Used by ``_marginal_curve_with_spread`` and ``_marginal_levels_with_spread``
    to build the per-grid-point uncertainty band. The line value (bin-mean)
    stays the existing ``nanmean(mu)``; the band is computed on per-trial
    means so it reflects between-trial heterogeneity in non-target
    covariates rather than within-trial bin noise.
    """
    if mu is None:
        nan = float("nan")
        return {"q05": nan, "q25": nan, "q75": nan, "q95": nan,
                "std": nan, "n_trials": 0}
    # Per-trial mean; drop NaN trials (a trial whose bins all NaN'd in the
    # rate eval contributes nothing).
    series = pd.Series(np.asarray(mu, dtype=np.float64))
    per_trial = series.groupby(np.asarray(trial_ids)).mean().to_numpy()
    per_trial = per_trial[np.isfinite(per_trial)]
    n = int(per_trial.size)
    if n == 0:
        nan = float("nan")
        return {"q05": nan, "q25": nan, "q75": nan, "q95": nan,
                "std": nan, "n_trials": 0}
    return {
        "q05": float(np.percentile(per_trial, 5)),
        "q25": float(np.percentile(per_trial, 25)),
        "q75": float(np.percentile(per_trial, 75)),
        "q95": float(np.percentile(per_trial, 95)),
        "std": float(np.std(per_trial)),
        "n_trials": n,
    }


def _marginal_curve(
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    cond_df: pd.DataFrame,
    *,
    sweep_var: str,          # "Speed" | "TF"
    grid: np.ndarray,
    fix_speed: float | None = None,
    fix_tf: float | None = None,
    fix_sf: float | None = None,
    fix_or: float | None = None,
) -> np.ndarray:
    """E_observed[μ | sweep_var=g] across the given condition's motion bins.

    At each grid point ``g``, build a design matrix of length ``n_bins``
    (the condition's observed motion rows) with:

    - ``sweep_var`` set to ``g`` (broadcast)
    - each other variable either pinned to its ``fix_*`` value (if provided)
      or drawn from the condition's observed per-bin value. ``time_since_onset``
      is always drawn from observed (that's the point — option B).

    The predicted rate is ``exp(X_aligned @ β)`` per bin; we return the
    mean across bins. A grid point whose design matrix is rank-deficient
    (``_predict_rate`` returns None) yields NaN.
    """
    n_bins = len(cond_df)
    if n_bins == 0:
        return np.full(len(grid), np.nan)
    obs_speed = cond_df["speed"].to_numpy(dtype=np.float64)
    obs_tf = cond_df["tf"].to_numpy(dtype=np.float64)
    obs_sf = cond_df["sf"].to_numpy(dtype=np.float64)
    obs_or = cond_df["orientation"].to_numpy(dtype=np.float64)
    onset_t = cond_df["time_since_onset"].to_numpy(dtype=np.float64)
    # Clip observed onset into the basis support; values ≤0 evaluate to 0
    # (the stationary pre-onset window) by construction in onset_kernel_basis.
    B_onset = onset_kernel_basis(onset_t, config.n_onset_bases, config.onset_range[1])

    rates = np.empty(len(grid))
    for i, g in enumerate(grid):
        if sweep_var == "Speed":
            speed = np.full(n_bins, g)
            tf = np.full(n_bins, fix_tf) if fix_tf is not None else obs_tf
        elif sweep_var == "TF":
            speed = np.full(n_bins, fix_speed) if fix_speed is not None else obs_speed
            tf = np.full(n_bins, g)
        else:
            raise ValueError(f"unsupported sweep_var {sweep_var!r}")
        sf = np.full(n_bins, fix_sf) if fix_sf is not None else obs_sf
        orientation = np.full(n_bins, fix_or) if fix_or is not None else obs_or
        mu = _predict_rate(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            speed, tf, B_onset, sf, orientation,
        )
        rates[i] = np.nan if mu is None else float(np.nanmean(mu))
    return rates


def _marginal_curve_with_spread(
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    cond_df: pd.DataFrame,
    *,
    sweep_var: str,
    grid: np.ndarray,
    fix_speed: float | None = None,
    fix_tf: float | None = None,
    fix_sf: float | None = None,
    fix_or: float | None = None,
) -> dict[str, np.ndarray]:
    """As ``_marginal_curve`` but also returns per-trial spread per grid point.

    Returns a dict with keys ``mean`` (matches ``_marginal_curve``'s return),
    ``q05`` / ``q25`` / ``q75`` / ``q95`` / ``std`` (per-trial), and
    ``n_trials``. The line on the plot stays at ``mean``; the band is
    drawn from the quantile arrays. ``mean`` is the bin-weighted mean (as
    before) for backwards compatibility; the per-trial mean is implicitly
    encoded in the spread arrays.
    """
    n_grid = len(grid)
    nan_arr = np.full(n_grid, np.nan)
    n_bins = len(cond_df)
    if n_bins == 0:
        return {
            "mean": nan_arr.copy(),
            "q05": nan_arr.copy(),
            "q25": nan_arr.copy(),
            "q75": nan_arr.copy(),
            "q95": nan_arr.copy(),
            "std": nan_arr.copy(),
            "n_trials": np.zeros(n_grid, dtype=int),
        }
    trial_ids = cond_df["trial_id"].to_numpy()
    obs_speed = cond_df["speed"].to_numpy(dtype=np.float64)
    obs_tf = cond_df["tf"].to_numpy(dtype=np.float64)
    obs_sf = cond_df["sf"].to_numpy(dtype=np.float64)
    obs_or = cond_df["orientation"].to_numpy(dtype=np.float64)
    onset_t = cond_df["time_since_onset"].to_numpy(dtype=np.float64)
    B_onset = onset_kernel_basis(onset_t, config.n_onset_bases, config.onset_range[1])

    out = {
        "mean": np.empty(n_grid),
        "q05": np.empty(n_grid),
        "q25": np.empty(n_grid),
        "q75": np.empty(n_grid),
        "q95": np.empty(n_grid),
        "std": np.empty(n_grid),
        "n_trials": np.zeros(n_grid, dtype=int),
    }
    for i, g in enumerate(grid):
        if sweep_var == "Speed":
            speed = np.full(n_bins, g)
            tf = np.full(n_bins, fix_tf) if fix_tf is not None else obs_tf
        elif sweep_var == "TF":
            speed = np.full(n_bins, fix_speed) if fix_speed is not None else obs_speed
            tf = np.full(n_bins, g)
        else:
            raise ValueError(f"unsupported sweep_var {sweep_var!r}")
        sf = np.full(n_bins, fix_sf) if fix_sf is not None else obs_sf
        orientation = np.full(n_bins, fix_or) if fix_or is not None else obs_or
        mu = _predict_rate(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            speed, tf, B_onset, sf, orientation,
        )
        out["mean"][i] = np.nan if mu is None else float(np.nanmean(mu))
        s = _per_trial_summary(mu, trial_ids)
        out["q05"][i] = s["q05"]
        out["q25"][i] = s["q25"]
        out["q75"][i] = s["q75"]
        out["q95"][i] = s["q95"]
        out["std"][i] = s["std"]
        out["n_trials"][i] = s["n_trials"]
    return out


def _simulated_curve_with_spread(
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    cond_df: pd.DataFrame,
    *,
    sweep_var: str,                # "Speed" | "TF"
    bin_edges: np.ndarray,         # display bin edges (cache's 5%-quantile)
    bin_centres: np.ndarray,       # display bin centres
    n_iterations: int = 100,
    rng_seed: int = 0,
) -> dict[str, np.ndarray]:
    """Parametric-bootstrap reconstruction of the tuning curve.

    Algorithm (matches Laura's spec, 2026-04-28):

    1. Predict ``λ_pred[t]`` per 100 ms time-bin using each trial's
       *actual* covariates (no forcing of the sweep variable). This is the
       training granularity — same as the GLM was fit on.
    2. Draw simulated spike counts ``y_sim[t] ~ Poisson(λ_pred[t] · Δt)``
       per 100 ms bin.
    3. Convert to firing rate ``fr_sim[t] = y_sim[t] / Δt``.
    4. Group those simulated per-bin firing rates by ``cond_df[sweep_var]``
       falling into the cache's display bins; per (trial × display-bin),
       compute the mean simulated firing rate.
    5. Per display bin: median + IQR ACROSS TRIALS (not across iterations).
    6. Average the per-display-bin median / IQR across ``n_iterations``
       iterations to dampen MC sampling noise.

    Returns a dict with the same shape as ``_marginal_curve_with_spread``
    so the existing ``_band_edges`` / ``_draw_band`` plumbing reuses
    untouched.
    """
    n_grid = bin_centres.size
    nan_arr = np.full(n_grid, np.nan)
    out_empty = {
        "mean": nan_arr.copy(),
        "q05": nan_arr.copy(),
        "q25": nan_arr.copy(),
        "q75": nan_arr.copy(),
        "q95": nan_arr.copy(),
        "std": nan_arr.copy(),
        "n_trials": np.zeros(n_grid, dtype=int),
    }

    n_bins = len(cond_df)
    if n_bins == 0:
        return out_empty

    bw = config.time_bin_width  # seconds per training bin (default 0.1)
    trial_ids = cond_df["trial_id"].to_numpy()
    obs_speed = cond_df["speed"].to_numpy(dtype=np.float64)
    obs_tf = cond_df["tf"].to_numpy(dtype=np.float64)
    obs_sf = cond_df["sf"].to_numpy(dtype=np.float64)
    obs_or = cond_df["orientation"].to_numpy(dtype=np.float64)
    onset_t = cond_df["time_since_onset"].to_numpy(dtype=np.float64)
    B_onset = onset_kernel_basis(onset_t, config.n_onset_bases, config.onset_range[1])

    # 1. Predict λ at each trial's actual covariates — one λ per 100 ms bin.
    mu = _predict_rate(
        beta, train_names, selected_vars, config, sf_ref, or_ref,
        obs_speed, obs_tf, B_onset, obs_sf, obs_or,
    )
    if mu is None:
        return out_empty
    mu = np.asarray(mu, dtype=np.float64)

    # Display-bin assignment of each 100ms bin based on its sweep_var value.
    sweep_vals = obs_speed if sweep_var == "Speed" else obs_tf
    # `np.digitize` returns 1..n_bins for in-range values; clamp anything
    # outside to (-1) so it gets dropped.
    digitised = np.digitize(sweep_vals, bin_edges) - 1
    valid_mask = (digitised >= 0) & (digitised < n_grid)

    # Per (trial, display-bin) we want the mean simulated FR. Set up
    # accumulators and a stable trial index ordering.
    unique_trials = np.unique(trial_ids)
    n_trials_total = unique_trials.size
    trial_to_idx = {int(t): i for i, t in enumerate(unique_trials)}
    bin_indices = np.array([trial_to_idx[int(t)] for t in trial_ids])

    # Pre-build (trial_idx, display_bin) flat indexer for valid bins.
    valid_idx = np.where(valid_mask)[0]
    flat_idx = bin_indices[valid_idx] * n_grid + digitised[valid_idx]
    n_cells = n_trials_total * n_grid

    rng = np.random.default_rng(rng_seed)

    # Accumulators: per-iteration per-display-bin median/quantiles, averaged
    # across iterations to reduce MC sampling noise.
    median_acc = np.zeros(n_grid)
    q05_acc = np.zeros(n_grid)
    q25_acc = np.zeros(n_grid)
    q75_acc = np.zeros(n_grid)
    q95_acc = np.zeros(n_grid)
    std_acc = np.zeros(n_grid)
    n_obs_acc = np.zeros(n_grid, dtype=np.int64)
    n_eff_iters = np.zeros(n_grid, dtype=np.int64)

    for _ in range(n_iterations):
        # Step 2: Poisson draws per 100ms bin.
        y_sim = rng.poisson(mu * bw)
        fr_sim = y_sim.astype(np.float64) / bw

        # Step 4: per (trial, display-bin) mean simulated FR. Use bincount
        # for vectorised aggregation.
        valid_fr = fr_sim[valid_idx]
        sums = np.bincount(flat_idx, weights=valid_fr, minlength=n_cells)
        counts = np.bincount(flat_idx, minlength=n_cells)
        with np.errstate(invalid="ignore", divide="ignore"):
            per_cell_mean = np.where(counts > 0, sums / counts, np.nan)
        per_cell_mean = per_cell_mean.reshape(n_trials_total, n_grid)

        # Step 5: per display bin, summarise across trials (NaN-tolerant).
        for j in range(n_grid):
            col = per_cell_mean[:, j]
            col = col[~np.isnan(col)]
            if col.size == 0:
                continue
            n_eff_iters[j] += 1
            n_obs_acc[j] += col.size
            median_acc[j] += np.median(col)
            q05_acc[j] += np.percentile(col, 5)
            q25_acc[j] += np.percentile(col, 25)
            q75_acc[j] += np.percentile(col, 75)
            q95_acc[j] += np.percentile(col, 95)
            std_acc[j] += np.std(col)

    # Step 6: average across iterations.
    out = {
        "mean":     np.where(n_eff_iters > 0, median_acc / np.maximum(n_eff_iters, 1), np.nan),
        "q05":      np.where(n_eff_iters > 0, q05_acc / np.maximum(n_eff_iters, 1), np.nan),
        "q25":      np.where(n_eff_iters > 0, q25_acc / np.maximum(n_eff_iters, 1), np.nan),
        "q75":      np.where(n_eff_iters > 0, q75_acc / np.maximum(n_eff_iters, 1), np.nan),
        "q95":      np.where(n_eff_iters > 0, q95_acc / np.maximum(n_eff_iters, 1), np.nan),
        "std":      np.where(n_eff_iters > 0, std_acc / np.maximum(n_eff_iters, 1), np.nan),
        # n_trials = average per-bin trial count across iterations (rounded).
        "n_trials": np.where(
            n_eff_iters > 0,
            (n_obs_acc / np.maximum(n_eff_iters, 1)).astype(int),
            0,
        ),
    }
    return out


def _marginal_levels(
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    cond_df: pd.DataFrame,
    *,
    level_kind: str,         # "SF" | "OR"
    levels: np.ndarray,
    fix_speed: float | None,
    fix_tf: float | None,
    fix_sf: float | None = None,
    fix_or: float | None = None,
) -> np.ndarray:
    """Trial-averaged version of _plot_levels — one rate per categorical level.

    At each level ``L``, the swept variable is pinned to ``L`` while
    Speed / TF / the other categorical can be fixed or marginalised over
    observed bins. Onset always uses the condition's observed distribution.
    """
    n_bins = len(cond_df)
    if n_bins == 0:
        return np.full(len(levels), np.nan)
    obs_speed = cond_df["speed"].to_numpy(dtype=np.float64)
    obs_tf = cond_df["tf"].to_numpy(dtype=np.float64)
    obs_sf = cond_df["sf"].to_numpy(dtype=np.float64)
    obs_or = cond_df["orientation"].to_numpy(dtype=np.float64)
    onset_t = cond_df["time_since_onset"].to_numpy(dtype=np.float64)
    B_onset = onset_kernel_basis(onset_t, config.n_onset_bases, config.onset_range[1])

    rates = np.empty(len(levels))
    for i, L in enumerate(levels):
        speed = np.full(n_bins, fix_speed) if fix_speed is not None else obs_speed
        tf = np.full(n_bins, fix_tf) if fix_tf is not None else obs_tf
        if level_kind == "SF":
            sf = np.full(n_bins, L)
            orientation = np.full(n_bins, fix_or) if fix_or is not None else obs_or
        else:
            sf = np.full(n_bins, fix_sf) if fix_sf is not None else obs_sf
            orientation = np.full(n_bins, L)
        mu = _predict_rate(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            speed, tf, B_onset, sf, orientation,
        )
        rates[i] = np.nan if mu is None else float(np.nanmean(mu))
    return rates


def _marginal_levels_with_spread(
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    cond_df: pd.DataFrame,
    *,
    level_kind: str,
    levels: np.ndarray,
    fix_speed: float | None,
    fix_tf: float | None,
    fix_sf: float | None = None,
    fix_or: float | None = None,
) -> dict[str, np.ndarray]:
    """Trial-spread analogue of ``_marginal_levels``."""
    n_levels = len(levels)
    nan_arr = np.full(n_levels, np.nan)
    n_bins = len(cond_df)
    if n_bins == 0:
        return {
            "mean": nan_arr.copy(),
            "q05": nan_arr.copy(),
            "q25": nan_arr.copy(),
            "q75": nan_arr.copy(),
            "q95": nan_arr.copy(),
            "std": nan_arr.copy(),
            "n_trials": np.zeros(n_levels, dtype=int),
        }
    trial_ids = cond_df["trial_id"].to_numpy()
    obs_speed = cond_df["speed"].to_numpy(dtype=np.float64)
    obs_tf = cond_df["tf"].to_numpy(dtype=np.float64)
    obs_sf = cond_df["sf"].to_numpy(dtype=np.float64)
    obs_or = cond_df["orientation"].to_numpy(dtype=np.float64)
    onset_t = cond_df["time_since_onset"].to_numpy(dtype=np.float64)
    B_onset = onset_kernel_basis(onset_t, config.n_onset_bases, config.onset_range[1])

    out = {
        "mean": np.empty(n_levels),
        "q05": np.empty(n_levels),
        "q25": np.empty(n_levels),
        "q75": np.empty(n_levels),
        "q95": np.empty(n_levels),
        "std": np.empty(n_levels),
        "n_trials": np.zeros(n_levels, dtype=int),
    }
    for i, L in enumerate(levels):
        speed = np.full(n_bins, fix_speed) if fix_speed is not None else obs_speed
        tf = np.full(n_bins, fix_tf) if fix_tf is not None else obs_tf
        if level_kind == "SF":
            sf = np.full(n_bins, L)
            orientation = np.full(n_bins, fix_or) if fix_or is not None else obs_or
        else:
            sf = np.full(n_bins, fix_sf) if fix_sf is not None else obs_sf
            orientation = np.full(n_bins, L)
        mu = _predict_rate(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            speed, tf, B_onset, sf, orientation,
        )
        out["mean"][i] = np.nan if mu is None else float(np.nanmean(mu))
        s = _per_trial_summary(mu, trial_ids)
        out["q05"][i] = s["q05"]
        out["q25"][i] = s["q25"]
        out["q75"][i] = s["q75"]
        out["q95"][i] = s["q95"]
        out["std"][i] = s["std"]
        out["n_trials"][i] = s["n_trials"]
    return out


def _band_edges(
    summary: dict[str, np.ndarray], mode: str,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Resolve (lo, hi) band edges from the per-trial summary dict.

    ``mode`` ∈ {"none", "covariate-spread", "simulated", "wide-quantile",
    "std"}. Returns ``None`` when no band should be drawn. ``"iqr"`` is a
    deprecated alias for ``"covariate-spread"`` (kept for back-compat with
    pre-2026-04-28 CLI invocations). NaN-padded grid points where
    ``n_trials < MIN_TRIALS_FOR_BAND`` are masked with NaN so
    ``fill_between`` skips them silently.
    """
    if mode == "none":
        return None
    if mode == "iqr":  # deprecated alias
        mode = "covariate-spread"
    n = summary.get("n_trials")
    sparse = None if n is None else (n < MIN_TRIALS_FOR_BAND)
    if mode in ("covariate-spread", "simulated"):
        # Both modes populate q25/q75 via the same dict shape; only the
        # underlying computation differs (per-trial covariate spread vs
        # parametric-bootstrap trial spread).
        lo = summary["q25"].copy()
        hi = summary["q75"].copy()
    elif mode == "wide-quantile":
        lo = summary["q05"].copy()
        hi = summary["q95"].copy()
    elif mode == "std":
        mean = summary["mean"]
        std = summary["std"]
        lo = mean - std
        hi = mean + std
    else:
        raise ValueError(f"unknown tuning_curve_uncertainty mode: {mode!r}")
    if sparse is not None and np.any(sparse):
        lo = lo.copy(); hi = hi.copy()
        lo[sparse] = np.nan
        hi[sparse] = np.nan
    return lo, hi


def _draw_band(ax, x: np.ndarray, summary: dict[str, np.ndarray],
               color, mode: str, *, alpha: float = 0.18) -> None:
    edges = _band_edges(summary, mode)
    if edges is None:
        return
    lo, hi = edges
    ax.fill_between(x, lo, hi, color=color, alpha=alpha, linewidth=0)


def _predict_model_row_trial_averaged(
    row_axes,
    beta: np.ndarray,
    train_names: list[str],
    selected_vars: list[str],
    config: GLMConfig,
    *,
    B_onset_stat: np.ndarray,
    sf_ref: np.ndarray,
    or_ref: np.ndarray,
    sf_plot: np.ndarray,
    or_plot: np.ndarray,
    mean_speed: float,
    mean_tf: float,
    mode_sf: float,
    mode_or: float,
    cond_bins: dict,
    spd_centres_by_cond: dict[str, np.ndarray] | None = None,
    tf_centres_by_cond: dict[str, np.ndarray] | None = None,
) -> None:
    """Trial-averaged analogue of ``_predict_model_row``.

    Condition semantics (same as the steady-state path, with the onset
    kernel marginalised over each condition's observed ``time_since_onset``
    distribution — Laura 2026-04-23 option B):

    - T_Vstatic Speed sweep: fix tf=0 (visual static), marginalise rest.
    - VT Speed sweep: marginalise everything over the VT condition's bins.
    - V TF sweep: fix speed=0 (mouse not translating), marginalise rest.
    - VT TF sweep: marginalise everything over the VT condition's bins.
    - SF / OR levels: marginalise within V and VT conditions, respectively.

    The stationary dot uses observed stationary bins (onset=0 by construction).
    Each condition line is evaluated at the same per-condition Speed/TF
    bin centres the Observed row uses (MATLAB 5%-quantile edges cached
    in ``csvs/{tuning,tf_tuning}_curves``) so x-axis alignment is exact.
    """
    ax_spd, ax_tf, ax_sf, ax_or = row_axes
    spd_by_cond = spd_centres_by_cond or {}
    tf_by_cond = tf_centres_by_cond or {}
    band_mode = getattr(config, "tuning_curve_uncertainty", "none")
    if band_mode == "iqr":  # back-compat alias
        band_mode = "covariate-spread"
    n_iter = int(getattr(config, "n_bootstrap_iterations", 100))

    def _summarise_speed(cond_df, grid, *, fix_tf=None):
        if band_mode == "simulated":
            edges = np.concatenate([
                [grid[0] - (grid[1] - grid[0]) / 2.0],
                0.5 * (grid[:-1] + grid[1:]),
                [grid[-1] + (grid[-1] - grid[-2]) / 2.0],
            ])
            return _simulated_curve_with_spread(
                beta, train_names, selected_vars, config, sf_ref, or_ref,
                cond_df, sweep_var="Speed", bin_edges=edges, bin_centres=grid,
                n_iterations=n_iter,
            )
        kwargs = {} if fix_tf is None else {"fix_tf": fix_tf}
        return _marginal_curve_with_spread(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            cond_df, sweep_var="Speed", grid=grid, **kwargs,
        )

    def _summarise_tf(cond_df, grid, *, fix_speed=None):
        if band_mode == "simulated":
            edges = np.concatenate([
                [grid[0] - (grid[1] - grid[0]) / 2.0],
                0.5 * (grid[:-1] + grid[1:]),
                [grid[-1] + (grid[-1] - grid[-2]) / 2.0],
            ])
            return _simulated_curve_with_spread(
                beta, train_names, selected_vars, config, sf_ref, or_ref,
                cond_df, sweep_var="TF", bin_edges=edges, bin_centres=grid,
                n_iterations=n_iter,
            )
        kwargs = {} if fix_speed is None else {"fix_speed": fix_speed}
        return _marginal_curve_with_spread(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            cond_df, sweep_var="TF", grid=grid, **kwargs,
        )

    # --- Speed sweep ---
    t_bins = cond_bins.get("T_Vstatic")
    spd_t = spd_by_cond.get("T_Vstatic")
    if t_bins is not None and spd_t is not None:
        summary = _summarise_speed(t_bins, spd_t, fix_tf=0.0)
        _draw_band(ax_spd, spd_t, summary, COND_COLORS["T_Vstatic"], band_mode)
        ax_spd.plot(spd_t, summary["mean"], "o-",
                    color=COND_COLORS["T_Vstatic"], markersize=4, linewidth=1.2)
    vt_bins = cond_bins.get("VT")
    spd_vt = spd_by_cond.get("VT")
    if vt_bins is not None and spd_vt is not None:
        summary = _summarise_speed(vt_bins, spd_vt)
        _draw_band(ax_spd, spd_vt, summary, COND_COLORS["VT"], band_mode)
        ax_spd.plot(spd_vt, summary["mean"], "o-",
                    color=COND_COLORS["VT"], markersize=4, linewidth=1.2)
    # Stationary dot at speed=0 — compute from observed stationary bins.
    _plot_point(
        ax_spd, beta, train_names, selected_vars, config, sf_ref, or_ref,
        speed=0.0, tf=0.0, onset_row=B_onset_stat[0],
        sf=np.nan, orientation=np.nan,
        color=COND_COLORS["stationary"], x=0.0,
    )

    # --- TF sweep ---
    v_bins = cond_bins.get("V")
    tf_v = tf_by_cond.get("V")
    if v_bins is not None and tf_v is not None:
        summary = _summarise_tf(v_bins, tf_v, fix_speed=0.0)
        _draw_band(ax_tf, tf_v, summary, COND_COLORS["V"], band_mode)
        ax_tf.plot(tf_v, summary["mean"], "o-",
                   color=COND_COLORS["V"], markersize=4, linewidth=1.2)
    tf_vt = tf_by_cond.get("VT")
    if vt_bins is not None and tf_vt is not None:
        summary = _summarise_tf(vt_bins, tf_vt)
        _draw_band(ax_tf, tf_vt, summary, COND_COLORS["VT"], band_mode)
        ax_tf.plot(tf_vt, summary["mean"], "o-",
                   color=COND_COLORS["VT"], markersize=4, linewidth=1.2)

    # --- SF levels (V and VT) ---
    sf_x = np.arange(sf_plot.size)
    if v_bins is not None:
        summary = _marginal_levels_with_spread(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            v_bins, level_kind="SF", levels=sf_plot,
            fix_speed=0.0, fix_tf=None, fix_or=mode_or,
        )
        _draw_band(ax_sf, sf_x, summary, COND_COLORS["V"], band_mode)
        ax_sf.plot(sf_x, summary["mean"], "o-",
                   color=COND_COLORS["V"], markersize=5, linewidth=1.2)
    if vt_bins is not None:
        summary = _marginal_levels_with_spread(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            vt_bins, level_kind="SF", levels=sf_plot,
            fix_speed=None, fix_tf=None, fix_or=mode_or,
        )
        _draw_band(ax_sf, sf_x, summary, COND_COLORS["VT"], band_mode)
        ax_sf.plot(sf_x, summary["mean"], "o-",
                   color=COND_COLORS["VT"], markersize=5, linewidth=1.2)
    ax_sf.set_xticks(sf_x)
    ax_sf.set_xticklabels([f"{v:.3f}" for v in sf_plot], fontsize=7)

    # --- OR levels (V and VT) ---
    or_x = np.arange(or_plot.size)
    if v_bins is not None:
        summary = _marginal_levels_with_spread(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            v_bins, level_kind="OR", levels=or_plot,
            fix_speed=0.0, fix_tf=None, fix_sf=mode_sf,
        )
        _draw_band(ax_or, or_x, summary, COND_COLORS["V"], band_mode)
        ax_or.plot(or_x, summary["mean"], "o-",
                   color=COND_COLORS["V"], markersize=5, linewidth=1.2)
    if vt_bins is not None:
        summary = _marginal_levels_with_spread(
            beta, train_names, selected_vars, config, sf_ref, or_ref,
            vt_bins, level_kind="OR", levels=or_plot,
            fix_speed=None, fix_tf=None, fix_sf=mode_sf,
        )
        _draw_band(ax_or, or_x, summary, COND_COLORS["VT"], band_mode)
        ax_or.plot(or_x, summary["mean"], "o-",
                   color=COND_COLORS["VT"], markersize=5, linewidth=1.2)
    ax_or.set_xticks(or_x)
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


def _trial_quantiles_by_bin(
    sub_df: pd.DataFrame,
    *,
    value_col: str,
    bin_edges: np.ndarray,
    bin_width: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-trial median / Q1 / Q3 rate inside each ``bin_edges`` bucket.

    Returns median, q1, q3 arrays with NaN where the cell has fewer than
    ``MIN_TRIALS_PER_CELL`` contributing trials (MATLAB required ≥2, we
    default to 1 so observed data is always visible — single-trial cells
    render with q1==q3==median, i.e. no visible error bar).

    Switched from mean±std to quartiles: firing-rate-per-trial distributions
    are right-skewed with occasional high-count trials, so std is a poor
    summary of where the middle of the distribution sits. Q1/Q3 is robust
    and renders asymmetric when the distribution is one-sided.
    """
    x = sub_df[value_col].to_numpy(dtype=np.float64)
    rate = sub_df["spike_count"].to_numpy(dtype=np.float64) / bin_width
    trial_ids = sub_df["trial_id"].to_numpy(dtype=np.int64)
    n_bins = bin_edges.size - 1
    median = np.full(n_bins, np.nan)
    q1 = np.full(n_bins, np.nan)
    q3 = np.full(n_bins, np.nan)
    good = ~np.isnan(x) & ~np.isnan(rate)
    x, rate, trial_ids = x[good], rate[good], trial_ids[good]
    if x.size == 0:
        return median, q1, q3
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
            median[b] = float(np.median(arr))
            if arr.size >= 2:
                q1[b] = float(np.quantile(arr, 0.25))
                q3[b] = float(np.quantile(arr, 0.75))
            else:
                q1[b] = median[b]
                q3[b] = median[b]
    return median, q1, q3


def _trial_quantiles_by_level(
    sub_df: pd.DataFrame,
    *,
    value_col: str,
    levels: np.ndarray,
    tol: float,
    bin_width: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-trial median / Q1 / Q3 rate at each discrete level (SF / OR)."""
    x = sub_df[value_col].to_numpy(dtype=np.float64)
    rate = sub_df["spike_count"].to_numpy(dtype=np.float64) / bin_width
    trial_ids = sub_df["trial_id"].to_numpy(dtype=np.int64)
    median = np.full(levels.size, np.nan)
    q1 = np.full(levels.size, np.nan)
    q3 = np.full(levels.size, np.nan)
    good = ~np.isnan(x) & ~np.isnan(rate)
    x, rate, trial_ids = x[good], rate[good], trial_ids[good]
    if x.size == 0:
        return median, q1, q3
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
            median[i] = float(np.median(arr))
            if arr.size >= 2:
                q1[i] = float(np.quantile(arr, 0.25))
                q3[i] = float(np.quantile(arr, 0.75))
            else:
                q1[i] = median[i]
                q3[i] = median[i]
    return median, q1, q3


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
