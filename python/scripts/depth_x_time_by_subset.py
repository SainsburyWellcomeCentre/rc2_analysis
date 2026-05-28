"""Depth × time z-FR heatmaps, faceted by GLM-selected-variables subset.

Companion to ``depth_variable_selection.py`` (the non-time-resolved
counterpart): for each mutually-exclusive cluster subset defined by the
GLM Selected-model boolean flags, plot the population-average z-FR as
a heatmap over (cortical depth, time relative to motion onset).

Layout: 7 rows (cluster subsets) × 3 columns (V / VT / T_Vstatic).
Subsets, per Laura's 2026-05-28 choice:

  All — every prefilter-passed cluster in the GLM cohort (n=80)
  only_S — Speed in Selected, ME and Visual absent
  only_M — ME in Selected, Speed and Visual absent
  S+M — Speed and ME in Selected, no Visual
  S+V — Speed and any of {TF, SF, OR} in Selected, no ME
  M+V — ME and any of {TF, SF, OR} in Selected, no Speed
  S+M+V — all three families in Selected

Inside each cell:
- y = cortical-layer label from the formatted .mat (VISp2/3 → VISp6b,
  ordered superficial → deep)
- x = time relative to motion onset, 100 ms bins from −2 s to +2 s
- colour = mean z-FR across the clusters belonging to that (subset,
  layer) cell. Z-score is per-cluster against its own baseline-window
  (−1 s to −0.1 s) statistics, matching the Lohuis 2024 convention.

Output: ``~/local_data/motion_clouds/figures/glm/exploration/depth_x_time_by_subset/{pdf,png}``.

Default cohort = GLM cohort (~80 clusters from
``glm_model_comparison.csv``).
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import FormattedDataReader, StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.io import load_probe_data


DATA_ROOT = Path("~/local_data/motion_clouds").expanduser()
MAT_DIR = DATA_ROOT / "formatted_data_3probe"
SEQ_PATH = DATA_ROOT / "motion_cloud_sequence_250414.mat"
FOL_PATH = DATA_ROOT / "image_folders.mat"
GLM_CSV = (
    DATA_ROOT / "figures" / "glm" / "current_with_ME_3_probes"
    / "glm_model_comparison.csv"
)
OUT_DIR = DATA_ROOT / "figures" / "glm" / "exploration" / "depth_x_time_by_subset"

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
CONDITIONS = ("V", "VT", "T_Vstatic")
TRIAL_CONDITION_MAP = {1: "VT", 2: "V", 3: "T_Vstatic"}

# Time grid relative to motion onset. 50 ms bins for the depth × time
# heatmaps (finer than the production GLM 100 ms but still smoothing
# enough for population averages); window ±4 s matches the fr_me_corr
# time axis so the figures can be cross-referenced.
BIN_WIDTH_S = 0.05
T_LO_S = -4.0
T_HI_S = 4.0
N_BINS = int(round((T_HI_S - T_LO_S) / BIN_WIDTH_S))   # 160
BIN_EDGES_REL = T_LO_S + BIN_WIDTH_S * np.arange(N_BINS + 1)
BIN_CENTRES_REL = 0.5 * (BIN_EDGES_REL[:-1] + BIN_EDGES_REL[1:])

# Baseline window for per-cluster z-score normalisation (Lohuis convention).
BASELINE_LO_S = -1.0
BASELINE_HI_S = -0.1

# Depth bin width (μm). Cohort depths span roughly 220–1160 μm; 50 μm
# bins → ~19 rows, fine enough to see layer-scale structure without
# leaving most bins empty.
DEPTH_BIN_UM = 50.0
DEPTH_LO_UM = 200.0
DEPTH_HI_UM = 1200.0

# Layer label order, superficial → deep. Used to derive data-driven
# layer-boundary depths (midpoints between consecutive layer means)
# in the same style as depth_variable_selection.py.
LAYER_ORDER = ("VISp2/3", "VISp4", "VISp5", "VISp6a", "VISp6b")

# Subset order — matches Laura's request, mutex categories.
SUBSET_ORDER = ("All", "only_S", "only_M", "S+M", "S+V", "M+V", "S+M+V")


# -- subset assignment ------------------------------------------------------
def _assign_subset(speed: bool, me: bool, visual: bool) -> str | None:
    """Map (Speed, ME, Visual) booleans to one of the mutex subset labels,
    or None if the cluster falls into a dropped category (V-only / none)."""
    s, m, v = bool(speed), bool(me), bool(visual)
    if not (s or m or v):
        return None
    if s and not m and not v:
        return "only_S"
    if m and not s and not v:
        return "only_M"
    if v and not s and not m:
        return None  # V-only is too sparse (n=2 on this cohort)
    if s and m and not v:
        return "S+M"
    if s and v and not m:
        return "S+V"
    if m and v and not s:
        return "M+V"
    if s and m and v:
        return "S+M+V"
    return None


def _load_cluster_table() -> pd.DataFrame:
    """One row per cohort cluster with its subset label."""
    df = pd.read_csv(GLM_CSV)
    s = df["time_is_speed_tuned"].astype(bool)
    m = df["time_is_me_face_tuned"].astype(bool)
    v = (
        df["time_is_tf_tuned"].astype(bool)
        | df["time_is_sf_tuned"].astype(bool)
        | df["time_is_or_tuned"].astype(bool)
    )
    df["subset"] = [
        _assign_subset(si, mi, vi) for si, mi, vi in zip(s, m, v)
    ]
    return df[["probe_id", "cluster_id", "subset"]].copy()


# -- per-probe binning ------------------------------------------------------
def _bin_one_probe(probe_stem: str, lookup: StimulusLookup) -> dict:
    """Bin each VISp cluster's spike train into 50 ms × ±4 s windows around
    motion onset, grouped by condition.

    Built on top of ``rc2_glm.io.load_probe_data`` — per the
    feedback-recycle-rc2-modules memory, trial wrappers
    (``trial.probe_t``, ``trial.motion_mask``, ``trial.condition``,
    ``trial.excluded``) are the canonical path. Don't re-derive any of
    that here. Cluster depth is the one piece not in ``ClusterData``;
    we read it from the reader for each cluster index.
    """
    mat = MAT_DIR / f"{probe_stem}.mat"
    probe = load_probe_data(
        mat, config=GLMConfig(), stimulus_lookup=lookup, visp_only=True,
    )
    # Per-cluster depth (not carried on ClusterData; cheap to fetch).
    with FormattedDataReader(mat) as r:
        cluster_depths = [r.cluster_depth(c.cluster_idx) for c in probe.clusters]

    # Per-trial motion-onset time in probe time, condition label.
    trial_t_onset: list[float] = []
    trial_cond: list[str | None] = []
    for trial in probe.trials:
        if trial.excluded or trial.condition not in CONDITIONS:
            trial_t_onset.append(float("nan"))
            trial_cond.append(None)
            continue
        motion_idx = np.flatnonzero(trial.motion_mask)
        if motion_idx.size == 0:
            trial_t_onset.append(float("nan"))
            trial_cond.append(None)
            continue
        trial_t_onset.append(float(trial.probe_t[int(motion_idx[0])]))
        trial_cond.append(trial.condition)

    # Per-cluster binning.
    out = {
        "probe_id": probe.probe_id,
        "cluster_ids": [],
        "cluster_regions": [],
        "cluster_depths": [],
        "fr_by_cluster_cond": [],  # per cluster: dict[condition] -> (n_trials, N_BINS)
    }
    for ci, cluster in enumerate(probe.clusters):
        per_cond_fr: dict[str, list[np.ndarray]] = {c: [] for c in CONDITIONS}
        for ti, (t_onset, cond) in enumerate(zip(trial_t_onset, trial_cond)):
            if not np.isfinite(t_onset) or cond is None:
                continue
            edges = t_onset + BIN_EDGES_REL
            counts, _ = np.histogram(cluster.spike_times, bins=edges)
            per_cond_fr[cond].append(counts.astype(np.float64) / BIN_WIDTH_S)
        stacked = {
            c: (np.vstack(lst) if lst else np.zeros((0, N_BINS), dtype=np.float64))
            for c, lst in per_cond_fr.items()
        }
        out["cluster_ids"].append(cluster.cluster_id)
        out["cluster_regions"].append(cluster.region)
        out["cluster_depths"].append(cluster_depths[ci])
        out["fr_by_cluster_cond"].append(stacked)
    return out


def _depth_bin_edges() -> np.ndarray:
    """Edges of the depth (μm) bins for the heatmap y-axis."""
    return np.arange(DEPTH_LO_UM, DEPTH_HI_UM + DEPTH_BIN_UM, DEPTH_BIN_UM)


def _layer_boundary_depths(probe_data: list[dict]) -> dict[str, float]:
    """Per-layer mean depth and the data-derived inter-layer boundaries.

    Same approach as ``depth_variable_selection.py``: mean depth of all
    clusters labelled with each ``VISp*`` region across all 3 probes,
    midpoint between consecutive layer means is the boundary line drawn
    on the figure. Layers absent from the cohort drop out silently.
    """
    pooled: dict[str, list[float]] = {L: [] for L in LAYER_ORDER}
    for pd_ in probe_data:
        for region, depth in zip(pd_["cluster_regions"], pd_["cluster_depths"]):
            if region in pooled:
                pooled[region].append(float(depth))
    layer_means: dict[str, float] = {}
    for L, depths in pooled.items():
        if depths:
            layer_means[L] = float(np.mean(depths))
    # Boundaries between consecutive layers (in superficial → deep order).
    ordered = [L for L in LAYER_ORDER if L in layer_means]
    boundaries: list[tuple[str, str, float]] = []
    for i in range(len(ordered) - 1):
        a, b = ordered[i], ordered[i + 1]
        boundaries.append((a, b, 0.5 * (layer_means[a] + layer_means[b])))
    return {"layer_means": layer_means, "boundaries": boundaries}


def _per_cluster_baseline_stats(fr_by_cond: dict[str, np.ndarray]) -> tuple[float, float]:
    """Pooled-across-conditions mean and std of FR in the baseline window
    [BASELINE_LO_S, BASELINE_HI_S]. Used for the z-scored variant of
    the depth × time figure (Lohuis convention).
    """
    baseline_bin_mask = (
        (BIN_CENTRES_REL >= BASELINE_LO_S)
        & (BIN_CENTRES_REL <= BASELINE_HI_S)
    )
    pool_chunks = []
    for c in CONDITIONS:
        a = fr_by_cond.get(c)
        if a is None or a.size == 0:
            continue
        pool_chunks.append(a[:, baseline_bin_mask].ravel())
    if not pool_chunks:
        return 0.0, 1.0
    pool = np.concatenate(pool_chunks)
    pool = pool[np.isfinite(pool)]
    if pool.size < 10:
        return 0.0, 1.0
    m = float(pool.mean())
    s = float(pool.std(ddof=0))
    return m, (s if s > 0 else 1.0)


def _aggregate_subset(
    probe_data: list[dict],
    cluster_table: pd.DataFrame,
    subset: str,
    condition: str,
    depth_edges: np.ndarray,
    *,
    z_score: bool = False,
) -> dict:
    """For one (subset, condition), aggregate mean FR over (depth-bin,
    time-bin). Per-cluster trial-mean FR is averaged across clusters
    whose depth falls in the bin.

    ``z_score=False`` (default): raw FR (Hz). Heterogeneous cluster
    rates mean high-rate clusters dominate the cell mean.

    ``z_score=True``: per-cluster z-score against the pooled-across-
    conditions baseline window stats (Lohuis convention). Each cluster
    contributes (FR − baseline_mean) / baseline_std.
    """
    n_depth = depth_edges.size - 1
    sums = np.zeros((n_depth, N_BINS), dtype=np.float64)
    counts = np.zeros((n_depth, N_BINS), dtype=np.int64)
    clusters_per_depth = np.zeros(n_depth, dtype=np.int64)
    # Cluster-weighted accumulators for the line plot — every cluster
    # contributes equally per time bin, matching how fr_me_correlation
    # Fig 2 averages FR across the pooled (cluster × trial) rows. Cluster
    # weighting differs from the depth-weighted nanmean over the heatmap
    # matrix when the cohort is unevenly distributed across depth bins.
    line_sums = np.zeros(N_BINS, dtype=np.float64)
    line_counts = np.zeros(N_BINS, dtype=np.int64)
    total = 0
    for pd_ in probe_data:
        probe_id = pd_["probe_id"]
        for ci, cid in enumerate(pd_["cluster_ids"]):
            row = cluster_table[
                (cluster_table.probe_id == probe_id)
                & (cluster_table.cluster_id == cid)
            ]
            if row.empty:
                continue
            cluster_subset = row.iloc[0]["subset"]
            if subset == "All":
                if cluster_subset is None:
                    continue
            elif cluster_subset != subset:
                continue
            depth = float(pd_["cluster_depths"][ci])
            bin_i = int(np.digitize(depth, depth_edges) - 1)
            if bin_i < 0 or bin_i >= n_depth:
                continue
            fr_trials = pd_["fr_by_cluster_cond"][ci][condition]
            if fr_trials.shape[0] == 0:
                continue
            with np.errstate(invalid="ignore"):
                fr_mean = np.nanmean(fr_trials, axis=0)
            if z_score:
                b_mean, b_std = _per_cluster_baseline_stats(
                    pd_["fr_by_cluster_cond"][ci]
                )
                fr_mean = (fr_mean - b_mean) / b_std
            finite = np.isfinite(fr_mean)
            sums[bin_i, finite] += fr_mean[finite]
            counts[bin_i, finite] += 1
            clusters_per_depth[bin_i] += 1
            line_sums[finite] += fr_mean[finite]
            line_counts[finite] += 1
            total += 1
    means = np.full((n_depth, N_BINS), np.nan, dtype=np.float64)
    ok = counts > 0
    means[ok] = sums[ok] / counts[ok]
    line_mean = np.full(N_BINS, np.nan, dtype=np.float64)
    line_ok = line_counts > 0
    line_mean[line_ok] = line_sums[line_ok] / line_counts[line_ok]
    return {
        "means": means,
        "counts": counts,
        "clusters_per_depth": clusters_per_depth,
        "line_mean": line_mean,
        "total": total,
    }


# -- plotting ---------------------------------------------------------------
def _plot_grid(
    aggregated: dict[tuple[str, str], dict],
    layer_info: dict,
    depth_edges: np.ndarray,
    out_pdf: Path,
    *,
    z_score: bool = False,
    xlim: tuple[float, float] = (T_LO_S, T_HI_S),
) -> None:
    """7 rows (subsets) × 3 cols (conditions). Each cell = depth-averaged
    line plot above a depth × time heatmap.

    Line plot (top sub-panel): mean z-FR across all clusters in the
    subset, per time bin — kept on a SHARED y-axis per row so V / VT / T
    can be compared at a glance (motion-onset peaks are visible here
    even when the heatmap colourmap flattens them).

    Heatmap (bottom sub-panel): mean z-FR per (depth-bin, time-bin),
    with the per-row colourmap scale (vmax = 95th-percentile of |z| in
    that subset). Dotted lines mark data-derived layer boundaries.
    """
    from matplotlib.gridspec import GridSpec

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    n_rows = len(SUBSET_ORDER)
    n_cols = len(CONDITIONS)
    # Each subset row is split into a (line, heatmap) pair via a sub-GridSpec.
    fig = plt.figure(figsize=(16, 3.2 * n_rows + 0.6))
    outer = GridSpec(
        n_rows, n_cols + 1,
        figure=fig,
        width_ratios=[*([1.0] * n_cols), 0.04],
        wspace=0.08, hspace=0.45,
    )
    depth_lo, depth_hi = float(depth_edges[0]), float(depth_edges[-1])
    extent = (T_LO_S, T_HI_S, depth_hi, depth_lo)
    for ri, subset in enumerate(SUBSET_ORDER):
        # Per-row heatmap colour range. Raw → 0..95th percentile (FR ≥ 0
        # so colourmap is sequential 'magma'); z-score → ±95th
        # percentile of |z| (diverging 'RdBu_r' around 0).
        row_vals = np.concatenate(
            [aggregated[(subset, c)]["means"].ravel() for c in CONDITIONS]
        )
        row_vals = row_vals[np.isfinite(row_vals)]
        if z_score:
            vmax = (
                max(0.3, float(np.nanpercentile(np.abs(row_vals), 95)))
                if row_vals.size else 1.0
            )
        else:
            vmax = (
                max(1.0, float(np.nanpercentile(row_vals, 95)))
                if row_vals.size else 1.0
            )
        # Per-row line-plot y range. Cluster-weighted mean FR
        # (each cluster contributes equally per time bin), parallel to
        # the fr_me_correlation Fig 2 trace aggregation. Computed in
        # _aggregate_subset alongside the depth-binned heatmap matrix.
        line_means = {
            c: aggregated[(subset, c)]["line_mean"] for c in CONDITIONS
        }
        line_vals = np.concatenate([v[np.isfinite(v)] for v in line_means.values()])
        if line_vals.size:
            line_ymax = max(1.0, float(np.nanpercentile(line_vals, 99)) * 1.1)
        else:
            line_ymax = 1.0
        line_ax_share = None
        im = None
        for ci, condition in enumerate(CONDITIONS):
            # Two sub-rows inside this (subset, condition) cell.
            inner = outer[ri, ci].subgridspec(
                2, 1, height_ratios=(0.35, 1.0), hspace=0.05,
            )
            ax_line = fig.add_subplot(inner[0], sharey=line_ax_share)
            ax = fig.add_subplot(inner[1], sharex=ax_line)
            line_ax_share = line_ax_share or ax_line

            # ─── line plot ───
            ax_line.plot(BIN_CENTRES_REL, line_means[condition],
                          color="black", lw=1.0)
            ax_line.axvline(0, color="black", lw=0.5, ls="--")
            ax_line.axvline(0.1, color="black", lw=0.5)
            if z_score:
                ax_line.axhline(0, color="0.7", lw=0.4)
                ax_line.set_ylim(-line_ymax, line_ymax)
            else:
                ax_line.set_ylim(0, line_ymax)
            ax_line.set_xlim(*xlim)
            ax_line.tick_params(labelbottom=False, labelsize=6)
            if ci != 0:
                ax_line.tick_params(labelleft=False)
            if ri == 0:
                ax_line.set_title(condition, fontsize=9)

            # ─── heatmap ───
            cell = aggregated[(subset, condition)]
            if z_score:
                im = ax.imshow(
                    cell["means"],
                    extent=extent, aspect="auto",
                    origin="upper",
                    cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    interpolation="nearest",
                )
            else:
                im = ax.imshow(
                    cell["means"],
                    extent=extent, aspect="auto",
                    origin="upper",
                    cmap="magma", vmin=0, vmax=vmax,
                    interpolation="nearest",
                )
            ax.axvline(0, color="black", lw=0.5, ls="--")
            ax.axvline(0.1, color="black", lw=0.5)
            for a_label, b_label, depth in layer_info["boundaries"]:
                ax.axhline(depth, color="0.25", lw=0.4, ls=":")
            ax.set_xlim(*xlim)
            ax.set_ylim(depth_hi, depth_lo)
            if ri == n_rows - 1:
                ax.set_xlabel("Time from motion onset (s)", fontsize=8)
            if ci == 0:
                ax_line.set_ylabel(
                    "z-FR" if z_score else "FR (Hz)", fontsize=6,
                )
                ax.set_ylabel(
                    f"{subset}\n(n={cell['total']} clusters)\n\ndepth (μm)",
                    fontsize=7,
                )
                for L, mean_d in layer_info["layer_means"].items():
                    ax.text(
                        xlim[0] - 0.05 * (xlim[1] - xlim[0]),
                        mean_d, L,
                        fontsize=5, ha="right", va="center", color="0.25",
                    )
            else:
                ax.tick_params(labelleft=False)
            ax.tick_params(labelsize=6)
        # Per-row colourbar in the rightmost outer column.
        cax = fig.add_subplot(outer[ri, -1])
        cb = fig.colorbar(im, cax=cax)
        if z_score:
            cb.set_label(f"z-FR (±{vmax:.1f})", fontsize=6)
        else:
            cb.set_label(f"FR (Hz)  vmax = {vmax:.1f}", fontsize=6)
        cb.ax.tick_params(labelsize=6)
    title_prefix = (
        f"z-FR per cluster (baseline {BASELINE_LO_S:+.1f} … {BASELINE_HI_S:+.1f} s)"
        if z_score else "Raw FR (Hz)"
    )
    fig.suptitle(
        f"{title_prefix} per cortical depth × time, faceted by GLM Selected-model subset.  "
        "Line plots: cluster-weighted mean (shared y per row).  "
        f"{int(BIN_WIDTH_S*1000)} ms time bins, {int(DEPTH_BIN_UM)} μm depth bins; "
        "per-row colour scale.  Dotted: data-derived layer boundaries.",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.04, top=0.94)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--probes", nargs="*", default=list(PROBES),
        help="Probe stems (default: 3-camera cohort).",
    )
    args = parser.parse_args()

    cluster_table = _load_cluster_table()
    n_total = int(cluster_table["subset"].notna().sum())
    print(f"GLM cohort = {n_total} / {len(cluster_table)} clusters "
          f"(dropped {len(cluster_table) - n_total} V-only / none)", flush=True)
    print(cluster_table["subset"].value_counts(dropna=False).to_string())

    lookup = StimulusLookup(SEQ_PATH, FOL_PATH)
    probe_data = []
    for stem in args.probes:
        print(f"[{stem}] binning …", flush=True)
        probe_data.append(_bin_one_probe(stem, lookup))

    layer_info = _layer_boundary_depths(probe_data)
    print("layer means (μm):", layer_info["layer_means"])
    print("boundaries (μm):", [(a, b, round(d, 1)) for a, b, d in layer_info["boundaries"]])

    depth_edges = _depth_bin_edges()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Three figures: raw FR full ±4 s, z-FR full ±4 s, and z-FR zoomed
    # to ±1 s around motion onset (Laura's 2026-05-28 ask — onset window
    # only, easier to compare layer structure within the response peak).
    variants = [
        ("raw_fr", False, (T_LO_S, T_HI_S)),
        ("z_score", True, (T_LO_S, T_HI_S)),
        ("z_score_zoom", True, (-1.0, 1.0)),
    ]
    for tag, z_flag, xlim in variants:
        aggregated: dict[tuple[str, str], dict] = {}
        for subset in SUBSET_ORDER:
            for condition in CONDITIONS:
                aggregated[(subset, condition)] = _aggregate_subset(
                    probe_data, cluster_table, subset, condition, depth_edges,
                    z_score=z_flag,
                )
        _plot_grid(
            aggregated, layer_info, depth_edges,
            OUT_DIR / f"depth_x_time_by_subset_{tag}.pdf",
            z_score=z_flag, xlim=xlim,
        )


if __name__ == "__main__":
    main()
