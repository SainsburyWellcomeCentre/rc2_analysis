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
from rc2_formatted_data_reader import masks
from rc2_glm.masks_helpers import acceleration_from_velocity


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
    """Return per-(cluster, condition) z-FR arrays of shape (n_trials_in_cond, N_BINS).

    Z-score is per-cluster against pooled (trial × baseline-bin) FR
    statistics. NaN where a trial doesn't cover a bin (short
    pre-motion-onset interval, etc.).
    """
    mat = MAT_DIR / f"{probe_stem}.mat"
    out = {
        "probe_id": probe_stem,
        "cluster_ids": [],
        "cluster_regions": [],
        "cluster_depths": [],
        "fr_by_cluster_cond": [],  # list per cluster: dict[condition] -> (n_trials, N_BINS) array
        "z_baseline_mean": [],
        "z_baseline_std": [],
    }
    with FormattedDataReader(mat) as r:
        fs = r.fs
        n_clusters = r.n_clusters
        n_trials = r.n_trials
        cluster_ids = r.cluster_ids()
        regions = r.cluster_regions()

        # Per trial: motion-onset time + condition
        trial_motion_time = np.full(n_trials, np.nan, dtype=np.float64)
        trial_cond = np.full(n_trials, None, dtype=object)
        for ti in range(n_trials):
            tid = r.trial_id(ti)
            protocol = r.trial_protocol(ti)
            if protocol not in ("StageOnly", "ReplayOnly"):
                continue
            condition = r.trial_condition(ti)
            if condition not in CONDITIONS:
                continue
            # excluded?
            params = lookup.stimulus_params(tid)
            if params.get("excluded", False):
                continue
            # Compute motion onset via the same logic as fr_me_corr:
            # trial_velocity → treadmill_motion_mask → analysis_mask → first True.
            v = r.trial_velocity(ti)
            accel = acceleration_from_velocity(v)
            m_mask = masks.treadmill_motion_mask(
                v, accel, fs, vel_thresh=1.0, acc_thresh=0.5, min_dur=0.2,
            )
            analysis_mask = r.trial_analysis_mask(ti)
            m_mask = m_mask & analysis_mask
            motion_idx = np.flatnonzero(m_mask)
            if motion_idx.size == 0:
                continue
            s_idx, _ = r.trial_bounds(ti)
            t_onset = float((s_idx + int(motion_idx[0])) / fs)
            trial_motion_time[ti] = t_onset
            trial_cond[ti] = condition

        # Per-cluster binning
        for ci in range(n_clusters):
            cid = int(cluster_ids[ci])
            region = regions[ci]
            depth = r.cluster_depth(ci)
            spikes = r.spike_times(ci)
            per_cond_fr: dict[str, list[np.ndarray]] = {c: [] for c in CONDITIONS}
            for ti in range(n_trials):
                t_onset = trial_motion_time[ti]
                cond = trial_cond[ti]
                if not np.isfinite(t_onset) or cond is None:
                    continue
                edges = t_onset + BIN_EDGES_REL
                counts, _ = np.histogram(spikes, bins=edges)
                fr = counts.astype(np.float64) / BIN_WIDTH_S  # Hz
                per_cond_fr[cond].append(fr)
            # Stack per condition
            stacked: dict[str, np.ndarray] = {}
            for c, lst in per_cond_fr.items():
                stacked[c] = (
                    np.vstack(lst) if lst else np.zeros((0, N_BINS), dtype=np.float64)
                )
            # Baseline stats — pooled across conditions, across (trial × baseline bin).
            baseline_mask = (
                (BIN_CENTRES_REL >= BASELINE_LO_S)
                & (BIN_CENTRES_REL <= BASELINE_HI_S)
            )
            pool = np.concatenate(
                [stacked[c][:, baseline_mask].ravel() for c in CONDITIONS]
                if any(stacked[c].size for c in CONDITIONS) else []
            ) if any(stacked[c].size for c in CONDITIONS) else np.zeros(0)
            if pool.size >= 10 and pool.std(ddof=0) > 0:
                b_mean = float(pool.mean())
                b_std = float(pool.std(ddof=0))
            else:
                b_mean, b_std = 0.0, 1.0
            out["cluster_ids"].append(cid)
            out["cluster_regions"].append(region)
            out["cluster_depths"].append(depth)
            out["fr_by_cluster_cond"].append(stacked)
            out["z_baseline_mean"].append(b_mean)
            out["z_baseline_std"].append(b_std)
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


def _aggregate_subset(
    probe_data: list[dict],
    cluster_table: pd.DataFrame,
    subset: str,
    condition: str,
    depth_edges: np.ndarray,
) -> dict:
    """For one (subset, condition), aggregate mean z-FR over (depth-bin, time-bin).

    Returns a dict with a 2D ``means`` array shaped (n_depth_bins,
    N_BINS), a parallel ``counts`` array (number of contributing
    clusters per cell), the total cluster count, and the per-bin
    cluster count along the depth axis (for caveat annotations).
    """
    n_depth = depth_edges.size - 1
    sums = np.zeros((n_depth, N_BINS), dtype=np.float64)
    counts = np.zeros((n_depth, N_BINS), dtype=np.int64)
    clusters_per_depth = np.zeros(n_depth, dtype=np.int64)
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
            b_mean = pd_["z_baseline_mean"][ci]
            b_std = pd_["z_baseline_std"][ci]
            fr_trials = pd_["fr_by_cluster_cond"][ci][condition]
            if fr_trials.shape[0] == 0:
                continue
            z = (fr_trials - b_mean) / b_std
            with np.errstate(invalid="ignore"):
                z_mean = np.nanmean(z, axis=0)
            finite = np.isfinite(z_mean)
            sums[bin_i, finite] += z_mean[finite]
            counts[bin_i, finite] += 1
            clusters_per_depth[bin_i] += 1
            total += 1
    means = np.full((n_depth, N_BINS), np.nan, dtype=np.float64)
    ok = counts > 0
    means[ok] = sums[ok] / counts[ok]
    return {
        "means": means,
        "counts": counts,
        "clusters_per_depth": clusters_per_depth,
        "total": total,
    }


# -- plotting ---------------------------------------------------------------
def _plot_grid(
    aggregated: dict[tuple[str, str], dict],
    layer_info: dict,
    depth_edges: np.ndarray,
    out_pdf: Path,
) -> None:
    """7 rows (subsets) × 3 cols (conditions) + per-row colourbar.

    Y axis: continuous depth in μm (50 μm bins). Dashed horizontal lines
    mark the data-derived layer boundaries from
    ``_layer_boundary_depths``. Each row uses its OWN colourmap scale
    (vmax = nanpercentile(|z|, 95) of that subset) so weaker subsets
    are not washed out by stronger ones — Laura's 2026-05-28 ask.
    """
    from matplotlib.gridspec import GridSpec

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    n_rows = len(SUBSET_ORDER)
    n_cols = len(CONDITIONS)
    # Width per condition column + extra for colourbar gutter + layer label gutter.
    fig = plt.figure(figsize=(16, 2.2 * n_rows + 0.6))
    gs = GridSpec(
        n_rows, n_cols + 1,
        figure=fig,
        width_ratios=[*([1.0] * n_cols), 0.04],
        wspace=0.08, hspace=0.35,
    )
    depth_lo, depth_hi = float(depth_edges[0]), float(depth_edges[-1])
    extent = (T_LO_S, T_HI_S, depth_hi, depth_lo)  # origin upper → invert y
    for ri, subset in enumerate(SUBSET_ORDER):
        # Per-row colour range: 95th-percentile of |z| across the three
        # conditions of this subset. Skip NaNs. Floor at 0.3 so a totally
        # flat row still has a tiny visible scale.
        row_vals = np.concatenate(
            [aggregated[(subset, c)]["means"].ravel() for c in CONDITIONS]
        )
        row_vals = row_vals[np.isfinite(row_vals)]
        vmax = (
            max(0.3, float(np.nanpercentile(np.abs(row_vals), 95)))
            if row_vals.size else 1.0
        )
        im = None
        for ci, condition in enumerate(CONDITIONS):
            ax = fig.add_subplot(gs[ri, ci])
            cell = aggregated[(subset, condition)]
            im = ax.imshow(
                cell["means"],
                extent=extent, aspect="auto",
                origin="upper",
                cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                interpolation="nearest",
            )
            ax.axvline(0, color="black", lw=0.5, ls="--")
            ax.axvline(0.1, color="black", lw=0.5)
            # Layer boundary lines + layer-name labels on the leftmost col.
            for a_label, b_label, depth in layer_info["boundaries"]:
                ax.axhline(depth, color="0.25", lw=0.4, ls=":")
            ax.set_xlim(T_LO_S, T_HI_S)
            ax.set_ylim(depth_hi, depth_lo)
            if ri == 0:
                ax.set_title(condition, fontsize=9)
            if ri == n_rows - 1:
                ax.set_xlabel("Time from motion onset (s)", fontsize=8)
            if ci == 0:
                ax.set_ylabel(
                    f"{subset}\n(n={cell['total']} clusters)\n\ndepth (μm)",
                    fontsize=7,
                )
                # Annotate layer names against their mean depth in the LH gutter.
                for L, mean_d in layer_info["layer_means"].items():
                    ax.text(
                        T_LO_S - 0.05 * (T_HI_S - T_LO_S),
                        mean_d, L,
                        fontsize=5, ha="right", va="center", color="0.25",
                    )
            else:
                ax.tick_params(labelleft=False)
            ax.tick_params(labelsize=6)
        # Per-row colourbar in the rightmost gridspec column.
        cax = fig.add_subplot(gs[ri, -1])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(f"z-FR (±{vmax:.1f})", fontsize=6)
        cb.ax.tick_params(labelsize=6)
    fig.suptitle(
        "z-FR per cortical depth × time, faceted by GLM Selected-model subset.  "
        f"z-score per cluster against baseline ({BASELINE_LO_S:+.1f} s, {BASELINE_HI_S:+.1f} s); "
        f"{int(BIN_WIDTH_S*1000)} ms time bins, {int(DEPTH_BIN_UM)} μm depth bins; "
        "per-row colour scale (vmax = 95th percentile of |z| in that subset).  "
        "Dotted lines: data-derived layer boundaries.",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.05, top=0.94)
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
    aggregated: dict[tuple[str, str], dict] = {}
    for subset in SUBSET_ORDER:
        for condition in CONDITIONS:
            aggregated[(subset, condition)] = _aggregate_subset(
                probe_data, cluster_table, subset, condition, depth_edges,
            )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _plot_grid(aggregated, layer_info, depth_edges, OUT_DIR / "depth_x_time_by_subset.pdf")


if __name__ == "__main__":
    main()
