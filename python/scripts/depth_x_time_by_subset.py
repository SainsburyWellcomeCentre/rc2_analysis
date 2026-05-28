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

# Time grid relative to motion onset. 100 ms bins matches the production GLM
# binning convention; ±2 s window covers a useful chunk of both the
# stationary prelude and the motion phase.
BIN_WIDTH_S = 0.10
T_LO_S = -2.0
T_HI_S = 2.0
N_BINS = int(round((T_HI_S - T_LO_S) / BIN_WIDTH_S))   # 40
BIN_EDGES_REL = T_LO_S + BIN_WIDTH_S * np.arange(N_BINS + 1)
BIN_CENTRES_REL = 0.5 * (BIN_EDGES_REL[:-1] + BIN_EDGES_REL[1:])

# Baseline window for per-cluster z-score normalisation (Lohuis convention).
BASELINE_LO_S = -1.0
BASELINE_HI_S = -0.1

# Layer order, superficial → deep. The 80-cluster cohort spans L2/3 down,
# so VISp1 is omitted.
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


def _aggregate_subset(
    probe_data: list[dict],
    cluster_table: pd.DataFrame,
    subset: str,
    condition: str,
) -> dict[str, np.ndarray]:
    """For one (subset, condition), aggregate mean z-FR over (depth-layer, time-bin).

    Returns a dict layer → (N_BINS,) array of mean z-FR (NaN where no
    cluster contributes), plus a layer-wise cluster-count dict.
    """
    sums = defaultdict(lambda: np.zeros(N_BINS, dtype=np.float64))
    counts = defaultdict(lambda: np.zeros(N_BINS, dtype=np.int64))
    cluster_counts: dict[str, int] = defaultdict(int)
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
            region = pd_["cluster_regions"][ci]
            if region not in LAYER_ORDER:
                continue
            b_mean = pd_["z_baseline_mean"][ci]
            b_std = pd_["z_baseline_std"][ci]
            fr_trials = pd_["fr_by_cluster_cond"][ci][condition]
            if fr_trials.shape[0] == 0:
                continue
            z = (fr_trials - b_mean) / b_std
            # Trial-mean per bin for this cluster
            with np.errstate(invalid="ignore"):
                z_mean = np.nanmean(z, axis=0)
            finite = np.isfinite(z_mean)
            sums[region][finite] += z_mean[finite]
            counts[region][finite] += 1
            cluster_counts[region] += 1
    means: dict[str, np.ndarray] = {}
    for layer in LAYER_ORDER:
        c = counts[layer]
        s = sums[layer]
        out = np.full(N_BINS, np.nan, dtype=np.float64)
        ok = c > 0
        out[ok] = s[ok] / c[ok]
        means[layer] = out
    return {"means": means, "cluster_counts": dict(cluster_counts)}


# -- plotting ---------------------------------------------------------------
def _plot_grid(
    aggregated: dict[tuple[str, str], dict],
    out_pdf: Path,
    *,
    vlim: float = 1.5,
) -> None:
    """7 rows (subsets) × 3 cols (conditions). Each cell = layer × time heatmap."""
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    n_rows = len(SUBSET_ORDER)
    n_cols = len(CONDITIONS)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(11, 1.6 * n_rows + 0.6), sharex=True,
    )
    extent = (T_LO_S, T_HI_S, len(LAYER_ORDER) - 0.5, -0.5)
    for ri, subset in enumerate(SUBSET_ORDER):
        for ci, condition in enumerate(CONDITIONS):
            ax = axes[ri, ci]
            cell = aggregated[(subset, condition)]
            mat = np.vstack([cell["means"][L] for L in LAYER_ORDER])
            im = ax.imshow(
                mat,
                extent=extent, aspect="auto",
                origin="upper",
                cmap="RdBu_r", vmin=-vlim, vmax=vlim,
                interpolation="nearest",
            )
            ax.axvline(0, color="black", lw=0.5, ls="--")
            ax.axvline(0.1, color="black", lw=0.5)
            ax.set_xlim(T_LO_S, T_HI_S)
            if ri == 0:
                ax.set_title(condition, fontsize=9)
            if ri == n_rows - 1:
                ax.set_xlabel("Time from motion onset (s)", fontsize=8)
            if ci == 0:
                counts = cell["cluster_counts"]
                total = sum(counts.get(L, 0) for L in LAYER_ORDER)
                ax.set_ylabel(f"{subset}\n(n={total})", fontsize=7)
                ax.set_yticks(np.arange(len(LAYER_ORDER)))
                ax.set_yticklabels(LAYER_ORDER, fontsize=6)
            else:
                ax.set_yticks([])
            ax.tick_params(labelsize=6)
            # Annotate per-layer cluster counts on the right edge
            if ci == n_cols - 1:
                counts = cell["cluster_counts"]
                for li, L in enumerate(LAYER_ORDER):
                    n = counts.get(L, 0)
                    ax.text(
                        T_HI_S + 0.05, li, f"n={n}",
                        fontsize=5, va="center", ha="left", color="0.3",
                    )
    fig.suptitle(
        "z-FR per cortical layer × time, faceted by GLM Selected-model subset.  "
        "z-score per cluster against baseline (−1 s, −0.1 s); 100 ms bins; "
        "diverging colour ±1.5 z.",
        fontsize=9,
    )
    fig.subplots_adjust(
        left=0.10, right=0.92, bottom=0.06, top=0.93, hspace=0.25, wspace=0.10,
    )
    # Single colourbar on the right.
    cbar_ax = fig.add_axes((0.94, 0.10, 0.012, 0.78))
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("z-FR", fontsize=8)
    cb.ax.tick_params(labelsize=7)
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

    aggregated: dict[tuple[str, str], dict] = {}
    for subset in SUBSET_ORDER:
        for condition in CONDITIONS:
            aggregated[(subset, condition)] = _aggregate_subset(
                probe_data, cluster_table, subset, condition,
            )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _plot_grid(aggregated, OUT_DIR / "depth_x_time_by_subset.pdf")


if __name__ == "__main__":
    main()
