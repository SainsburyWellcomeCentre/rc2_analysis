"""Lohuis Fig 2F-style correlation of firing rate to face motion energy.

Per the 2026-05-28 stated plan (CLAUDE.md "stated-findings-before-action"
trail): across-trial Pearson r at each time bin between per-cluster
firing rate and face motion energy, computed on stationary→motion
aligned trials. Two figures:

- **Figure 1** (one per cluster, multi-page PDF per probe):
  3 cols (V / VT / T_Vstatic) × 2 rows (passive replay profile 1 / 2).
  Each cell: thick line = across-trial Pearson r at each bin, shaded =
  bootstrap 95 % CI over trials. Bins where < 70 % of the cell's trials
  cover are greyed.

- **Figure 2** (single PDF, all 3 probes pooled):
  cols = V / VT / T_Vstatic, rows = animals (243, 244, 466), plus a
  bottom row pooling all 3 animals. Per cell: pool trials across both
  profiles and all clusters of that animal; z-score FR per cluster, ME
  per probe; pooled Pearson r per bin + bootstrap 95 % CI.

Conventions (from the stated plan):
- 30 ms bins (matched to 60 Hz video).
- Spike trains smoothed with a causal half-Gaussian, 50 ms SD (Lohuis).
- ME re-binned at 30 ms from raw camera0 (face) signal, z-scored over
  pooled motion-row support per probe.
- Motion onset = first True in ``trial.motion_mask`` (protocol-aware
  via ``rc2_glm.io._load_trial``); stationary block = last contiguous
  stationary run before motion onset (matches existing
  ``time_binning._bin_stationary_prelude``).
- **Cluster cohort = the production GLM set** (~80 clusters across
  3 probes, from ``glm_model_comparison.csv``), NOT raw VISp /
  Kilosort output. Motion-clouds per-cluster analyses always run on
  this prefiltered cohort. ``--me-filter`` narrows further to the
  ME-selected subset (~55 clusters); ``--no-cohort-filter`` opts back
  into all VISp (debug only).
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.io import ProbeData, load_probe_data


# --- Constants -------------------------------------------------------------
PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
BIN_WIDTH = 0.030  # 30 ms
FR_SMOOTH_SD_S = 0.050  # 50 ms causal half-Gaussian
FR_SMOOTH_TRUNC_SD = 3.0  # truncate kernel at 3 SD
TRIAL_COVERAGE_THRESH = 0.70  # grey bins where < 70 % of trials cover
BOOTSTRAP_N = 500
BOOTSTRAP_RNG_SEED = 20260528
# Fixed x-axis limits so per-cluster panels are homogeneous and Figure 2
# panels match Figure 1 frame-for-frame. Per Laura's 2026-05-28 request.
XLIM_S = (-4.0, 4.0)
# Y-axis ranges for the Pearson-r panels:
# - Figure 2 (pooled) sits at r ~ 0.05–0.15 and needs a tight axis to
#   show structure. ±0.2 per Laura's 2026-05-28 follow-up.
# - Figure 1 (per-cluster) has 36-trial per-bin SD ≈ 0.17; a tight axis
#   visibly clips the natural scatter, so ±1.0 stays honest.
YLIM_R_POOLED = (-0.2, 0.2)
YLIM_R_PER_CLUSTER = (-1.0, 1.0)
# Paradigm-specific landmark to draw on every panel (solid vertical bar).
T100_MS_S = 0.10
CONDITIONS = ("V", "VT", "T_Vstatic")
PROFILES = (1, 2)

DATA_ROOT = Path("~/local_data/motion_clouds").expanduser()
MAT_DIR = DATA_ROOT / "formatted_data_3probe"
SEQ_PATH = DATA_ROOT / "motion_cloud_sequence_250414.mat"
FOL_PATH = DATA_ROOT / "image_folders.mat"
OUT_ROOT = DATA_ROOT / "figures" / "glm" / "exploration" / "fr_me_corr"
OUT_ROOT_ME_FILTER = DATA_ROOT / "figures" / "glm" / "exploration" / "fr_me_corr_me_selected"
GLM_CSV_DEFAULT = (
    DATA_ROOT / "figures" / "glm" / "current_with_ME_3_probes"
    / "glm_model_comparison.csv"
)
CACHE_ROOT = DATA_ROOT / "fr_me_corr_cache"


# --- Per-trial binning -----------------------------------------------------
@dataclass
class TrialBins:
    """Per-trial 30 ms bin grid centred on motion onset."""

    trial_id: int
    condition: str
    profile_id: int
    # Integer bin offsets relative to motion onset (centre of bin 0 = t=0).
    bin_offsets: np.ndarray   # shape (n_bins,) dtype int64; e.g. [-78, ..., 118]
    # Bin edges in absolute probe_t seconds (n_bins + 1).
    bin_edges: np.ndarray     # shape (n_bins+1,)
    me_raw: np.ndarray        # shape (n_bins,) face ME bin-mean, NaN where no camera samples


def _causal_half_gaussian_kernel(sd_s: float, bin_width: float, trunc_sd: float) -> np.ndarray:
    """Causal half-Gaussian kernel for convolution at ``bin_width`` resolution.

    Only past samples contribute. Lohuis 2024 methods: "convolved with a
    causal half-Gaussian window with 50-ms standard deviation".
    """
    n = max(1, int(np.ceil(trunc_sd * sd_s / bin_width)))
    x = np.arange(n) * bin_width  # 0, dt, 2dt, ...
    k = np.exp(-0.5 * (x / sd_s) ** 2)
    k /= k.sum()
    return k


def _bin_continuous_mean(
    values: np.ndarray, sample_t: np.ndarray, bin_edges: np.ndarray
) -> np.ndarray:
    """Bin-mean of ``values`` (sampled at ``sample_t``) onto ``bin_edges``.

    Returns NaN for bins with no contributing samples. Mirrors
    ``rc2_glm.time_binning._bin_continuous_to_edges`` (kept inline so the
    script is self-contained and so we can adjust the NaN handling later
    if needed).
    """
    n_bins = bin_edges.size - 1
    if n_bins <= 0:
        return np.zeros(0, dtype=np.float64)
    sample_bin = np.digitize(sample_t, bin_edges, right=False) - 1
    valid = (sample_bin >= 0) & (sample_bin < n_bins)
    bins_vec = sample_bin[valid]
    vals_vec = values[valid].astype(np.float64)
    n_per_bin = np.bincount(bins_vec, minlength=n_bins)
    sum_per_bin = np.bincount(bins_vec, weights=vals_vec, minlength=n_bins)
    out = np.full(n_bins, np.nan, dtype=np.float64)
    has_data = n_per_bin > 0
    out[has_data] = sum_per_bin[has_data] / n_per_bin[has_data]
    return out


def _build_trial_bins(trial, probe: ProbeData) -> TrialBins | None:
    """Build the 30 ms bin grid for one trial, anchored at motion onset.

    Returns None when the trial has no motion or no usable stationary
    prelude. The stationary block is the *last* contiguous stationary
    run before motion onset (matches ``_bin_stationary_prelude``).
    """
    if trial.excluded:
        return None

    motion_idx = np.flatnonzero(trial.motion_mask)
    if motion_idx.size == 0:
        return None
    motion_start_idx = int(motion_idx[0])
    motion_end_idx = int(motion_idx[-1])
    t_motion_start = float(trial.probe_t[motion_start_idx])
    t_motion_end = float(trial.probe_t[motion_end_idx])

    # Last contiguous stationary run before motion onset.
    pre = trial.stationary_mask.copy()
    pre[motion_start_idx:] = False
    stat_idx = np.flatnonzero(pre)
    if stat_idx.size == 0:
        return None
    gaps = np.flatnonzero(np.diff(stat_idx) > 1)
    if gaps.size > 0:
        stat_start_idx = int(stat_idx[gaps[-1] + 1])
    else:
        stat_start_idx = int(stat_idx[0])
    t_stat_start = float(trial.probe_t[stat_start_idx])

    # Integer bin offsets from motion onset. Bin k spans
    # [t_motion_start + (k - 0.5) * dt, t_motion_start + (k + 0.5) * dt],
    # i.e. centre at t_motion_start + k * dt.
    dt = BIN_WIDTH
    k_lo = int(np.ceil((t_stat_start - t_motion_start) / dt - 0.5))
    k_hi = int(np.floor((t_motion_end - t_motion_start) / dt + 0.5))
    if k_lo >= k_hi:
        return None

    bin_offsets = np.arange(k_lo, k_hi + 1, dtype=np.int64)
    n_bins = bin_offsets.size
    bin_edges = t_motion_start + (bin_offsets.astype(np.float64) - 0.5) * dt
    bin_edges = np.concatenate([bin_edges, [bin_edges[-1] + dt]])

    # ME signal (camera0) binned to these edges.
    if probe.camera0 is None or probe.camera_t is None:
        me_raw = np.full(n_bins, np.nan)
    else:
        cam_vals = np.asarray(probe.camera0, dtype=np.float64).ravel()
        cam_t = np.asarray(probe.camera_t, dtype=np.float64).ravel()
        n = min(cam_vals.size, cam_t.size)
        me_raw = _bin_continuous_mean(cam_vals[:n], cam_t[:n], bin_edges)

    return TrialBins(
        trial_id=trial.trial_id,
        condition=trial.condition,
        profile_id=trial.profile_id,
        bin_offsets=bin_offsets,
        bin_edges=bin_edges,
        me_raw=me_raw,
    )


def _bin_cluster_fr(
    cluster_spike_times: np.ndarray,
    bin_edges: np.ndarray,
    kernel: np.ndarray,
) -> np.ndarray:
    """Histogram spikes into ``bin_edges`` (Hz), then causal-smooth."""
    counts, _ = np.histogram(cluster_spike_times, bins=bin_edges)
    rate = counts.astype(np.float64) / BIN_WIDTH  # Hz
    # Causal convolution: y[t] = sum_{i=0..K-1} k[i] * x[t-i]. np.convolve
    # with mode='full' then truncate to the original length, taking the
    # first n_bins samples (so x[0] only sees the present, x[1] sees x[1]
    # and x[0], etc.).
    if kernel.size > 1:
        rate = np.convolve(rate, kernel, mode="full")[: bin_edges.size - 1]
    return rate


# --- Per-probe aggregation -------------------------------------------------
@dataclass
class ProbeBins:
    probe_id: str
    cluster_ids: np.ndarray            # (n_clusters,)
    trials: list[TrialBins]            # one per usable trial
    # Per-cluster smoothed FR per trial, parallel to ``trials``.
    # fr[cluster_idx][trial_idx] = np.ndarray of length matching trial.bin_offsets
    fr: list[list[np.ndarray]]
    # Per-probe ME z-score stats (computed over motion bins, pooled).
    me_mean: float
    me_std: float
    # Per-cluster z-score stats (computed over all bins).
    fr_mean: np.ndarray                # (n_clusters,)
    fr_std: np.ndarray                 # (n_clusters,)


def _aggregate_probe(probe: ProbeData) -> ProbeBins:
    kernel = _causal_half_gaussian_kernel(FR_SMOOTH_SD_S, BIN_WIDTH, FR_SMOOTH_TRUNC_SD)
    trial_bins: list[TrialBins] = []
    for trial in probe.trials:
        tb = _build_trial_bins(trial, probe)
        if tb is not None:
            trial_bins.append(tb)

    n_clusters = len(probe.clusters)
    cluster_ids = np.array([c.cluster_id for c in probe.clusters], dtype=np.int64)

    # Per-cluster smoothed FR per trial.
    fr: list[list[np.ndarray]] = []
    for ci, cluster in enumerate(probe.clusters):
        spikes = cluster.spike_times
        rates_per_trial = [
            _bin_cluster_fr(spikes, tb.bin_edges, kernel) for tb in trial_bins
        ]
        fr.append(rates_per_trial)
        if (ci + 1) % 50 == 0 or (ci + 1) == n_clusters:
            print(f"  [{probe.probe_id}] binned cluster {ci+1}/{n_clusters}", flush=True)

    # ME z-score stats over motion bins pooled across all trials.
    me_motion_values: list[np.ndarray] = []
    for tb in trial_bins:
        motion_mask_bins = tb.bin_offsets >= 0
        vals = tb.me_raw[motion_mask_bins]
        me_motion_values.append(vals[np.isfinite(vals)])
    me_all = np.concatenate(me_motion_values) if me_motion_values else np.array([])
    if me_all.size >= 10:
        me_mean = float(me_all.mean())
        me_std = float(me_all.std(ddof=0)) or 1.0
    else:
        me_mean, me_std = 0.0, 1.0

    # Per-cluster FR z-score stats over all of that cluster's bins.
    fr_mean = np.zeros(n_clusters, dtype=np.float64)
    fr_std = np.ones(n_clusters, dtype=np.float64)
    for ci in range(n_clusters):
        all_rates = np.concatenate(fr[ci]) if fr[ci] else np.array([])
        if all_rates.size > 10:
            fr_mean[ci] = float(all_rates.mean())
            s = float(all_rates.std(ddof=0))
            fr_std[ci] = s if s > 0 else 1.0

    return ProbeBins(
        probe_id=probe.probe_id,
        cluster_ids=cluster_ids,
        trials=trial_bins,
        fr=fr,
        me_mean=me_mean,
        me_std=me_std,
        fr_mean=fr_mean,
        fr_std=fr_std,
    )


# --- Cell-level correlation ------------------------------------------------
def _stack_cell_arrays(
    pb: ProbeBins,
    cluster_idx: int,
    cell_trial_indices: list[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stack FR and ME into (n_trials, n_bins) arrays aligned on motion onset.

    Bin offsets are unified across the cell as the min/max of all trial
    bin_offsets. Bins where a given trial has no data are NaN.

    Returns ``(fr, me, bin_offsets)``.
    """
    if not cell_trial_indices:
        return (
            np.zeros((0, 0)),
            np.zeros((0, 0)),
            np.zeros(0, dtype=np.int64),
        )
    lo = min(int(pb.trials[ti].bin_offsets[0]) for ti in cell_trial_indices)
    hi = max(int(pb.trials[ti].bin_offsets[-1]) for ti in cell_trial_indices)
    bin_offsets = np.arange(lo, hi + 1, dtype=np.int64)
    n_bins = bin_offsets.size
    n_trials = len(cell_trial_indices)

    fr = np.full((n_trials, n_bins), np.nan, dtype=np.float64)
    me = np.full((n_trials, n_bins), np.nan, dtype=np.float64)
    for row, ti in enumerate(cell_trial_indices):
        tb = pb.trials[ti]
        start = int(tb.bin_offsets[0]) - lo
        end = start + tb.bin_offsets.size
        fr[row, start:end] = pb.fr[cluster_idx][ti]
        me[row, start:end] = tb.me_raw
    return fr, me, bin_offsets


def _pearson_r_per_bin(fr: np.ndarray, me: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pearson r across rows (trials) at each column (bin), NaN-aware.

    Returns ``(r, n_used)`` where ``n_used[bin]`` is the number of rows
    used at that bin (both fr and me finite).
    """
    if fr.size == 0:
        return np.zeros(0), np.zeros(0, dtype=np.int64)
    valid = np.isfinite(fr) & np.isfinite(me)
    n_used = valid.sum(axis=0)
    # Replace invalid with 0 so vectorised sums are safe; then mask out
    # bins with insufficient trials.
    fr_v = np.where(valid, fr, 0.0)
    me_v = np.where(valid, me, 0.0)
    n = n_used.astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        sum_fr = fr_v.sum(axis=0)
        sum_me = me_v.sum(axis=0)
        sum_fr2 = (fr_v * fr_v).sum(axis=0)
        sum_me2 = (me_v * me_v).sum(axis=0)
        sum_frme = (fr_v * me_v).sum(axis=0)
        num = n * sum_frme - sum_fr * sum_me
        den = np.sqrt(
            (n * sum_fr2 - sum_fr * sum_fr) * (n * sum_me2 - sum_me * sum_me)
        )
        r = np.where((den > 0) & (n_used >= 3), num / den, np.nan)
    return r, n_used


def _bootstrap_ci(
    fr: np.ndarray,
    me: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
    group_ids: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Trial-resampling bootstrap 95 % CI of the per-bin Pearson r.

    For Figure 1 (one cluster per cell), rows of ``fr`` / ``me`` are
    already one-per-trial — call without ``group_ids`` and the bootstrap
    resamples rows directly.

    For Figure 2 (pooled across clusters × trials), the same trial is
    replicated across all clusters of that animal. The independent unit
    is the trial, not the (cluster, trial) row. Pass ``group_ids`` of
    shape ``(n_rows,)`` carrying a trial identifier per row; the
    bootstrap then resamples *unique trial IDs* and gathers all rows
    sharing each sampled ID. This is the cluster-aware variant; without
    it the CI would be artificially tight by a factor of ~sqrt(n_clusters).
    """
    n_rows = fr.shape[0]
    if n_rows < 3:
        nb = fr.shape[1]
        return np.full(nb, np.nan), np.full(nb, np.nan)
    boots = np.empty((n_boot, fr.shape[1]), dtype=np.float64)

    if group_ids is None:
        for b in range(n_boot):
            idx = rng.integers(0, n_rows, size=n_rows)
            r_b, _ = _pearson_r_per_bin(fr[idx], me[idx])
            boots[b] = r_b
    else:
        # Build group → row-index list once; resample groups, gather rows.
        group_ids = np.asarray(group_ids)
        unique_groups, inverse = np.unique(group_ids, return_inverse=True)
        n_groups = unique_groups.size
        # rows_by_group[g] = np.ndarray of row indices for group g
        rows_by_group: list[np.ndarray] = [
            np.flatnonzero(inverse == g) for g in range(n_groups)
        ]
        if n_groups < 3:
            nb = fr.shape[1]
            return np.full(nb, np.nan), np.full(nb, np.nan)
        for b in range(n_boot):
            sampled = rng.integers(0, n_groups, size=n_groups)
            # Concatenate the row indices of all sampled groups (with
            # repeats — that's how the bootstrap is supposed to behave).
            idx = np.concatenate([rows_by_group[g] for g in sampled])
            r_b, _ = _pearson_r_per_bin(fr[idx], me[idx])
            boots[b] = r_b
    lo = np.nanpercentile(boots, 2.5, axis=0)
    hi = np.nanpercentile(boots, 97.5, axis=0)
    return lo, hi


# --- Figure 1: per-cluster grid --------------------------------------------
def _cell_trial_indices(pb: ProbeBins, condition: str, profile: int) -> list[int]:
    return [
        i for i, tb in enumerate(pb.trials)
        if tb.condition == condition and tb.profile_id == profile
    ]


def _mean_trace_r(fr: np.ndarray, me: np.ndarray) -> float:
    """Pearson r between the across-trial mean-FR and mean-ME curves.

    Captures "how similar are the trial-averaged shapes" — the quantity
    the eye reads off the trace plot. Distinct from the trial-residual
    per-bin r returned by ``_pearson_r_per_bin``.
    """
    with np.errstate(invalid="ignore"):
        m_fr = np.nanmean(fr, axis=0)
        m_me = np.nanmean(me, axis=0)
    ok = np.isfinite(m_fr) & np.isfinite(m_me)
    if int(ok.sum()) < 3:
        return float("nan")
    x = m_fr[ok]
    y = m_me[ok]
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt((x * x).sum() * (y * y).sum())
    if denom <= 0:
        return float("nan")
    return float((x * y).sum() / denom)


def _draw_corr_cell(
    ax,
    fr: np.ndarray,
    me: np.ndarray,
    bin_offsets: np.ndarray,
    rng: np.random.Generator,
    *,
    title: str,
    group_ids: np.ndarray | None = None,
    n_label: int | None = None,
    ylim: tuple[float, float] = YLIM_R_POOLED,
) -> None:
    """Plot one cell: per-bin Pearson r curve + bootstrap CI.

    Annotates ``mean-trace r`` (Pearson r of the trial-averaged FR / ME
    curves) in the upper-right corner — this is the "do the shapes look
    correlated" number, which is much higher than the trial-residual r
    sitting in the curve. Two distinct quantities, both informative.

    ``group_ids`` (optional): trial identifier per row, for the
    cluster-aware trial-grouped bootstrap on the Figure 2 pooled cells.
    ``n_label``: number to report in the title — for Figure 2 this is
    the unique trial count, since the row count is trial × cluster.
    """
    if fr.shape[0] < 3:
        ax.set_title(f"{title}\n(n={fr.shape[0]} — skipped)", fontsize=8)
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="black", lw=0.5, ls="--")
        ax.axvline(T100_MS_S, color="black", lw=0.7, alpha=0.7)
        ax.set_ylim(*ylim)
        ax.set_xlim(*XLIM_S)
        return

    r, n_used = _pearson_r_per_bin(fr, me)
    lo, hi = _bootstrap_ci(fr, me, BOOTSTRAP_N, rng, group_ids=group_ids)
    t = bin_offsets.astype(np.float64) * BIN_WIDTH

    coverage = n_used.astype(np.float64) / max(fr.shape[0], 1)
    well_covered = coverage >= TRIAL_COVERAGE_THRESH

    # Light-grey shading on poorly-covered tail bins.
    if (~well_covered).any():
        ax.fill_between(
            t, ylim[0] * 1.05, ylim[1] * 1.05, where=~well_covered,
            color="0.92", linewidth=0, zorder=0,
        )

    ax.fill_between(t, lo, hi, color="C0", alpha=0.25, linewidth=0)
    ax.plot(t, r, color="C0", lw=1.2)
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="black", lw=0.5, ls="--")
    # Solid +100 ms landmark.
    ax.axvline(T100_MS_S, color="black", lw=0.7, alpha=0.7)
    ax.set_ylim(*ylim)
    ax.set_xlim(*XLIM_S)
    n_for_title = n_label if n_label is not None else fr.shape[0]
    ax.set_title(f"{title}  (n={n_for_title})", fontsize=8)
    # Mean-trace r annotation (the "what your eye sees" number).
    mtr = _mean_trace_r(fr, me)
    if np.isfinite(mtr):
        ax.text(
            0.98, 0.95, f"mean-trace r = {mtr:+.2f}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=7, color="0.25",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1.5),
        )


def _figure_for_cluster(
    pb: ProbeBins,
    cluster_idx: int,
    rng: np.random.Generator,
):
    """Per-cluster figure: 4 rows × 3 cols.

    Rows alternate (raw FR/ME traces, Pearson r curve) for each of the
    two passive replay profiles. Cols = V / VT / T_Vstatic.
    """
    cid = int(pb.cluster_ids[cluster_idx])
    n_rows = 2 * len(PROFILES)
    n_cols = len(CONDITIONS)
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(11, 1.85 * n_rows))
    height_ratios: list[float] = []
    for _ in range(len(PROFILES)):
        height_ratios.extend([0.7, 1.0])
    gs = GridSpec(
        n_rows, n_cols, figure=fig,
        height_ratios=height_ratios, hspace=0.45, wspace=0.22,
    )

    for pi, profile in enumerate(PROFILES):
        for ci, condition in enumerate(CONDITIONS):
            idxs = _cell_trial_indices(pb, condition, profile)
            fr, me, bin_offsets = _stack_cell_arrays(pb, cluster_idx, idxs)
            ax_tr = fig.add_subplot(gs[2 * pi, ci])
            ax_r = fig.add_subplot(gs[2 * pi + 1, ci], sharex=ax_tr)
            _draw_traces_cell(
                ax_tr, fr, me, bin_offsets,
                title=f"{condition}  profile {profile}",
            )
            _draw_corr_cell(
                ax_r, fr, me, bin_offsets, rng,
                title=f"{condition}  profile {profile}",
                ylim=YLIM_R_PER_CLUSTER,
            )
            ax_tr.tick_params(labelbottom=False)
            if ci == 0:
                ax_tr.set_ylabel(f"profile {profile}\nFR (Hz)", fontsize=7)
                ax_r.set_ylabel("Pearson r", fontsize=8)
            if pi == len(PROFILES) - 1:
                ax_r.set_xlabel("Time from motion onset (s)")

    fig.suptitle(
        f"{pb.probe_id}  cluster {cid}  "
        f"(30 ms bins, FR causal half-Gauss SD 50 ms; "
        f"CI: {BOOTSTRAP_N}-resample trial bootstrap)",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def _emit_figure1(pb: ProbeBins, out_pdf: Path, rng: np.random.Generator) -> None:
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    with PdfPages(out_pdf) as pdf:
        for ci in range(pb.cluster_ids.size):
            fig = _figure_for_cluster(pb, ci, rng)
            pdf.savefig(fig)
            plt.close(fig)
            if (ci + 1) % 25 == 0 or (ci + 1) == pb.cluster_ids.size:
                elapsed = time.perf_counter() - t0
                print(
                    f"  [{pb.probe_id}] fig1 page {ci+1}/{pb.cluster_ids.size}"
                    f"  ({elapsed:.1f}s)",
                    flush=True,
                )
    print(f"  [{pb.probe_id}] wrote {out_pdf} ({pb.cluster_ids.size} pages)", flush=True)


# --- Figure 2: animal-level pooled grid ------------------------------------
def _stack_animal_condition(
    pb: ProbeBins, condition: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pooled blocks for one animal × condition.

    Returns ``(fr_z, me_z, fr_raw, me_raw, bin_offsets, trial_row_id)``.

    - ``fr_z`` / ``me_z`` feed the correlation: FR z-scored per cluster,
      ME z-scored per probe. (cluster × trial, n_bins).
    - ``fr_raw`` / ``me_raw`` feed the trace plot: FR in Hz per cluster,
      ME in raw camera units. Same shape; mean over rows gives the
      population-mean trace at each bin.
    - ``trial_row_id`` (n_rows,) carries the trial-table index per row
      so the cluster-aware bootstrap can resample trials (each trial →
      all its cluster rows). IDs are local to this probe — the
      multi-animal stacker offsets them per probe to keep groups
      disjoint.
    """
    trial_idxs = [
        i for i, tb in enumerate(pb.trials) if tb.condition == condition
    ]
    if not trial_idxs:
        return (
            np.zeros((0, 0)), np.zeros((0, 0)),
            np.zeros((0, 0)), np.zeros((0, 0)),
            np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64),
        )

    lo = min(int(pb.trials[ti].bin_offsets[0]) for ti in trial_idxs)
    hi = max(int(pb.trials[ti].bin_offsets[-1]) for ti in trial_idxs)
    bin_offsets = np.arange(lo, hi + 1, dtype=np.int64)
    n_bins = bin_offsets.size
    n_clusters = pb.cluster_ids.size
    n_trials = len(trial_idxs)

    # Pre-compute ME row per trial (raw and z-scored, common grid).
    me_block_raw = np.full((n_trials, n_bins), np.nan, dtype=np.float64)
    me_block_z = np.full((n_trials, n_bins), np.nan, dtype=np.float64)
    for row, ti in enumerate(trial_idxs):
        tb = pb.trials[ti]
        s = int(tb.bin_offsets[0]) - lo
        e = s + tb.bin_offsets.size
        me_block_raw[row, s:e] = tb.me_raw
        me_block_z[row, s:e] = (tb.me_raw - pb.me_mean) / pb.me_std

    fr_z_rows: list[np.ndarray] = []
    fr_raw_rows: list[np.ndarray] = []
    me_z_rows: list[np.ndarray] = []
    me_raw_rows: list[np.ndarray] = []
    group_rows: list[np.ndarray] = []
    trial_ids_per_row_block = np.array(trial_idxs, dtype=np.int64)
    for ci in range(n_clusters):
        c_mean = pb.fr_mean[ci]
        c_std = pb.fr_std[ci]
        if c_std <= 0:
            continue
        fr_z_block = np.full((n_trials, n_bins), np.nan, dtype=np.float64)
        fr_raw_block = np.full((n_trials, n_bins), np.nan, dtype=np.float64)
        for row, ti in enumerate(trial_idxs):
            tb = pb.trials[ti]
            s = int(tb.bin_offsets[0]) - lo
            e = s + tb.bin_offsets.size
            fr_raw_block[row, s:e] = pb.fr[ci][ti]
            fr_z_block[row, s:e] = (pb.fr[ci][ti] - c_mean) / c_std
        fr_z_rows.append(fr_z_block)
        fr_raw_rows.append(fr_raw_block)
        me_z_rows.append(me_block_z)
        me_raw_rows.append(me_block_raw)
        group_rows.append(trial_ids_per_row_block)

    if not fr_z_rows:
        return (
            np.zeros((0, 0)), np.zeros((0, 0)),
            np.zeros((0, 0)), np.zeros((0, 0)),
            bin_offsets, np.zeros(0, dtype=np.int64),
        )
    return (
        np.vstack(fr_z_rows),
        np.vstack(me_z_rows),
        np.vstack(fr_raw_rows),
        np.vstack(me_raw_rows),
        bin_offsets,
        np.concatenate(group_rows),
    )


def _stack_multi_animal(
    pbs: list[ProbeBins], condition: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Concatenate the pooled blocks of multiple animals at a common bin grid.

    Trial-row IDs are offset per probe so the bootstrap treats each
    (probe, trial) as a distinct group. Returns the same 6-tuple as
    ``_stack_animal_condition``: ``(fr_z, me_z, fr_raw, me_raw,
    bin_offsets, trial_row_id)``.
    """
    blocks = [_stack_animal_condition(pb, condition) for pb in pbs]
    blocks = [b for b in blocks if b[0].shape[0] > 0]
    if not blocks:
        return (
            np.zeros((0, 0)), np.zeros((0, 0)),
            np.zeros((0, 0)), np.zeros((0, 0)),
            np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64),
        )
    lo = min(int(b[4][0]) for b in blocks)
    hi = max(int(b[4][-1]) for b in blocks)
    bin_offsets = np.arange(lo, hi + 1, dtype=np.int64)
    n_bins = bin_offsets.size
    fr_z_full: list[np.ndarray] = []
    me_z_full: list[np.ndarray] = []
    fr_raw_full: list[np.ndarray] = []
    me_raw_full: list[np.ndarray] = []
    grp_full: list[np.ndarray] = []
    offset = 0
    for fr_z, me_z, fr_raw, me_raw, off, grp in blocks:
        s = int(off[0]) - lo
        e = s + off.size
        n = fr_z.shape[0]
        for src, dst in (
            (fr_z, fr_z_full), (me_z, me_z_full),
            (fr_raw, fr_raw_full), (me_raw, me_raw_full),
        ):
            pad = np.full((n, n_bins), np.nan)
            pad[:, s:e] = src
            dst.append(pad)
        grp_full.append(grp + offset)
        offset += int(grp.max()) + 1 if grp.size else 0
    return (
        np.vstack(fr_z_full), np.vstack(me_z_full),
        np.vstack(fr_raw_full), np.vstack(me_raw_full),
        bin_offsets, np.concatenate(grp_full),
    )


def _draw_traces_cell(
    ax,
    fr: np.ndarray,
    me: np.ndarray,
    bin_offsets: np.ndarray,
    *,
    title: str,
    fr_units: str = "Hz",
    me_units: str = "ME (raw)",
) -> None:
    """Plot mean FR (black, left y) + mean ME (red, right y) in raw units.

    ``fr`` and ``me`` are raw (no z-scoring) (rows, n_bins) matrices.
    Rows are trials for Figure 1, cluster × trial for Figure 2 — the
    nan-mean along the row axis collapses to the trace.

    Twin y-axes since FR (Hz) and ME (camera units) are not directly
    comparable in scale. Matches Laura's 2026-05-28 ask: see raw values
    instead of the motion-row-z-scored signal whose baseline drifts
    negative.
    """
    if fr.shape[0] < 3:
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="black", lw=0.5, ls="--")
        ax.axvline(T100_MS_S, color="black", lw=0.7, alpha=0.7)
        ax.set_xlim(*XLIM_S)
        ax.set_title(f"{title}\n(insufficient data)", fontsize=8)
        return
    t = bin_offsets.astype(np.float64) * BIN_WIDTH
    with np.errstate(invalid="ignore"):
        mean_fr = np.nanmean(fr, axis=0)
        mean_me = np.nanmean(me, axis=0)
    ax.plot(t, mean_fr, color="black", lw=1.2, label=f"mean FR ({fr_units})")
    ax.set_xlim(*XLIM_S)
    ax.tick_params(axis="y", labelcolor="black", labelsize=7)
    ax.axvline(0, color="black", lw=0.5, ls="--")
    ax.axvline(T100_MS_S, color="black", lw=0.7, alpha=0.7)
    ax.set_title(title, fontsize=8)
    ax_me = ax.twinx()
    ax_me.plot(t, mean_me, color="C3", lw=1.2, label=f"mean {me_units}")
    ax_me.tick_params(axis="y", labelcolor="C3", labelsize=7)


def _emit_figure2(pbs: list[ProbeBins], out_pdf: Path, rng: np.random.Generator) -> None:
    """Emit the pooled per-animal grid.

    Each animal contributes a *pair* of rows (traces + correlation); a
    final pair pools across all 3 animals. Columns are V / VT / T.
    """
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    n_blocks = len(pbs) + 1   # 3 animals + pooled
    n_rows = 2 * n_blocks      # (trace, r) per block
    n_cols = len(CONDITIONS)
    # Use a GridSpec so the trace row can be a bit shorter than the r row.
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(11, 1.8 * n_rows))
    height_ratios: list[float] = []
    for _ in range(n_blocks):
        height_ratios.extend([0.7, 1.0])  # trace, r
    gs = GridSpec(
        n_rows, n_cols, figure=fig,
        height_ratios=height_ratios, hspace=0.45, wspace=0.18,
    )

    def _block(block_idx: int, label: str, stack_fn):
        for ci, condition in enumerate(CONDITIONS):
            fr_z, me_z, fr_raw, me_raw, bin_offsets, groups = stack_fn(condition)
            n_unique = int(np.unique(groups).size) if groups.size else 0
            ax_tr = fig.add_subplot(gs[2 * block_idx, ci])
            ax_r = fig.add_subplot(gs[2 * block_idx + 1, ci], sharex=ax_tr)
            _draw_traces_cell(
                ax_tr, fr_raw, me_raw, bin_offsets,
                title=f"{label}  {condition}  (n_trials={n_unique})",
            )
            _draw_corr_cell(
                ax_r, fr_z, me_z, bin_offsets, rng,
                title=f"{label}  {condition}",
                group_ids=groups,
                n_label=n_unique,
            )
            ax_tr.tick_params(labelbottom=False)
            if ci == 0:
                ax_tr.set_ylabel(f"{label}\nFR (Hz)", fontsize=7, color="black")
                ax_r.set_ylabel("Pearson r", fontsize=8)
            if block_idx == n_blocks - 1:
                ax_r.set_xlabel("Time from motion onset (s)")

    for i, pb in enumerate(pbs):
        _block(
            i, pb.probe_id,
            lambda condition, _pb=pb: _stack_animal_condition(_pb, condition),
        )
    _block(
        n_blocks - 1, "all 3 animals",
        lambda condition: _stack_multi_animal(pbs, condition),
    )

    fig.suptitle(
        "Pooled FR-vs-ME Pearson r per bin "
        f"(z-FR per cluster, z-ME per probe; cluster × trial pooled; "
        f"{BOOTSTRAP_N}-resample trial-grouped bootstrap CI — "
        "each resampled trial pulls all its cluster rows together)",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"wrote {out_pdf}", flush=True)


# --- Cache I/O -------------------------------------------------------------
def _cache_path(probe_id: str) -> Path:
    return CACHE_ROOT / f"{probe_id}.npz"


def _save_cache(pb: ProbeBins, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Variable-length per-trial arrays → pack as object arrays.
    n_trials = len(pb.trials)
    n_clusters = pb.cluster_ids.size
    trial_ids = np.array([tb.trial_id for tb in pb.trials], dtype=np.int64)
    conditions = np.array([tb.condition for tb in pb.trials], dtype=object)
    profiles = np.array([tb.profile_id for tb in pb.trials], dtype=np.int64)
    bin_offsets = np.empty(n_trials, dtype=object)
    me_raw = np.empty(n_trials, dtype=object)
    for ti, tb in enumerate(pb.trials):
        bin_offsets[ti] = tb.bin_offsets.astype(np.int32)
        me_raw[ti] = tb.me_raw.astype(np.float32)
    fr = np.empty((n_clusters, n_trials), dtype=object)
    for ci in range(n_clusters):
        for ti in range(n_trials):
            fr[ci, ti] = pb.fr[ci][ti].astype(np.float32)
    np.savez_compressed(
        path,
        probe_id=pb.probe_id,
        cluster_ids=pb.cluster_ids,
        trial_ids=trial_ids,
        conditions=conditions,
        profiles=profiles,
        bin_offsets=bin_offsets,
        me_raw=me_raw,
        fr=fr,
        me_mean=pb.me_mean,
        me_std=pb.me_std,
        fr_mean=pb.fr_mean,
        fr_std=pb.fr_std,
    )


def _load_cache(path: Path) -> ProbeBins:
    z = np.load(path, allow_pickle=True)
    n_trials = int(z["trial_ids"].size)
    trials = []
    for ti in range(n_trials):
        trials.append(TrialBins(
            trial_id=int(z["trial_ids"][ti]),
            condition=str(z["conditions"][ti]),
            profile_id=int(z["profiles"][ti]),
            bin_offsets=np.asarray(z["bin_offsets"][ti], dtype=np.int64),
            bin_edges=np.zeros(0),  # not cached; not used downstream
            me_raw=np.asarray(z["me_raw"][ti], dtype=np.float64),
        ))
    fr_arr = z["fr"]
    n_clusters = fr_arr.shape[0]
    fr = [[np.asarray(fr_arr[ci, ti], dtype=np.float64) for ti in range(n_trials)]
          for ci in range(n_clusters)]
    return ProbeBins(
        probe_id=str(z["probe_id"]),
        cluster_ids=np.asarray(z["cluster_ids"], dtype=np.int64),
        trials=trials,
        fr=fr,
        me_mean=float(z["me_mean"]),
        me_std=float(z["me_std"]),
        fr_mean=np.asarray(z["fr_mean"], dtype=np.float64),
        fr_std=np.asarray(z["fr_std"], dtype=np.float64),
    )


# --- Figure 3: pooled std diagnostic ---------------------------------------
def _emit_figure3_pooled_std(pbs: list[ProbeBins], out_pdf: Path) -> None:
    """Pooled trial-to-trial SD of FR (Hz) and ME (raw) across all 3 animals.

    Added 2026-05-28 as a diagnostic for the trial-residual r drop seen in
    VT / T motion in Figure 2. If ME variance collapses post motion-onset
    (the stage motion is stereotyped), Pearson r becomes noise-dominated
    and the drop is partly a denominator effect rather than a real
    decoupling. This figure plots the SD curves directly so the cause is
    visible.

    Per cell (3 conditions): black = SD of FR across the pooled
    (cluster × trial) rows at each bin, red = SD of ME across trials at
    each bin. Twin y-axes since FR (Hz) and raw ME aren't comparable.
    Same x-axis convention as Figure 2.
    """
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(CONDITIONS), figsize=(11, 2.8), sharex=True)
    for ci, condition in enumerate(CONDITIONS):
        fr_z, me_z, fr_raw, me_raw, bin_offsets, groups = _stack_multi_animal(
            pbs, condition,
        )
        ax = axes[ci]
        if fr_raw.shape[0] < 3:
            ax.set_title(f"{condition}\n(insufficient data)", fontsize=8)
            ax.set_xlim(*XLIM_S)
            continue
        t = bin_offsets.astype(np.float64) * BIN_WIDTH
        # Trial-axis SD: for FR pool every (cluster × trial) row; for ME
        # de-replicate to one row per unique trial first (the matrix has
        # ME repeated across clusters within a trial, which doesn't
        # change np.nanstd but is conceptually cleaner to report).
        with np.errstate(invalid="ignore"):
            fr_sd = np.nanstd(fr_raw, axis=0, ddof=1)
            unique_groups, first_idx = np.unique(groups, return_index=True)
            me_unique = me_raw[first_idx]
            me_sd = np.nanstd(me_unique, axis=0, ddof=1)
        n_trials = int(unique_groups.size)
        ax.plot(t, fr_sd, color="black", lw=1.2, label="FR SD (Hz)")
        ax.set_xlim(*XLIM_S)
        ax.set_xlabel("Time from motion onset (s)")
        ax.tick_params(axis="y", labelcolor="black", labelsize=7)
        ax.axvline(0, color="black", lw=0.5, ls="--")
        ax.axvline(T100_MS_S, color="black", lw=0.7, alpha=0.7)
        ax.set_title(
            f"{condition}  (n_trials={n_trials}, n_clusters_x_trials={fr_raw.shape[0]})",
            fontsize=8,
        )
        if ci == 0:
            ax.set_ylabel("Trial-to-trial SD of FR (Hz)", fontsize=8)
        ax_me = ax.twinx()
        ax_me.plot(t, me_sd, color="C3", lw=1.2, label="ME SD")
        ax_me.tick_params(axis="y", labelcolor="C3", labelsize=7)
        if ci == len(CONDITIONS) - 1:
            ax_me.set_ylabel("Trial-to-trial SD of ME (raw)", fontsize=8, color="C3")
    fig.suptitle(
        "Trial-to-trial SD of FR (black, left) and ME (red, right) pooled across all 3 animals — "
        "diagnostic for the trial-residual r drop in Figure 2",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"wrote {out_pdf}", flush=True)


# --- Main ------------------------------------------------------------------
def _load_or_build(probe_stem: str, *, use_cache: bool) -> ProbeBins:
    cache = _cache_path(probe_stem)
    if use_cache and cache.exists():
        print(f"[{probe_stem}] loading cache {cache.name}", flush=True)
        return _load_cache(cache)
    print(f"[{probe_stem}] building from .mat …", flush=True)
    lookup = StimulusLookup(SEQ_PATH, FOL_PATH)
    probe = load_probe_data(
        MAT_DIR / f"{probe_stem}.mat",
        config=GLMConfig(),
        stimulus_lookup=lookup,
        visp_only=True,
    )
    pb = _aggregate_probe(probe)
    _save_cache(pb, cache)
    print(f"[{probe_stem}] cached {cache.name}", flush=True)
    return pb


def _glm_cohort_cluster_ids(glm_csv: Path) -> dict[str, set[int]]:
    """Per-probe set of all cluster IDs in the production GLM cohort.

    This is the *default* cohort for motion-clouds per-cluster analyses
    — the ~80 prefilter-passed clusters across the 3 probes
    (CAA-1123243 / 244 / 466). Reading the cluster_id column of
    ``glm_model_comparison.csv`` gives the full set without any
    additional filtering.
    """
    import pandas as pd
    df = pd.read_csv(glm_csv)
    out: dict[str, set[int]] = {}
    for probe_id, ids in df.groupby("probe_id")["cluster_id"]:
        out[str(probe_id)] = set(int(x) for x in ids)
    return out


def _me_selected_cluster_ids(glm_csv: Path) -> dict[str, set[int]]:
    """Per-probe set of cluster IDs whose GLM Selected model includes any
    ME_face component (main effect or ``ME_face × Speed`` interaction).

    Reads ``glm_model_comparison.csv`` (the production output of
    ``rc2_glm.pipeline``) and uses the existing boolean columns
    ``time_is_me_face_tuned`` and ``time_has_me_face_x_speed`` to define
    the cohort.
    """
    import pandas as pd
    df = pd.read_csv(glm_csv)
    mask = df["time_is_me_face_tuned"] | df["time_has_me_face_x_speed"]
    out: dict[str, set[int]] = {}
    for probe_id, ids in df.loc[mask].groupby("probe_id")["cluster_id"]:
        out[str(probe_id)] = set(int(x) for x in ids)
    return out


def _restrict_to_clusters(pb: ProbeBins, keep_ids: set[int]) -> ProbeBins:
    """Return a copy of ``pb`` keeping only clusters whose ID is in ``keep_ids``."""
    keep_idx = [i for i, cid in enumerate(pb.cluster_ids) if int(cid) in keep_ids]
    return ProbeBins(
        probe_id=pb.probe_id,
        cluster_ids=pb.cluster_ids[keep_idx],
        trials=pb.trials,
        fr=[pb.fr[i] for i in keep_idx],
        me_mean=pb.me_mean,
        me_std=pb.me_std,
        fr_mean=pb.fr_mean[keep_idx],
        fr_std=pb.fr_std[keep_idx],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probes", nargs="*", default=list(PROBES))
    parser.add_argument(
        "--skip-fig1", action="store_true",
        help="Skip per-cluster figures (Figure 1)."
    )
    parser.add_argument(
        "--skip-fig2", action="store_true",
        help="Skip the per-animal pooled grid (Figure 2)."
    )
    parser.add_argument(
        "--skip-fig3", action="store_true",
        help="Skip the pooled trial-to-trial SD diagnostic (Figure 3)."
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Always re-bin from .mat, ignore any existing cache."
    )
    parser.add_argument(
        "--cluster-limit", type=int, default=None,
        help="Only emit the first N clusters of each probe (debug)."
    )
    parser.add_argument(
        "--me-filter", action="store_true",
        help="Restrict to clusters whose GLM Selected model includes any "
             "ME_face component (main effect or ME_face × Speed). Output goes to "
             "a separate exploration folder so the GLM-cohort run isn't overwritten.",
    )
    parser.add_argument(
        "--no-cohort-filter", action="store_true",
        help="Disable the default ~80-cluster GLM cohort filter and run on "
             "all VISp clusters instead. Off by default — motion-clouds "
             "per-cluster analyses use the GLM cohort, never raw Kilosort/VISp.",
    )
    parser.add_argument(
        "--glm-csv", type=Path, default=GLM_CSV_DEFAULT,
        help=f"Path to glm_model_comparison.csv defining the cluster cohort (default: {GLM_CSV_DEFAULT}).",
    )
    args = parser.parse_args()

    rng = np.random.default_rng(BOOTSTRAP_RNG_SEED)
    out_root = OUT_ROOT_ME_FILTER if args.me_filter else OUT_ROOT
    out_root.mkdir(parents=True, exist_ok=True)

    # Default cohort = the ~80 prefilter-passed clusters in the production
    # GLM CSV. ME-filter narrows further to ME-selected. --no-cohort-filter
    # opts back into all VISp (debug only).
    cohort_ids: dict[str, set[int]] | None
    if args.no_cohort_filter:
        cohort_ids = None
        print("--no-cohort-filter: running on all VISp clusters", flush=True)
    elif args.me_filter:
        cohort_ids = _me_selected_cluster_ids(args.glm_csv)
        total = sum(len(s) for s in cohort_ids.values())
        print(
            f"ME-filter: keeping {total} clusters across "
            f"{len(cohort_ids)} probes (from {args.glm_csv})",
            flush=True,
        )
    else:
        cohort_ids = _glm_cohort_cluster_ids(args.glm_csv)
        total = sum(len(s) for s in cohort_ids.values())
        print(
            f"GLM cohort (default): keeping {total} clusters across "
            f"{len(cohort_ids)} probes (from {args.glm_csv})",
            flush=True,
        )

    pbs: list[ProbeBins] = []
    for stem in args.probes:
        pb = _load_or_build(stem, use_cache=not args.no_cache)
        if args.cluster_limit is not None:
            keep = slice(0, args.cluster_limit)
            pb = ProbeBins(
                probe_id=pb.probe_id,
                cluster_ids=pb.cluster_ids[keep],
                trials=pb.trials,
                fr=pb.fr[keep],
                me_mean=pb.me_mean,
                me_std=pb.me_std,
                fr_mean=pb.fr_mean[keep],
                fr_std=pb.fr_std[keep],
            )
        if cohort_ids is not None:
            keep_ids = cohort_ids.get(pb.probe_id, set())
            before = pb.cluster_ids.size
            pb = _restrict_to_clusters(pb, keep_ids)
            print(
                f"[{pb.probe_id}] cohort filter: {before} → {pb.cluster_ids.size} clusters",
                flush=True,
            )
        pbs.append(pb)

    if not args.skip_fig1:
        for pb in pbs:
            out_pdf = out_root / "per_cluster" / f"{pb.probe_id}.pdf"
            _emit_figure1(pb, out_pdf, rng)

    if not args.skip_fig2:
        out_pdf = out_root / "per_animal_grid.pdf"
        _emit_figure2(pbs, out_pdf, rng)

    if not args.skip_fig3:
        out_pdf = out_root / "pooled_trial_to_trial_sd.pdf"
        _emit_figure3_pooled_std(pbs, out_pdf)


if __name__ == "__main__":
    main()
