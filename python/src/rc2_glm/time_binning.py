"""Time-bin trials into a per-cluster long-format table.

Replicates Section 3b of scripts/glm_single_cluster_analysis.m
(lines ~830-1107). For each trial:

- Identify the **motion period** as first..last True sample in
  motion_mask. Bin from t_motion_start..t_motion_end with a fixed bin
  width (default 0.1 s).
- Per bin: count spikes, compute mean |velocity| over the motion samples.
  Keep only bins where the motion-sample fraction >= 0.5.
- Per condition (T_Vstatic / V / VT) populate `speed`, `tf`, `sf`,
  `orientation`, `batch_gain` per the MATLAB switch (lines 982-1007).
- For each trial also collect the *last* contiguous stationary period
  before motion onset, with condition='stationary'.

Output columns match MATLAB T_master_time:
    probe_id, cluster_id, trial_id, condition, speed, tf, sf,
    orientation, batch_gain, spike_count, time_in_trial, time_since_onset
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from rc2_glm.config import GLMConfig
from rc2_glm.io import ClusterData, ProbeData, TrialData

COLUMNS: list[str] = [
    "probe_id", "cluster_id", "trial_id", "profile_id", "condition",
    "speed", "tf", "sf", "orientation", "batch_gain",
    "spike_count", "time_in_trial", "time_since_onset",
]


def bin_probe(probe: ProbeData) -> pd.DataFrame:
    """Bin every (cluster, trial) pair in `probe` and return one DataFrame."""
    frames = []
    for cluster in probe.clusters:
        df = bin_cluster(probe, cluster)
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame(columns=COLUMNS)
    return pd.concat(frames, ignore_index=True, copy=False)


def bin_cluster(probe: ProbeData, cluster: ClusterData) -> pd.DataFrame:
    config = probe.config
    rows: list[pd.DataFrame] = []
    for trial in probe.trials:
        if trial.excluded:
            continue
        df = _bin_trial_for_cluster(probe, trial, cluster, config)
        if not df.empty:
            rows.append(df)
    if not rows:
        return pd.DataFrame(columns=COLUMNS)
    return pd.concat(rows, ignore_index=True, copy=False)


def _bin_trial_for_cluster(
    probe: ProbeData,
    trial: TrialData,
    cluster: ClusterData,
    config: GLMConfig,
) -> pd.DataFrame:
    motion_idx = np.flatnonzero(trial.motion_mask)
    if motion_idx.size == 0:
        return pd.DataFrame(columns=COLUMNS)

    motion_start_idx = int(motion_idx[0])
    motion_end_idx = int(motion_idx[-1])
    t_motion_start = float(trial.probe_t[motion_start_idx])
    t_motion_end = float(trial.probe_t[motion_end_idx])
    if (t_motion_end - t_motion_start) < config.time_bin_width:
        return pd.DataFrame(columns=COLUMNS)

    motion_df = _bin_motion_period(
        trial, cluster, t_motion_start, t_motion_end, config
    )

    stat_df = _bin_stationary_prelude(
        trial, cluster, motion_start_idx, t_motion_start, config
    )

    parts = [df for df in (motion_df, stat_df) if not df.empty]
    if not parts:
        return pd.DataFrame(columns=COLUMNS)

    out = pd.concat(parts, ignore_index=True, copy=False)
    out.insert(0, "probe_id", probe.probe_id)
    out["cluster_id"] = cluster.cluster_id
    out["trial_id"] = trial.trial_id
    out["profile_id"] = trial.profile_id
    return out[COLUMNS]


def _bin_motion_period(
    trial: TrialData,
    cluster: ClusterData,
    t_motion_start: float,
    t_motion_end: float,
    config: GLMConfig,
) -> pd.DataFrame:
    bin_edges = np.arange(
        t_motion_start, t_motion_end + 1e-12, config.time_bin_width, dtype=np.float64
    )
    if bin_edges.size < 2:
        return pd.DataFrame(columns=COLUMNS[3:])
    n_bins = bin_edges.size - 1

    sample_bin = np.digitize(trial.probe_t, bin_edges, right=False) - 1
    sample_bin[sample_bin < 0] = -1
    sample_bin[sample_bin >= n_bins] = -1
    valid = sample_bin >= 0
    bins_vec = sample_bin[valid]
    mmask_vec = trial.motion_mask[valid].astype(np.int64)
    absvel_vec = np.abs(trial.velocity[valid]) * mmask_vec

    n_samp_per_bin = np.bincount(bins_vec, minlength=n_bins)
    n_motion_per_bin = np.bincount(bins_vec, weights=mmask_vec, minlength=n_bins)
    sum_speed_per_bin = np.bincount(bins_vec, weights=absvel_vec, minlength=n_bins)

    good = (n_samp_per_bin > 0) & (
        n_motion_per_bin >= config.motion_fraction_threshold * n_samp_per_bin
    )
    if not good.any():
        return pd.DataFrame(columns=COLUMNS[3:])

    good_idx = np.flatnonzero(good)
    mean_speed = sum_speed_per_bin[good_idx] / n_motion_per_bin[good_idx]
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    time_in_trial = bin_centres[good_idx] - t_motion_start

    spike_times = cluster.spike_times
    counts_all, _ = np.histogram(spike_times, bins=bin_edges)
    spike_counts = counts_all[good_idx].astype(np.int64)

    n_good = good_idx.size
    speed_v, tf_v, sf_v, or_v, gain_v = _condition_vectors(trial, mean_speed)

    return pd.DataFrame({
        "condition": np.full(n_good, trial.condition, dtype=object),
        "speed": speed_v,
        "tf": tf_v,
        "sf": sf_v,
        "orientation": or_v,
        "batch_gain": gain_v,
        "spike_count": spike_counts,
        "time_in_trial": time_in_trial,
        "time_since_onset": time_in_trial.copy(),
    })


def _condition_vectors(trial: TrialData, mean_speed: np.ndarray):
    n = mean_speed.size
    nan = np.full(n, np.nan)
    if trial.condition == "T_Vstatic":
        return mean_speed.copy(), np.zeros(n), nan.copy(), nan.copy(), nan.copy()
    if trial.condition == "V":
        return (
            np.zeros(n),
            trial.batch_gain * mean_speed,
            np.full(n, trial.sf),
            np.full(n, trial.orientation),
            np.full(n, trial.batch_gain),
        )
    if trial.condition == "VT":
        return (
            mean_speed.copy(),
            trial.batch_gain * mean_speed,
            np.full(n, trial.sf),
            np.full(n, trial.orientation),
            np.full(n, trial.batch_gain),
        )
    # Unknown condition: degrade to T_Vstatic-style row (matches MATLAB otherwise)
    return mean_speed.copy(), np.zeros(n), nan.copy(), nan.copy(), nan.copy()


def _bin_stationary_prelude(
    trial: TrialData,
    cluster: ClusterData,
    motion_start_idx: int,
    t_motion_start: float,
    config: GLMConfig,
) -> pd.DataFrame:
    pre = trial.stationary_mask.copy()
    pre[motion_start_idx:] = False
    stat_idx = np.flatnonzero(pre)
    if stat_idx.size == 0:
        return pd.DataFrame(columns=COLUMNS[3:])

    # Last contiguous stationary run before motion onset (gap > 1 sample)
    gaps = np.flatnonzero(np.diff(stat_idx) > 1)
    if gaps.size > 0:
        stat_start_idx = int(stat_idx[gaps[-1] + 1])
    else:
        stat_start_idx = int(stat_idx[0])
    stat_end_idx = int(stat_idx[-1])

    t_stat_start = float(trial.probe_t[stat_start_idx])
    t_stat_end = float(trial.probe_t[stat_end_idx])
    if (t_stat_end - t_stat_start) < config.time_bin_width:
        return pd.DataFrame(columns=COLUMNS[3:])

    bin_edges = np.arange(
        t_stat_start, t_stat_end + 1e-12, config.time_bin_width, dtype=np.float64
    )
    if bin_edges.size < 2:
        return pd.DataFrame(columns=COLUMNS[3:])

    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    time_since_onset = bin_centres - t_motion_start  # negative

    counts, _ = np.histogram(cluster.spike_times, bins=bin_edges)
    n = counts.size
    nan = np.full(n, np.nan)

    return pd.DataFrame({
        "condition": np.full(n, "stationary", dtype=object),
        "speed": np.zeros(n),
        "tf": np.zeros(n),
        "sf": nan.copy(),
        "orientation": nan.copy(),
        "batch_gain": nan.copy(),
        "spike_count": counts.astype(np.int64),
        "time_in_trial": np.zeros(n),
        "time_since_onset": time_since_onset,
    })
