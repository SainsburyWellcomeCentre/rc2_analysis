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
    # Face/body camera motion energy bin-mean (raw, NOT z-scored).
    # Per prompt 06 (2026-04-30). NaN when the configured camera is
    # absent from the session, or when a bin has no camera samples.
    # Z-scoring happens downstream in pipeline.py over the motion-row
    # support so stationary preludes don't dominate the std.
    "me_face_raw",
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
    me_signal = _resolve_me_signal(probe)
    rows: list[pd.DataFrame] = []
    for trial in probe.trials:
        if trial.excluded:
            continue
        df = _bin_trial_for_cluster(probe, trial, cluster, config, me_signal=me_signal)
        if not df.empty:
            rows.append(df)
    if not rows:
        return pd.DataFrame(columns=COLUMNS)
    return pd.concat(rows, ignore_index=True, copy=False)


# Module-level memo so the truncation / placeholder messages print
# once per (probe_id, cam_attr) instead of once per cluster.
_ME_SIGNAL_CACHE: dict[tuple[str, str], tuple[np.ndarray, np.ndarray] | None] = {}
_ME_SIGNAL_LOGGED: set[tuple[str, str]] = set()


def _resolve_me_signal(probe: ProbeData) -> tuple[np.ndarray, np.ndarray] | None:
    """Pick the configured camera trace (camera0 or camera1) + its timebase.

    Returns None when the camera is absent from the session — downstream
    sees NaN ``me_face_raw`` and the cluster is excluded with
    ``excluded_reason="no_camera"`` at the pipeline level.

    Defensive: truncates camera_vals and camera_t to the common length.
    Some sessions have a few extra samples on one or the other (observed
    on CAA-1123244_rec1 where camera_t had 377147 samples and camera0
    had 375887 — ~1260-sample drift, likely from the MATLAB RC2Format
    write order). Truncation aligns them; the dropped tail is unanalysed.
    """
    cam_attr = getattr(probe.config, "motion_energy_camera", "camera0")
    key = (probe.probe_id, cam_attr)
    if key in _ME_SIGNAL_CACHE:
        return _ME_SIGNAL_CACHE[key]
    log_once = key not in _ME_SIGNAL_LOGGED
    _ME_SIGNAL_LOGGED.add(key)

    cam_vals = getattr(probe, cam_attr, None)
    if cam_vals is None or probe.camera_t is None:
        _ME_SIGNAL_CACHE[key] = None
        return None
    cam_vals = np.asarray(cam_vals, dtype=np.float64).ravel()
    cam_t = np.asarray(probe.camera_t, dtype=np.float64).ravel()

    # Sentinel detection: some sessions write a (2,) all-zero placeholder
    # for camera0/camera1 instead of leaving the dataset absent (observed
    # 2026-04-30 on CAA-1123244_rec1: shape=(2,), dtype=uint64, values
    # [0, 0]). Camera_t is still ~60Hz over the full session in those
    # cases, so the None check above doesn't fire. Heuristic:
    #   - cam_vals.size < 1000 (< ~17s at 60Hz) → almost certainly placeholder
    #   - or all-zero → no real motion energy signal
    # Either condition: treat as "camera not recorded" and bail.
    if cam_vals.size < 1000:
        if log_once:
            print(
                f"[time_binning] probe {probe.probe_id}: {cam_attr} has only "
                f"{cam_vals.size} samples (camera_t has {cam_t.size}); "
                f"treating as 'no camera recorded'.",
                flush=True,
            )
        _ME_SIGNAL_CACHE[key] = None
        return None
    if np.all(cam_vals == 0.0):
        if log_once:
            print(
                f"[time_binning] probe {probe.probe_id}: {cam_attr} is all-zero "
                f"(placeholder); treating as 'no camera recorded'.",
                flush=True,
            )
        _ME_SIGNAL_CACHE[key] = None
        return None

    n = min(cam_vals.size, cam_t.size)
    if cam_vals.size != cam_t.size and log_once:
        print(
            f"[time_binning] probe {probe.probe_id}: {cam_attr} length "
            f"{cam_vals.size} != camera_t length {cam_t.size}; truncating "
            f"both to {n}. Trials beyond camera_t[{n-1}] = {cam_t[n-1]:.1f}s "
            "will have NaN me_face_raw and be dropped from ME-related fits.",
            flush=True,
        )
    out = (cam_vals[:n], cam_t[:n])
    _ME_SIGNAL_CACHE[key] = out
    return out


def _bin_continuous_to_edges(
    values: np.ndarray, sample_t: np.ndarray, bin_edges: np.ndarray
) -> np.ndarray:
    """Per-bin mean of a continuous signal over arbitrary bin edges.

    Aggregates ``values`` (sampled at ``sample_t``) onto ``bin_edges``
    via bin-mean: ``out[i] = mean(values[sample_t in [edges[i], edges[i+1])])``.
    Returns NaN for bins with no contributing samples.
    """
    n_bins = bin_edges.size - 1
    if n_bins <= 0:
        return np.zeros(0, dtype=np.float64)
    sample_bin = np.digitize(sample_t, bin_edges, right=False) - 1
    valid = (sample_bin >= 0) & (sample_bin < n_bins)
    bins_vec = sample_bin[valid]
    vals_vec = values[valid]
    n_per_bin = np.bincount(bins_vec, minlength=n_bins)
    sum_per_bin = np.bincount(
        bins_vec, weights=vals_vec.astype(np.float64), minlength=n_bins
    )
    out = np.full(n_bins, np.nan, dtype=np.float64)
    has_data = n_per_bin > 0
    out[has_data] = sum_per_bin[has_data] / n_per_bin[has_data]
    return out


def _bin_trial_for_cluster(
    probe: ProbeData,
    trial: TrialData,
    cluster: ClusterData,
    config: GLMConfig,
    me_signal: tuple[np.ndarray, np.ndarray] | None = None,
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
        trial, cluster, t_motion_start, t_motion_end, config, me_signal=me_signal
    )

    stat_df = _bin_stationary_prelude(
        trial, cluster, motion_start_idx, t_motion_start, config, me_signal=me_signal
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
    me_signal: tuple[np.ndarray, np.ndarray] | None = None,
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

    if me_signal is not None:
        me_per_bin = _bin_continuous_to_edges(me_signal[0], me_signal[1], bin_edges)[good_idx]
    else:
        me_per_bin = np.full(n_good, np.nan)

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
        "me_face_raw": me_per_bin,
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
    me_signal: tuple[np.ndarray, np.ndarray] | None = None,
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

    if me_signal is not None:
        me_per_bin = _bin_continuous_to_edges(me_signal[0], me_signal[1], bin_edges)
    else:
        me_per_bin = nan.copy()

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
        "me_face_raw": me_per_bin,
    })
