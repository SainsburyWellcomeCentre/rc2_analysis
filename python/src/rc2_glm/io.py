"""Per-trial data access wrapper around rc2_formatted_data_reader.

Provides structured per-trial records with the velocity-derived motion
mask and (optionally) the per-trial stimulus parameters from a
StimulusLookup.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from rc2_formatted_data_reader import FormattedDataReader, StimulusLookup, masks
from rc2_glm.config import GLMConfig
from rc2_glm.masks_helpers import acceleration_from_velocity, filter_velocity


@dataclass
class TrialData:
    trial_idx: int
    trial_id: int
    condition: str            # 'VT' | 'V' | 'T_Vstatic'
    protocol: str
    probe_t: np.ndarray
    velocity: np.ndarray      # filtered velocity (cm/s)
    motion_mask: np.ndarray   # treadmill_motion_mask logic
    stationary_mask: np.ndarray
    sf: float = float("nan")
    orientation: float = float("nan")
    batch_gain: float = float("nan")
    excluded: bool = False


@dataclass
class ClusterData:
    cluster_idx: int
    cluster_id: int
    region: str
    spike_times: np.ndarray   # full session, seconds


@dataclass
class ProbeData:
    """All trial / cluster data for a single formatted .mat file."""
    probe_id: str
    mat_path: Path
    fs: float
    config: GLMConfig
    trials: list[TrialData] = field(default_factory=list)
    clusters: list[ClusterData] = field(default_factory=list)


def load_probe_data(
    mat_path: str | Path,
    config: GLMConfig | None = None,
    stimulus_lookup: StimulusLookup | None = None,
    cluster_indices: Iterable[int] | None = None,
    trial_indices: Iterable[int] | None = None,
    visp_only: bool = True,
) -> ProbeData:
    """Load a probe .mat file into structured ProbeData."""
    config = config or GLMConfig()
    with FormattedDataReader(mat_path) as reader:
        fs = reader.fs

        if cluster_indices is None:
            cluster_indices = (
                reader.visp_cluster_indices().tolist()
                if visp_only
                else list(range(reader.n_clusters))
            )
        cluster_indices = list(cluster_indices)
        cluster_ids = reader.cluster_ids()

        clusters = [
            ClusterData(
                cluster_idx=int(ci),
                cluster_id=int(cluster_ids[ci]),
                region=reader.cluster_region(int(ci)),
                spike_times=reader.spike_times(int(ci)),
            )
            for ci in cluster_indices
        ]

        if trial_indices is None:
            trial_indices = list(range(reader.n_trials))
        trial_indices = list(trial_indices)

        trials = [
            _load_trial(reader, int(ti), config, stimulus_lookup)
            for ti in trial_indices
        ]

        return ProbeData(
            probe_id=reader.probe_id,
            mat_path=Path(mat_path),
            fs=fs,
            config=config,
            trials=trials,
            clusters=clusters,
        )


def _load_trial(
    reader: FormattedDataReader,
    trial_idx: int,
    config: GLMConfig,
    stimulus_lookup: StimulusLookup | None,
) -> TrialData:
    s, e = reader.trial_bounds(trial_idx)
    trial_id = reader.trial_id(trial_idx)
    condition = reader.trial_condition(trial_idx)
    protocol = reader.trial_protocol(trial_idx)

    raw_v = np.asarray(reader.velocity(s, e), dtype=np.float64)
    v = (
        filter_velocity(raw_v, reader.fs, config.filter_cutoff_hz, config.filter_order)
        if config.apply_velocity_filter
        else raw_v
    )

    accel = acceleration_from_velocity(v)
    m_mask = masks.treadmill_motion_mask(
        v,
        accel,
        reader.fs,
        vel_thresh=config.velocity_threshold,
        acc_thresh=config.acceleration_threshold,
        min_dur=config.min_stationary_duration,
    )
    analysis_mask = reader.trial_analysis_mask(trial_idx)
    m_mask = m_mask & analysis_mask
    s_mask = (~m_mask) & analysis_mask

    sf = orient = gain = float("nan")
    excluded = False
    if stimulus_lookup is not None:
        params = stimulus_lookup.stimulus_params(trial_id)
        sf = float(params["sf"])
        orient = float(params["orientation"])
        gain = float(params["batch_gain"])
        excluded = bool(params["excluded"])

    return TrialData(
        trial_idx=trial_idx,
        trial_id=trial_id,
        condition=condition,
        protocol=protocol,
        probe_t=np.asarray(reader.probe_t(s, e), dtype=np.float64),
        velocity=v,
        motion_mask=m_mask,
        stationary_mask=s_mask,
        sf=sf,
        orientation=orient,
        batch_gain=gain,
        excluded=excluded,
    )


def iter_cluster_spike_slices(
    cluster: ClusterData, trial: TrialData
) -> np.ndarray:
    """Spike times falling within the trial's probe_t window."""
    if trial.probe_t.size == 0:
        return np.array([], dtype=np.float64)
    t0, t1 = trial.probe_t[0], trial.probe_t[-1]
    st = cluster.spike_times
    return st[(st >= t0) & (st <= t1)]
