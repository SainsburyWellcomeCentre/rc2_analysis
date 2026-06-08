"""Per-trial data access wrapper around rc2_formatted_data_reader.

Provides structured per-trial records with the velocity-derived motion
mask and (optionally) the per-trial stimulus parameters from a
StimulusLookup.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from rc2_formatted_data_reader import FormattedDataReader, StimulusLookup, masks

logger = logging.getLogger("rc2_glm.io")
from rc2_glm.config import GLMConfig
from rc2_glm.masks_helpers import acceleration_from_velocity


@dataclass
class TrialData:
    trial_idx: int
    trial_id: int
    condition: str            # 'VT' | 'V' | 'T_Vstatic'
    protocol: str
    velocity_channel: str     # 'filtered_teensy' | 'stage' | 'multiplexer_output' | ...
    probe_t: np.ndarray
    velocity: np.ndarray      # filtered velocity (cm/s)
    motion_mask: np.ndarray   # treadmill_motion_mask logic
    stationary_mask: np.ndarray
    sf: float = float("nan")
    orientation: float = float("nan")
    batch_gain: float = float("nan")
    excluded: bool = False
    # Speed-profile id (1 or 2) — which of the two reproduced velocity
    # trajectories this trial belongs to. Mirrors MATLAB
    # ``sp_fold`` (glm_single_cluster_analysis.m:2291-2293). Populated
    # from StimulusLookup.trial_profile_id when a lookup is provided;
    # 0 means "unknown" (no lookup available). Consumed by the
    # speed-profile CV strategy in rc2_glm.cross_validation.
    profile_id: int = 0


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
    # Camera motion energy (session-level, prompt 06, 2026-04-30).
    # Pre-computed pixel-variance traces for the face (camera0) and body
    # (camera1) cameras at camera fs. ``camera_t`` carries the camera
    # sample times in seconds (probe-aligned). Any of the three is
    # ``None`` for sessions where the camera was not recorded; cluster
    # exclusion with ``excluded_reason="no_camera"`` is handled by
    # downstream code that consumes these.
    camera0: np.ndarray | None = None
    camera1: np.ndarray | None = None
    camera_t: np.ndarray | None = None


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

        # Re-label profile_id from the actual recorded velocity trajectory.
        # The stimulus trial-order halving does not track the two reproduced
        # trajectories (they're interleaved across trial_id), so cluster the
        # onset-aligned velocity shape instead. Default on (config flag).
        if getattr(config, "profile_from_velocity", True):
            _assign_profile_ids_by_velocity(trials)

        # Session-level camera motion-energy traces (face + body) and
        # their shared timebase. Reader returns None when the camera
        # group is absent from the .mat. Cheap to load whole (~30 Hz
        # over a session ≈ 54k samples per camera, a few MB at most),
        # so fetch both and let downstream pick via
        # config.motion_energy_camera.
        camera0 = reader.camera0()
        camera1 = reader.camera1()
        camera_t = reader.camera_t()

        return ProbeData(
            probe_id=reader.probe_id,
            mat_path=Path(mat_path),
            fs=fs,
            config=config,
            trials=trials,
            clusters=clusters,
            camera0=camera0,
            camera1=camera1,
            camera_t=camera_t,
        )


def _assign_profile_ids_by_velocity(
    trials: list["TrialData"],
    grid_max_s: float = 3.5,
    n_grid: int = 70,
    dip_window_s: tuple[float, float] = (1.6, 2.1),
) -> None:
    """Re-label ``profile_id`` (1/2) from the recorded velocity trajectory.

    The two reproduced velocity profiles are interleaved across ``trial_id``,
    so the stimulus trial-order halving mislabels them. Each trial's
    onset-aligned ``|velocity|`` is resampled to a common grid; trials are
    clustered (KMeans, k=2) on the **shape** (per-trace z-scored, so the
    V-replay vs stage amplitude difference doesn't dominate). Profile 1 is
    assigned to the cluster with the lower mean velocity in ``dip_window_s``
    (the dip trajectory), profile 2 to the other — matching the
    speed_profiles_reference convention. Mutates ``trials`` in place.

    No-op (leaves existing profile_id) when fewer than 4 trials have a usable
    motion trace or scikit-learn is unavailable.
    """
    grid = np.linspace(0.0, grid_max_s, n_grid)
    idx_usable: list[int] = []
    raw: list[np.ndarray] = []
    for i, t in enumerate(trials):
        mi = np.flatnonzero(t.motion_mask)
        if mi.size == 0:
            continue
        t_rel = t.probe_t - float(t.probe_t[int(mi[0])])
        if float(t_rel.max()) < grid_max_s * 0.7:
            continue  # motion period too short to characterise the trajectory
        vel = np.abs(np.asarray(t.velocity, dtype=np.float64))
        raw.append(np.interp(grid, t_rel, vel))
        idx_usable.append(i)

    if len(idx_usable) < 4:
        logger.warning(
            "profile-from-velocity: only %d usable trials; keeping stimulus "
            "profile_id", len(idx_usable),
        )
        return

    X_raw = np.asarray(raw)
    mu = X_raw.mean(axis=1, keepdims=True)
    sd = X_raw.std(axis=1, keepdims=True)
    X_z = (X_raw - mu) / np.where(sd > 0, sd, 1.0)   # cluster on shape

    try:
        from sklearn.cluster import KMeans
    except Exception as exc:  # pragma: no cover
        logger.warning("profile-from-velocity: sklearn unavailable (%s); "
                       "keeping stimulus profile_id", exc)
        return

    labels = KMeans(n_clusters=2, n_init=10, random_state=0).fit_predict(X_z)

    # Deterministic mapping: profile 1 = the dip trajectory (lower mean
    # |velocity| in the dip window), profile 2 = the other.
    dip = (grid >= dip_window_s[0]) & (grid <= dip_window_s[1])
    mean_dip = {c: X_raw[labels == c][:, dip].mean() for c in (0, 1)}
    dip_cluster = min(mean_dip, key=mean_dip.get)
    to_profile = {dip_cluster: 1, 1 - dip_cluster: 2}

    for i, lab in zip(idx_usable, labels):
        trials[i].profile_id = int(to_profile[int(lab)])

    n1 = sum(1 for c in labels if to_profile[int(c)] == 1)
    logger.info(
        "profile-from-velocity: relabelled %d trials by trajectory shape "
        "(profile1/dip=%d, profile2=%d)", len(idx_usable), n1, len(labels) - n1,
    )


def _load_trial(
    reader: FormattedDataReader,
    trial_idx: int,
    config: GLMConfig,
    stimulus_lookup: StimulusLookup | None,
) -> TrialData:
    """Build one TrialData record with the condition-aware motion mask.

    The motion mask is applied to ``trial_velocity(trial_idx)``, which
    returns the protocol-appropriate channel per the MATLAB ``Trial.m``
    switch (see PROTOCOL_VELOCITY_CHANNEL in the reader). This gives a
    **condition-appropriate motion mask**:

    - T_Vstatic / VT (Coupled / EncoderOnly / StageOnly trials): motion =
      active translation of treadmill or stage.
    - V (ReplayOnly trials): motion = visual flow active (multiplexer_output
      carries the replayed visual velocity; mouse does not translate).

    The stationary mask is the ``analysis_mask \\ motion_mask`` complement,
    which for V trials is the pre-stimulus analysis window (Laura's
    2026-04-23 clarification: "the stationary period [...] is an analysis
    window before"). Regression test:
    ``tests/test_protocol_velocity.py::test_motion_mask_reflects_condition_appropriate_signal``.
    """
    s, e = reader.trial_bounds(trial_idx)
    trial_id = reader.trial_id(trial_idx)
    condition = reader.trial_condition(trial_idx)
    protocol = reader.trial_protocol(trial_idx)
    velocity_channel = reader.trial_velocity_channel(trial_idx)

    v = reader.trial_velocity(
        trial_idx,
        apply_filter=config.apply_velocity_filter,
        cutoff_hz=config.filter_cutoff_hz,
        filter_order=config.filter_order,
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
    profile_id = 0
    if stimulus_lookup is not None:
        params = stimulus_lookup.stimulus_params(trial_id)
        sf = float(params["sf"])
        orient = float(params["orientation"])
        gain = float(params["batch_gain"])
        excluded = bool(params["excluded"])
        profile_id = int(stimulus_lookup.trial_profile_id(trial_id))

    return TrialData(
        trial_idx=trial_idx,
        trial_id=trial_id,
        condition=condition,
        protocol=protocol,
        velocity_channel=velocity_channel,
        probe_t=np.asarray(reader.probe_t(s, e), dtype=np.float64),
        velocity=v,
        motion_mask=m_mask,
        stationary_mask=s_mask,
        sf=sf,
        orientation=orient,
        batch_gain=gain,
        excluded=excluded,
        profile_id=profile_id,
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
