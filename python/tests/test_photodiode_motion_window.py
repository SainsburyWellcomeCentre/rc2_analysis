"""Photodiode-defined motion window for the goggles cohort.

The goggles V/VT trials have a real stationary↔motion GAP: the cloud only
displays ~0.6–1.0 s after the velocity command starts (command→display
latency), and lasts ~4 s. ``motion_window_source="photodiode"`` defines motion
as the photodiode visual-stimulus window for the visual conditions and excludes
the gap; ``"velocity"`` (default, screens-safe) keeps the velocity threshold
mask. ``T_Vstatic`` has no visual stimulus → flat photodiode → velocity
fallback. See ``reference_motion_clouds_goggles_trial_structure``.

Skips if the goggles formatted-data file is not present locally.
"""
from __future__ import annotations

import dataclasses
import os
from pathlib import Path

import numpy as np
import pytest

from rc2_formatted_data_reader import FormattedDataReader
from rc2_glm.config import GLMConfig
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

GOGGLES_MAT_NAME = "CAA-1124370_rec1_rec2_rec3.mat"


def _find_goggles_mat() -> Path | None:
    env_dir = os.environ.get("RC2_GOGGLES_DATA_DIR")
    candidates: list[Path] = []
    if env_dir:
        candidates.append(Path(env_dir) / GOGGLES_MAT_NAME)
    candidates.append(
        Path.home() / "local_data" / "motion_clouds" / "formatted_data_goggles"
        / GOGGLES_MAT_NAME
    )
    for c in candidates:
        if c.is_file():
            return c
    return None


@pytest.fixture(scope="module")
def goggles_mat() -> Path:
    p = _find_goggles_mat()
    if p is None:
        pytest.skip(
            f"{GOGGLES_MAT_NAME} not found — set RC2_GOGGLES_DATA_DIR or drop it "
            f"under ~/local_data/motion_clouds/formatted_data_goggles/."
        )
    return p


def test_screens_default_is_velocity() -> None:
    """The default config must stay on the velocity definition (screens-safe)."""
    assert GLMConfig().motion_window_source == "velocity"


def test_visual_window_per_condition(goggles_mat: Path) -> None:
    """V/VT (visual) trials yield a ~4 s photodiode window; T_Vstatic yields None."""
    with FormattedDataReader(goggles_mat) as r:
        durs_visual, n_tvstatic_none = [], 0
        for ti in range(r.n_trials):
            cond = r.trial_condition(ti)
            win = r.trial_visual_window(ti)
            if cond in ("V", "VT"):
                assert win is not None, f"trial {ti} ({cond}) has no visual window"
                on, off = win
                durs_visual.append((off - on) / r.fs)
            elif cond == "T_Vstatic":
                if win is None:
                    n_tvstatic_none += 1
        assert durs_visual, "no V/VT trials found"
        # visual stimulus is ~4 s (3.5–4.3 s across trials)
        assert 3.4 < np.median(durs_visual) < 4.4
        # T_Vstatic has no visual stimulus → flat photodiode → None
        assert n_tvstatic_none > 0


def test_photodiode_motion_mask_and_gap(goggles_mat: Path) -> None:
    """Photodiode mode: motion = visual window (~4 s) with a gap; velocity mode:
    longer motion, no gap fields. A V trial is used."""
    base = GLMConfig()
    cfg_pd = dataclasses.replace(base, motion_window_source="photodiode")
    probe_pd = load_probe_data(goggles_mat, config=cfg_pd, cluster_set="selected")
    v_trials = [t for t in probe_pd.trials if t.condition == "V"]
    assert v_trials
    tr = v_trials[0]
    assert tr.command_onset_idx is not None
    assert tr.visual_onset_idx is not None and tr.visual_offset_idx is not None
    fs = 1.0 / float(np.median(np.diff(tr.probe_t)))
    gap_s = (tr.visual_onset_idx - tr.command_onset_idx) / fs
    assert 0.3 < gap_s < 1.3, f"gap {gap_s:.3f}s out of expected range"
    motion_s = int(np.asarray(tr.motion_mask, bool).sum()) / fs
    assert 3.4 < motion_s < 4.4, f"photodiode motion {motion_s:.3f}s not ~4 s"

    # velocity mode: same trial, no gap fields, motion runs longer (to vel end)
    cfg_vel = dataclasses.replace(base, motion_window_source="velocity")
    probe_vel = load_probe_data(goggles_mat, config=cfg_vel, cluster_set="selected")
    tr_vel = next(t for t in probe_vel.trials if t.trial_id == tr.trial_id)
    assert tr_vel.command_onset_idx is None
    motion_vel_s = int(np.asarray(tr_vel.motion_mask, bool).sum()) / fs
    assert motion_vel_s > motion_s, "velocity motion should exceed photodiode window"


def test_binned_gap_excluded(goggles_mat: Path) -> None:
    """Under photodiode mode the binned df has stationary (negative tso) and
    motion (~0..4 s) with the gap unrepresented (no bins in the gap)."""
    cfg = dataclasses.replace(GLMConfig(), motion_window_source="photodiode")
    probe = load_probe_data(goggles_mat, config=cfg, cluster_set="selected")
    cl = probe.clusters[0]
    df = bin_cluster(probe, cl)
    v_tids = {t.trial_id for t in probe.trials if t.condition == "V"}
    sub = df[df["trial_id"].isin(v_tids)]
    motion = sub[sub["condition"] != "stationary"]
    assert not motion.empty
    # motion bins start at ~0 (visual onset) and span ~4 s
    assert motion["time_since_onset"].min() >= -0.05
    assert 3.0 < motion["time_since_onset"].max() < 4.5
