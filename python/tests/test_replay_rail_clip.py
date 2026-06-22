"""Rail-clip for ReplayOnly (V) trials — the ``replay_rail_clip`` flag, and the
no-gap velocity motion mask.

V (ReplayOnly) trials store ``forward_limit=NaN``, so the analysis-window rail
clip no-ops and V motion runs to the velocity end (~4.75 s) instead of the rail
(~3.8 s, matching MATLAB ``to_aligned``). ``replay_rail_clip`` backfills the rail
clip from a StageOnly trial so V and VT motion masks are defined the same way.
Default ``False`` = screens byte-identical. Motion/stationary are the velocity
threshold mask and its complement — they abut, no gap (the cloud displays in sync
with motion). See ``project_motion_clouds_goggles_motion_mask_fix``.

Skips if the goggles formatted-data file is not present locally.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from rc2_glm.config import GLMConfig
from rc2_glm.io import load_probe_data

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


def _motion_lengths(mat: Path, replay: bool) -> dict[str, np.ndarray]:
    probe = load_probe_data(mat, config=GLMConfig(replay_rail_clip=replay),
                            cluster_set="selected")
    out: dict[str, list[float]] = {"V": [], "VT": []}
    for tr in probe.trials:
        if tr.condition not in out:
            continue
        fs = 1.0 / np.median(np.diff(np.asarray(tr.probe_t, float)))
        out[tr.condition].append(np.asarray(tr.motion_mask, bool).sum() / fs)
    return {k: np.asarray(v) for k, v in out.items()}


def test_default_off_screens_safe() -> None:
    """replay_rail_clip is off by default so screens runs stay byte-identical."""
    assert GLMConfig().replay_rail_clip is False


def test_replay_rail_clip_shortens_V(goggles_mat: Path) -> None:
    """ON: V motion rail-clips to ~3.8 s (MATLAB to_aligned); VT unchanged."""
    off = _motion_lengths(goggles_mat, False)
    on = _motion_lengths(goggles_mat, True)
    # V clips well below the un-clipped velocity-end (~4.75 s), to the rail.
    assert np.median(on["V"]) < np.median(off["V"]) - 0.5
    assert 3.6 < np.median(on["V"]) < 4.0
    # VT (StageOnly, finite forward_limit) is untouched by the backfill.
    np.testing.assert_allclose(np.median(on["VT"]), np.median(off["VT"]), atol=1e-6)


def test_masks_abut_no_gap(goggles_mat: Path) -> None:
    """Velocity-mode motion and stationary don't overlap, and the pre-motion
    stationary abuts the motion onset (gap = 1 sample) — no gap."""
    probe = load_probe_data(goggles_mat, config=GLMConfig(replay_rail_clip=True),
                            cluster_set="selected")
    checked = 0
    for tr in probe.trials:
        if tr.condition not in ("V", "VT"):
            continue
        m = np.flatnonzero(np.asarray(tr.motion_mask, bool))
        s = np.flatnonzero(np.asarray(tr.stationary_mask, bool))
        if m.size == 0 or s.size == 0:
            continue
        assert not np.intersect1d(m, s).size  # no overlap
        pre = s[s < m[0]]
        if pre.size:
            assert m[0] - pre.max() == 1   # stationary abuts motion onset
        checked += 1
    assert checked > 0
