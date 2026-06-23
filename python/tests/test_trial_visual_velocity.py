"""Regression tests for the VF/T separation behind the VT trial-structure figure.

The trial-structure poster figure draws BOTH a visual-flow (VF) row and a stage
(T, translation) row. The two come from DIFFERENT session channels, and which
one ``TrialData.velocity`` carries depends on the protocol:

  - V (ReplayOnly): protocol velocity == ``multiplexer_output`` (the visual
    command); the stage is ~0 (mouse does not translate).
  - VT (StageOnly): protocol velocity == ``stage`` (the TRANSLATION); the visual
    command lives on ``multiplexer_output``, a separate channel.

So for VT the figure cannot read VF from ``velocity`` — it must read the
multiplexer explicitly. ``TrialData.visual_velocity`` (populated from
``FormattedDataReader.trial_visual_velocity``) is that channel, filtered like
``velocity``. These tests pin: (1) the reader helper reads the multiplexer, not
the protocol channel, and guards a missing channel; (2) on real data VT's
visual_velocity is the multiplexer while its velocity is the stage, and the two
coincide for V (same channel).
"""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rc2_formatted_data_reader import FormattedDataReader
from rc2_glm.config import GLMConfig
from rc2_glm.io import _load_trial


def test_trial_visual_velocity_returns_empty_when_multiplexer_absent() -> None:
    """The helper must guard a missing channel with an empty array (not raise),
    so probes/cohorts without ``multiplexer_output`` leave visual_velocity empty
    rather than crashing the loader."""
    fake = SimpleNamespace(_sess_group={"stage": object()})  # no multiplexer_output
    out = FormattedDataReader.trial_visual_velocity(fake, 0)
    assert isinstance(out, np.ndarray)
    assert out.size == 0


# --- Data-dependent: the channel separation on a real goggles recording. -------

GOGGLES_MAT = "CAA-1124371_rec1_rec2_rec3.mat"


def _find_goggles_mat() -> Path | None:
    env_dir = os.environ.get("RC2_FORMATTED_DATA_GOGGLES_DIR")
    candidates: list[Path] = []
    if env_dir:
        candidates.append(Path(env_dir) / GOGGLES_MAT)
    candidates.append(
        Path.home() / "local_data" / "motion_clouds"
        / "formatted_data_goggles" / GOGGLES_MAT
    )
    for c in candidates:
        if c.is_file():
            return c
    return None


@pytest.fixture(scope="module")
def goggles_mat_path() -> Path:
    p = _find_goggles_mat()
    if p is None:
        pytest.skip(
            f"{GOGGLES_MAT} not found — set RC2_FORMATTED_DATA_GOGGLES_DIR or drop "
            f"the file under ~/local_data/motion_clouds/formatted_data_goggles/."
        )
    return p


def test_vt_velocity_is_stage_while_visual_velocity_is_multiplexer(
    goggles_mat_path: Path,
) -> None:
    """VT (StageOnly): ``velocity`` is the stage TRANSLATION and ``visual_velocity``
    is the (distinct) visual command; V (ReplayOnly): both are the multiplexer."""
    config = GLMConfig()
    with FormattedDataReader(goggles_mat_path) as r:
        vt_i = next(i for i in range(r.n_trials) if r.trial_condition(i) == "VT")
        v_i = next(i for i in range(r.n_trials) if r.trial_condition(i) == "V")
        vt = _load_trial(r, vt_i, config, None)
        v = _load_trial(r, v_i, config, None)

        mux_vt = r.trial_visual_velocity(
            vt_i,
            apply_filter=config.apply_velocity_filter,
            cutoff_hz=config.filter_cutoff_hz,
            filter_order=config.filter_order,
        )

    # VT: protocol velocity is the stage channel (the translation drive).
    assert vt.velocity_channel == "stage"
    # visual_velocity is populated, same length as velocity, and IS the multiplexer.
    assert vt.visual_velocity.size == vt.velocity.size > 0
    assert np.allclose(vt.visual_velocity, mux_vt)
    # The two channels are genuinely different for VT (stage has the negative
    # return; the visual command is clamped) — VF must NOT be read from velocity.
    assert not np.allclose(vt.visual_velocity, vt.velocity)

    # V: protocol velocity already IS the multiplexer, so the two coincide.
    assert v.velocity_channel == "multiplexer_output"
    assert v.visual_velocity.size == v.velocity.size > 0
    assert np.allclose(v.visual_velocity, v.velocity)
