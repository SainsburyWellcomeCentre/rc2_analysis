"""End-to-end regression test for the protocol-dependent velocity fix.

Before prompt 02.5, running the Python pipeline on the all-passive
``CAA-1123243_rec1.mat`` recording produced empty ``glm_model_comparison.csv``
and ``glm_coefficients.csv`` because ``FormattedDataReader.velocity`` read
``filtered_teensy`` unconditionally, and that channel is ~0 on StageOnly /
ReplayOnly trials. This test asserts the visible consequences of the fix on
real data:

1. Motion fraction on the chosen velocity channel is > 10% for both
   StageOnly and ReplayOnly trials.
2. ``bin_cluster`` returns at least one row for an arbitrary VISp cluster —
   i.e. the pipeline has something to fit.

The test skips if the file is not present in the user's ``local_data``. The
asserts are deliberately loose: exact per-trial motion fractions depend on
treadmill/stage geometry, but > 10% is nowhere near the near-zero numbers we
see if the bug regresses.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from rc2_formatted_data_reader import FormattedDataReader
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster


# The passive-paradigm recording this prompt is scoped to.
PASSIVE_MAT_NAME = "CAA-1123243_rec1.mat"


def _find_passive_mat() -> Path | None:
    """Locate CAA-1123243_rec1.mat across the machines we run on."""
    env_dir = os.environ.get("RC2_FORMATTED_DATA_DIR")
    candidates: list[Path] = []
    if env_dir:
        candidates.append(Path(env_dir) / PASSIVE_MAT_NAME)
    env_path = os.environ.get("RC2_FORMATTED_DATA_PATH")
    if env_path and Path(env_path).name == PASSIVE_MAT_NAME:
        candidates.append(Path(env_path))
    candidates.append(
        Path.home() / "local_data" / "motion_clouds" / "formatted_data" / PASSIVE_MAT_NAME
    )
    for c in candidates:
        if c.is_file():
            return c
    return None


@pytest.fixture(scope="module")
def passive_mat_path() -> Path:
    p = _find_passive_mat()
    if p is None:
        pytest.skip(
            f"{PASSIVE_MAT_NAME} not found — set RC2_FORMATTED_DATA_DIR "
            f"or drop the file under ~/local_data/motion_clouds/formatted_data/."
        )
    return p


def test_chosen_channel_motion_fraction_stage_only(passive_mat_path: Path) -> None:
    """StageOnly trials must have > 10% motion on the protocol-appropriate channel."""
    fracs: list[float] = []
    with FormattedDataReader(passive_mat_path) as r:
        for ti in range(r.n_trials):
            if r.trial_protocol(ti) != "StageOnly":
                continue
            v = r.trial_velocity(ti)
            fracs.append(float((np.abs(v) >= 1.0).mean()))
    assert fracs, "no StageOnly trials found in recording"
    mean_frac = float(np.mean(fracs))
    assert mean_frac > 0.10, (
        f"StageOnly chosen-channel motion fraction {mean_frac*100:.2f}% is "
        "suspiciously low — velocity fix has regressed."
    )


def test_chosen_channel_motion_fraction_replay_only(passive_mat_path: Path) -> None:
    """ReplayOnly trials must have > 10% motion on multiplexer_output."""
    fracs: list[float] = []
    with FormattedDataReader(passive_mat_path) as r:
        for ti in range(r.n_trials):
            if r.trial_protocol(ti) != "ReplayOnly":
                continue
            assert r.trial_velocity_channel(ti) == "multiplexer_output"
            v = r.trial_velocity(ti)
            fracs.append(float((np.abs(v) >= 1.0).mean()))
    assert fracs, "no ReplayOnly trials found in recording"
    mean_frac = float(np.mean(fracs))
    assert mean_frac > 0.10, (
        f"ReplayOnly chosen-channel motion fraction {mean_frac*100:.2f}% is "
        "suspiciously low — multiplexer_output is empty."
    )


def test_bin_cluster_returns_rows_on_passive_data(passive_mat_path: Path) -> None:
    """``bin_cluster`` must produce > 0 rows on a VISp cluster of passive data.

    Before the fix, ``motion_mask`` was all False on these trials so every
    cluster returned an empty DataFrame; ``forward_select`` then skipped all
    of them and the model CSVs were empty.
    """
    probe = load_probe_data(passive_mat_path, visp_only=True)
    assert probe.clusters, "no VISp clusters in passive recording"

    # Pick the first cluster with a non-trivial spike count so the test does
    # not depend on cluster ordering.
    cluster = next(
        (c for c in probe.clusters if c.spike_times.size >= 100),
        probe.clusters[0],
    )
    df = bin_cluster(probe, cluster)
    assert not df.empty, (
        f"bin_cluster returned 0 rows for cluster {cluster.cluster_id} — "
        "velocity-mask regression suspected."
    )
    # Motion bins must exist too (not just stationary bins), since motion is
    # what the GLM fits against.
    n_motion = int((df["condition"] != "stationary").sum())
    assert n_motion > 0, (
        f"cluster {cluster.cluster_id} produced only stationary bins — the "
        "chosen velocity channel is not driving any motion mask."
    )
