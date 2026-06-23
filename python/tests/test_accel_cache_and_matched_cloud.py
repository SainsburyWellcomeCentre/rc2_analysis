"""Tests for the per-condition-tuning change: the acceleration MATLAB cache
loader and the matched-cloud trial selection.

Acceleration cache: MATLAB writes ``csvs/acceleration_tuning_curves/<probe>.mat``
with the SAME layout as the Speed/TF caches, EXCEPT it has only the StageOnly
trial groups (VT, T_Vstatic) — acceleration is undefined for the replay (V)
condition. So the generic ``_read_cache_file`` (which defaults to requiring all
three groups) must reject it, and the acceleration path must pass the relaxed
``required_groups=("VT","T_Vstatic")``.

Matched-cloud selection: the figure pairs a V and a VT trial that show the SAME
motion cloud; the picker chooses the cloud where the cluster is most active.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from rc2_glm.precomputed_bins import (
    _TRIAL_GROUPS,
    _read_cache_file,
    load_precomputed_bin_edges,
)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from make_fens_poster_figures import select_matched_cloud_trials  # noqa: E402


# --- matched-cloud selection (pure) ------------------------------------------

def test_select_matched_cloud_trials_picks_most_active_shared_cloud() -> None:
    df = pd.DataFrame({
        "trial_id":   [1,   2,   3,   4],
        "condition":  ["V", "VT", "V", "VT"],
        "spike_count":[10,  20,  100, 200],   # cloud B (3,4) far more active
    })
    trials = {
        1: SimpleNamespace(cloud_name="A"), 2: SimpleNamespace(cloud_name="A"),
        3: SimpleNamespace(cloud_name="B"), 4: SimpleNamespace(cloud_name="B"),
    }
    cloud, v_tid, vt_tid, v_spk, vt_spk = select_matched_cloud_trials(df, trials)
    assert cloud == "B"
    assert (v_tid, vt_tid) == (3, 4)        # the V and VT trial of the SAME cloud
    assert (v_spk, vt_spk) == (100, 200)


def test_select_matched_cloud_trials_none_when_no_shared_cloud() -> None:
    # cloud A only has a V trial, cloud B only a VT trial → nothing pairs.
    df = pd.DataFrame({
        "trial_id": [1, 2], "condition": ["V", "VT"], "spike_count": [10, 20],
    })
    trials = {1: SimpleNamespace(cloud_name="A"), 2: SimpleNamespace(cloud_name="B")}
    assert select_matched_cloud_trials(df, trials) is None


# --- acceleration cache (data-dependent) -------------------------------------

GOGGLES_MAT = "CAA-1124371_rec1_rec2_rec3.mat"


def _goggles_mat() -> Path | None:
    env = os.environ.get("RC2_FORMATTED_DATA_GOGGLES_DIR")
    cands = ([Path(env) / GOGGLES_MAT] if env else []) + [
        Path.home() / "local_data" / "motion_clouds"
        / "formatted_data_goggles" / GOGGLES_MAT
    ]
    return next((c for c in cands if c.is_file()), None)


@pytest.fixture(scope="module")
def goggles_mat() -> Path:
    p = _goggles_mat()
    if p is None or not (p.parent / "csvs" / "acceleration_tuning_curves"
                         / GOGGLES_MAT).is_file():
        pytest.skip("goggles formatted .mat + acceleration_tuning_curves cache "
                    "not present")
    return p


def test_accel_cache_is_two_group_and_rejected_by_default_required_groups(
    goggles_mat: Path,
) -> None:
    accel_path = (goggles_mat.parent / "csvs" / "acceleration_tuning_curves"
                  / goggles_mat.name)
    # The default required set is all 3 groups → the accel cache (no 'V') is rejected.
    with pytest.raises(KeyError):
        _read_cache_file(accel_path)            # defaults to _TRIAL_GROUPS
    assert "V" in _TRIAL_GROUPS                  # guard: 'V' really is in the default
    # Relaxed to the StageOnly groups → it loads.
    groups = _read_cache_file(accel_path, required_groups=("VT", "T_Vstatic"))
    assert set(groups) == {"VT", "T_Vstatic"}


def test_load_precomputed_populates_accel_for_stageonly_only(
    goggles_mat: Path,
) -> None:
    pc = load_precomputed_bin_edges(goggles_mat)
    assert pc is not None
    assert set(pc.accel_by_group) == {"VT", "T_Vstatic"}
    # VT/T_Vstatic have acceleration tuning; V does not (cache has no V group).
    assert pc.accel_centres("VT") is not None
    assert pc.accel_tuning("V", 14) is None
    # Speed/TF caches still load (additive change didn't break them).
    assert pc.speed_centres("V") is not None and pc.tf_centres("VT") is not None
