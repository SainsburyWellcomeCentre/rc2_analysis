"""Unit + integration tests for ``rc2_glm.precomputed_bins``.

The integration test is gated on the presence of a real formatted .mat beside
``csvs/tuning_curves/`` (local-only — skipped in CI). The unit tests cover
the soft-fail path (missing files) and the prc_per_bin sanity check so the
pipeline doesn't silently diverge if MATLAB changes the bin convention.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from rc2_glm.precomputed_bins import (
    PrecomputedBinEdges,
    _read_edges_file,
    load_precomputed_bin_edges,
)


def test_missing_csvs_dir_returns_none(tmp_path):
    """Pipeline keeps running on recordings where edges haven't been precomputed."""
    fake_mat = tmp_path / "PROBE_X_rec1.mat"
    fake_mat.write_bytes(b"")  # parent exists, csvs/ subtree does not
    out = load_precomputed_bin_edges(fake_mat)
    assert out is None


def test_partial_csvs_dir_returns_none(tmp_path):
    """Soft-fail when only one of the two edge files exists."""
    probe = "PROBE_Y_rec1"
    (tmp_path / f"{probe}.mat").write_bytes(b"")
    speed_dir = tmp_path / "csvs" / "tuning_curves"
    speed_dir.mkdir(parents=True)
    (speed_dir / f"{probe}.mat").write_bytes(b"")  # only speed file present
    out = load_precomputed_bin_edges(tmp_path / f"{probe}.mat")
    assert out is None


def test_prc_per_bin_sanity_check(tmp_path):
    """If MATLAB changes prc_per_bin away from 5.0 we must fail loud.

    Rationale: Python downstream plots and tuning-curve bin counts assume
    5% bins. Silently using 10% or 1% edges would produce tuning curves
    that look like they match MATLAB but actually don't.
    """
    path = tmp_path / "fake.mat"
    _write_fake_edges_file(path, prc_per_bin=10.0)
    with pytest.raises(ValueError, match="prc_per_bin=10"):
        _read_edges_file(path)


def test_read_edges_file_parses_all_three_trial_groups(tmp_path):
    """Round-trip: write a fake HDF5 edges file, read it back with expected keys."""
    path = tmp_path / "fake.mat"
    _write_fake_edges_file(path, prc_per_bin=5.0)
    out = _read_edges_file(path)
    assert set(out.keys()) == {"VT", "V", "T_Vstatic"}
    for edges in out.values():
        assert edges.dtype == np.float64
        assert edges.size == 21  # 5% × 20 bins + 1 edge


def test_precomputed_bin_edges_accessors():
    """The two accessors fall back cleanly on unknown conditions."""
    pb = PrecomputedBinEdges(
        speed_by_group={"VT": np.array([0.0, 1.0])},
        tf_by_group={"V": np.array([0.0, 2.0])},
    )
    assert pb.speed_edges("VT") is not None
    assert pb.speed_edges("T_Vstatic") is None
    assert pb.tf_edges("V") is not None
    assert pb.tf_edges("VT") is None


@pytest.mark.skipif(
    not Path(
        "/Users/lauraporta/local_data/motion_clouds/formatted_data"
        "/csvs/tuning_curves/CAA-1123243_rec1.mat"
    ).exists(),
    reason="local formatted_data not present; run on metonymy to exercise",
)
def test_loader_on_real_formatted_data():
    """End-to-end: parse real MATLAB-written edges and check shapes + range."""
    path = Path(
        "/Users/lauraporta/local_data/motion_clouds/formatted_data/CAA-1123243_rec1.mat"
    )
    edges = load_precomputed_bin_edges(path)
    assert edges is not None
    assert set(edges.speed_by_group.keys()) == {"VT", "V", "T_Vstatic"}
    assert set(edges.tf_by_group.keys()) == {"VT", "V", "T_Vstatic"}
    for cond in ("VT", "V", "T_Vstatic"):
        sp = edges.speed_by_group[cond]
        tf = edges.tf_by_group[cond]
        # 5% bins → 20 bins → 21 edges.
        assert sp.size == 21
        assert tf.size == 21
        # Speed edges span roughly 0..55 cm/s for motion clouds.
        assert 0.0 - 0.1 < sp[0] < 0.5
        assert 50.0 < sp[-1] < 60.0
        # TF edges span roughly 0..8 Hz.
        assert 0.0 - 0.01 < tf[0] < 0.1
        assert 5.0 < tf[-1] < 10.0


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _write_fake_edges_file(path: Path, *, prc_per_bin: float) -> None:
    """Write a minimal v7.3-style HDF5 that mimics the MATLAB layout."""
    with h5py.File(path, "w") as f:
        refs = f.create_group("#refs#")
        tg_refs = []
        tc_refs = []
        for name in ("VT", "V", "T_Vstatic"):
            name_arr = np.array([ord(c) for c in name], dtype=np.uint16)
            name_ds = refs.create_dataset(f"name_{name}", data=name_arr)
            tg_refs.append(name_ds.ref)

            # One trial group → one field group with cluster-wise edges/prc
            grp = refs.create_group(f"group_{name}")
            edges = np.linspace(0.0, 10.0, 21, dtype=np.float64)
            prc = np.array([[prc_per_bin]], dtype=np.float64)
            # 1-column ragged-style dataset of refs
            edges_ds = refs.create_dataset(
                f"edges_{name}", data=edges.reshape(-1, 1),
            )
            prc_ds = refs.create_dataset(f"prc_{name}", data=prc)

            # Build (96, 1) object-ref arrays pointing at the same underlying arrays.
            edges_refs = np.array(
                [[edges_ds.ref] for _ in range(96)], dtype=h5py.ref_dtype,
            )
            prc_refs = np.array(
                [[prc_ds.ref] for _ in range(96)], dtype=h5py.ref_dtype,
            )
            grp.create_dataset("bin_edges", data=edges_refs)
            grp.create_dataset("prc_per_bin", data=prc_refs)
            tc_refs.append(grp.ref)

        f.create_dataset(
            "trial_groups",
            data=np.array([[r] for r in tg_refs], dtype=h5py.ref_dtype),
        )
        f.create_dataset(
            "tuning_curves",
            data=np.array([[r] for r in tc_refs], dtype=h5py.ref_dtype),
        )
