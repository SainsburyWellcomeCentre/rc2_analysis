"""Unit + integration tests for ``rc2_glm.precomputed_bins``.

The integration test is gated on the presence of a real formatted .mat beside
``csvs/tuning_curves/`` (local-only — skipped in CI). The unit tests cover
the soft-fail path (missing files) and the prc_per_bin sanity check so the
pipeline doesn't silently diverge if MATLAB changes the bin convention.

Updated 2026-04-28 (prompt-12 follow-up): the cache exposes per-trial-per-bin
firing rates (Laura's call to use the .mat-cached values directly for the
Observed-row plot rather than recomputing from raw spikes). Tests now cover
the new ``tuning`` / ``trial_ids`` / ``stationary_fr`` / ``bin_centres``
fields too.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from rc2_glm.precomputed_bins import (
    PerConditionTuning,
    PrecomputedBinEdges,
    _read_cache_file,
    load_precomputed_bin_edges,
)


def test_missing_csvs_dir_returns_none(tmp_path):
    """Pipeline keeps running on recordings where caches haven't been precomputed."""
    fake_mat = tmp_path / "PROBE_X_rec1.mat"
    fake_mat.write_bytes(b"")
    out = load_precomputed_bin_edges(fake_mat)
    assert out is None


def test_partial_csvs_dir_returns_none(tmp_path):
    """Soft-fail when only one of the two cache files exists."""
    probe = "PROBE_Y_rec1"
    (tmp_path / f"{probe}.mat").write_bytes(b"")
    speed_dir = tmp_path / "csvs" / "tuning_curves"
    speed_dir.mkdir(parents=True)
    (speed_dir / f"{probe}.mat").write_bytes(b"")  # only speed file present
    out = load_precomputed_bin_edges(tmp_path / f"{probe}.mat")
    assert out is None


def test_prc_per_bin_sanity_check(tmp_path):
    """If MATLAB changes prc_per_bin away from 5.0 we must fail loud."""
    path = tmp_path / "fake.mat"
    _write_fake_cache_file(path, prc_per_bin=10.0)
    with pytest.raises(ValueError, match="prc_per_bin=10"):
        _read_cache_file(path)


def test_read_cache_file_parses_all_three_trial_groups(tmp_path):
    """Round-trip: write a fake HDF5 cache, read it back with expected keys + content."""
    path = tmp_path / "fake.mat"
    _write_fake_cache_file(path, prc_per_bin=5.0)
    out = _read_cache_file(path)
    assert set(out.keys()) == {"VT", "V", "T_Vstatic"}
    for cache in out.values():
        assert isinstance(cache, PerConditionTuning)
        assert cache.bin_edges.dtype == np.float64
        assert cache.bin_edges.size == 21  # 5% × 20 bins + 1 edge
        assert cache.bin_centres.size == 20
        assert cache.cluster_ids.size == 3  # fixture writes 3 clusters
        # Per-cluster fields populated for each cluster_id
        for cid in cache.cluster_ids:
            assert cache.tuning[cid].shape == (5, 20)  # 5 trials × 20 bins
            assert cache.trial_ids[cid].shape == (5,)
            assert cache.stationary_fr[cid].shape == (5,)
            assert cache.stationary_time[cid].shape == (5,)


def test_precomputed_bin_edges_accessors():
    """Edge / centres / per-cluster accessors fall back cleanly on unknown keys."""
    cache_vt = PerConditionTuning(
        bin_edges=np.array([0.0, 1.0]),
        bin_centres=np.array([0.5]),
        cluster_ids=np.array([42]),
        tuning={42: np.array([[3.0]])},
        trial_ids={42: np.array([1])},
        stationary_fr={42: np.array([0.5])},
        stationary_time={42: np.array([4.0])},
    )
    pb = PrecomputedBinEdges(
        speed_by_group={"VT": cache_vt},
        tf_by_group={},
    )
    assert pb.speed_edges("VT") is not None
    assert pb.speed_edges("T_Vstatic") is None
    assert pb.tf_edges("V") is None
    assert pb.speed_centres("VT").size == 1
    assert pb.speed_tuning("VT", 42).shape == (1, 1)
    assert pb.speed_tuning("VT", 999) is None  # unknown cluster → None
    assert pb.tf_tuning("VT", 42) is None       # not in tf store
    assert pb.stationary_fr("speed", "VT", 42)[0] == 0.5
    assert pb.stationary_fr("speed", "VT", 999) is None


@pytest.mark.skipif(
    not Path(
        "/Users/lauraporta/local_data/motion_clouds/formatted_data"
        "/csvs/tuning_curves/CAA-1123243_rec1.mat"
    ).exists(),
    reason="local formatted_data not present; run on metonymy to exercise",
)
def test_loader_on_real_formatted_data():
    """End-to-end: parse real MATLAB-written cache and check shapes + range."""
    path = Path(
        "/Users/lauraporta/local_data/motion_clouds/formatted_data/CAA-1123243_rec1.mat"
    )
    cache = load_precomputed_bin_edges(path)
    assert cache is not None
    assert set(cache.speed_by_group.keys()) == {"VT", "V", "T_Vstatic"}
    assert set(cache.tf_by_group.keys()) == {"VT", "V", "T_Vstatic"}
    for cond in ("VT", "V", "T_Vstatic"):
        sp = cache.speed_by_group[cond]
        tf = cache.tf_by_group[cond]
        assert sp.bin_edges.size == 21
        assert tf.bin_edges.size == 21
        assert sp.bin_centres.size == 20
        assert tf.bin_centres.size == 20
        # Speed edges span roughly 0..55 cm/s for motion clouds.
        assert 0.0 - 0.1 < sp.bin_edges[0] < 0.5
        assert 50.0 < sp.bin_edges[-1] < 60.0
        # TF edges span roughly 0..8 Hz.
        assert 0.0 - 0.01 < tf.bin_edges[0] < 0.1
        assert 5.0 < tf.bin_edges[-1] < 10.0
        # At least one cluster present, with consistent shapes.
        assert sp.cluster_ids.size > 0
        first_cid = int(sp.cluster_ids[0])
        n_trials = sp.trial_ids[first_cid].size
        assert sp.tuning[first_cid].shape == (n_trials, 20)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _write_fake_cache_file(path: Path, *, prc_per_bin: float) -> None:
    """Write a minimal v7.3-style HDF5 that mimics the MATLAB cache layout.

    Three trial groups, three clusters per group, five trials per cluster,
    twenty bins per trial. All per-cluster arrays share the same content
    (this is just a parsing test, not a numerical test).
    """
    n_clusters = 3
    n_trials = 5
    n_bins = 20

    edges = np.linspace(0.0, 10.0, n_bins + 1, dtype=np.float64)
    centres = 0.5 * (edges[:-1] + edges[1:])
    tuning = (np.arange(n_trials * n_bins, dtype=np.float64) % 7.0).reshape(n_trials, n_bins)
    trial_ids = np.arange(n_trials, dtype=np.float64)
    stat_fr = np.full(n_trials, 1.5, dtype=np.float64)
    stat_time = np.full(n_trials, 4.0, dtype=np.float64)
    prc = np.array([[prc_per_bin]], dtype=np.float64)

    with h5py.File(path, "w") as f:
        refs = f.create_group("#refs#")
        tg_refs: list = []
        tc_refs: list = []
        for name in ("VT", "V", "T_Vstatic"):
            name_arr = np.array([ord(c) for c in name], dtype=np.uint16)
            name_ds = refs.create_dataset(f"name_{name}", data=name_arr)
            tg_refs.append(name_ds.ref)

            grp = refs.create_group(f"group_{name}")

            # Per-cluster underlying datasets — ref arrays of shape (n_clusters, 1)
            def _per_cluster(field_name: str, payload: np.ndarray):
                ds = refs.create_dataset(
                    f"{field_name}_{name}", data=payload,
                )
                ref_arr = np.array(
                    [[ds.ref] for _ in range(n_clusters)],
                    dtype=h5py.ref_dtype,
                )
                grp.create_dataset(field_name, data=ref_arr)

            _per_cluster("bin_edges", edges.reshape(-1, 1))
            _per_cluster("bin_centers", centres.reshape(-1, 1))  # MATLAB spelling
            _per_cluster("tuning", tuning.flatten().reshape(-1, 1))
            _per_cluster("trial_ids", trial_ids.reshape(-1, 1))
            _per_cluster("stationary_fr", stat_fr.reshape(-1, 1))
            _per_cluster("stationary_time", stat_time.reshape(-1, 1))
            _per_cluster("prc_per_bin", prc)

            # cluster_id needs distinct values per row (so the cache builds a
            # valid per-cluster dict). Write per-cluster scalar refs.
            cid_per_cluster_refs = []
            for j in range(n_clusters):
                cid_ds = refs.create_dataset(
                    f"cid_{name}_{j}", data=np.array([[100 + j]], dtype=np.float64),
                )
                cid_per_cluster_refs.append(cid_ds.ref)
            grp.create_dataset(
                "cluster_id",
                data=np.array(
                    [[r] for r in cid_per_cluster_refs], dtype=h5py.ref_dtype,
                ),
            )
            tc_refs.append(grp.ref)

        f.create_dataset(
            "trial_groups",
            data=np.array([[r] for r in tg_refs], dtype=h5py.ref_dtype),
        )
        f.create_dataset(
            "tuning_curves",
            data=np.array([[r] for r in tc_refs], dtype=h5py.ref_dtype),
        )
