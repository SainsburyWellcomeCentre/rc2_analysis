"""Integration tests for FormattedDataReader against a real .mat file."""

from __future__ import annotations

import re

import h5py
import numpy as np
import pytest

from rc2_formatted_data_reader import (
    CONDITION_MAP,
    FormattedDataReader,
    StimulusLookup,
    masks,
    parse_cloud_name,
)


# --- Reader: top-level metadata ---


def test_open_and_close(formatted_mat_path):
    reader = FormattedDataReader(formatted_mat_path)
    assert reader.probe_id
    reader.close()


def test_context_manager(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        assert r.fs > 0
        assert r.n_clusters > 0
        assert r.n_trials > 0


def test_probe_id_matches_filename(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        assert r.probe_id == formatted_mat_path.stem


def test_n_clusters_matches_dataset_length(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r, \
         h5py.File(formatted_mat_path, "r") as f:
        assert r.n_clusters == f["clusters"]["id"].shape[0]


# --- Reader: clusters ---


def test_cluster_ids_shape(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        ids = r.cluster_ids()
        assert ids.shape == (r.n_clusters,)
        assert ids.dtype.kind in ("i", "u")


def test_visp_filtering(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        idx = r.visp_cluster_indices()
        regions = r.cluster_regions()
        for i in idx:
            assert re.match(r"VISp", regions[i])
        for i, region in enumerate(regions):
            if i not in idx:
                assert not re.match(r"VISp", region)


def test_spike_times_for_one_cluster(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        st = r.spike_times(0)
        assert st.ndim == 1
        assert st.dtype == np.float64
        assert st.size > 0
        assert (st >= 0).all()
        assert np.all(np.diff(st) >= 0)  # sorted


def test_visp_clusters_have_visp_region(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        idx = r.visp_cluster_indices()
        if idx.size == 0:
            pytest.skip("No VISp clusters in this file")
        for i in idx[: min(5, len(idx))]:
            assert re.match(r"VISp", r.cluster_region(int(i)))


# --- Reader: trials and session arrays ---


def test_trial_bounds_and_slicing(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        s, e = r.trial_bounds(0)
        assert e > s
        pt = r.probe_t(s, e)
        v = r.velocity(s, e)
        assert pt.shape == v.shape == (e - s,)


def test_trial_condition_mapping(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        # Sample first / middle / last trial
        n = r.n_trials
        for ti in [0, n // 2, n - 1]:
            cond = r.trial_condition(ti)
            assert cond in CONDITION_MAP.values()


def test_trial_protocol_is_string(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        proto = r.trial_protocol(0)
        assert isinstance(proto, str)
        assert len(proto) > 0


def test_session_arrays_lazy(formatted_mat_path):
    """Slicing session arrays must not load the full array into memory."""
    with FormattedDataReader(formatted_mat_path) as r:
        s, e = r.trial_bounds(0)
        v_slice = r.velocity(s, e)
        full_n = r.n_samples
        assert v_slice.size == e - s
        assert v_slice.size < full_n  # we never materialised the full array


# --- Reader: motion / stationary masks ---


def test_motion_and_stationary_masks_are_complementary(formatted_mat_path):
    with FormattedDataReader(formatted_mat_path) as r:
        s, e = r.trial_bounds(0)
        m = r.motion_mask(s, e, threshold=1.0)
        st = r.stationary_mask(s, e, threshold=1.0)
        assert m.shape == st.shape == (e - s,)
        assert (m == ~st).all()
        assert m.dtype == bool


# --- masks module pure functions ---


def test_treadmill_motion_mask_short_stationary_promoted_to_motion():
    fs = 1000.0
    # 1 s motion, 0.05 s stationary (< min_dur 0.2 s), 1 s motion
    velocity = np.concatenate([
        np.full(1000, 5.0),
        np.zeros(50),
        np.full(1000, 5.0),
    ])
    accel = np.zeros_like(velocity)
    mask = masks.treadmill_motion_mask(velocity, accel, fs)
    # The short stationary bout is promoted to motion → all True
    assert mask.all()


def test_treadmill_motion_mask_long_stationary_kept_stationary():
    fs = 1000.0
    velocity = np.concatenate([
        np.full(1000, 5.0),
        np.zeros(500),  # 0.5 s > min_dur 0.2 s
        np.full(1000, 5.0),
    ])
    accel = np.zeros_like(velocity)
    mask = masks.treadmill_motion_mask(velocity, accel, fs)
    assert mask[:1000].all()
    assert (~mask[1000:1500]).all()
    assert mask[1500:].all()


# --- Stimulus lookup ---


def test_parse_cloud_name_excluded():
    name = "theta0p000_Btheta3p142_sf00p006_Bsf0p004_VX0p000_BV2p000"
    p = parse_cloud_name(name)
    assert p["excluded"] is True
    assert np.isnan(p["sf"])


def test_parse_cloud_name_orientation_at_start():
    # 'Btheta0p785' is bandwidth, not orientation — orientation must come
    # from the leading 'theta-0p785' token only.
    name = "theta-0p785_Btheta0p785_sf00p003_Bsf0p002_VX1p002_BV0p200"
    p = parse_cloud_name(name)
    assert p["sf"] == 0.003
    assert np.isclose(p["orientation"], -np.pi / 4)
    assert np.isclose(p["vx"], 1.002)
    assert np.isclose(p["batch_gain"], 1 / 30)
    assert p["excluded"] is False


def test_stimulus_lookup_on_real_files(mc_sequence_path, mc_folders_path):
    lookup = StimulusLookup(mc_sequence_path, mc_folders_path)
    assert lookup.n_clouds > 0
    assert lookup.n_trials > 0
    name = lookup.cloud_name(1)
    assert name is not None
    assert name == lookup.cloud_name(1)  # idempotent
    params = lookup.stimulus_params(1)
    assert "sf" in params
    assert "orientation" in params
    assert "batch_gain" in params
    assert "excluded" in params
    assert params["cloud_name"] == name


def test_stimulus_lookup_aligns_with_reader_trials(
    formatted_mat_path, mc_sequence_path, mc_folders_path
):
    lookup = StimulusLookup(mc_sequence_path, mc_folders_path)
    with FormattedDataReader(formatted_mat_path) as r:
        # Each reader trial_id must be a valid 1-based index into the
        # presentation sequence.
        for ti in [0, r.n_trials // 2, r.n_trials - 1]:
            tid = r.trial_id(ti)
            assert 1 <= tid <= lookup.n_trials
