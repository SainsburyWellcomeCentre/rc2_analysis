"""Protocol-dependent velocity selection and zero-phase filter tests.

Covers ``FormattedDataReader.trial_velocity`` and the private
``_filter_trace`` helper introduced to fix the silent-empty-bins bug on
passive (StageOnly / ReplayOnly) recordings — see prompt 02.5.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from rc2_formatted_data_reader import FormattedDataReader
from rc2_formatted_data_reader.reader import (
    PROTOCOL_VELOCITY_CHANNEL,
    _filter_trace,
)


# --------------------------------------------------------------------------- #
# Pure protocol-to-channel mapping
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "protocol, expected_channel",
    [
        ("Coupled", "filtered_teensy"),
        ("EncoderOnly", "filtered_teensy"),
        ("CoupledMismatch", "filtered_teensy_2"),
        ("EncoderOnlyMismatch", "filtered_teensy_2"),
        ("StageOnly", "stage"),
        ("ReplayOnly", "multiplexer_output"),
    ],
)
def test_protocol_velocity_channel_mapping(protocol: str, expected_channel: str) -> None:
    assert PROTOCOL_VELOCITY_CHANNEL[protocol] == expected_channel


# --------------------------------------------------------------------------- #
# Zero-phase filter (matches lib/fcn/general/filter_trace.m)
# --------------------------------------------------------------------------- #


def test_filter_trace_preserves_length() -> None:
    fs = 10_000.0
    x = np.random.default_rng(0).standard_normal(5000)
    y = _filter_trace(x, fs)
    assert y.shape == x.shape
    assert y.dtype == np.float64


def test_filter_trace_is_zero_phase_on_in_band_sine() -> None:
    """A 5 Hz sine with cutoff 50 Hz should come through with no phase shift.

    ``filtfilt`` passes the signal forward then backward; any linear phase
    introduced on the forward pass is cancelled on the backward pass. A phase
    shift here would break motion-mask edges (the mask is thresholded on
    |v|, so a phase shift of even a few samples would shift onset times).
    """
    fs = 10_000.0
    t = np.arange(1.0, step=1.0 / fs)  # 1 s
    x = np.sin(2.0 * np.pi * 5.0 * t)
    y = _filter_trace(x, fs, cutoff_hz=50.0, order=3)
    # In-band amplitude should come through ~unattenuated, with matching phase.
    # Skip a handful of edge samples to avoid filtfilt's padding artefacts.
    # Tolerance 5e-3 accommodates the ~0.5% ripple that order-3 Butterworth
    # introduces well below its cutoff; a phase shift would show up as a
    # larger sinusoidal-envelope error.
    edge = 200
    assert np.allclose(x[edge:-edge], y[edge:-edge], atol=5e-3)


def test_filter_trace_attenuates_out_of_band_noise() -> None:
    fs = 10_000.0
    t = np.arange(1.0, step=1.0 / fs)
    signal = np.sin(2.0 * np.pi * 5.0 * t)               # 5 Hz — in-band
    noise = 0.5 * np.sin(2.0 * np.pi * 1000.0 * t)       # 1 kHz — well above cutoff
    y = _filter_trace(signal + noise, fs, cutoff_hz=50.0, order=3)
    # Residual high-frequency content should be essentially gone.
    edge = 200
    residual = y[edge:-edge] - signal[edge:-edge]
    assert np.max(np.abs(residual)) < 0.05


def test_filter_trace_returns_unfiltered_for_tiny_arrays() -> None:
    """``filtfilt`` errors on arrays shorter than its padding; we fall back."""
    fs = 10_000.0
    tiny = np.array([1.0, 2.0, 3.0])
    out = _filter_trace(tiny, fs, order=3)
    assert np.array_equal(out, tiny)
    assert out.dtype == np.float64


# --------------------------------------------------------------------------- #
# Synthetic-HDF5 round-trip: trial_velocity picks the right channel
# --------------------------------------------------------------------------- #


def _write_synthetic_mat(path: Path, protocol: str, fs: float = 10_000.0) -> None:
    """Create a minimal HDF5 with one trial and distinct signals per channel.

    Each candidate velocity channel (``filtered_teensy``, ``filtered_teensy_2``,
    ``stage``, ``multiplexer_output``) gets a unique low-frequency sine so the
    test can identify which one ``trial_velocity`` returned.
    """
    n = int(fs)  # 1 s trial
    t = np.arange(n) / fs

    per_channel = {
        "filtered_teensy": 1.0 * np.sin(2.0 * np.pi * 2.0 * t),
        "filtered_teensy_2": 2.0 * np.sin(2.0 * np.pi * 3.0 * t),
        "stage": 3.0 * np.sin(2.0 * np.pi * 4.0 * t),
        "multiplexer_output": 4.0 * np.sin(2.0 * np.pi * 5.0 * t),
    }

    with h5py.File(path, "w") as f:
        # MATLAB v7.3 stores strings as character codes; the reader uses
        # ``_deref_string`` / ``_hdf5_string`` which maps ``chr(int(c))``.
        f.create_dataset(
            "probe_id",
            data=np.array([ord(c) for c in "synthetic"], dtype=np.uint16).reshape(-1, 1),
        )

        sess = f.create_group("sessions")
        sess.create_dataset("fs", data=np.array([[fs]], dtype=np.float64))
        sess.create_dataset("n_samples", data=np.array([[n]], dtype=np.int64))
        sess.create_dataset("n_trials", data=np.array([[1]], dtype=np.int64))
        for name, trace in per_channel.items():
            sess.create_dataset(name, data=trace.reshape(1, -1))
        sess.create_dataset("probe_t", data=t.reshape(1, -1))
        sess.create_dataset("solenoid", data=np.zeros((1, n)))

        # MATLAB stores per-trial scalars and strings as object refs. Recreate
        # that structure so the reader's deref path exercises the real code.
        scratch = f.create_group("__scratch")
        start_d = scratch.create_dataset("start", data=np.array([[0]], dtype=np.int64))
        end_d = scratch.create_dataset("end", data=np.array([[n]], dtype=np.int64))
        tid_d = scratch.create_dataset("tid", data=np.array([[1]], dtype=np.int64))
        proto_codes = np.array([ord(c) for c in protocol], dtype=np.uint16).reshape(-1, 1)
        proto_d = scratch.create_dataset("protocol", data=proto_codes)
        cfg_grp = scratch.create_group("cfg")
        cfg_grp.create_dataset("trial_sequence", data=np.array([[1]], dtype=np.int64))
        cfg_grp.create_dataset("start_pos", data=np.array([[0.0]], dtype=np.float64))
        cfg_grp.create_dataset("forward_limit", data=np.array([[0.0]], dtype=np.float64))

        ref_dtype = h5py.ref_dtype
        trials = sess.create_group("trials")
        trials.create_dataset("start_idx", data=np.array([[start_d.ref]], dtype=ref_dtype))
        trials.create_dataset("end_idx", data=np.array([[end_d.ref]], dtype=ref_dtype))
        trials.create_dataset("trial_id", data=np.array([[tid_d.ref]], dtype=ref_dtype))
        trials.create_dataset("protocol", data=np.array([[proto_d.ref]], dtype=ref_dtype))
        trials.create_dataset("config", data=np.array([[cfg_grp.ref]], dtype=ref_dtype))

        # The cluster dataset is never read by trial_velocity, but reader
        # consumers sometimes touch ``n_clusters`` — provide an empty stub.
        clusters = f.create_group("clusters")
        clusters.create_dataset("id", data=np.zeros((0, 1), dtype=ref_dtype))


@pytest.mark.parametrize(
    "protocol, expected_channel",
    [
        ("Coupled", "filtered_teensy"),
        ("CoupledMismatch", "filtered_teensy_2"),
        ("StageOnly", "stage"),
        ("ReplayOnly", "multiplexer_output"),
    ],
)
def test_trial_velocity_returns_protocol_channel(
    tmp_path: Path, protocol: str, expected_channel: str,
) -> None:
    mat_path = tmp_path / f"{protocol.lower()}.mat"
    _write_synthetic_mat(mat_path, protocol)

    with FormattedDataReader(mat_path) as r:
        assert r.trial_velocity_channel(0) == expected_channel

        # Unfiltered path should exactly reproduce the channel we wrote.
        raw = r.trial_velocity(0, apply_filter=False)
        with h5py.File(mat_path, "r") as f:
            expected = np.asarray(f["sessions"][expected_channel][0, :])
        assert raw.shape == expected.shape
        assert np.allclose(raw, expected)

        # Filtered path should not time-shift the in-band signal. Order-3
        # Butterworth introduces ~0.5% ripple on in-band sinusoids — hence
        # the rtol; a phase shift of even a millisecond would blow this up.
        filt = r.trial_velocity(0, apply_filter=True, cutoff_hz=50.0, filter_order=3)
        edge = 200
        assert np.allclose(
            filt[edge:-edge], expected[edge:-edge], rtol=1e-2, atol=5e-3
        )


def test_unknown_protocol_raises(tmp_path: Path) -> None:
    mat_path = tmp_path / "unknown.mat"
    _write_synthetic_mat(mat_path, protocol="NotAProtocol")
    with FormattedDataReader(mat_path) as r:
        with pytest.raises(ValueError, match="Unknown protocol"):
            r.trial_velocity_channel(0)
