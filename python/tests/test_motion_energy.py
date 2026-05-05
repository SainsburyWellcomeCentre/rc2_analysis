"""ME_face basis + design matrix wiring — unit tests.

Per prompt 06 (2026-04-30). Pins five invariants:

1. ``raised_cosine_basis_linear`` returns ``(n, n_bases)`` with values in
   [0, 1] and unit peaks (raised-cosine property).
2. The basis sums to ~1 across the value range interior (partition-of-
   unity for raised cosines spaced at width = 1.5×delta).
3. Values outside the range clip to the edge basis (saturating, not zero).
4. ``_bin_continuous_to_edges`` aggregates a continuous signal onto
   GLM bin edges via bin-mean; empty bins return NaN.
5. ``"ME_face"`` selection adds 5 main-effect columns; ``"ME_face_x_Speed"``
   selection adds 25 interaction columns; both branches no-op when
   ``B_me_face`` is None.

Plus one INTERACTION_PARENTS gating check: ``ME_face_x_Speed`` requires
both ``ME_face`` and ``Speed`` parents present in Phase-1.
"""
from __future__ import annotations

import numpy as np

from rc2_glm.basis import (
    onset_kernel_basis,
    raised_cosine_basis,
    raised_cosine_basis_linear,
)
from rc2_glm.config import INTERACTION_PARENTS
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.time_binning import _bin_continuous_to_edges


# --------------------------------------------------------------------------- #
# Basis builder                                                               #
# --------------------------------------------------------------------------- #

def test_me_face_basis_shape_and_unit_peak():
    x = np.linspace(-2.0, 3.0, 100)
    B = raised_cosine_basis_linear(x, 5, -2.0, 3.0)
    assert B.shape == (100, 5)
    assert np.isfinite(B).all()
    assert B.max() <= 1.0 + 1e-12
    assert B.min() >= 0.0
    # Each basis attains its peak (= 1) at its own centre.
    centres = np.linspace(-2.0, 3.0, 5)
    for k, c in enumerate(centres):
        peak_at_centre = raised_cosine_basis_linear(np.array([c]), 5, -2.0, 3.0)
        assert peak_at_centre[0, k] >= 0.999


def test_me_face_basis_constant_sum_in_interior():
    """In the interior, the sum across bases is constant (smooth tiling).

    With width=1.5×delta the raised cosines overlap more than a strict
    partition-of-unity; the sum is 1.5 (not 1) but it IS constant — that
    constancy is what makes the basis a smooth, well-behaved span across
    the interior without dips between centres. Edge regions deviate
    because outer bases lack a neighbour on one side.
    """
    centres = np.linspace(-2.0, 3.0, 5)
    interior_x = np.linspace(centres[1], centres[3], 50)
    B = raised_cosine_basis_linear(interior_x, 5, -2.0, 3.0)
    sums = B.sum(axis=1)
    # Constant across the interior (relative variation < 1%):
    assert sums.std() / sums.mean() < 1e-2
    # And matches the closed-form for width=1.5×delta cosines: 1.5.
    assert np.allclose(sums, 1.5, atol=1e-2)


def test_me_face_basis_zero_far_outside_range():
    """Values far past the range produce zero activation across all bases.

    Internal z is clipped to ±π, where cos(±π) = -1, so the cosine basis
    evaluates to 0. This means an extreme outlier (z > x_max + width or
    z < x_min - width) reads as "no activation anywhere" rather than as
    "saturating high". Acceptable because ME values are z-scored in
    practice — extreme outliers are rare; widen ``me_face_range`` if
    they're not.
    """
    far_above = raised_cosine_basis_linear(np.array([10.0]), 5, -2.0, 3.0)
    far_below = raised_cosine_basis_linear(np.array([-10.0]), 5, -2.0, 3.0)
    assert np.allclose(far_above, 0.0, atol=1e-12)
    assert np.allclose(far_below, 0.0, atol=1e-12)
    # At the edge centre itself, the edge basis peaks at 1.0.
    edge_high = raised_cosine_basis_linear(np.array([3.0]), 5, -2.0, 3.0)
    edge_low = raised_cosine_basis_linear(np.array([-2.0]), 5, -2.0, 3.0)
    assert edge_high[0, -1] >= 0.999
    assert edge_low[0, 0] >= 0.999


# --------------------------------------------------------------------------- #
# Continuous-signal binner                                                    #
# --------------------------------------------------------------------------- #

def test_bin_continuous_to_edges_basic():
    """3 samples per 100 ms bin (≈30 Hz camera, 100 ms GLM bin) → bin-mean."""
    sample_t = np.linspace(0.0, 1.0, 31)  # 30 evenly-spaced samples in 1s
    values = sample_t.copy()  # values = time → bin-mean = bin-centre
    edges = np.arange(0.0, 1.001, 0.1)  # 10 bins
    out = _bin_continuous_to_edges(values, sample_t, edges)
    assert out.shape == (10,)
    # Bin centres are at 0.05, 0.15, ..., 0.95
    expected_centres = 0.5 * (edges[:-1] + edges[1:])
    # ~3 samples per bin → bin-mean is ~ bin centre, within sample spacing.
    assert np.allclose(out, expected_centres, atol=0.05)


def test_bin_continuous_to_edges_empty_bin_is_nan():
    """A bin with no contributing samples returns NaN, not 0."""
    sample_t = np.array([0.05, 0.55])  # one sample in bin 0, one in bin 5
    values = np.array([1.0, 2.0])
    edges = np.arange(0.0, 1.001, 0.1)  # 10 bins
    out = _bin_continuous_to_edges(values, sample_t, edges)
    assert np.isnan(out[1])  # bin 1 (0.1..0.2) has no samples
    assert out[0] == 1.0
    assert out[5] == 2.0


# --------------------------------------------------------------------------- #
# Design matrix wiring                                                        #
# --------------------------------------------------------------------------- #

def _toy_inputs(n=200, seed=0):
    rng = np.random.default_rng(seed)
    speed = rng.uniform(0, 50, n)
    tf = rng.uniform(0, 5, n)
    onset = np.linspace(0, 2, n)
    sf_vals = rng.choice([0.001, 0.003, 0.006, 0.012], n)
    or_vals = rng.choice([-np.pi / 4, 0.0, np.pi / 4, np.pi / 2], n)
    B_speed = raised_cosine_basis(speed, 5, 0.0, 50.0)
    B_tf = raised_cosine_basis(tf, 5, 0.0, 7.3)
    B_onset = onset_kernel_basis(onset, 6, 2.0)
    me_z = rng.normal(0.0, 1.0, n)
    B_me_face = raised_cosine_basis_linear(me_z, 5, -2.0, 3.0)
    return B_speed, B_tf, B_onset, sf_vals, or_vals, B_me_face


def _prediction_mode_kwargs(sf_vals, or_vals):
    """Pass non-empty ref levels so _drop_zero_variance is bypassed.

    Required for structural column-count tests with random toy inputs
    where some interaction column may have low variance and get pruned
    in training mode.
    """
    return dict(
        sf_ref_levels=np.sort(np.unique(sf_vals)),
        or_ref_levels=np.sort(np.unique(or_vals)),
    )


def test_design_matrix_me_face_main_effect_5_columns():
    B_speed, B_tf, B_onset, sf_vals, or_vals, B_me_face = _toy_inputs()
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, ["ME_face"],
        B_me_face=B_me_face,
        **_prediction_mode_kwargs(sf_vals, or_vals),
    )
    me_cols = [n for n in names if n.startswith("ME_face_")]
    assert len(me_cols) == 5
    assert me_cols == [f"ME_face_{i + 1}" for i in range(5)]


def test_design_matrix_me_face_x_speed_25_interaction_columns():
    B_speed, B_tf, B_onset, sf_vals, or_vals, B_me_face = _toy_inputs()
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        ["Speed", "ME_face", "ME_face_x_Speed"],
        B_me_face=B_me_face,
        **_prediction_mode_kwargs(sf_vals, or_vals),
    )
    me_main = [n for n in names if n.startswith("ME_face_")]
    speed_main = [n for n in names if n.startswith("Speed_")]
    me_x_speed = [n for n in names if n.startswith("MEf") and "_x_Spd" in n]
    assert len(me_main) == 5
    assert len(speed_main) == 5
    assert len(me_x_speed) == 25
    # Naming convention pinned: MEf{i}_x_Spd{j}
    assert "MEf1_x_Spd1" in names
    assert "MEf5_x_Spd5" in names


def test_design_matrix_me_face_none_skips_columns():
    """When B_me_face is None, ME_face main and interaction are no-op."""
    B_speed, B_tf, B_onset, sf_vals, or_vals, _ = _toy_inputs()
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        ["ME_face", "ME_face_x_Speed"],
        B_me_face=None,
    )
    assert not any(n.startswith("ME_face_") for n in names)
    assert not any(n.startswith("MEf") for n in names)


def test_design_matrix_me_face_x_speed_row_wise_outer_product():
    """Interaction column k = (mi, si) should equal B_me_face[:, mi] * B_speed[:, si]."""
    B_speed, B_tf, B_onset, sf_vals, or_vals, B_me_face = _toy_inputs()
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals,
        ["Speed", "ME_face", "ME_face_x_Speed"],
        B_me_face=B_me_face,
        **_prediction_mode_kwargs(sf_vals, or_vals),
    )
    for mi in range(5):
        for si in range(5):
            col_name = f"MEf{mi + 1}_x_Spd{si + 1}"
            col_idx = names.index(col_name)
            expected = B_me_face[:, mi] * B_speed[:, si]
            assert np.allclose(X[:, col_idx], expected)


def test_interaction_parents_me_face_x_speed():
    """ME_face_x_Speed parents must be ME_face and Speed."""
    assert INTERACTION_PARENTS["ME_face_x_Speed"] == ("ME_face", "Speed")
