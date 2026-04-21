"""Design matrix structural tests."""

from __future__ import annotations

import numpy as np

from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis
from rc2_glm.design_matrix import (
    assemble_design_matrix,
    assemble_design_matrix_selected,
)


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
    return B_speed, B_tf, B_onset, sf_vals, or_vals


def test_null_model_has_intercept_and_onset_only():
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    X, names = assemble_design_matrix(
        B_speed, B_tf, B_onset, sf_vals, or_vals, "Null"
    )
    assert names[0] == "Intercept"
    assert (X[:, 0] == 1.0).all()
    assert sum(1 for n in names if n.startswith("Onset_")) == B_onset.shape[1]


def test_additive_model_includes_all_main_effects():
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    X, names = assemble_design_matrix(
        B_speed, B_tf, B_onset, sf_vals, or_vals, "Additive"
    )
    assert any(n.startswith("Speed_") for n in names)
    assert any(n.startswith("TF_") for n in names)
    assert any(n.startswith("SF_") for n in names)
    assert any(n.startswith("OR_") for n in names)


def test_full_interaction_model_has_interactions():
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    X, names = assemble_design_matrix(
        B_speed, B_tf, B_onset, sf_vals, or_vals, "FullInteraction"
    )
    assert any("_x_" in n or "Spd" in n for n in names)


def test_sf_reference_coding_drops_first_level():
    """SF dummies should be (n_unique - 1)."""
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, ["SF"]
    )
    n_sf_dummies = sum(1 for n in names if n.startswith("SF_"))
    assert n_sf_dummies == len(np.unique(sf_vals)) - 1


def test_zero_variance_columns_dropped_during_training():
    """Columns of all-zero (other than intercept) should be removed."""
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    sf_vals_one_level = np.full_like(sf_vals, 0.003)
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals_one_level, or_vals, ["SF"]
    )
    assert not any(n.startswith("SF_") for n in names)


def test_prediction_mode_keeps_zero_variance_columns():
    """If sf_ref_levels is supplied, no zero-variance pruning."""
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    sf_vals_one_level = np.full_like(sf_vals, 0.003)
    sf_ref = np.array([0.001, 0.003, 0.006, 0.012])
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals_one_level, or_vals, ["SF"],
        sf_ref_levels=sf_ref,
    )
    assert sum(1 for n in names if n.startswith("SF_")) == sf_ref.size - 1
