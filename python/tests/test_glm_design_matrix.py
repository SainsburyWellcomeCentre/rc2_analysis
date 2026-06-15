"""Design matrix structural tests."""

from __future__ import annotations

import numpy as np

from rc2_glm.basis import (
    onset_kernel_basis,
    raised_cosine_basis,
    raised_cosine_basis_linear,
)
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


def test_me_face_main_effect_in_selected():
    """ME_face main-effect adds 5 raised-cosine columns when basis is supplied."""
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    rng = np.random.default_rng(42)
    me_z = rng.normal(0.0, 1.0, B_speed.shape[0])
    B_me_face = raised_cosine_basis_linear(me_z, 5, -2.0, 3.0)
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, ["ME_face"],
        B_me_face=B_me_face,
    )
    assert sum(1 for n in names if n.startswith("ME_face_")) == 5


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


# --------------------------------------------------------------------------- #
# Acceleration / ME_face full-pairwise interaction set (2026-06-15)
# --------------------------------------------------------------------------- #

_NEW_PAIRWISE = (
    "Speed_x_Acceleration", "TF_x_ME_face", "TF_x_Acceleration",
    "SF_x_ME_face", "SF_x_Acceleration", "OR_x_ME_face",
    "OR_x_Acceleration", "ME_face_x_Acceleration",
)


def _me_accel_bases(n, seed=7):
    rng = np.random.default_rng(seed)
    me_z = rng.normal(0.0, 1.0, n)
    acc_z = rng.normal(0.0, 1.0, n)
    B_me_face = raised_cosine_basis_linear(me_z, 5, -2.0, 3.0)
    B_accel = raised_cosine_basis_linear(acc_z, 5, -3.0, 3.0)
    return B_me_face, B_accel


def test_new_pairwise_interactions_built_when_bases_present():
    """Each of the 8 new interactions produces its own correctly-named
    product columns when the value bases are supplied (prediction mode so no
    zero-variance pruning hides a column)."""
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    n = B_speed.shape[0]
    B_me_face, B_accel = _me_accel_bases(n)
    sf_ref = np.array([0.001, 0.003, 0.006, 0.012])
    or_ref = np.array([-np.pi / 4, 0.0, np.pi / 4, np.pi / 2])
    n_sf = sf_ref.size - 1   # 3 dummies
    n_or = or_ref.size - 1   # 3 dummies
    _, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, list(_NEW_PAIRWISE),
        sf_ref_levels=sf_ref, or_ref_levels=or_ref,
        B_me_face=B_me_face, B_accel=B_accel,
    )
    counts = {
        "Spd_x_Acc": sum(n.startswith("Spd") and "_x_Acc" in n for n in names),
        "TF_x_MEf": sum(n.startswith("TF") and "_x_MEf" in n for n in names),
        "TF_x_Acc": sum(n.startswith("TF") and "_x_Acc" in n for n in names),
        "SF_x_MEf": sum(n.startswith("SF") and "_x_MEf" in n for n in names),
        "SF_x_Acc": sum(n.startswith("SF") and "_x_Acc" in n for n in names),
        "OR_x_MEf": sum(n.startswith("OR") and "_x_MEf" in n for n in names),
        "OR_x_Acc": sum(n.startswith("OR") and "_x_Acc" in n for n in names),
        "MEf_x_Acc": sum(n.startswith("MEf") and "_x_Acc" in n for n in names),
    }
    assert counts["Spd_x_Acc"] == 5 * 5
    assert counts["TF_x_MEf"] == 5 * 5
    assert counts["TF_x_Acc"] == 5 * 5
    assert counts["MEf_x_Acc"] == 5 * 5
    assert counts["SF_x_MEf"] == n_sf * 5
    assert counts["SF_x_Acc"] == n_sf * 5
    assert counts["OR_x_MEf"] == n_or * 5
    assert counts["OR_x_Acc"] == n_or * 5


def test_new_interactions_inert_without_value_bases():
    """Selecting the new interactions with no ME / accel basis (the token /
    screens path) must add NO columns — keeps those runs byte-identical."""
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    _, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf_vals, or_vals, list(_NEW_PAIRWISE),
        B_me_face=None, B_accel=None,
    )
    assert not any("_x_Acc" in n for n in names)
    assert not any("_x_MEf" in n for n in names)


def test_full_interaction_label_inert_without_me_accel():
    """The FullInteraction ceiling lists the new interactions, but with no
    ME/accel basis it must reduce to the legacy stimulus-only design."""
    B_speed, B_tf, B_onset, sf_vals, or_vals = _toy_inputs()
    _, names = assemble_design_matrix(
        B_speed, B_tf, B_onset, sf_vals, or_vals, "FullInteraction",
    )
    assert not any("_x_Acc" in n for n in names)
    assert not any("_x_MEf" in n for n in names)
    assert not any(n.startswith("Acceleration_") for n in names)
