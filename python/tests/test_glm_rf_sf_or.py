"""RF-local SF/OR path (config sf_or_source="rf_local").

Test names read as the lessons (Pattern 8):
- the rf_local plumbing must NOT perturb the default token design matrix;
- OR is encoded π-periodically (0° ≡ 180°);
- continuous SF/OR produce SF_n / OR_n columns, and the main-effect name
  filter must exclude the SF_i_x_OR_j interaction columns;
- the trial→cloud join keys on (theta, sf, VX), so it is robust to the BV
  velocity-bandwidth token differing between the rendered-frame folders
  (BV0p100) and the presentation metadata (BV0p200).
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from rc2_glm.basis import circular_basis, raised_cosine_basis_linear
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.rf_sf_or import cloud_key, load_rf_sf_or

_COHORT = os.path.expanduser(
    "~/local_data/motion_clouds/saved_goggles/_extract/cohort"
)


def _toy_bases(n=300, seed=0):
    rng = np.random.default_rng(seed)
    return (
        rng.random((n, 5)),  # B_speed
        rng.random((n, 5)),  # B_tf
        rng.random((n, 6)),  # B_onset
    )


def test_rf_local_kwargs_do_not_perturb_token_design():
    """B_sf=B_or=None (token mode) is byte-identical to omitting them."""
    B_speed, B_tf, B_onset = _toy_bases()
    n = B_speed.shape[0]
    rng = np.random.default_rng(1)
    sf = rng.choice([np.nan, 0.003, 0.006, 0.012], n)
    orr = rng.choice([np.nan, -0.785, 0.0, 0.785, 1.571], n)
    sel = ["Speed", "TF", "SF", "OR", "Speed_x_SF", "SF_x_OR"]
    refs = dict(sf_ref_levels=[0.003, 0.006, 0.012],
                or_ref_levels=[-0.785, 0.0, 0.785, 1.571])
    X0, n0 = assemble_design_matrix_selected(B_speed, B_tf, B_onset, sf, orr, sel, **refs)
    X1, n1 = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf, orr, sel, **refs, B_sf=None, B_or=None
    )
    assert n0 == n1
    assert np.array_equal(X0, X1)
    # Categorical token names preserved exactly.
    assert "SF_0.0060" in n0 and "Spd1_x_SF0.0060" in n0


def test_or_circular_basis_is_pi_periodic():
    """Orientation basis treats 0° and 180° as identical (π-periodic)."""
    B0 = circular_basis(np.array([0.0]), 6)
    B180 = circular_basis(np.array([180.0]), 6)
    assert np.allclose(B0, B180)
    assert np.allclose(circular_basis(np.array([20.0]), 6),
                       circular_basis(np.array([200.0]), 6))
    B = circular_basis(np.linspace(0, 180, 50), 6)
    assert B.shape == (50, 6)
    assert B.max() <= 1.0 + 1e-9 and B.min() >= 0.0


def test_continuous_sf_or_columns_and_interaction_filtering():
    """rf_local design names SF_n / OR_n and keeps SF_i_x_OR_j separate."""
    B_speed, B_tf, B_onset = _toy_bases()
    n = B_speed.shape[0]
    rng = np.random.default_rng(2)
    sf = np.clip(rng.normal(0.05, 0.02, n), 0.02, 0.14)
    orr = rng.uniform(0, 180, n)
    B_sf = raised_cosine_basis_linear(sf, 4, 0.02, 0.15)
    B_or = circular_basis(orr, 6)
    sel = ["SF", "OR", "Speed_x_SF", "SF_x_OR"]
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf, orr, sel,
        sf_ref_levels=[], or_ref_levels=[], B_sf=B_sf, B_or=B_or,
    )
    main_sf = [c for c in names if c in {f"SF_{i}" for i in range(1, 5)}]
    main_or = [c for c in names if c in {f"OR_{i}" for i in range(1, 7)}]
    assert len(main_sf) == 4 and len(main_or) == 6
    assert "SF_1_x_OR_1" in names and "Spd1_x_SF_1" in names
    # The SF_i_x_OR_j columns must NOT be mistaken for SF main effects.
    assert not any(c.startswith("SF_") and "_x_" in c for c in main_sf)
    assert np.isfinite(X).all()


def test_sf_basis_zero_where_undefined():
    """SF/OR bases are zero on rows where the value is NaN (grey screen)."""
    sf = np.array([0.03, np.nan, 0.05])
    fin = np.isfinite(sf)
    B = np.zeros((3, 4))
    B[fin] = raised_cosine_basis_linear(sf[fin], 4, 0.02, 0.15)
    assert np.allclose(B[1], 0.0)
    assert B[0].sum() > 0 and B[2].sum() > 0


@pytest.mark.skipif(not os.path.isdir(_COHORT), reason="goggle cohort parquet not local")
def test_trial_cloud_join_is_bv_agnostic():
    """A BV0p200 metadata cloud name retrieves the BV0p100 rendered-frame data."""
    lookup = load_rf_sf_or(_COHORT, "CAA-1124370_rec1_rec2_rec3")
    assert len(lookup.clusters) > 0
    cid = sorted(lookup.clusters)[0]
    cloud_bv200 = "theta0p000_Btheta0p785_sf00p016_Bsf0p005_VX0p191_BV0p200"
    cloud_bv100 = cloud_bv200.replace("BV0p200", "BV0p100")
    got200 = lookup.get(cid, cloud_bv200)
    got100 = lookup.get(cid, cloud_bv100)
    assert got200 is not None
    assert got100 is not None
    assert np.array_equal(got200[0], got100[0])  # same SF(frame) array
    # The join key ignores BV.
    assert cloud_key(cloud_bv200) == cloud_key(cloud_bv100)
