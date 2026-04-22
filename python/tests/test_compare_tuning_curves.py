"""Prediction-space tuning-curve parity is basis-rotation invariant.

Two β vectors can produce the same tuning curve (B @ β) while having
different individual signs — this is exactly the identifiability the
Speed/TF raw-β sign check was tripping on. This test synthesises such a
pair and asserts the tuning-curve Pearson check passes while the raw-β
sign-agreement is low on the same input.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from rc2_glm.basis import raised_cosine_basis
from rc2_glm.compare import (
    _tuning_curve_pearson,
    _variable_betas_per_cluster,
)
from rc2_glm.config import GLMConfig


def _coef_df_for(cluster_id: int, variable: str, beta: np.ndarray) -> pd.DataFrame:
    """Build a glm_coefficients-shaped DataFrame for one cluster × variable."""
    rows = [
        {
            "probe_id": "TEST_probe",
            "cluster_id": cluster_id,
            "glm_type": "time",
            "model": "Selected",
            "coefficient": f"{variable}_{i + 1}",
            "estimate": float(beta[i]),
            "se": 0.1,
        }
        for i in range(beta.size)
    ]
    return pd.DataFrame(rows)


def test_tuning_curve_pearson_sees_through_raw_beta_disagreement():
    """Curve Pearson is high even when raw-β sign agreement is below threshold.

    The β pair is taken from the real Python/MATLAB fits of cluster 116 on
    ``CAA-1123243_rec1`` (Selected model, Speed bases). These are the
    coefficients the 02.7 prompt is scoped around — they give the fitter
    near-identical predictions while landing on different raw-β rotations
    because the 5 raised-cosine bases are correlated on the actual design
    matrix (Speed rows co-occur with Onset rows, share SF/OR dummies, etc).
    """
    config = GLMConfig()

    # Real Python Selected-model Speed β for cluster 116
    beta_py = np.array([-2.365849, -0.630662, -1.137915, -0.776475, -1.238769])
    # Real MATLAB Selected-model Speed β for the same cluster
    beta_ref = np.array([-1.375586,  0.225434, -0.275601, -0.139787, -0.296917])

    # Raw-β sign agreement is below the 0.85 gating threshold the old check
    # was tripping on — only index 1 flips, but that's enough for 4/5 = 0.8.
    sign_agree = float((np.sign(beta_py) == np.sign(beta_ref)).mean())
    assert sign_agree < 0.85, sign_agree

    coef_py = _coef_df_for(cluster_id=116, variable="Speed", beta=beta_py)
    coef_ref = _coef_df_for(cluster_id=116, variable="Speed", beta=beta_ref)

    row = _tuning_curve_pearson(coef_py, coef_ref, "Speed", config)
    assert row["check"] == "tuning_curve_pearson_Speed"
    # Curve Pearson r clears the per-cluster threshold (0.85) — the two
    # pipelines agree on shape even though raw-β sign does not.
    assert row["value"] > 0.85, row
    # And clears the mean-threshold since this is the only cluster in the
    # synthetic table.
    assert row["pass"] is True, row


def test_tuning_curve_pearson_fails_on_genuinely_different_curves():
    """Negative control: unrelated β vectors should yield low Pearson r."""
    config = GLMConfig()
    rng = np.random.default_rng(42)

    beta_py = np.array([1.0, 0.5, 0.2, -0.1, -0.4])
    # Pick βs that give a meaningfully different curve shape, not just a
    # different amplitude (which Pearson is invariant to).
    beta_ref = np.array([-0.4, -0.1, 0.2, 0.5, 1.0])   # reversed

    coef_py = _coef_df_for(cluster_id=1, variable="Speed", beta=beta_py)
    coef_ref = _coef_df_for(cluster_id=1, variable="Speed", beta=beta_ref)

    row = _tuning_curve_pearson(coef_py, coef_ref, "Speed", config)
    assert row["value"] < 0.5, row
    assert row["pass"] is False, row


def test_variable_betas_per_cluster_respects_expected_ordering():
    """Ensure β are returned in Speed_1..Speed_N order regardless of CSV row order."""
    beta = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    df = _coef_df_for(cluster_id=7, variable="Speed", beta=beta)
    # Shuffle the rows — _variable_betas_per_cluster must re-sort by name.
    df = df.sample(frac=1.0, random_state=0).reset_index(drop=True)

    out = _variable_betas_per_cluster(df, "Speed", "Selected", n_bases=5)
    key = ("TEST_probe", 7)
    assert key in out
    assert np.allclose(out[key], beta)
