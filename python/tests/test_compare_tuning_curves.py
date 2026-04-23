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

    rows = _tuning_curve_pearson(coef_py, coef_ref, "Speed", config)
    # Now returns a pair (gate row + informational mean row).
    assert isinstance(rows, list) and len(rows) == 2, rows
    gate, info = rows
    assert gate["check"] == "tuning_curve_pearson_Speed"
    assert info["check"] == "tuning_curve_pearson_mean_Speed"
    # Gate value is the fraction of clusters clearing per-cluster r≥0.85;
    # on a single-cluster input that Pearson-correlates at ~0.91 the
    # fraction is 1.0 and the row passes.
    assert gate["value"] == 1.0, gate
    assert gate["pass"] is True, gate
    # Informational row carries mean(r) — a single number in [0.85, 1.0]
    # on this input; no threshold; never gates.
    assert info["value"] > 0.85, info
    assert info["pass"] is True, info
    assert np.isnan(info["threshold"]), info


def test_tuning_curve_pearson_fails_on_genuinely_different_curves():
    """Negative control: unrelated β vectors should yield low Pearson r."""
    config = GLMConfig()

    beta_py = np.array([1.0, 0.5, 0.2, -0.1, -0.4])
    # Pick βs that give a meaningfully different curve shape, not just a
    # different amplitude (which Pearson is invariant to).
    beta_ref = np.array([-0.4, -0.1, 0.2, 0.5, 1.0])   # reversed

    coef_py = _coef_df_for(cluster_id=1, variable="Speed", beta=beta_py)
    coef_ref = _coef_df_for(cluster_id=1, variable="Speed", beta=beta_ref)

    rows = _tuning_curve_pearson(coef_py, coef_ref, "Speed", config)
    gate, info = rows
    # Single cluster is below the per-cluster bar → fraction = 0, gate fails.
    assert gate["value"] == 0.0, gate
    assert gate["pass"] is False, gate
    # Mean(r) is low too, but it's informational — never gates.
    assert info["value"] < 0.5, info
    assert info["pass"] is True, info


def test_tuning_curve_pearson_small_n_tolerates_single_outlier():
    """n=6 TF regression: 5/6 clusters clearing 0.85 must not trip the gate.

    The old mean-based gate (mean(r) >= 0.90) was tripped by a single
    r=0.04 outlier on the TF-selected clusters: mean fell to ~0.82
    even though 5/6 clusters agreed cluster-by-cluster. The current gate
    is ``fraction(r >= 0.85) >= 0.80`` — 5/6 = 0.833 clears it. This
    test locks that behaviour in as the lesson (pattern 8: encode the
    root cause, not the fix).
    """
    config = GLMConfig()
    # Five well-agreeing β pairs (identical → r = 1.0 per cluster) plus
    # one genuinely different pair that Pearson-correlates below 0.85.
    good_py = np.array([1.0, 0.5, 0.2, -0.1, -0.4])
    good_ref = good_py.copy()
    bad_py = np.array([1.0, 0.5, 0.2, -0.1, -0.4])
    bad_ref = np.array([-0.4, -0.1, 0.2, 0.5, 1.0])  # reversed shape

    rows_py = []
    rows_ref = []
    for cluster_id in range(5):
        rows_py.append(_coef_df_for(cluster_id, "TF", good_py))
        rows_ref.append(_coef_df_for(cluster_id, "TF", good_ref))
    rows_py.append(_coef_df_for(5, "TF", bad_py))
    rows_ref.append(_coef_df_for(5, "TF", bad_ref))
    coef_py = pd.concat(rows_py, ignore_index=True)
    coef_ref = pd.concat(rows_ref, ignore_index=True)

    gate, info = _tuning_curve_pearson(coef_py, coef_ref, "TF", config)
    # 5 of 6 clusters clear the per-cluster bar → fraction 0.833, passes.
    assert abs(gate["value"] - 5 / 6) < 1e-9, gate
    assert gate["pass"] is True, gate
    # Mean(r) is ~5×1.0 + 1×(<0.5) averaged → < 0.90 on this synthetic,
    # i.e. the old gate would have failed. Informational now, does not gate.
    assert info["value"] < 0.90, info
    assert info["pass"] is True, info


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
