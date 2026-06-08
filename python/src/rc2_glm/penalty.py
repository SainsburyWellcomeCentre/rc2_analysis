"""Structured penalty (prior precision) matrices for the Poisson GLM.

The default penalty is isotropic ridge: ``lambda * I`` with the intercept
left unpenalised — exactly what ``_fit_irls`` built inline before this module
existed. ``build_penalty_matrix`` reproduces that when no smoothness prior is
requested, so passing its output is a no-op relative to the historical fit.

Spike-history smoothness prior (opt-in, ``history_smooth_lambda`` set)
---------------------------------------------------------------------
At fine bin widths the 5 raised-cosine history weights can take large
opposite-sign values on overlapping bases, so the *reconstructed lag-space
filter* ``f = B @ beta_hist`` oscillates wildly (the ±5 swings seen in the
2026-04-29 20 ms sweep). Isotropic ridge does not stop this — it penalises
coefficient magnitude, not filter roughness.

This module instead penalises the squared second difference of the
reconstructed filter::

    penalty(beta_hist) = lambda_smooth * || D2 @ (B @ beta_hist) ||^2
                       = beta_hist.T @ (lambda_smooth * B.T D2.T D2 B) @ beta_hist

so smooth (constant / linear) filters are free and curvature is penalised.
This is the fixed-hyperparameter form of the ASD / smooth-Gaussian prior
(Park & Pillow 2011). A small ``floor * I`` keeps the block positive-definite
(the curvature term alone is rank-deficient — constant and linear filters lie
in its null space, which is the intended behaviour).

At 100 ms / 200 ms window there are only 2 lag bins, so ``D2`` is empty and
the history block falls back to the ridge floor — smoothing only bites once
the grid is fine enough to resolve >=3 lag bins (e.g. 20 ms).
"""

from __future__ import annotations

import numpy as np


def second_difference_matrix(n: int) -> np.ndarray:
    """Discrete second-difference operator, shape ``(n - 2, n)``.

    Row ``i`` encodes ``f[i] - 2 f[i+1] + f[i+2]``. Returns a ``(0, n)``
    array when ``n < 3`` (no curvature is defined on fewer than 3 points).
    """
    if n < 3:
        return np.zeros((0, n), dtype=np.float64)
    D = np.zeros((n - 2, n), dtype=np.float64)
    rows = np.arange(n - 2)
    D[rows, rows] = 1.0
    D[rows, rows + 1] = -2.0
    D[rows, rows + 2] = 1.0
    return D


def build_penalty_matrix(
    col_names: list[str],
    lambda_ridge: float,
    *,
    lambda_min: float = 1e-6,
    history_smooth_lambda: float | None = None,
    history_basis_mat: np.ndarray | None = None,
    history_floor: float | None = None,
) -> np.ndarray:
    """Build the ``(p, p)`` penalty matrix added to ``X.T W X`` in IRLS.

    Non-history coefficients get isotropic ridge ``max(lambda_ridge,
    lambda_min)``; any column named ``"Intercept"`` is unpenalised. When
    ``history_smooth_lambda`` is set and ``History_*`` columns are present,
    that block is replaced by ``history_smooth_lambda * B.T D2.T D2 B +
    floor * I`` (see module docstring). ``history_floor`` defaults to the
    same ridge level used elsewhere.

    ``history_basis_mat`` must be the ``(n_lag, n_bases)`` array from
    ``basis.history_basis`` for the *same* config (n_bases / window /
    bin_width) that produced the ``History_*`` design columns.
    """
    p = len(col_names)
    lambda_eff = max(float(lambda_ridge), float(lambda_min))
    P = lambda_eff * np.eye(p)
    for i, nm in enumerate(col_names):
        if nm == "Intercept":
            P[i, i] = 0.0

    if history_smooth_lambda is None or history_basis_mat is None:
        return P

    hist_idx = [i for i, nm in enumerate(col_names) if nm.startswith("History_")]
    if not hist_idx:
        return P

    B = np.asarray(history_basis_mat, dtype=np.float64)
    n_lag, n_bases = B.shape
    if n_bases != len(hist_idx):
        raise ValueError(
            f"history basis has {n_bases} columns but design matrix has "
            f"{len(hist_idx)} History_* columns; config mismatch"
        )

    D2 = second_difference_matrix(n_lag)
    if D2.shape[0] > 0:
        S = B.T @ D2.T @ D2 @ B          # (n_bases, n_bases), curvature
    else:
        S = np.zeros((n_bases, n_bases))  # <3 lag bins: nothing to smooth

    floor = lambda_eff if history_floor is None else float(history_floor)
    block = history_smooth_lambda * S + floor * np.eye(n_bases)

    idx = np.array(hist_idx)
    P[np.ix_(idx, idx)] = block
    return P
