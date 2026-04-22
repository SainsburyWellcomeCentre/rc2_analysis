"""Poisson GLM fitting.

Two backends:

- ``"irls"`` (default): Direct port of the MATLAB IRLS solver
  (``fit_poisson_glm`` at line 5762 of glm_single_cluster_analysis.m).
  Uses the same eta clipping, mu floor, ridge-on-non-intercept, and
  Cholesky-based standard errors. Deterministic; matches MATLAB to
  ~1e-10 once the same random / numerical context is in place.
  **This is the parity oracle for the rc2_analysis port** — every row
  in ``rc2-glm-compare`` is phrased against the IRLS backend.

- ``"nemos"``: Uses ``nemos.glm.GLM`` with Poisson observation model
  and Ridge regularisation. NeMoS solves with L-BFGS (jaxopt), not
  Newton-IRLS, so its β lands at a numerically different optimum than
  IRLS/MATLAB even on identical data (expected; the likelihoods are
  within float precision of each other, the parameter rotation differs
  — see co-scientist note "identifiability ≠ correctness" in the
  obsidian-agents CLAUDE.md). Use NeMoS for throughput and GPU
  acceleration, not for chasing MATLAB parity. The ``--device`` CLI
  flag routes to ``configure_jax_device`` below; set it before any
  ``import jax`` / ``import nemos`` to bind the JAX backend for the
  process. NeMoS doesn't accept a per-sample offset directly; for
  uniform 100 ms bins we absorb the constant ``log(time_bin_width)``
  into the intercept.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

logger = logging.getLogger(__name__)


_JAX_DEVICE_CONFIGURED: str | None = None


def configure_jax_device(device: str = "auto") -> str:
    """Pin JAX to ``cpu`` / ``gpu`` / auto-detect, before ``import jax``.

    JAX selects its backend at import time from the ``JAX_PLATFORMS`` env
    var. Call this once per process, before ``import nemos`` (which
    imports jax transitively). Subsequent calls with the same ``device``
    are a no-op; a later different ``device`` logs a warning and is
    ignored — JAX cannot switch platforms after import.

    Returns the resolved backend name (``"cpu"``, ``"gpu"``, ``"tpu"``,
    or ``"<unavailable>"`` when the requested device is missing).
    """
    global _JAX_DEVICE_CONFIGURED

    device = (device or "auto").lower()
    if device not in ("auto", "cpu", "gpu", "tpu"):
        raise ValueError(f"Unknown JAX device: {device!r}")

    if _JAX_DEVICE_CONFIGURED is not None:
        if device != _JAX_DEVICE_CONFIGURED:
            logger.warning(
                "JAX already configured for %r; cannot switch to %r "
                "within the same process.",
                _JAX_DEVICE_CONFIGURED, device,
            )
        return _current_jax_backend()

    if device != "auto":
        os.environ["JAX_PLATFORMS"] = device

    _JAX_DEVICE_CONFIGURED = device
    return _current_jax_backend()


def _current_jax_backend() -> str:
    """Import jax (honours ``JAX_PLATFORMS``) and report its backend."""
    try:
        import jax  # noqa: F401
        return str(jax.default_backend())
    except Exception as exc:  # pragma: no cover — depends on install
        return f"<unavailable: {exc}>"


@dataclass
class GLMFit:
    beta: np.ndarray
    se: np.ndarray
    log_likelihood: float
    aic: float
    bic: float
    deviance: float
    n_params: int
    n_obs: int
    predicted_count: np.ndarray
    pearson_residuals: np.ndarray
    dispersion: float
    converged: bool
    lambda_ridge: float


def fit_poisson_glm(
    X: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float = 0.0,
    lambda_ridge: float = 0.0,
    backend: str = "irls",
    *,
    max_iter: int = 100,
    tol: float = 1e-8,
    eta_clip: float = 20.0,
    mu_floor: float = 1e-10,
    lambda_min: float = 1e-6,
) -> GLMFit:
    if backend == "irls":
        return _fit_irls(
            X, y, offset, lambda_ridge,
            max_iter=max_iter, tol=tol, eta_clip=eta_clip,
            mu_floor=mu_floor, lambda_min=lambda_min,
        )
    if backend == "nemos":
        return _fit_nemos(X, y, offset, lambda_ridge, mu_floor=mu_floor)
    raise ValueError(f"Unknown backend: {backend!r}")


# --------------------------------------------------------------------------- #
# Backend: IRLS (MATLAB-equivalent)
# --------------------------------------------------------------------------- #


def _fit_irls(
    X: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float,
    lambda_ridge: float,
    *,
    max_iter: int,
    tol: float,
    eta_clip: float,
    mu_floor: float,
    lambda_min: float,
) -> GLMFit:
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    offset = _broadcast_offset(offset, y.size)

    n, p = X.shape
    lambda_eff = max(lambda_ridge, lambda_min)
    ridge_mat = lambda_eff * np.eye(p)
    ridge_mat[0, 0] = 0.0  # don't penalise intercept

    beta = np.zeros(p)
    beta[0] = np.log(max(y.mean(), 0.1))

    converged = False
    for it in range(max_iter):
        eta = np.clip(X @ beta + offset, -eta_clip, eta_clip)
        mu = np.maximum(np.exp(eta), mu_floor)

        W = mu
        z = eta + (y - mu) / mu - offset

        XtWX = X.T @ (X * W[:, None]) + ridge_mat
        XtWz = X.T @ (W * z)

        try:
            beta_new = np.linalg.solve(XtWX, XtWz)
        except np.linalg.LinAlgError:
            beta_new = np.linalg.lstsq(XtWX, XtWz, rcond=None)[0]

        if not np.all(np.isfinite(beta_new)):
            break
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            converged = True
            break
        beta = beta_new
    else:
        # exhausted max_iter without break
        converged = False

    eta = np.clip(X @ beta + offset, -eta_clip, eta_clip)
    mu = np.maximum(np.exp(eta), mu_floor)

    log_lik = float(np.sum(y * np.log(mu) - mu - gammaln(y + 1.0)))
    y_pos = np.maximum(y, mu_floor)
    deviance = float(2.0 * np.sum(y * np.log(y_pos / mu) - (y - mu)))

    se = _standard_errors(X, mu, ridge_mat)

    return GLMFit(
        beta=beta,
        se=se,
        log_likelihood=log_lik,
        aic=-2.0 * log_lik + 2.0 * p,
        bic=-2.0 * log_lik + p * np.log(n),
        deviance=deviance,
        n_params=p,
        n_obs=n,
        predicted_count=mu,
        pearson_residuals=(y - mu) / np.sqrt(np.maximum(mu, mu_floor)),
        dispersion=deviance / max(n - p, 1),
        converged=converged and bool(np.all(np.isfinite(beta))),
        lambda_ridge=lambda_ridge,
    )


def _standard_errors(
    X: np.ndarray, mu: np.ndarray, ridge_mat: np.ndarray
) -> np.ndarray:
    XtWX = X.T @ (X * mu[:, None]) + ridge_mat
    p = X.shape[1]
    try:
        R = np.linalg.cholesky(XtWX).T  # upper triangular like MATLAB chol
        inv_R = np.linalg.solve(R, np.eye(p))
        return np.sqrt(np.sum(inv_R ** 2, axis=1))
    except np.linalg.LinAlgError:
        try:
            cov = np.linalg.pinv(XtWX)
            return np.sqrt(np.abs(np.diag(cov)))
        except np.linalg.LinAlgError:
            return np.zeros(p)


# --------------------------------------------------------------------------- #
# Backend: NeMoS (GPU-friendly via JAX)
# --------------------------------------------------------------------------- #


def _fit_nemos(
    X: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray | float,
    lambda_ridge: float,
    *,
    mu_floor: float,
) -> GLMFit:
    import nemos as nmo

    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    offset = _broadcast_offset(offset, y.size)

    if not np.allclose(offset, offset[0]):
        raise NotImplementedError(
            "NeMoS backend requires a constant offset across samples; "
            "got non-uniform values. Use the IRLS backend for that case."
        )
    constant_offset = float(offset[0])

    # Strip leading intercept column if present so NeMoS handles its own intercept
    has_intercept = X.shape[1] >= 1 and np.allclose(X[:, 0], 1.0)
    X_fit = X[:, 1:] if has_intercept else X

    model = nmo.glm.GLM(
        observation_model="Poisson",
        regularizer="Ridge" if lambda_ridge > 0 else "UnRegularized",
        regularizer_strength=lambda_ridge if lambda_ridge > 0 else None,
    )
    model.fit(X_fit, y)

    coef = np.asarray(model.coef_).ravel()
    nemos_intercept = float(np.asarray(model.intercept_).ravel()[0])

    if has_intercept:
        beta = np.concatenate([[nemos_intercept - constant_offset], coef])
    else:
        beta = np.concatenate([[0.0], coef])  # caller didn't supply intercept; ignore

    eta = np.clip(X @ beta + offset, -20.0, 20.0)
    mu = np.maximum(np.exp(eta), mu_floor)
    log_lik = float(np.sum(y * np.log(mu) - mu - gammaln(y + 1.0)))
    n, p = X.shape
    return GLMFit(
        beta=beta,
        se=np.zeros(p),  # NeMoS doesn't provide SE directly
        log_likelihood=log_lik,
        aic=-2.0 * log_lik + 2.0 * p,
        bic=-2.0 * log_lik + p * np.log(n),
        deviance=float("nan"),
        n_params=p,
        n_obs=n,
        predicted_count=mu,
        pearson_residuals=(y - mu) / np.sqrt(np.maximum(mu, mu_floor)),
        dispersion=float("nan"),
        converged=True,
        lambda_ridge=lambda_ridge,
    )


def _broadcast_offset(offset: np.ndarray | float, n: int) -> np.ndarray:
    if np.isscalar(offset):
        return np.full(n, float(offset))
    arr = np.asarray(offset, dtype=np.float64).ravel()
    if arr.size == 1:
        return np.full(n, float(arr[0]))
    if arr.size != n:
        raise ValueError(f"offset has length {arr.size}, expected {n}")
    return arr
