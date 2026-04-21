"""IRLS Poisson fit tests on synthetic data."""

from __future__ import annotations

import numpy as np

from rc2_glm.fitting import fit_poisson_glm


def _synth_poisson(n=2000, seed=0, intercept=-1.0, beta=(0.5, -0.3, 0.2)):
    rng = np.random.default_rng(seed)
    p = len(beta)
    X = np.column_stack([np.ones(n)] + [rng.normal(size=n) for _ in range(p)])
    eta = X @ np.r_[intercept, beta]
    mu = np.exp(eta)
    y = rng.poisson(mu)
    return X, y, np.r_[intercept, beta]


def test_irls_recovers_known_coefficients():
    X, y, beta_true = _synth_poisson(n=5000, seed=42)
    fit = fit_poisson_glm(X, y, offset=0.0, lambda_ridge=0.0, backend="irls")
    assert fit.converged
    np.testing.assert_allclose(fit.beta, beta_true, atol=0.1)


def test_irls_intercept_is_not_ridge_penalised():
    """With huge ridge, the intercept must still pick up the mean rate."""
    X, y, _ = _synth_poisson(n=2000, seed=1, intercept=2.0, beta=(0.0, 0.0, 0.0))
    fit = fit_poisson_glm(X, y, offset=0.0, lambda_ridge=1e6, backend="irls")
    # Slope coefficients shrunk to ~0
    np.testing.assert_allclose(fit.beta[1:], 0.0, atol=1e-3)
    # Intercept ≈ log(mean(y))
    np.testing.assert_allclose(fit.beta[0], np.log(y.mean()), atol=0.05)


def test_offset_shifts_intercept():
    """A constant offset c should shift the recovered intercept by -c."""
    X, y, beta_true = _synth_poisson(n=3000, seed=3)
    fit_no_off = fit_poisson_glm(X, y, offset=0.0, lambda_ridge=0.0)
    fit_w_off = fit_poisson_glm(X, y, offset=1.5, lambda_ridge=0.0)
    np.testing.assert_allclose(
        fit_w_off.beta[0], fit_no_off.beta[0] - 1.5, atol=1e-3
    )
    np.testing.assert_allclose(fit_w_off.beta[1:], fit_no_off.beta[1:], atol=1e-3)


def test_standard_errors_positive_finite():
    X, y, _ = _synth_poisson(n=1500, seed=5)
    fit = fit_poisson_glm(X, y, offset=0.0)
    assert (fit.se >= 0).all()
    assert np.all(np.isfinite(fit.se))
