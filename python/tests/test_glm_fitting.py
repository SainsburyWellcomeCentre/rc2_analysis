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


def test_default_lambda_ridge_is_nonzero_for_correlated_basis_identifiability():
    """GLMConfig().lambda_ridge must be > 0.

    The raised-cosine Speed/TF bases are correlated by construction (Park
    et al. 2014) so the Poisson likelihood has flat directions — the
    unregularised fit falls onto a different β rotation than MATLAB's
    glmnet did, making raw-β sign agreement unreliable.

    Keeping a small non-zero ridge stabilises the fit without materially
    biasing the tuning curve (B @ β is insensitive to the rotation that
    ridge prefers). If someone sets this back to 0.0 to "match MATLAB
    exactly", they should rename this test rather than silently flipping
    back.
    """
    from rc2_glm.config import GLMConfig
    cfg = GLMConfig()
    assert cfg.lambda_ridge > 0.0, cfg.lambda_ridge


def test_nemos_backend_recovers_known_coefficients():
    """NeMoS fit is within a looser tolerance of the IRLS ground truth.

    NeMoS uses L-BFGS (jaxopt) rather than Newton-IRLS, and L-BFGS's
    default stopping tolerance is looser than our IRLS ``tol=1e-8`` —
    so β won't match IRLS bit-for-bit, but on well-conditioned synthetic
    data it should land within ~0.1 of ground truth. The test exists to
    catch silent breakage of the nemos path (e.g. a future NeMoS API
    change that leaves model.coef_ un-fit), not to enforce numerical
    parity with IRLS.
    """
    pytest = __import__("pytest")
    nmo = pytest.importorskip("nemos")  # type: ignore[attr-defined]
    # nemos is declared in pyproject dependencies — skip only in environments
    # where jax fails to initialise a backend at all.
    del nmo

    X, y, beta_true = _synth_poisson(n=5000, seed=42)
    fit = fit_poisson_glm(X, y, offset=0.0, lambda_ridge=0.0, backend="nemos")
    np.testing.assert_allclose(fit.beta, beta_true, atol=0.1)
    assert fit.converged  # NeMoS reports True after model.fit returns
    assert np.all(np.isfinite(fit.beta))


def test_configure_jax_device_accepts_auto_and_cpu():
    """configure_jax_device is idempotent and returns a backend name.

    We can't switch JAX backends mid-process (JAX binds at import time),
    so the test only asserts: (a) the helper runs without raising for
    valid inputs, (b) it returns a non-empty backend name string, and
    (c) re-invoking with the same device is a no-op.
    """
    from rc2_glm.fitting import configure_jax_device

    backend = configure_jax_device("auto")
    assert isinstance(backend, str) and backend
    again = configure_jax_device("auto")
    assert again == backend
