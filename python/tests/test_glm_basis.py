"""Unit tests for raised-cosine and onset-kernel bases.

The MATLAB references are ``make_raised_cosine_basis`` (line 5506) and
``make_onset_kernel_basis`` (line 5531). Tests check shape, partition-of-
unity-style invariants, and direct analytic values to ~1e-12.
"""

from __future__ import annotations

import numpy as np

from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis


def test_raised_cosine_shape_and_range():
    x = np.linspace(0.0, 50.0, 200)
    B = raised_cosine_basis(x, 5, 0.0, 50.0)
    assert B.shape == (200, 5)
    assert (B >= 0.0).all() and (B <= 1.0 + 1e-12).all()


def test_raised_cosine_centers_peak_at_one():
    """At x = center, the cosine argument is 0 → basis value is 1."""
    n_bases = 5
    x_min, x_max = 0.0, 50.0
    epsilon = 0.5
    log_min = np.log(x_min + epsilon)
    log_max = np.log(x_max + epsilon)
    centers_log = np.linspace(log_min, log_max, n_bases)
    centers_x = np.exp(centers_log) - epsilon
    B = raised_cosine_basis(centers_x, n_bases, x_min, x_max)
    diag = np.diag(B)
    np.testing.assert_allclose(diag, np.ones(n_bases), atol=1e-12)


def test_raised_cosine_clipped_outside_support():
    """Far outside the [-pi, pi] window the basis must be 0."""
    x = np.array([1e6])
    B = raised_cosine_basis(x, 5, 0.0, 50.0)
    np.testing.assert_allclose(B[0, 0], 0.0, atol=1e-12)


def test_onset_kernel_zero_for_negative_t():
    t = np.array([-1.0, -0.5, -0.01])
    B = onset_kernel_basis(t, 6, 2.0)
    assert np.all(B == 0.0)


def test_onset_kernel_first_basis_at_t_zero():
    """At t=0 (= centers[0]), basis 0 should equal 1 exactly."""
    t = np.array([0.0])
    B = onset_kernel_basis(t, 6, 2.0)
    np.testing.assert_allclose(B[0, 0], 1.0, atol=1e-12)


def test_onset_kernel_shape():
    t = np.linspace(0.0, 2.0, 100)
    B = onset_kernel_basis(t, 6, 2.0)
    assert B.shape == (100, 6)
    assert (B >= 0.0).all()
