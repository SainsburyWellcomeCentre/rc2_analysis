"""Prefilter tests: Wilcoxon signed-rank stationary vs motion."""

from __future__ import annotations

import numpy as np

from rc2_glm.prefilter import stationary_vs_motion_test


def test_wilcoxon_detects_motion_increase():
    rng = np.random.default_rng(0)
    stat = rng.normal(2.0, 0.3, 30)
    mot = stat + 1.5  # motion rate clearly higher
    res = stationary_vs_motion_test(stat, mot)
    assert res is not None
    assert res.significant
    assert res.direction == 1


def test_wilcoxon_detects_motion_decrease():
    rng = np.random.default_rng(1)
    stat = rng.normal(5.0, 0.3, 30)
    mot = stat - 2.0
    res = stationary_vs_motion_test(stat, mot)
    assert res is not None
    assert res.significant
    assert res.direction == -1


def test_wilcoxon_no_difference_is_ns():
    rng = np.random.default_rng(2)
    vals = rng.normal(3.0, 0.3, 40)
    res = stationary_vs_motion_test(vals.copy(), vals + rng.normal(0, 1e-6, 40))
    # With near-identical paired data the test p-value is high
    assert res is not None
    assert not res.significant
    assert res.direction == 0


def test_wilcoxon_nan_handling():
    stat = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    mot = np.array([2.0, 3.0, 5.0, np.nan, 7.0])
    res = stationary_vs_motion_test(stat, mot)
    assert res is not None
    assert res.n_trials == 3   # two pairs dropped for NaN


def test_wilcoxon_too_few_returns_none():
    res = stationary_vs_motion_test(np.array([1.0]), np.array([2.0]))
    assert res is None
