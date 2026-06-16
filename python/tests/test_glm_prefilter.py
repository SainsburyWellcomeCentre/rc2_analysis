"""Prefilter tests: Wilcoxon signed-rank stationary vs motion."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from rc2_glm.config import GLMConfig
from rc2_glm.pipeline import _log_prefilter_summary
from rc2_glm.prefilter import stationary_vs_motion_test


def _toy_prefilter_df() -> pd.DataFrame:
    """5 clusters: 2 motion-responsive (should_run_glm), 3 not."""
    return pd.DataFrame(
        {
            "cluster_id": [10, 11, 12, 13, 14],
            "should_run_glm": [True, True, False, False, False],
            "category": [
                "all_three_significant", "VT_only_significant",
                "none_significant", "none_significant", "T_only_significant",
            ],
        }
    )


def test_apply_prefilter_defaults_off():
    """Policy: the whole selected cohort is fit by default — the stationary-
    vs-motion prefilter is a diagnostic, not the selection gate."""
    assert GLMConfig().apply_prefilter is False


def test_prefilter_summary_not_gating_keeps_whole_cohort(caplog):
    df = _toy_prefilter_df()
    keep_ids = set(df.loc[df["should_run_glm"], "cluster_id"])  # {10, 11}
    with caplog.at_level(logging.INFO, logger="rc2_glm"):
        _log_prefilter_summary(df, keep_ids, gating=False)
    text = "\n".join(r.getMessage() for r in caplog.records)
    # All 5 clusters are kept, not just the 2 motion-responsive ones.
    assert "whole cohort kept: 5" in text
    # The prefilter input total must never be reported as "the cohort".
    assert "total clusters" not in text


def test_prefilter_summary_gating_reports_kept_cohort(caplog):
    df = _toy_prefilter_df()
    keep_ids = set(df.loc[df["should_run_glm"], "cluster_id"])
    with caplog.at_level(logging.INFO, logger="rc2_glm"):
        _log_prefilter_summary(df, keep_ids, gating=True)
    text = "\n".join(r.getMessage() for r in caplog.records)
    assert "selected cohort (prefilter-gated): 2" in text


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
