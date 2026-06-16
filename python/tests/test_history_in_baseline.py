"""history_in_baseline mode (2026-06-16).

When ``config.history_in_baseline`` is True the spike-History term is an
always-on baseline/nuisance regressor — in the null model and every fitted
model, like the onset kernel — but is NEVER a forward-selection candidate and
NEVER displayed. Tests:

1. Baseline Null design CONTAINS History columns; legacy Null does not.
2. History is never selected under a strong history signal in baseline mode;
   in legacy include_history mode it IS offered as a candidate.
3. ``SelectionResult.history_mode`` reports off / candidate / baseline.
4. ``forward_select`` fails loud if history_in_baseline but B_history is None.
5. Display suppression: ``_plot_beta_swarm(hide_history=True)`` drops History;
   ``plot_cluster_kernels`` renders no History panel when the flag is on.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from rc2_glm.config import GLMConfig
from rc2_glm.cross_validation import make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.forward_selection import forward_select


def _toy(n, freq, rng):
    t = np.linspace(0, 1, n)
    return np.stack([np.cos(freq * np.pi * t + 2 * np.pi * k / 5) for k in range(5)], 1) \
        + rng.normal(0, 0.05, (n, 5))


def _data(seed=0, history_strength=1.0, n_folds=10):
    rng = np.random.default_rng(seed)
    nt, nb = 30, 30
    n = nt * nb
    trial_ids = np.repeat(np.arange(nt), nb)
    cond = np.array(["A" if t % 2 == 0 else "B" for t in range(nt)])[trial_ids]
    B_speed, B_tf = _toy(n, 2.0, rng), _toy(n, 3.0, rng)
    B_history = rng.normal(0, 1, (n, 4))
    y = rng.poisson(
        np.exp(-1.0 + 1.5 * B_speed[:, 0] + history_strength * B_history[:, 0]) * 0.1
    )
    folds = make_trial_folds(trial_ids, n_folds=n_folds, seed=0,
                             condition_labels_per_bin=cond)
    return dict(
        B_speed=B_speed, B_tf=B_tf, B_onset=np.zeros((n, 0)),
        sf=np.zeros(n), or_=np.zeros(n), y=y, offset=float(np.log(0.1)),
        folds=folds, B_history=B_history, n=n,
    )


def _cfg(**kw):
    return GLMConfig(
        include_onset_kernel=False, include_me_face=False,
        main_effects=("Speed", "TF", "SF", "OR"), interactions=(), n_folds=10,
        **kw,
    )


def _select(d, cfg):
    return forward_select(
        d["B_speed"], d["B_tf"], d["B_onset"], d["sf"], d["or_"],
        d["y"], d["offset"], d["folds"], config=cfg, backend="irls",
        B_history=d["B_history"],
    )


def test_baseline_null_contains_history_legacy_does_not():
    d = _data()
    _, names_base = assemble_design_matrix_selected(
        d["B_speed"], d["B_tf"], d["B_onset"], d["sf"], d["or_"], [],
        B_history=d["B_history"], include_onset_kernel=False,
        history_in_baseline=True,
    )
    _, names_legacy = assemble_design_matrix_selected(
        d["B_speed"], d["B_tf"], d["B_onset"], d["sf"], d["or_"], [],
        B_history=d["B_history"], include_onset_kernel=False,
        history_in_baseline=False,
    )
    assert any(nm.startswith("History_") for nm in names_base)
    assert not any(nm.startswith("History_") for nm in names_legacy)


def test_history_never_selected_in_baseline_mode():
    d = _data(history_strength=2.0)
    res = _select(d, _cfg(history_in_baseline=True))
    assert "History" not in res.selected_vars
    assert res.history_mode == "baseline"
    # History is not even a tested candidate in any round.
    for r in res.history:
        assert "History" not in r.tested


def test_history_is_candidate_in_legacy_mode():
    d = _data(history_strength=2.0)
    res = _select(d, _cfg(include_history=True))
    assert res.history_mode == "candidate"
    # History appears as a tested candidate in round 1.
    assert "History" in res.history[0].tested


def test_history_mode_off_when_no_history():
    d = _data()
    res = forward_select(
        d["B_speed"], d["B_tf"], d["B_onset"], d["sf"], d["or_"],
        d["y"], d["offset"], d["folds"], config=_cfg(), backend="irls",
        B_history=None,
    )
    assert res.history_mode == "off"


def test_baseline_requires_b_history():
    d = _data()
    with pytest.raises(ValueError, match="history_in_baseline"):
        forward_select(
            d["B_speed"], d["B_tf"], d["B_onset"], d["sf"], d["or_"],
            d["y"], d["offset"], d["folds"],
            config=_cfg(history_in_baseline=True), backend="irls",
            B_history=None,
        )


def test_beta_swarm_hides_history():
    from rc2_glm.plots import _plot_beta_swarm
    import matplotlib.pyplot as plt

    beta = np.array([0.1, 0.2, 0.3, 0.4])
    cols = ["Intercept", "History_1", "History_2", "Speed_1"]
    fig, ax = plt.subplots()
    _plot_beta_swarm(ax, beta, cols, hide_history=True)
    tick_labels = [t.get_text() for t in ax.get_xticklabels()]
    assert "History" not in tick_labels
    assert "Speed" in tick_labels
    plt.close(fig)


def test_cluster_kernels_render_no_history_panel_when_baseline():
    from rc2_glm.plots import plot_cluster_kernels
    import matplotlib.pyplot as plt

    # Minimal coefficient table with a Selected model carrying History_* +
    # Speed_*; in baseline mode the History rows must be dropped before plotting.
    cfg = _cfg(history_in_baseline=True)
    rows = []
    for model in ("Null", "Selected", "Additive", "FullInteraction"):
        rows.append({"model": model, "coefficient": "Intercept", "estimate": -1.0, "se": 0.1})
        if model != "Null":
            for i in range(1, 4):  # History basis count is irrelevant — filtered out
                rows.append({"model": model, "coefficient": f"History_{i}",
                             "estimate": 0.2, "se": 0.05})
            for i in range(1, cfg.n_speed_bases + 1):  # Speed needs the full basis set
                rows.append({"model": model, "coefficient": f"Speed_{i}",
                             "estimate": 0.3, "se": 0.05})
    coef_df = pd.DataFrame(rows)
    fig = plot_cluster_kernels("P", 1, coef_df, cfg)
    # No axis anywhere should be titled for History.
    titles = [ax.get_title() for ax in fig.axes]
    assert not any("History" in (t or "") for t in titles), titles
    plt.close(fig)


def _fold_fixture():
    """A small cluster_df + beta/train_names for _fold_history_into_intercept,
    using an identity history basis with 2 lag bins (window 0.04 / bin 0.02)."""
    n = 80
    trial_id = np.repeat(np.arange(8), 10)
    spike_count = (np.arange(n) % 4).astype(float)  # nonzero, varied
    condition = np.where(np.arange(n) % 2 == 0, "motion", "stationary")
    cluster_df = pd.DataFrame(
        {"spike_count": spike_count, "trial_id": trial_id, "condition": condition}
    )
    train_names = ["Intercept", "History_1", "History_2"]
    beta = np.array([-1.0, 0.5, 0.5])
    return beta, train_names, cluster_df


def test_history_fold_fires_in_baseline_mode():
    """The marginal History fold (mean contribution → Intercept) must fire in
    history_in_baseline mode, not only legacy include_history — else the
    History contribution is silently dropped from every predicted curve."""
    from rc2_glm.plots import _fold_history_into_intercept

    beta, names, df = _fold_fixture()
    cfg = _cfg(history_in_baseline=True, history_basis_kind="identity",
               history_window_s=0.04, time_bin_width=0.02)
    out = _fold_history_into_intercept(beta, names, df, cfg)
    # Intercept shifted by the (nonzero) History mean contribution; History
    # coefficients themselves unchanged.
    assert out[0] != beta[0]
    np.testing.assert_array_equal(out[1:], beta[1:])


def test_history_fold_noop_when_both_history_modes_off():
    from rc2_glm.plots import _fold_history_into_intercept

    beta, names, df = _fold_fixture()
    cfg = _cfg(include_history=False, history_in_baseline=False,
               history_basis_kind="identity", history_window_s=0.04,
               time_bin_width=0.02)
    out = _fold_history_into_intercept(beta, names, df, cfg)
    np.testing.assert_array_equal(out, beta)
