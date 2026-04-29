"""Per-cluster pre-filtering: stationary vs motion firing-rate test.

Mirrors MATLAB ``is_stationary_vs_motion_significant`` (FormattedData.m
line 648) and the decision tree in ``glm_single_cluster_analysis.m``
(lines 364-458). Despite the MATLAB ``svm_table`` naming, the
significance test is a **Wilcoxon signed-rank test** on paired
per-trial firing rates, not an SVM.

For each (cluster, condition) pair:

1. Compute per-trial mean firing rate during the stationary mask and
   during the motion mask of that trial.
2. Drop trials where either rate is NaN (e.g. zero-duration period).
3. Run ``scipy.stats.wilcoxon(stat, mot)`` and resolve direction the
   same way the MATLAB helper does (sign of the difference of medians,
   tie-broken by signed-rank statistic of the swapped test).

Then categorise each cluster by which conditions are significant:

- ``VT`` significant alone or in combination with T or V → run GLM
- ``T + V`` both significant (no VT) → run GLM
- otherwise → not selected
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from rc2_glm.config import GLMConfig
from rc2_glm.io import ClusterData, ProbeData, TrialData, iter_cluster_spike_slices


CONDITIONS = ("T_Vstatic", "V", "VT")


@dataclass
class ConditionTest:
    condition: str
    significant: bool
    p_value: float
    direction: int           # +1 motion>stationary, -1 motion<stationary, 0 ns
    stationary_median: float
    motion_median: float
    n_trials: int


@dataclass
class PrefilterResult:
    cluster_idx: int
    cluster_id: int
    region: str
    tests: dict[str, ConditionTest]
    is_speed_tuned: bool
    should_run_glm: bool
    category: str


def per_trial_firing_rate(
    cluster: ClusterData, trial: TrialData, fs: float
) -> tuple[float, float]:
    """Return (stationary_rate, motion_rate) in spikes/s for one trial.

    Rates are NaN if the corresponding mask covers no samples.
    """
    spikes = iter_cluster_spike_slices(cluster, trial)
    t = trial.probe_t
    if t.size == 0:
        return float("nan"), float("nan")

    # Bin spikes back onto the per-sample timebase (cheap because
    # spikes-per-trial is small relative to samples).
    # Defensive: when probe_t is float32 and trial spans many seconds,
    # adjacent samples can collide after the +1/fs offset → np.histogram
    # rejects "bins must increase monotonically". Cast to float64 +
    # enforce monotonicity by clipping rather than crashing.
    edges = np.concatenate([t.astype(np.float64), [float(t[-1]) + 1.0 / fs]])
    if not np.all(np.diff(edges) > 0):
        # Force strict monotonicity by adding ε between any non-increasing
        # adjacent pair. Preserves total span; correct counts to within
        # one sample of which spike-time tiebreaks where.
        eps = (1.0 / fs) * 1e-3
        diffs = np.diff(edges)
        for i in np.where(diffs <= 0)[0]:
            edges[i + 1] = edges[i] + eps
    counts, _ = np.histogram(spikes, bins=edges)

    stat_n = int(trial.stationary_mask.sum())
    mot_n = int(trial.motion_mask.sum())

    stat_rate = (
        float(counts[trial.stationary_mask].sum()) / (stat_n / fs)
        if stat_n > 0 else float("nan")
    )
    mot_rate = (
        float(counts[trial.motion_mask].sum()) / (mot_n / fs)
        if mot_n > 0 else float("nan")
    )
    return stat_rate, mot_rate


def stationary_vs_motion_test(
    stationary_rates: np.ndarray,
    motion_rates: np.ndarray,
    alpha: float = 0.05,
) -> ConditionTest | None:
    """Paired Wilcoxon signed-rank test, returning ConditionTest or None.

    Returns ``None`` if fewer than two valid (non-NaN) paired trials
    survive — signrank is undefined in that case.
    """
    s = np.asarray(stationary_rates, dtype=np.float64)
    m = np.asarray(motion_rates, dtype=np.float64)
    keep = ~(np.isnan(s) | np.isnan(m))
    s, m = s[keep], m[keep]
    n = s.size
    if n < 2 or np.all(s == m):
        return None

    try:
        result = wilcoxon(s, m, zero_method="wilcox", alternative="two-sided")
        p_value = float(result.pvalue)
    except ValueError:
        return None

    significant = p_value < alpha
    s_med = float(np.median(s))
    m_med = float(np.median(m))

    if not significant:
        direction = 0
    elif m_med > s_med:
        direction = 1
    elif m_med < s_med:
        direction = -1
    else:
        # MATLAB tie-break: compare swapped signed-rank statistic.
        stat_sm = float(getattr(wilcoxon(s, m, zero_method="wilcox"), "statistic", 0.0))
        stat_ms = float(getattr(wilcoxon(m, s, zero_method="wilcox"), "statistic", 0.0))
        direction = 1 if stat_ms > stat_sm else -1

    return ConditionTest(
        condition="",   # filled in by caller
        significant=significant,
        p_value=p_value,
        direction=direction,
        stationary_median=s_med,
        motion_median=m_med,
        n_trials=int(n),
    )


def prefilter_cluster(
    probe: ProbeData,
    cluster: ClusterData,
    alpha: float = 0.05,
    svm_lookup: dict[tuple[int, int], tuple[float, float]] | None = None,
) -> PrefilterResult:
    """Run the stationary/motion test for T_Vstatic, V, VT and categorise."""
    by_cond: dict[str, list[TrialData]] = {c: [] for c in CONDITIONS}
    for trial in probe.trials:
        if trial.condition in by_cond:
            by_cond[trial.condition].append(trial)

    tests: dict[str, ConditionTest] = {}
    for cond in CONDITIONS:
        trials = by_cond[cond]
        if not trials:
            continue
        stat_rates = np.empty(len(trials))
        mot_rates = np.empty(len(trials))
        for i, trial in enumerate(trials):
            if svm_lookup is not None:
                key = (cluster.cluster_id, trial.trial_id)
                if key in svm_lookup:
                    s, m = svm_lookup[key]
                else:
                    s, m = float("nan"), float("nan")
            else:
                s, m = per_trial_firing_rate(cluster, trial, probe.fs)
            stat_rates[i] = s
            mot_rates[i] = m
        result = stationary_vs_motion_test(stat_rates, mot_rates, alpha=alpha)
        if result is not None:
            result.condition = cond
            tests[cond] = result

    sig_T = tests.get("T_Vstatic", None)
    sig_V = tests.get("V", None)
    sig_VT = tests.get("VT", None)

    is_T = bool(sig_T and sig_T.significant)
    is_V = bool(sig_V and sig_V.significant)
    is_VT = bool(sig_VT and sig_VT.significant)
    n_sig = int(is_T) + int(is_V) + int(is_VT)

    is_speed_tuned = is_T or is_VT

    if is_VT and is_T and is_V:
        category, should_run = "all_three_significant", True
    elif is_VT and is_T:
        category, should_run = "T_and_VT_significant", True
    elif is_VT and is_V:
        category, should_run = "V_and_VT_significant", True
    elif is_VT:
        category, should_run = "VT_only_significant", True
    elif is_T and is_V:
        category, should_run = "T_and_V_significant", True
    elif n_sig == 0:
        category, should_run = "none_significant", False
    elif is_T:
        category, should_run = "T_only_significant", False
    else:
        category, should_run = "V_only_significant", False

    return PrefilterResult(
        cluster_idx=cluster.cluster_idx,
        cluster_id=cluster.cluster_id,
        region=cluster.region,
        tests=tests,
        is_speed_tuned=is_speed_tuned,
        should_run_glm=should_run,
        category=category,
    )


def prefilter_probe(
    probe: ProbeData,
    config: GLMConfig | None = None,
    alpha: float = 0.05,
    use_precomputed: bool = True,
) -> pd.DataFrame:
    """Run the prefilter for every cluster in ``probe`` and return a table."""
    import logging
    logger = logging.getLogger(__name__)

    config = config or probe.config

    svm_lookup = None
    if use_precomputed:
        svm_csv = probe.mat_path.parent / "csvs" / "stationary_vs_motion_fr" / probe.mat_path.with_suffix(".csv").name
        if svm_csv.exists():
            logger.info("Loading precomputed stationary vs motion FR from %s", svm_csv.name)
            try:
                df = pd.read_csv(svm_csv)
                # Dictionary mapping (cluster_id, trial_id) -> (stationary_fr, motion_fr)
                svm_lookup = {
                    (int(r.cluster_id), int(r.trial_id)): (float(r.stationary_fr), float(r.motion_fr))
                    for r in df.itertuples()
                }
            except Exception as e:
                logger.warning("Failed to load %s: %s", svm_csv.name, e)
        else:
            logger.info("No precomputed stationary vs motion FR found at %s. Recalculating.", svm_csv.name)

    rows: list[dict] = []
    for cluster in probe.clusters:
        res = prefilter_cluster(probe, cluster, alpha=alpha, svm_lookup=svm_lookup)
        row = {
            "probe_id": probe.probe_id,
            "cluster_idx": res.cluster_idx,
            "cluster_id": res.cluster_id,
            "region": res.region,
            "is_speed_tuned": res.is_speed_tuned,
            "should_run_glm": res.should_run_glm,
            "category": res.category,
        }
        for cond in CONDITIONS:
            t = res.tests.get(cond)
            row[f"{cond}_p"] = t.p_value if t else float("nan")
            row[f"{cond}_sig"] = bool(t.significant) if t else False
            row[f"{cond}_dir"] = t.direction if t else 0
            row[f"{cond}_stat_med"] = t.stationary_median if t else float("nan")
            row[f"{cond}_mot_med"] = t.motion_median if t else float("nan")
            row[f"{cond}_n"] = t.n_trials if t else 0
        rows.append(row)
    return pd.DataFrame(rows)
