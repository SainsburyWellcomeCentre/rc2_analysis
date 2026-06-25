"""FENS poster figure 1 — LEFT panel (Hardcastle 2017 fig 1C/D style).

TEST DRIVER (iteration stage). Builds, for a chosen cluster, a cumulative
model-buildup on one example trial: starting from a fixed baseline
(``Intercept + Onset + History``), components are added one at a time and the
predicted firing-rate trace is overlaid on the observed FR — annotating the
gain in cross-validated bits/spike (cv-bps) at each addition. This mirrors
Hardcastle et al. 2017 (Neuron) panel D (example-cell prediction improving with
model complexity); the right panel (forward-selection process, top of their
panel C) is deferred.

Status / honesty note: this is intentionally a *standalone* script, NOT yet
folded into ``rc2_glm.plots`` / ``rc2_glm.pipeline``. The per-cluster design
prep below is a faithful MIRROR of ``pipeline._fit_one_cluster`` (lines ~489-642
on the glm-improvements branch). To guard against silent divergence, the script
oracle-checks itself: the buildup's cv-bps for each cluster's stored *Selected*
variable set must reproduce ``time_Selected_cv_bps`` in the run's
``glm_model_comparison.csv``. Once Laura picks a cluster and locks the figure
design, the fit helper migrates into the package (a shared cluster-prep helper
reused by both ``_fit_one_cluster`` and this buildup) per the display contract.

Data source: the goggles History-in-baseline run config (``make_config_histbase_all``
from ``run_glm_goggles_rf_sfor_20ms``), read against the
``figures/glm/current_rf_sfor_20ms_histbase_goggles_all`` outputs. History is an
always-on baseline nuisance there (never a forward-selection candidate), and the
selection rule is the condition-stratified 10-fold one-sided Wilcoxon signed-rank
test. The fig2c bottom panel plots each candidate at the CUMULATIVE-over-baseline
bits/spike its model would reach (running cumulative + its paired fold-mean Δ,
error bar = the per-fold Δ SD over the 10 folds, the signed-rank input), with a
green vertical segment marking the accepted step from one column to the next.

Usage:
    python scripts/make_fens_poster_figures.py \
        --probe CAA-1124370_rec1_rec2_rec3 --clusters 194 125 20 29
"""

from __future__ import annotations

# Pin BLAS/OMP to a single thread BEFORE numpy imports so the IRLS cv-bps is
# reproducible: unpinned, thread-order float summation makes cv-bps wander by
# ~few×1e-3 across process launches (a property of the production pipeline's
# solver, not this script). setdefault lets the caller override.
import os as _os

for _v in (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    _os.environ.setdefault(_v, "1")

import argparse
import logging
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from rc2_glm.basis import (
    circular_basis,
    convolve_history,
    history_basis,
    onset_kernel_basis,
    raised_cosine_basis_linear,
    value_basis,
)
from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.fitting import fit_poisson_glm
from rc2_glm.io import load_probe_data
from rc2_glm.tuning_significance import (
    LINEAR_FAMILIES as LINEAR_FAMILIES_TS,
    MIN_TRIALS,
    TuningSignificance,
    evaluate as eval_tuning_fit,
    fit_tuning,
    per_trial_bin_matrix,
    rsq_against_mean,
    select_best,
    tuning_significance,
)
from rc2_glm.precomputed_bins import load_precomputed_bin_edges
from rc2_glm.rf_sf_or import load_rf_sf_or
from rc2_glm.time_binning import bin_cluster

# The _all run config + paths come straight from the production driver so this
# stays byte-faithful to how the cohort was fit.
from run_glm_goggles_rf_sfor_20ms import (
    FORMATTED_DIR,
    OUT_HISTBASE_ALL,
    RF_PARQUET_DIR,
    _lookup,
    make_config_histbase_all,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("fens_poster")

HOME = Path.home()
# Figures land in the FENS_figures_poster ROOT (the /test iteration subdir is
# retired now the figure design is locked on the histbase fits).
OUT_DIR = HOME / "local_data" / "motion_clouds" / "figures" / "glm" / "FENS_figures_poster"
# The histbase "_all" run outputs (oracle check + forward-selection history; the
# fits themselves come from the formatted .mat). History is a baseline nuisance
# here (history_in_baseline) — never selected, always in every model's baseline.
_GLM_DIR = HOME / "local_data" / "motion_clouds" / "figures" / "glm"
RUN_ALL_DIR = _GLM_DIR / "current_rf_sfor_20ms_histbase_goggles_all"

# Global rf_local SF range (cpd) across the whole goggles cohort extraction
# (parquet p0.5–p99.5 ≈ 0.027–0.130, full 0.021–0.151) — fixed SF y-axis so SF
# panels are comparable across trials/clusters.
SF_YLIM_CPD = (0.02, 0.15)

# The buildup follows the cluster's OWN forward-selection result: rows = the
# regressors forward selection chose, added one at a time in selection order
# (so row count varies per cluster). History, when selected, is folded into the
# baseline (Laura's "baseline = Intercept + Onset + History"); a cluster that did
# NOT select History gets a baseline of Intercept + Onset only. Onset is
# auto-added by the design assembler (include_onset_kernel=True).


def load_selection_order(probe_id: str, cluster_id: int) -> list[str] | None:
    """Ordered list of forward-selected terms for a cluster (selection order),
    read from the run's glm_selection_history.csv. None if unavailable."""
    csv = RUN_ALL_DIR / "_runs" / probe_id / "glm_selection_history.csv"
    if not csv.exists():
        return None
    sh = pd.read_csv(csv)
    sub = sh[sh["cluster_id"] == cluster_id].copy()
    if sub.empty:
        return None
    added = sub["added"].astype(str).str.lower() == "true"
    sub = sub[added].sort_values("round")
    return [str(c) for c in sub["best_candidate"].tolist() if str(c) not in ("", "nan")]


def buildup_steps_for(selected_ordered: list[str]) -> list[tuple[str, list[str]]]:
    """(label, cumulative-vars) steps. History (if selected) lives in the
    baseline; every other selected term becomes its own cumulative row, in
    selection order."""
    history = "History" in selected_ordered
    base_vars = ["History"] if history else []
    base_label = "baseline (Intercept+Onset" + ("+History)" if history else ")")
    steps = [(base_label, list(base_vars))]
    cum = list(base_vars)
    for v in selected_ordered:
        if v == "History":
            continue
        cum = cum + [v]
        steps.append((f"+ {v}", list(cum)))
    return steps


# --------------------------------------------------------------------------- #
# Per-cluster design prep — MIRROR of pipeline._fit_one_cluster.
# --------------------------------------------------------------------------- #
def prepare_cluster_design(df: pd.DataFrame, config) -> dict:
    """Build the bases / y / offset / folds for one binned cluster.

    Faithful mirror of pipeline._fit_one_cluster (the rf_local + history + ME +
    accel branches as configured by make_config_histbase_all). Returns a dict of
    everything assemble_design_matrix_selected / cross_validate_glm need.
    """
    motion_mask = (df["condition"] != "stationary").to_numpy(dtype=bool)

    speed = df["speed"].to_numpy(dtype=np.float64)
    tf = df["tf"].to_numpy(dtype=np.float64)
    onset = df["time_since_onset"].to_numpy(dtype=np.float64)
    sf_vals = df["sf"].to_numpy(dtype=np.float64)
    or_vals = df["orientation"].to_numpy(dtype=np.float64)
    y = df["spike_count"].to_numpy(dtype=np.float64)
    trial_ids = df["trial_id"].to_numpy(dtype=np.int64)
    condition_labels = df["condition"].to_numpy(dtype=object)
    profile_ids = (
        df["profile_id"].to_numpy(dtype=np.int64) if "profile_id" in df.columns else None
    )

    spacing = getattr(config, "speed_tf_basis_spacing", "log")
    B_speed = value_basis(speed, config.n_speed_bases, *config.speed_range, spacing=spacing)
    B_tf = value_basis(tf, config.n_tf_bases, *config.tf_range, spacing=spacing)
    B_onset = onset_kernel_basis(onset, config.n_onset_bases, config.onset_range[1])

    # History (identity 2-lag). Built when EITHER include_history (legacy
    # candidate mode) OR history_in_baseline (the histbase config: History is an
    # always-on baseline nuisance, never a candidate) — mirrors the pipeline,
    # which builds B_history under either flag. Missing this would silently drop
    # History from every model's baseline in the histbase fits.
    if getattr(config, "include_history", False) or getattr(config, "history_in_baseline", False):
        h_basis = history_basis(
            n_bases=config.n_history_bases,
            t_max_s=config.history_window_s,
            bin_width_s=config.time_bin_width,
            kind=getattr(config, "history_basis_kind", "raised_cosine"),
        )
        B_history = convolve_history(y, trial_ids, h_basis)
    else:
        B_history = None

    # Face ME (z-scored over finite motion rows).
    B_me_face = None
    if "me_face_raw" in df.columns and getattr(config, "include_me_face", True):
        me_raw = df["me_face_raw"].to_numpy(dtype=np.float64)
        me_motion_finite = np.isfinite(me_raw) & motion_mask
        if int(me_motion_finite.sum()) >= 10:
            me_mean = float(me_raw[me_motion_finite].mean())
            me_std = float(me_raw[me_motion_finite].std(ddof=0)) or 1.0
            me_z = np.where(np.isfinite(me_raw), (me_raw - me_mean) / me_std, 0.0)
            B_me_face = raised_cosine_basis_linear(
                me_z, config.n_me_face_bases, config.me_face_range[0], config.me_face_range[1]
            )

    # Acceleration (z-scored signed accel, clipped to accel_range).
    B_accel = None
    if "acceleration" in df.columns and getattr(config, "include_acceleration", False):
        acc_raw = df["acceleration"].to_numpy(dtype=np.float64)
        acc_motion = np.isfinite(acc_raw) & motion_mask
        if int(acc_motion.sum()) >= 10:
            a_mean = float(acc_raw[acc_motion].mean())
            a_std = float(acc_raw[acc_motion].std(ddof=0)) or 1.0
            acc_z = np.where(np.isfinite(acc_raw), (acc_raw - a_mean) / a_std, 0.0)
            acc_z = np.clip(acc_z, config.accel_range[0], config.accel_range[1])
            B_accel = raised_cosine_basis_linear(
                acc_z, config.n_accel_bases, config.accel_range[0], config.accel_range[1]
            )

    # RF-local SF/OR value bases.
    B_sf = B_or = None
    if getattr(config, "sf_or_source", "tokens") == "rf_local":
        n_rows = sf_vals.size
        sf_fin = np.isfinite(sf_vals)
        B_sf = np.zeros((n_rows, config.n_sf_bases), dtype=np.float64)
        if sf_fin.any():
            B_sf[sf_fin] = raised_cosine_basis_linear(
                sf_vals[sf_fin], config.n_sf_bases, *config.sf_cpd_range
            )
        or_fin = np.isfinite(or_vals)
        B_or = np.zeros((n_rows, config.n_or_bases), dtype=np.float64)
        if or_fin.any():
            B_or[or_fin] = circular_basis(or_vals[or_fin], config.n_or_bases)

    offset = float(np.log(config.time_bin_width))
    fold_ids = make_trial_folds(
        trial_ids,
        config.n_folds,
        config.cv_seed,
        condition_labels_per_bin=condition_labels,
        strategy=config.cv_strategy,
        profile_ids_per_bin=profile_ids,
    )

    sf_valid = sf_vals[(sf_vals != 0.0) & ~np.isnan(sf_vals)]
    sf_ref_levels = np.sort(np.unique(sf_valid)).tolist() if sf_valid.size > 0 else []
    or_valid = or_vals[(or_vals != 0.0) & ~np.isnan(or_vals)]
    or_ref_levels = np.sort(np.unique(or_valid)).tolist() if or_valid.size > 0 else []

    return dict(
        B_speed=B_speed, B_tf=B_tf, B_onset=B_onset, sf_vals=sf_vals, or_vals=or_vals,
        y=y, offset=offset, fold_ids=fold_ids, motion_mask=motion_mask,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
        B_history=B_history, B_me_face=B_me_face, B_accel=B_accel, B_sf=B_sf, B_or=B_or,
    )


def _design(prep: dict, vars_: list[str], config):
    return assemble_design_matrix_selected(
        prep["B_speed"], prep["B_tf"], prep["B_onset"], prep["sf_vals"], prep["or_vals"],
        vars_,
        sf_ref_levels=prep["sf_ref_levels"], or_ref_levels=prep["or_ref_levels"],
        B_history=prep["B_history"], B_me_face=prep["B_me_face"], B_accel=prep["B_accel"],
        B_sf=prep["B_sf"], B_or=prep["B_or"],
        include_onset_kernel=getattr(config, "include_onset_kernel", True),
        history_in_baseline=getattr(config, "history_in_baseline", False),
    )


def cv_bps_for(prep: dict, vars_: list[str], config, backend: str = "irls") -> float:
    X, _ = _design(prep, vars_, config)
    if X.shape[1] == 0 or X.shape[1] >= prep["y"].size:
        return float("nan")
    cv = cross_validate_glm(
        X, prep["y"], prep["offset"], prep["fold_ids"],
        lambda_ridge=config.lambda_ridge, backend=backend,
    )
    return float(cv.cv_bits_per_spike)


def cv_folds_for(prep: dict, vars_: list[str], config, backend: str = "irls") -> np.ndarray | None:
    """Per-fold cv bits/spike for a model (the 10 condition-stratified folds),
    on the SAME fold partition every model sees — so two models' fold arrays are
    paired and can be differenced fold-by-fold (exactly the signed-rank input).
    None when the model is degenerate (no columns / rank-deficient)."""
    X, _ = _design(prep, vars_, config)
    if X.shape[1] == 0 or X.shape[1] >= prep["y"].size:
        return None
    cv = cross_validate_glm(
        X, prep["y"], prep["offset"], prep["fold_ids"],
        lambda_ridge=config.lambda_ridge, backend=backend,
    )
    return np.asarray(cv.fold_bits_per_spike, dtype=np.float64)


def _paired_delta(f_model: np.ndarray | None, f_cur: np.ndarray | None) -> np.ndarray | None:
    """Fold-by-fold paired Δ bits/spike (model+candidate minus current model).
    None if either model is degenerate or the fold arrays don't align."""
    if f_model is None or f_cur is None or f_model.shape != f_cur.shape:
        return None
    return np.asarray(f_model, np.float64) - np.asarray(f_cur, np.float64)


def _delta_stats(d: np.ndarray | None) -> tuple[float, float, float]:
    """(mean, SD, one-sided Wilcoxon signed-rank p) of a paired per-fold Δ.
    SD is the across-fold spread (Laura's choice); the p-value is the same
    one-sided 'candidate improves bits/spike' test the production selection uses.
    """
    if d is None or d.size == 0 or not np.all(np.isfinite(d)):
        return float("nan"), float("nan"), float("nan")
    mean = float(np.mean(d))
    sd = float(np.std(d, ddof=1)) if d.size > 1 else 0.0
    try:
        from scipy.stats import wilcoxon
        p = float(wilcoxon(d, alternative="greater", zero_method="wilcox").pvalue)
    except Exception:  # all-zero/degenerate differences → no test
        p = float("nan")
    return mean, sd, p


def insample_rate_for(prep: dict, vars_: list[str], config, backend: str = "irls") -> np.ndarray:
    """Per-bin predicted firing rate (Hz) for an in-sample refit — ALL rows
    (stationary + motion), so the full trial window can be plotted."""
    X, _ = _design(prep, vars_, config)
    fit = fit_poisson_glm(
        X, prep["y"], prep["offset"], lambda_ridge=config.lambda_ridge, backend=backend
    )
    return np.exp(np.clip(X @ fit.beta, -20.0, 20.0))  # exp(X@beta) = count/bin = Hz


# --------------------------------------------------------------------------- #
# MATLAB FiringRate.get_convolution port (Gaussian, sigma=20 ms) — the
# "smoothed firing rate at a good timescale" property of the formatted-data
# class (lib/classes/analysis/FiringRate.m). Builds a spike train at the
# timebase fs, convolves fs*train with a normalised Gaussian (width=sigma,
# length=window), interpolates onto T. prepad/postpad avoid edge effects.
# --------------------------------------------------------------------------- #
def fr_convolution(
    spike_times: np.ndarray, T: np.ndarray,
    width: float = 0.02, length: float = 1.0, pad: float = 1.0,
) -> np.ndarray:
    fs = 1.0 / (T[1] - T[0])
    t0, t1 = T[0] - pad, T[-1] + pad
    st = spike_times[(spike_times > t0) & (spike_times < t1)]
    n = int(round((t1 - t0) * fs)) + 1
    train = np.zeros(n, dtype=np.float64)
    idx = np.round((st - t0) * fs).astype(int)
    idx = idx[(idx >= 0) & (idx < n)]
    np.add.at(train, idx, 1.0)
    sigma = width * fs
    sz = int(round(length * fs))
    x = np.linspace(-sz / 2, sz / 2, sz)
    g = np.exp(-x**2 / (2 * sigma**2))
    g /= g.sum()
    conv = np.convolve(fs * train, g, mode="same")
    sr_t = t0 + np.arange(n) / fs
    return np.interp(T, sr_t, conv)


def select_vt_trials(
    df: pd.DataFrame, n: int, baseline_min: int, baseline_max: int,
) -> list[tuple[int, int, int]]:
    """Top-n VT trials with 'a few' (baseline_min..baseline_max) stationary
    spikes, ranked by motion-spike count. Returns (trial_id, motion, baseline)."""
    stat = df[df["condition"] == "stationary"].groupby("trial_id")["spike_count"].sum()
    mot = df[df["condition"] == "VT"].groupby("trial_id")["spike_count"].sum()
    cand = pd.DataFrame({"motion": mot})
    cand["baseline"] = stat.reindex(cand.index).fillna(0).astype(int)
    ok = cand[(cand["baseline"] >= baseline_min) & (cand["baseline"] <= baseline_max)]
    ok = ok.sort_values("motion", ascending=False).head(n)
    return [(int(t), int(r.motion), int(r.baseline)) for t, r in ok.iterrows()]


def select_v_trial(
    df: pd.DataFrame, baseline_max: int = 10, condition: str = "V",
) -> tuple[int, int, int] | None:
    """Pick a representative trial of ``condition`` (V or VT) for the
    trial-structure schematic: rank that condition's trials by motion-spike count
    (so the pink FR is illustrative), preferring a sane handful of baseline spikes
    (≤baseline_max) when available. Returns (trial_id, motion_spk, baseline_spk)
    or None."""
    mot = df[df["condition"] == condition].groupby("trial_id")["spike_count"].sum()
    if mot.empty:
        return None
    stat = df[df["condition"] == "stationary"].groupby("trial_id")["spike_count"].sum()
    cand = pd.DataFrame({"motion": mot})
    cand["baseline"] = stat.reindex(cand.index).fillna(0).astype(int)
    cand = cand.sort_values("motion", ascending=False)
    pref = cand[cand["baseline"] <= baseline_max]
    pick = pref if not pref.empty else cand
    t = pick.index[0]
    return int(t), int(pick.iloc[0]["motion"]), int(pick.iloc[0]["baseline"])


def select_matched_cloud_trials(df, trials_by_id):
    """Pick a V trial and a VT trial that show the SAME motion cloud, so the two
    figures are directly comparable. Each goggles cloud is presented in both V and
    VT (2 repeats each); we pick the cloud where this cluster is most active across
    V+VT, then the most-active V and VT trial of that cloud. Returns
    (cloud_name, v_tid, vt_tid, v_motion_spk, vt_motion_spk) or None when no cloud
    has both a V and a VT trial."""
    cloud = {int(tid): getattr(tr, "cloud_name", None)
             for tid, tr in trials_by_id.items()}
    spk = (df[df["condition"].isin(("V", "VT"))]
           .groupby(["trial_id", "condition"])["spike_count"].sum())
    best: dict[str, dict[str, tuple[int, int]]] = {}
    for (tid, cond), s in spk.items():
        c = cloud.get(int(tid))
        if c is None:
            continue
        d = best.setdefault(c, {})
        if cond not in d or int(s) > d[cond][1]:
            d[cond] = (int(tid), int(s))
    cands = [(c, d["V"], d["VT"]) for c, d in best.items() if "V" in d and "VT" in d]
    if not cands:
        return None
    cands.sort(key=lambda x: x[1][1] + x[2][1], reverse=True)
    c, (v_tid, v_spk), (vt_tid, vt_spk) = cands[0]
    return c, v_tid, vt_tid, v_spk, vt_spk


def rank_clusters_by_firing_rate(probe, by_id, rf_lookup, config):
    """Rank the probe's selected cohort by overall mean firing rate, for browsing.

    Ranking by tuning R² surfaced near-silent cells whose 'tuning' was a few
    spikes (cluster 360: 3 Hz peak, p=0.08) — so instead pick the most ACTIVE
    cells, which have meaningful tuning curves to read. Score = mean FR
    (spike_count / bin_width) over all the cluster's motion bins. Returns a list
    of (cluster_id, mean_fr_hz) sorted descending."""
    bw = float(config.time_bin_width)
    scored = []
    for cid, cluster in by_id.items():
        df = bin_cluster(probe, cluster, rf_lookup=rf_lookup)
        if df.empty:
            continue
        mot = df[df["condition"] != "stationary"]
        if mot.empty:
            continue
        mean_fr = float(mot["spike_count"].mean()) / bw
        scored.append((cid, mean_fr))
    scored.sort(key=lambda kv: kv[1], reverse=True)
    return scored


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Observed FR-vs-value tuning (right block of the trial-structure figure).
# Delegates the binning to rc2_glm.tuning_significance.per_trial_bin_matrix —
# the ONE shared path (20 equal-count 5%-quantile bins, per-trial mean FR) the
# significance test also uses, parametrised by CONDITION (V or VT). From that
# (n_trials × n_bins) matrix we derive the per-bin spread (mean — the model's
# fit target — plus median and Q1/Q3 IQR); the model + p come from
# tuning_significance() on the same matrix. Returns None when there is nothing
# to bin (no rows in the condition, missing column, degenerate value range).
# --------------------------------------------------------------------------- #
def _tuning_stats(matrix, centres, source):
    """Per-bin mean/SD/median/IQR across trials from a (n_trials × n_bins) matrix."""
    # Pooled/empty bins → all-NaN columns; silence the expected empty-slice
    # warning, matching plots.py.
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(matrix, axis=0)
        sd = np.nanstd(matrix, axis=0)
        n_per_bin = np.sum(np.isfinite(matrix), axis=0)
        sem = sd / np.sqrt(np.where(n_per_bin > 0, n_per_bin, np.nan))
        median = np.nanmedian(matrix, axis=0)
        q1 = np.nanquantile(matrix, 0.25, axis=0)
        q3 = np.nanquantile(matrix, 0.75, axis=0)
    return dict(centres=np.asarray(centres, float), matrix=matrix, mean=mean,
                sd=sd, sem=sem, median=median, q1=q1, q3=q3,
                n_trials=matrix.shape[0], source=source)


def _observed_value_tuning(df, value_col, key, cond, bw, *, cluster_id=None,
                           pc=None, n_bins=20):
    """Per-trial-per-bin tuning stats for (value, condition).

    For **TF** (and Speed), prefer the MATLAB precomputed cache (``pc``) — its
    20 quantile bins + per-trial FR (own per-bin duration denominator) are what
    MATLAB's tuning-curve PDFs plot, so they MATCH (the convention,
    reference_motion_clouds_tuning_curve_20bins). Recomputing quantile bins from
    20 ms spike counts does NOT match. For **SF/OR** there is no MATLAB cache
    (rf_local Gabor values) → recompute via the shared per_trial_bin_matrix."""
    if pc is not None and cluster_id is not None and key in ("tf", "speed", "accel"):
        get_tun, get_cen = {
            "tf": (pc.tf_tuning, pc.tf_centres),
            "speed": (pc.speed_tuning, pc.speed_centres),
            "accel": (pc.accel_tuning, pc.accel_centres),
        }[key]
        matrix, centres = get_tun(cond, cluster_id), get_cen(cond)
        if matrix is not None and centres is not None:
            return _tuning_stats(matrix, centres, source="matlab_cache")
    # SF/OR (rf_local, no MATLAB cache): recompute, but with edges POOLED across
    # the visual conditions (V+VT) so bin k is the same interval in both columns
    # — the project convention (_pooled_quantile_edges / _plot_sf_or_scatter).
    from rc2_glm.plots import _pooled_quantile_edges
    pooled = np.array([])
    if value_col in df.columns:
        pooled = df.loc[df["condition"].isin(("V", "VT")), value_col].to_numpy(float)
        pooled = pooled[np.isfinite(pooled)]
    edges = None
    if pooled.size >= n_bins:
        edges, _ = _pooled_quantile_edges(pooled, n_bins=n_bins)
    matrix, centres = per_trial_bin_matrix(
        df, value_col, bw, condition=cond, n_bins=n_bins, edges=edges)
    if matrix is None:
        return None
    return _tuning_stats(matrix, centres, source="recomputed")


# Tuning rendering shared by the trial-structure column and the tuning-only grid.
TUN_XLABEL = {"tf": "TF (Hz)", "sf": "SF (cpd)", "or": "orientation (deg)",
              "speed": "speed (cm/s)", "accel": "acceleration (cm/s²)"}
ASYM_ONLY = ("asym_gaussian",)   # TF/SF: asymmetric Gaussian only
LINEAR_ONLY = ("linear",)        # TF/SF: straight-line fit only
# Default TF/SF family set = the full MATLAB ModelSelectionTuning classes
# (linear/quadratic/cubic/Gaussian/asymmetric-Gaussian/sigmoid), BIC-selected.
TUN_LINEAR_FAMILIES = LINEAR_FAMILIES_TS  # imported from rc2_glm.tuning_significance


def _render_tuning_panel(ax, df, value_col, key, cond, bw, col, *,
                         probe_id=None, cluster_id=None, pc=None, aggregate="flat",
                         select_criterion="bic", display="median_iqr",
                         fits_lookup=None, tuning_n_reps=1000,
                         linear_families=TUN_LINEAR_FAMILIES):
    """Draw one observed-tuning + best-fit-model panel for (value, condition).

    Observed points (TF from the MATLAB cache when ``pc`` is given; SF/OR
    recomputed) are connected by a line with error bars chosen by ``display``:
    ``"median_iqr"`` (median + asymmetric Q1/Q3 IQR), ``"mean_sd"`` (mean ± SD),
    or ``"mean_sem"`` (mean ± SEM = SD/√n per bin). The model (black) is the
    BIC-best of ``linear_families`` for TF/SF, or the von Mises for OR, with its
    bootstrap-null p. When ``fits_lookup`` (a ``{(cluster_id, value, condition):
    row}`` dict from a prior ``tuning_fits.csv``) supplies a row, the fit + p are
    REUSED from it — no re-fit / re-bootstrap (the fits are deterministic). Only
    the observed binning is redone (cheap). Returns the TuningSignificance, or
    None when there is no data."""
    kind = "circular" if key == "or" else "linear"
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(TUN_XLABEL[key], fontsize=8.5)
    t = _observed_value_tuning(df, value_col, key, cond, bw,
                               cluster_id=cluster_id, pc=pc)
    if t is None:
        ax.text(0.5, 0.5, f"no {cond}\ntuning data", ha="center", va="center",
                transform=ax.transAxes, color="#888", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        return None
    cen = t["centres"]
    if display in ("mean_sd", "mean_sem"):
        centre, gm = t["mean"], np.isfinite(t["mean"])
        yerr = (t["sem"] if display == "mean_sem" else t["sd"])[gm]  # symmetric
    else:
        centre, gm = t["median"], np.isfinite(t["median"])
        yerr = np.vstack([t["median"][gm] - t["q1"][gm],
                          t["q3"][gm] - t["median"][gm]])  # asymmetric IQR
    ax.errorbar(cen[gm], centre[gm], yerr=yerr,
                fmt="o-", color=col, ms=3.5, lw=1.1, capsize=2.5,
                elinewidth=0.8, zorder=2)
    cached = (None if fits_lookup is None
              else fits_lookup.get((probe_id, cluster_id, key, cond)))
    if cached is not None:
        sig = _sig_from_cache(cached)  # reuse stored fit + p, skip the bootstrap
        if np.isnan(sig.rsq_mean) and sig.params is not None:
            # old CSV without the column → recompute the (cheap) mean-curve R².
            sig.rsq_mean = rsq_against_mean(t["matrix"], cen, sig.best_model, sig.params)
    else:
        sig = tuning_significance(t["matrix"], cen, value=key, condition=cond,
                                  kind=kind, aggregate=aggregate,
                                  select_criterion=select_criterion,
                                  n_reps=tuning_n_reps, linear_families=linear_families)
        sig.data_source = t["source"]  # matlab_cache (TF) | recomputed (SF/OR)
    if sig.best_model is not None:
        xs = np.linspace(float(np.nanmin(cen)), float(np.nanmax(cen)), 200)
        ax.plot(xs, eval_tuning_fit(
                    dict(name=sig.best_model, params=sig.params), xs),
                color="black", lw=1.6, zorder=3)
        if np.isfinite(sig.p):
            pstr = "p<0.001" if sig.p < 1e-3 else f"p={sig.p:.3f}"
        else:
            pstr = "p=n/a"
        # R²(mean) = fit-vs-mean-curve; R²(trials) = across all per-trial points.
        ax.text(0.04, 0.96,
                f"{sig.best_model}\nR²mean={sig.rsq_mean:.2f} · "
                f"R²trials={sig.rsq:.3f}\n{pstr}",
                transform=ax.transAxes, ha="left", va="top", fontsize=6.0,
                color="black", zorder=5,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", alpha=0.85))
    return sig


# --------------------------------------------------------------------------- #
# Tuning-only browsing figure: rows TF/SF/OR × cols V|VT, observed tuning + the
# best-fit model. No time-series block — for picking which cells carry tuning.
# --------------------------------------------------------------------------- #
def _sig_from_cache(row):
    """Rebuild a TuningSignificance from a cached tuning_fits.csv row (for reuse
    of the deterministic fit + bootstrap p, skipping recomputation)."""
    bm = row.get("best_model")
    if bm is None or (isinstance(bm, float) and np.isnan(bm)) or str(bm) == "nan":
        bm, params = None, None
    else:
        bm = str(bm)
        ps = str(row.get("params", "") or "")
        params = (np.array([float(x) for x in ps.split(";")], float) if ps else None)
    return TuningSignificance(
        value=str(row["value"]), condition=str(row["condition"]),
        kind=str(row["kind"]), best_model=bm, params=params,
        rsq=float(row["rsq"]), rsq_mean=float(row.get("rsq_mean", np.nan)),
        bic=float(row["bic"]), p=float(row["p"]),
        aggregate=str(row.get("aggregate", "")),
        select_criterion=str(row.get("select_criterion", "bic")),
        data_source=str(row.get("data_source", "")),
        n_trials=int(row.get("n_trials", 0)), n_bins=int(row.get("n_bins", 0)),
        n_reps=int(row.get("n_reps", 0)), seed=int(row.get("seed", 1)),
        null_scheme=str(row.get("null_scheme", "")),
    )


def _load_fits_lookup(csv_path):
    """``{(probe_id, cluster_id, value, condition): row}`` from a tuning_fits.csv.

    Keyed on probe_id too: cluster IDs are NOT unique across the two goggles
    probes (6 collide), so omitting probe_id would cross-wire those cells."""
    fits = pd.read_csv(csv_path)
    return {(str(r["probe_id"]), int(r["cluster_id"]),
             str(r["value"]), str(r["condition"])): r
            for r in fits.to_dict("records")}


def _sig_row(probe_id, cluster_id, sig):
    """One flat record (the BIC-selected fit) for the tuning-fits CSV."""
    params = ("" if sig.params is None
              else ";".join(f"{v:.6g}" for v in np.asarray(sig.params, float)))
    return dict(
        probe_id=probe_id, cluster_id=cluster_id, value=sig.value,
        condition=sig.condition, kind=sig.kind, aggregate=sig.aggregate,
        select_criterion=sig.select_criterion,
        data_source=sig.data_source, best_model=sig.best_model,
        n_params=(0 if sig.params is None else int(np.asarray(sig.params).size)),
        params=params, rsq=sig.rsq, rsq_mean=sig.rsq_mean, bic=sig.bic, p=sig.p,
        n_trials=sig.n_trials, n_bins=sig.n_bins,
        null_scheme=sig.null_scheme, n_reps=sig.n_reps, seed=sig.seed,
    )


# ---------------------------------------------------------------------------
# Population tuning-significance summaries (read tuning_fits.csv; no re-fit).
# ---------------------------------------------------------------------------
# Cohort-level views over the per-cluster tuning-significance fits: one strip
# per variable (V tuned vs not-tuned), and the V→VT transition of that
# significance call. Both consume an existing tuning_fits.csv — they launch no
# compute. "tuned" is the same p<0.05 call used per cluster.
TUNING_VARS = (("tf", "TF"), ("sf", "SF"), ("or", "OR"))
TUNED_ALPHA = 0.05
_C_TUNED, _C_NOT = "#d62728", "#4c72b0"   # tuned (red) / not-tuned (blue)
_C_V, _C_VT = "#4c72b0", "#dd8452"        # V (blue) / VT (orange) in the pair


def tuning_transition_category(tuned_v, tuned_vt, model_v, model_vt):
    """V→VT transition class for one (cluster, variable) significance pair.

    1 = tuned in V, not in VT;  2 = not in V, tuned in VT;
    3 = tuned in both, *different* selected best_model;
    4 = tuned in both, *same* selected best_model;
    0 = not tuned in either (excluded from the transition figure).
    Same/different is exact best_model string identity. Pure — unit-tested."""
    tuned_v, tuned_vt = bool(tuned_v), bool(tuned_vt)
    if tuned_v and not tuned_vt:
        return 1
    if not tuned_v and tuned_vt:
        return 2
    if tuned_v and tuned_vt:
        return 4 if str(model_v) == str(model_vt) else 3
    return 0


def _p_floor(values):
    """Half the smallest positive p — a display floor so exact-zero bootstrap
    p's (0 of n_reps beat the real fit) are visible on a log axis. Display
    only; never used for the p<0.05 classification."""
    pos = np.asarray(values, float)
    pos = pos[pos > 0]
    return (pos.min() / 2.0) if pos.size else 1e-3


def plot_tuning_significance_summary(fits_df, out_path, *, condition="V"):
    """2×3 strip plot for one condition: tuned vs not-tuned p per TF/SF/OR.

    Top row = all clusters; bottom row = only rsq_mean above that variable's
    median (good-fit subset). y is the bootstrap tuning p on a log axis."""
    v = fits_df[(fits_df.condition == condition) & fits_df.p.notna()].copy()
    v["tuned"] = v.p < TUNED_ALPHA
    floor = _p_floor(v.p.values)
    v["p_plot"] = v.p.clip(lower=floor)
    med = {k: v[v.value == k].rsq_mean.median() for k, _ in TUNING_VARS}
    cats = [("tuned\n(p<%.2g)" % TUNED_ALPHA, True, _C_TUNED),
            ("not tuned", False, _C_NOT)]
    rng = np.random.default_rng(0)

    def panel(ax, s, title):
        for i, (_lab, tuned, col) in enumerate(cats):
            ss = s[s.tuned == tuned]
            x = i + rng.uniform(-0.18, 0.18, len(ss))
            ax.scatter(x, ss.p_plot, s=24, c=col, alpha=0.7, edgecolors="none")
            ax.text(i, 1.5, f"n={len(ss)}", ha="center", fontsize=9)
        ax.set_yscale("log")
        ax.axhline(TUNED_ALPHA, ls="--", lw=1, color="k", alpha=0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([c[0] for c in cats])
        ax.set_xlim(-0.5, 1.5)
        ax.set_title(title, pad=22)

    fig, axes = plt.subplots(2, 3, figsize=(9.5, 8.8), sharey=True)
    for j, (key, lab) in enumerate(TUNING_VARS):
        s = v[v.value == key]
        panel(axes[0, j], s, f"{lab}  (all, n={len(s)})")
        sf = s[s.rsq_mean > med[key]]
        panel(axes[1, j], sf, f"{lab}  (rsq>{med[key]:.2f}, n={len(sf)})")
    for r in range(2):
        axes[r, 0].set_ylabel("tuning p-value (log)")
    fig.suptitle(
        f"{condition} condition — tuning significance  "
        f"(top: all  |  bottom: rsq_mean > per-variable median)", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)


def plot_tuning_transitions(fits_df, out_path):
    """Paired V→VT significance transitions, one subplot per TF/SF/OR.

    Each cluster pairs its V (left) and VT (right) tuning-p, joined by a line,
    grouped into the four transition categories. Subplot title carries the
    VT-tuned count. Excludes cells not tuned in either condition (cat 0)."""
    d = fits_df[fits_df.p.notna()].copy()
    d["tuned"] = d.p < TUNED_ALPHA
    key = ["probe_id", "cluster_id", "value"]
    piv = d.pivot_table(index=key, columns="condition",
                        values=["p", "tuned", "best_model"], aggfunc="first")
    piv.columns = [f"{a}_{b}" for a, b in piv.columns]
    piv = piv.reset_index().dropna(subset=["p_V", "p_VT"])
    piv["tuned_V"] = piv["tuned_V"].astype(bool)
    piv["tuned_VT"] = piv["tuned_VT"].astype(bool)
    piv["cat"] = [tuning_transition_category(r.tuned_V, r.tuned_VT,
                                             r.best_model_V, r.best_model_VT)
                  for r in piv.itertuples()]
    floor = _p_floor(np.r_[piv.p_V.values, piv.p_VT.values])
    clip = lambda a: np.clip(a, floor, None)  # noqa: E731
    catlabels = {1: "tuned→\nnot tuned", 2: "not tuned\n→tuned",
                 3: "tuned→tuned\ndiff model", 4: "tuned→tuned\nsame model"}
    gap, step = 0.32, 1.4
    rng = np.random.default_rng(0)

    fig, axes = plt.subplots(1, 3, figsize=(13, 5.2), sharey=True)
    for ax, (key_, lab) in zip(axes, TUNING_VARS):
        s = piv[piv.value == key_]
        n_vt = int(s.tuned_VT.sum())
        centers = []
        for k, c in enumerate((1, 2, 3, 4)):
            sc = s[s.cat == c]
            xV, xVT = k * step - gap / 2, k * step + gap / 2
            centers.append(k * step)
            pv, pvt = clip(sc.p_V.values), clip(sc.p_VT.values)
            jV = rng.uniform(-0.05, 0.05, len(sc))
            jVT = rng.uniform(-0.05, 0.05, len(sc))
            for a, b, ya, yb in zip(xV + jV, xVT + jVT, pv, pvt):
                ax.plot([a, b], [ya, yb], color="0.7", lw=0.5, alpha=0.6, zorder=1)
            ax.scatter(xV + jV, pv, s=20, c=_C_V, alpha=0.8, edgecolors="none", zorder=2)
            ax.scatter(xVT + jVT, pvt, s=20, c=_C_VT, alpha=0.8, edgecolors="none", zorder=2)
            ax.text(k * step, 1.6, f"n={len(sc)}", ha="center", fontsize=9)
        ax.set_yscale("log")
        ax.axhline(TUNED_ALPHA, ls="--", lw=1, color="k", alpha=0.6)
        ax.set_xticks(centers)
        ax.set_xticklabels([catlabels[c] for c in (1, 2, 3, 4)], fontsize=8)
        ax.set_title(f"{lab}  ({n_vt}/{len(s)} tuned in VT)", pad=20)
        ax.set_xlim(-step / 2, 3 * step + step / 2)
    axes[0].set_ylabel("tuning p-value (log)")
    axes[-1].legend(
        handles=[Line2D([], [], marker="o", ls="", color=_C_V, label="V"),
                 Line2D([], [], marker="o", ls="", color=_C_VT, label="VT")],
        frameon=False, loc="lower right")
    fig.suptitle("V→VT tuning-significance transitions "
                 "(paired: V left, VT right)", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)


ERR_LABEL = {"median_iqr": "median + IQR", "mean_sd": "mean ± SD",
             "mean_sem": "mean ± SEM"}


def plot_tuning_grid(probe_id, cluster_id, df, config, out_path, *,
                     pc=None, display="median_iqr", select_criterion="bic",
                     fits_lookup=None, tuning_n_reps=1000,
                     linear_families=TUN_LINEAR_FAMILIES):
    """Render the TF/SF/OR × V|VT grid; return the per-panel fit rows (for CSV)."""
    bw = float(config.time_bin_width)
    rows = [("tf", "tf", "#d1701a"), ("sf", "sf", "#6b8e23"),
            ("or", "orientation", "#c0392b")]
    conds = ("V", "VT")
    fig, axes = plt.subplots(len(rows), len(conds), figsize=(7.5, 8.6),
                             constrained_layout=True)
    axes = np.atleast_2d(axes)
    fit_rows = []
    for r, (key, vcol, col) in enumerate(rows):
        for c, cond in enumerate(conds):
            ax = axes[r, c]
            sig = _render_tuning_panel(
                ax, df, vcol, key, cond, bw, col,
                probe_id=probe_id, cluster_id=cluster_id, pc=pc, display=display,
                select_criterion=select_criterion, fits_lookup=fits_lookup,
                tuning_n_reps=tuning_n_reps, linear_families=linear_families)
            if sig is not None:
                fit_rows.append(_sig_row(probe_id, cluster_id, sig))
            if r == 0:
                ax.set_title(cond, fontsize=11, fontweight="bold")
            if c == 0:
                ax.set_ylabel("FR (Hz)", fontsize=8)
    probe_short = probe_id.split("_rec")[0]
    fam = (linear_families[0] if len(linear_families) == 1 else "BIC-best")
    err_lbl = ERR_LABEL.get(display, display)
    fig.suptitle(
        f"{probe_short} cl {cluster_id} · observed FR tuning ({err_lbl} "
        f"error bars, 20 equal-count bins) + {fam} fit / von Mises (OR) · "
        f"R² across all trials (MATLAB way) · bootstrap-p",
        fontsize=8.5, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)
    return fit_rows


def run_tuning_browse(probe_id, config, out_dir, *, all_clusters=False,
                      top_fr=None, display="median_iqr", select_criterion="bic",
                      fits_lookup=None, tuning_n_reps=1000):
    """Load ONE goggles probe and emit per-cluster tuning grids into ``out_dir``.

    Cohort = the formatted .mat ``selected_clusters`` (``cluster_set="selected"``,
    the 49/57 goggles cohort — NOT Kilosort/VISp). With ``all_clusters`` every
    selected cluster is rendered; otherwise the top-``top_fr`` by mean firing
    rate. ``display`` picks the error bars; ``fits_lookup`` (from a prior
    tuning_fits.csv) reuses the cached fits instead of recomputing."""
    log.info("loading %s ...", probe_id)
    probe = load_probe_data(
        FORMATTED_DIR / f"{probe_id}.mat", config=config,
        stimulus_lookup=_lookup(), cluster_set="selected",
    )
    rf_lookup = load_rf_sf_or(
        config.rf_sf_or_parquet_dir, probe.probe_id,
        min_concentration=getattr(config, "rf_min_concentration", 0.0),
    )
    by_id = {c.cluster_id: c for c in probe.clusters}
    pc = load_precomputed_bin_edges(FORMATTED_DIR / f"{probe_id}.mat")

    ranked = rank_clusters_by_firing_rate(probe, by_id, rf_lookup, config)
    if all_clusters:
        cids = [cid for cid, _fr in ranked]  # every selected cluster (FR order)
        log.info("%s: all %d selected clusters (FR-ordered)", probe_id, len(cids))
    else:
        cids = [cid for cid, _fr in ranked[:top_fr]]
        log.info("%s: top %d/%d clusters by mean FR: %s", probe_id, len(cids),
                 len(ranked), ", ".join(f"{c}({fr:.1f}Hz)"
                                        for c, fr in ranked[:top_fr]))
    fit_rows = []
    for cid in cids:
        df = bin_cluster(probe, by_id[cid], rf_lookup=rf_lookup)
        if df.empty:
            log.warning("cluster %d: empty binned df — skipping", cid)
            continue
        fit_rows += plot_tuning_grid(
            probe.probe_id, cid, df, config,
            out_dir / f"tuning_{probe.probe_id}_cluster_{cid}",
            pc=pc, display=display, select_criterion=select_criterion,
            fits_lookup=fits_lookup, tuning_n_reps=tuning_n_reps,
        )
    return fit_rows


# Poster column 1 — TRIAL-STRUCTURE schematic (mirrors the MATLAB trial
# layout). For one V (ReplayOnly) trial of a defined-RF cluster, a stacked
# column on a SHARED onset-relative time axis (x = bin_centre − motion_onset,
# i.e. df.time_since_onset): baseline (negative x) → the small registered gap →
# motion (positive x), spanning >8 s. Rows top→bottom:
#   VF    visual-flow velocity (trial.velocity = the replayed multiplexer
#         velocity in a V trial), cm/s;
#   T     treadmill speed — flat 0 by design in a ReplayOnly trial (schematic);
#   TF    per-bin temporal frequency (Hz), = gain × replay-speed;
#   SF    per-bin RF-local spatial frequency (cpd, Gabor extraction — varies
#         over the trial, hence a defined-RF cluster);
#   OR    per-bin RF-local orientation (deg, circular);
#   FR    the pink MATLAB Gaussian-σ20 ms convolved firing rate (fr_convolution).
# Stationary vs motion are shaded from the binned mask extents; the unshaded
# sliver between them is the real stationary↔motion mask gap (NOT faked as
# contiguous, unlike the fig2c synthetic axis). Cluster is fully parametrised
# (--clusters / --trial) so swapping the example is a one-flag change.
# --------------------------------------------------------------------------- #
def plot_trial_structure(
    probe_id, cluster_id, df, config, tid, trial, spike_times, out_path,
    trials_by_id, tuning_n_reps=1000, pc=None, condition="V",
    tuning_clusters=None, save=True, fr_stat="median",
):
    # ``tuning_clusters`` = [(cid, cluster_df), ...] for the RIGHT tuning grid
    # (one column per cluster). The LEFT trial-structure block (incl. the pink FR)
    # is always the PRIMARY cluster (``cluster_id``/``df``/``spike_times``). Default
    # = a single column for the primary. ``save=False`` returns ``(fig, tun_axes)``
    # without writing, so the caller can harmonise FR y-limits across the V and VT
    # figures (the per-regressor tuning y-axis is shared across cluster columns and
    # between conditions) before saving. ``tun_axes`` maps regressor-key → [axes].
    tuning_clusters = list(tuning_clusters) if tuning_clusters else [(cluster_id, df)]
    # SOLENOID-aligned continuous axis: x = 0 at the velocity-command (solenoid)
    # onset. Stationary < 0; the REAL gap [0, gap] command→photodiode-onset drawn
    # TO SCALE; motion [gap, x_mend]. The command→motion-end interval is fixed
    # (~4.74 s, ±0.025) so command-alignment also lines up the motion END across
    # trials; only the photodiode visual onset scatters (the velocity-profile-
    # dependent gap), which is now shown honestly. SF/OR are drawn ONLY inside the
    # photodiode visual window (where the cloud actually drives the cell).
    sub = df[df["trial_id"] == tid].copy()
    if sub.empty:
        log.warning("cluster %d trial %d: no binned rows — skipping", cluster_id, tid)
        return
    sub = sub.iloc[np.argsort(sub["time_since_onset"].to_numpy(float))]
    tso = sub["time_since_onset"].to_numpy(float)  # 0 == photodiode visual onset
    tf = sub["tf"].to_numpy(float)
    sf = sub["sf"].to_numpy(float)
    # rf_local orientation is already stored in DEGREES (0–180), same stimulus
    # convention as the token (verified: token 135° vs rf_local circ-mean 145°).
    # Render it cl90-style: TOKEN-CENTRED (OR − token on a ±90° axis) with the
    # line SEAM-BROKEN at the 0/180 wrap so circular jumps don't draw as streaks.
    or_raw = sub["orientation"].to_numpy(float)
    or_token = float(np.degrees(trial.orientation)) % 180.0
    or_real = or_raw % 180.0  # REAL orientation degrees (0–180); token = ref line
    sf_token_cpd = float(trial.sf) * 9.77  # cpp token × screen calib (units note)

    def _seam_break(x, y, thr=90.0):
        """cl90 break_wrap: NaN-insert where the line jumps the circular seam."""
        y = np.asarray(y, float)
        j = np.where(np.abs(np.diff(y)) > thr)[0]
        return np.insert(np.asarray(x, float), j + 1, np.nan), np.insert(y, j + 1, np.nan)

    def _periods(tr):
        """Command(solenoid)-aligned anchors: time from the velocity-command
        onset (x=0). gap = command→photodiode-onset (real); x_vis = gap;
        x_mend = motion end ≈ 4.74 s (stable); xs_start = stationary start (<0)."""
        pt = np.asarray(tr.probe_t, np.float64)
        midx = np.flatnonzero(np.asarray(tr.motion_mask, bool))
        t_vis = float(pt[int(midx[0])])
        has_g = (getattr(tr, "command_onset_idx", None) is not None
                 and getattr(tr, "visual_onset_idx", None) is not None)
        t_cmd = float(pt[tr.command_onset_idx]) if has_g else t_vis
        cut = int(tr.command_onset_idx) if has_g else int(midx[0])
        gap = t_vis - t_cmd
        mend_vis = (float(pt[tr.visual_offset_idx]) - t_vis) if has_g \
            else float(pt[int(midx[-1])] - t_vis)
        pre = np.flatnonzero(np.asarray(tr.stationary_mask, bool))
        pre = pre[pre < cut]
        if pre.size:
            brk = np.flatnonzero(np.diff(pre) > 1)
            s_idx = int(pre[brk[-1] + 1]) if brk.size else int(pre[0])
            xs_start = float(pt[s_idx]) - t_cmd
        else:
            xs_start = 0.0
        return dict(t_cmd=t_cmd, gap=gap, x_vis=gap, x_mend=gap + mend_vis,
                    xs_start=xs_start, has_gap=has_g)

    sel = _periods(trial)

    # --- The cluster's trials of THIS condition (V or VT), for the across-trial
    # overlay (col 2). ---
    ov_ids = pd.unique(df.loc[df["condition"] == condition, "trial_id"])
    ov_trials = [(int(vt), trials_by_id.get(int(vt))) for vt in ov_ids]
    ov_trials = [
        (vt, tr) for vt, tr in ov_trials
        if tr is not None and np.flatnonzero(np.asarray(tr.motion_mask, bool)).size
    ]
    ov_sub = {vt: sub_v.iloc[np.argsort(sub_v["time_since_onset"].to_numpy(float))]
              for vt, _ in ov_trials
              for sub_v in [df[df["trial_id"] == vt]]}
    n_ov = len(ov_trials)
    per = {vt: _periods(tr) for vt, tr in ov_trials}

    # Command-aligned window (0 = solenoid command onset).
    x_lo = sel["xs_start"] - 0.2
    x_hi = sel["x_mend"] + 0.3
    log.info("cluster %d %s trial %d (solenoid-aligned): stationary %.2fs | "
             "gap(real) %.2fs | motion %.2fs | x_mend(motion end rel command) %.2fs",
             cluster_id, condition, tid, -sel["xs_start"], sel["gap"],
             sel["x_mend"] - sel["gap"], sel["x_mend"])

    def _vel_source(key, tr):
        """Probe_t-sampled source trace for the two velocity rows, or None when
        the row is the flat ``≡ 0`` line. VF = the visual command
        (``multiplexer_output``); T = the stage TRANSLATION (the ``stage``
        channel, which IS ``velocity`` for a StageOnly/VT trial). For V
        (ReplayOnly) the stage is ~0 so the T row is drawn flat, and VF reads the
        protocol ``velocity`` (= the multiplexer command)."""
        if key == "vf":
            src = tr.visual_velocity if condition == "VT" else tr.velocity
            return np.asarray(src, float)
        if key == "t":
            if condition == "V":
                return None                          # replay: no translation
            return np.asarray(tr.velocity, float)    # VT: stage translation
        return None

    def _xy(key, vt, tr):
        """One trial's command-aligned (x, y) for a row. SF/OR/TF are restricted
        to the visual window (tso ≥ 0); the velocity rows span probe_t."""
        vs = ov_sub[vt]
        p = per[vt]
        vtso = vs["time_since_onset"].to_numpy(float)
        x = p["gap"] + vtso                       # bins → command-aligned
        vis = vtso >= -1e-9                        # only where the visual cloud is on
        if key == "tf":
            return x[vis], vs["tf"].to_numpy(float)[vis]
        if key in ("sf", "or"):
            # SF/OR are defined only inside the visual window, BUT during the
            # stationary baseline the cloud sits at its first frame (static), so
            # its SF/OR are well-defined and constant: prepend that first-frame
            # value across the stationary window (xs_start → motion onset).
            if key == "sf":
                yv = vs["sf"].to_numpy(float)[vis]
                xv = x[vis]
                if len(yv):
                    xv = np.concatenate([[p["xs_start"]], xv])
                    yv = np.concatenate([[yv[0]], yv])
                return xv, yv
            yor = vs["orientation"].to_numpy(float)[vis] % 180.0  # real degrees
            xor = x[vis]
            if len(yor):
                xor = np.concatenate([[p["xs_start"]], xor])
                yor = np.concatenate([[yor[0]], yor])
            return _seam_break(xor, yor)
        if key in ("vf", "t"):
            src = _vel_source(key, tr)
            if src is None:
                return None, None
            vpt = np.asarray(tr.probe_t, float)
            xx = vpt - p["t_cmd"]
            m = (xx >= x_lo) & (xx <= x_hi)
            return xx[m], src[m]
        return None, None

    # --- FR (Gaussian σ20 ms): convolved per trial on the command-aligned grid;
    # continuous through the (realistic) gap. col 2 = the across-trial central
    # trace, ``fr_stat``: "median" (the typical trial; no band) or "mean" (with a
    # ± SEM band). On sparse cells the median collapses toward 0 because >half the
    # trials are spike-empty at any instant, so "mean" reads the PSTH-like bump. ---
    x_grid = np.arange(x_lo, x_hi + 1e-9, 0.001)
    fr_pink = fr_convolution(spike_times, sel["t_cmd"] + x_grid)
    fr_stack = [fr_convolution(spike_times, per[vt]["t_cmd"] + x_grid)
                for vt, _ in ov_trials]
    with np.errstate(all="ignore"):
        if fr_stack:
            _st = np.vstack(fr_stack)
            if fr_stat == "mean":
                fr_central = np.nanmean(_st, 0)
                _n = np.sum(np.isfinite(_st), 0)
                fr_sem = np.nanstd(_st, 0, ddof=1) / np.sqrt(np.maximum(_n, 1))
            else:
                fr_central = np.nanmedian(_st, 0)
                fr_sem = None
        else:
            fr_central, fr_sem = fr_pink, None

    t_ylab = "T tread.\n(cm/s)" if condition == "V" else "T stage\n(cm/s)"
    rows = [
        ("VF\n(cm/s)", "vf", "#1f4e79"),
        (t_ylab, "t", "#7a4a12"),
        ("TF\n(Hz)", "tf", "#d1701a"),
        ("SF\n(cpd)", "sf", "#6b8e23"),
        ("OR\n(deg)", "or", "#c0392b"),
        ("FR\n(Hz)", "fr", "#ff4da6"),
    ]
    # Layout: LEFT time-series block (cols 0,1: shared time x-axis, shared y per
    # row) + a RIGHT tuning column for THIS figure's OWN condition only, with a
    # condition-appropriate regressor set — V is visual {TF,SF,OR}; VT also carries
    # the stage-translation regressors {Speed, Acceleration}. (Speed≡0 / no stage
    # accel in V replay, so V omits them — "only fits relevant for their
    # condition".) The tuning panels are a standalone vertical stack, decoupled
    # from the 6 time-series rows since the regressor count differs by condition.
    # Speed/Accel/TF observed tuning come from the MATLAB cache (pc); SF/OR are
    # recomputed (rf_local, no cache).
    TUN_REG = {
        "V":  [("tf", "tf"), ("sf", "sf"), ("or", "orientation")],
        "VT": [("speed", "speed"), ("accel", "acceleration"),
               ("tf", "tf"), ("sf", "sf"), ("or", "orientation")],
    }[condition]
    TUN_COL = {"speed": "#7a4a12", "accel": "#8e44ad", "tf": "#d1701a",
               "sf": "#6b8e23", "or": "#c0392b"}
    n_rows = len(rows)
    n_tcl = len(tuning_clusters)
    fig = plt.figure(figsize=(13.5 + 3.2 * n_tcl, 13.5), constrained_layout=True)
    gs = fig.add_gridspec(n_rows, 3, width_ratios=[1.2, 1.2, 1.25 * n_tcl])
    axes = np.empty((n_rows, 2), dtype=object)
    for r in range(n_rows):
        for c in range(2):
            kw = {}
            if not (r == 0 and c == 0):
                kw["sharex"] = axes[0, 0]   # all time-series share the time axis
            if c == 1:
                kw["sharey"] = axes[r, 0]   # col1 shares y with col0 (per row)
            axes[r, c] = fig.add_subplot(gs[r, c], **kw)
    # Right tuning GRID: K regressor-rows × n_tcl cluster-columns. FR (y) is shared
    # across the cluster columns of a regressor row (sharey); cross-condition (V↔VT)
    # harmonisation is done by the caller. tun_axes[key] = [ax per cluster column].
    tun_gs = gs[:, 2].subgridspec(len(TUN_REG), n_tcl, hspace=0.6, wspace=0.3)
    tun_axes = {}
    for i, (key, _vc) in enumerate(TUN_REG):
        row_axes = []
        for j in range(n_tcl):
            shared = {"sharey": row_axes[0]} if j else {}
            row_axes.append(fig.add_subplot(tun_gs[i, j], **shared))
        tun_axes[key] = row_axes

    def _decorate(ax, key):
        """Shared chrome: period bands (real gap to scale), boundaries, refs."""
        ax.axvspan(sel["xs_start"], 0.0, color="#d6e4f0", alpha=0.7, lw=0, zorder=0)
        ax.axvspan(0.0, sel["gap"], facecolor="0.92", edgecolor="0.6",
                   hatch="///", lw=0.0, zorder=0)
        ax.axvspan(sel["gap"], sel["x_mend"], color="#fbe6d4", alpha=0.85, lw=0, zorder=0)
        # dashed: stationary start · command onset (0) · visual onset · motion end
        for xb in (sel["xs_start"], 0.0, sel["gap"], sel["x_mend"]):
            ax.axvline(xb, color="0.4", ls="--", lw=0.9, zorder=2)
        if key == "t" and condition == "V":
            ax.set_ylim(-1, 1)             # V: flat ≡ 0; VT autoscales to stage
        elif key == "sf":
            ax.axhline(sf_token_cpd, ls="--", color="0.35", lw=0.8)
            ax.set_ylim(*SF_YLIM_CPD)  # fixed to the global possible range
        elif key == "or":
            ax.axhline(or_token, ls="--", color="r", lw=0.8)  # token reference (deg)
            ax.set_ylim(0, 180)
            ax.set_yticks([0, 45, 90, 135, 180])
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(x_lo, x_hi)

    for r, (ylab, key, col) in enumerate(rows):
        ax0, ax1 = axes[r, 0], axes[r, 1]
        _decorate(ax0, key)
        _decorate(ax1, key)
        vis = tso >= -1e-9  # photodiode-on mask for SF/OR (selected trial)

        # ---- column 1: the selected trial (single-trial detail) ----
        if key in ("vf", "t"):
            src = _vel_source(key, trial)
            if src is None:                       # V: stage ≡ 0 (no translation)
                ax0.axhline(0.0, color=col, lw=1.4)
                ax0.text(0.985, 0.78, "≡ 0 (replay trial)", transform=ax0.transAxes,
                         ha="right", va="top", fontsize=6.5, color="0.4")
            else:
                vpt = np.asarray(trial.probe_t, float)
                xx = vpt - sel["t_cmd"]
                m = (xx >= x_lo) & (xx <= x_hi)
                ax0.plot(xx[m], src[m], color=col, lw=1.0)
        elif key == "fr":
            ax0.plot(x_grid, fr_pink, color=col, lw=1.0)
        elif key == "tf":
            ax0.plot(sel["gap"] + tso[vis], tf[vis], color=col, lw=1.0, marker="o", ms=2)
        elif key == "sf":
            # baseline: cloud at its first frame → constant SF across stationary.
            yv, xv = sf[vis], sel["gap"] + tso[vis]
            if len(yv):
                xv = np.concatenate([[sel["xs_start"]], xv])
                yv = np.concatenate([[yv[0]], yv])
            ax0.plot(xv, yv, color=col, lw=1.0, marker="o", ms=2)
            ax0.text(0.985, 0.04, f"token {sf_token_cpd:.3f}", transform=ax0.transAxes,
                     ha="right", va="bottom", fontsize=6, color="0.4")
        elif key == "or":
            yo, xv = or_real[vis], sel["gap"] + tso[vis]
            if len(yo):
                xv = np.concatenate([[sel["xs_start"]], xv])
                yo = np.concatenate([[yo[0]], yo])
            ax0.plot(*_seam_break(xv, yo), color=col, lw=1.0, marker="o", ms=2)
            ax0.text(0.985, 0.04, f"token {or_token:.0f}°", transform=ax0.transAxes,
                     ha="right", va="bottom", fontsize=6, color="0.4")
        ax0.set_ylabel(ylab, fontsize=10, rotation=0, ha="right", va="center")

        # ---- column 2: all this-condition trials overlaid (selected solid);
        # FR = median ----
        if key == "t" and condition == "V":
            ax1.axhline(0.0, color=col, lw=1.4)   # replay: no stage translation
        elif key == "fr":
            if fr_sem is not None:
                ax1.fill_between(x_grid, fr_central - fr_sem, fr_central + fr_sem,
                                 color=col, alpha=0.25, lw=0, zorder=1)
            ax1.plot(x_grid, fr_central, color=col, lw=1.4, zorder=2)
            _lbl = ("mean ± SEM" if fr_stat == "mean" else "median") + f" of {n_ov} trials"
            ax1.text(0.985, 0.92, _lbl, transform=ax1.transAxes,
                     ha="right", va="top", fontsize=6.5, color="0.4")
        else:
            for vt, tr in ov_trials:
                x, y = _xy(key, vt, tr)
                if x is None or not len(x):
                    continue
                s = vt == tid
                ax1.plot(x, y, color=col, alpha=1.0 if s else 0.2,
                         lw=1.1 if s else 0.7, zorder=4 if s else 1,
                         **(dict(marker="o", ms=2) if s and key not in ("vf", "t") else {}))

    # ---- right block: observed FR-vs-value tuning (median line + IQR band) +
    # the best-fit model, for THIS condition only, over its own regressor set.
    # Shared renderer (_render_tuning_panel): asym-Gaussian/linear families for
    # Speed/Accel/TF/SF (BIC-best), von Mises for OR. Speed/Accel/TF observed come
    # from the MATLAB cache (pc.{speed,accel,tf}_tuning); SF/OR are recomputed.
    bw = float(config.time_bin_width)
    for i, (key, vcol) in enumerate(TUN_REG):
        for j, (cid, cdf) in enumerate(tuning_clusters):
            ax = tun_axes[key][j]
            _render_tuning_panel(ax, cdf, vcol, key, condition, bw, TUN_COL[key],
                                 cluster_id=cid, pc=pc, tuning_n_reps=tuning_n_reps)
            ax.set_ylabel("FR (Hz)" if j == 0 else "", fontsize=9)
            if i == 0:                       # cluster id labels the column (top row)
                ax.set_title(f"cl {cid}", fontsize=11, fontweight="bold")

    # Period labels above both top axes.
    for ax in (axes[0, 0], axes[0, 1]):
        tr_ax = ax.get_xaxis_transform()
        ax.text(0.5 * sel["xs_start"], 1.04, "stationary\n(baseline)", transform=tr_ax,
                ha="center", va="bottom", fontsize=7, color="#2b5d8a")
        ax.text(0.5 * sel["gap"], 1.06, "gap", transform=tr_ax,
                ha="center", va="bottom", fontsize=6.5, color="0.35")
        motion_lab = ("motion\n(visual stimulus)" if condition == "V"
                      else "motion\n(stage + visual)")
        ax.text(0.5 * (sel["gap"] + sel["x_mend"]), 1.04, motion_lab,
                transform=tr_ax, ha="center", va="bottom", fontsize=7, color="#7a4a12")
    for c in (0, 1):
        axes[-1, c].set_xlabel(
            "time from command onset / solenoid (s)   "
            "(dashed: stat start · command · visual onset · motion end)",
            fontsize=6.6,
        )
    n_spk = int(sub["spike_count"].sum())
    dur_note = (f"stat {-sel['xs_start']:.1f}s · gap {sel['gap']:.2f}s · "
                f"motion {sel['x_mend'] - sel['gap']:.1f}s")
    probe_short = probe_id.split("_rec")[0]
    tcl_ids = ", ".join(str(c) for c, _ in tuning_clusters)
    reg_lab = "TF·SF·OR" if condition == "V" else "Speed·Accel·TF·SF·OR"
    left_note = ("OR real deg; SF/OR held at first frame through baseline"
                 if condition == "V" else
                 "VF = visual command, T = stage translation; OR real deg; "
                 "SF/OR held at first frame through baseline")
    head = (
        f"{probe_short} · {condition} trials · solenoid-aligned ({left_note})\n"
        f"LEFT (FR = cl {cluster_id}): trial {tid} ({n_spk} spk) + {n_ov}-trial overlay "
        f"(FR = {'mean ± SEM' if fr_stat == 'mean' else 'median'}) · "
        f"RIGHT: {condition} tuning [{reg_lab}] per cluster "
        f"[{tcl_ids}] — FR y shared across clusters & V↔VT; 20 bins"
        + ("; Speed/Accel/TF from MATLAB cache" if condition == "VT" else "")
        + f" + BIC-best & bootstrap-p · {dur_note}"
    )
    fig.suptitle(head, fontsize=11, fontweight="bold")
    if not save:
        return fig, tun_axes
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)
    return fig, tun_axes


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def plot_buildup(
    probe_id, cluster_id, df, config, step_results, tid,
    spike_times, t_motion_start, out_path, baseline_spk=None,
):
    """Stacked single-trial buildup: 1 column, one row per cumulative addition.

    Per row: raw 20 ms observed FR (gray line, NO smoothing), the MATLAB
    Gaussian FR (sigma=20 ms; pink), and the current cumulative model
    prediction (black). Window = full trial: stationary prelude THEN motion
    (~4 + ~4 s), motion onset marked. Rows/predictions are taken in df row
    order (the order the GLM saw → history stays self-consistent) on a synthetic
    axis (stationary is time_in_trial=0 in the binning, so given bin-width
    spacing before motion). The pink FR is computed in ABSOLUTE probe time over
    [t_motion_start - stat_dur, t_motion_start + motion_dur] and mapped onto the
    same synthetic axis (x = abs_t - t_lo).
    """
    bw = float(config.time_bin_width)
    df = df.reset_index(drop=True)
    cond = df["condition"].to_numpy(dtype=object)
    tit = df["time_in_trial"].to_numpy(dtype=float)
    trial_ids = df["trial_id"].to_numpy()

    pos = np.where(trial_ids == tid)[0]
    stat = pos[cond[pos] == "stationary"]
    mot = pos[cond[pos] != "stationary"]
    mot = mot[np.argsort(tit[mot])]  # temporal order within motion
    n_stat = stat.size
    onset_t = n_stat * bw
    motion_dur = float(tit[mot].max()) if mot.size else 0.0
    order = np.concatenate([stat, mot])
    t = np.concatenate([np.arange(n_stat) * bw, onset_t + tit[mot]])
    n_spk = int(df.iloc[order]["spike_count"].sum())
    cond_label = str(cond[mot[0]]) if mot.size else "?"

    # Pink: MATLAB Gaussian FR over the absolute trial window, mapped to x.
    t_lo = t_motion_start - onset_t
    t_hi = t_motion_start + motion_dur
    T = np.arange(t_lo, t_hi + 1e-9, 0.001)  # 1 ms timebase
    fr_pink = fr_convolution(spike_times, T)
    x_pink = T - t_lo

    n = len(step_results)
    fig, axes = plt.subplots(
        n, 1, figsize=(8.0, 1.9 * n), sharex=True, sharey=True, constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    baseline_cv = step_results[0]["cv"]
    for k, (ax, sr) in enumerate(zip(axes, step_results)):
        pred = sr["rate"][order]
        ax.plot(x_pink, fr_pink, color="#ff4da6", lw=1.0, zorder=1,
                label="observed FR (Gauss σ=20 ms)" if k == 0 else None)
        ax.plot(t, pred, color="black", lw=1.3, zorder=3,
                label="model prediction" if k == 0 else None)
        ax.axvline(onset_t, color="0.8", ls="--", lw=0.8, zorder=2)
        delta = sr["cv"] - (step_results[k - 1]["cv"] if k > 0 else sr["cv"])
        label = sr["label"].replace("\n", " ")
        if k == 0:
            txt = label
        else:
            txt = f"{label}   +{delta:.3f} bits/spike   (cum +{sr['cv'] - baseline_cv:.3f})"
        ax.set_title(txt, fontsize=8, loc="left")
        ax.set_ylabel("FR (Hz)", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        if k == 0:
            ax.legend(fontsize=7, loc="upper left", framealpha=0.9, ncol=2)

    axes[-1].set_xlabel(
        "time from stationary onset (s)   —   motion onset at dashed line", fontsize=9
    )
    base_note = f", {baseline_spk} baseline spk" if baseline_spk is not None else ""
    fig.suptitle(
        f"{probe_id}  cluster {cluster_id}  |  {cond_label} trial {tid} "
        f"({n_spk} spk{base_note}, full trial)  |  total gain "
        f"{step_results[-1]['cv'] - baseline_cv:+.3f} bps",
        fontsize=10,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)


# --------------------------------------------------------------------------- #
# Hardcastle 2017 Fig-2C-style forward-selection panel for ONE cluster + trial.
# Top row: cumulative-model FR reconstructions (pink observed Gaussian FR +
# black prediction). Bottom row: per-step candidate Δ cv-bps (the bits/spike
# increase each candidate would add at that step), the chosen term highlighted,
# a rejected term in red, and the 0.005 selection threshold as a dotted line.
# --------------------------------------------------------------------------- #
MAIN_POOL = ["Acceleration", "Speed", "TF", "SF", "OR", "ME_face"]
ABBR = {"Acceleration": "A", "Speed": "S", "TF": "TF", "SF": "SF",
        "OR": "OR", "ME_face": "ME"}


def _trial_window(df, config, tid, t_motion_start, spike_times):
    """Shared trial-window machinery: row order, synthetic time axis, motion
    onset, and the pink Gaussian-FR trace in mapped coordinates."""
    bw = float(config.time_bin_width)
    cond = df["condition"].to_numpy(dtype=object)
    tit = df["time_in_trial"].to_numpy(dtype=float)
    pos = np.where(df["trial_id"].to_numpy() == tid)[0]
    stat = pos[cond[pos] == "stationary"]
    mot = pos[cond[pos] != "stationary"]
    mot = mot[np.argsort(tit[mot])]
    n_stat = stat.size
    onset_t = n_stat * bw
    motion_dur = float(tit[mot].max()) if mot.size else 0.0
    order = np.concatenate([stat, mot])
    t = np.concatenate([np.arange(n_stat) * bw, onset_t + tit[mot]])
    t_lo = t_motion_start - onset_t
    T = np.arange(t_lo, t_motion_start + motion_dur + 1e-9, 0.001)
    fr_pink = fr_convolution(spike_times, T)
    x_pink = T - t_lo
    return order, t, onset_t, x_pink, fr_pink


def plot_forward_panel(
    probe_id, cluster_id, df, prep, config, tid, spike_times, t_motion_start,
    baseline_vars, accepted, rejected, out_path, backend="irls",
):
    order, t, onset_t, x_pink, fr_pink = _trial_window(
        df, config, tid, t_motion_start, spike_times)

    # Columns: baseline, then one per accepted term, then the rejected example.
    top_models = [("baseline", list(baseline_vars))]
    cur = list(baseline_vars)
    for v in accepted:
        cur = cur + [v]
        top_models.append((f"+ {v}", list(cur)))
    top_models.append((f"+ {rejected} (rejected)", list(cur) + [rejected]))
    n_col = len(top_models)
    base_cv = cv_bps_for(prep, baseline_vars, config, backend)

    # Bottom panels: one per transition (under columns 1..n_col-1). For every
    # still-available main effect, `_candidates` returns the PER-FOLD paired Δ
    # bits/spike over the current model (mean ± SD across the 10 condition-
    # stratified folds — the same paired Δ the signed-rank selection consumes);
    # the plotting loop below anchors each onto the cumulative axis. The chosen
    # term is highlighted with its one-sided Wilcoxon signed-rank p-value.
    def _candidates(cur, chosen, chosen_color):
        f_cur = cv_folds_for(prep, cur, config, backend)
        out = []
        for m in (mm for mm in MAIN_POOL if mm not in cur):
            mean, sd, p = _delta_stats(
                _paired_delta(cv_folds_for(prep, cur + [m], config, backend), f_cur))
            hl = m == chosen
            out.append((m, mean, sd, hl, chosen_color if hl else "black",
                        p if hl else None))
        return out

    # Bottom panels are anchored on a CUMULATIVE-over-baseline bits/spike axis:
    # each candidate sits at (running cumulative) + (its paired fold-mean Δ), so
    # the panels climb left→right and the chosen term's height carries into the
    # next column. `carried` = the cumulative of the current model (0 = baseline);
    # `is_accepted` gates the green accepted-step segment (the rejected column
    # shows where its term would land but takes no step, so cum does not advance).
    bottom = []  # (col_index, carried_cum, is_accepted, cands)
    cur = list(baseline_vars)
    cum = 0.0
    for i, v in enumerate(accepted):
        cands = _candidates(cur, v, "black")
        bottom.append((i + 1, cum, True, cands))
        chosen_mean = next((mn for m, mn, *_ in cands if m == v), float("nan"))
        cum = cum + (chosen_mean if np.isfinite(chosen_mean) else 0.0)
        cur = cur + [v]
    # Rejected example: the term that did NOT clear the signed-rank test — it is
    # plotted against the final accepted cumulative but does not advance it.
    bottom.append((n_col - 1, cum, False, _candidates(cur, rejected, "#d62728")))

    fig, axes = plt.subplots(
        2, n_col, figsize=(2.5 * n_col, 4.6),
        gridspec_kw={"height_ratios": [1.1, 1.0]}, constrained_layout=True,
    )

    # Top row: reconstructions (shared y so the prediction visibly catches up).
    rates = [insample_rate_for(prep, v, config, backend)[order] for _, v in top_models]
    cvs = [cv_bps_for(prep, v, config, backend) for _, v in top_models]
    top_ymax = max(float(fr_pink.max()), max(float(r.max()) for r in rates))
    x_end = float(t[-1])
    for k, (label, vars_) in enumerate(top_models):
        ax = axes[0, k]
        # Mask-derived periods: stationary prelude (shaded) vs motion.
        ax.axvspan(0, onset_t, color="#d6e4f0", alpha=0.8, lw=0, zorder=0)
        ax.plot(x_pink, fr_pink, color="#ff4da6", lw=0.9,
                label="Rec. FR (Gauss 20 ms)" if k == 0 else None)
        ax.plot(t, rates[k], color="black", lw=1.1,
                label="Pred. FR" if k == 0 else None)
        ax.axvline(onset_t, color="0.45", ls="--", lw=0.9)
        ax.set_xlim(0, x_end)
        ax.set_ylim(0, 1.05 * top_ymax)
        delta = cvs[k] - cvs[k - 1] if k > 0 else 0.0
        cum = cvs[k] - cvs[0]
        ttl = label if k == 0 else f"{label}\nΔ+{delta:.3f}  (cum +{cum:.3f})"
        ax.set_title(ttl, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        if k == 0:
            ax.set_ylabel("firing rate (Hz)", fontsize=8)
            ax.legend(fontsize=6, loc="lower left", framealpha=0.9)
            tr = ax.get_xaxis_transform()
            ax.text(onset_t * 0.5, 0.94, "stationary\n(baseline)", transform=tr,
                    ha="center", va="top", fontsize=6.5, color="#2b5d8a")
            ax.text(onset_t + (x_end - onset_t) * 0.5, 0.94, "motion", transform=tr,
                    ha="center", va="top", fontsize=6.5, color="#7a4a12")
        else:
            ax.set_yticklabels([])

    # Bottom row: each candidate plotted at the CUMULATIVE bits/spike its model
    # would reach (carried-in cumulative + its paired fold-mean Δ), error bar =
    # the paired per-fold Δ SD (the signed-rank input), shared y across columns.
    # A green vertical segment marks the accepted step (carried → chosen height).
    axes[1, 0].axis("off")  # no transition produces the baseline column
    pts = [(cc + mn, sd) for _, cc, _, cands in bottom for _, mn, sd, _, _, _ in cands
           if np.isfinite(mn)]
    ymax = max((y + sd for y, sd in pts), default=0.1)
    ymin = min([y - sd for y, sd in pts] + [0.0], default=-0.02)
    for col, carried, is_accepted, cands in bottom:
        ax = axes[1, col]
        # faint anchor at the cumulative carried in from the current model
        ax.axhline(carried, color="0.8", lw=0.7, ls=":")
        for xi, (m, mn, sd, hl, color, p) in enumerate(cands):
            if not np.isfinite(mn):
                continue
            y = carried + mn
            # green segment = the accepted increment (Δ) this step adds to the model
            if hl and is_accepted:
                ax.plot([xi, xi], [carried, y], color="#2ca02c", lw=2.4,
                        solid_capstyle="round", zorder=2.5)
            ax.errorbar([xi], [y], yerr=[sd], fmt="o",
                        ms=7 if hl else 4, color=color if hl else "0.6",
                        ecolor=color if hl else "0.7",
                        elinewidth=1.4 if hl else 0.9, capsize=3,
                        markeredgecolor="black" if hl else "none",
                        zorder=3 if hl else 2)
            if hl and p is not None and np.isfinite(p):
                ax.annotate(f"signed-rank\np={p:.3f}", (xi, y + sd),
                            textcoords="offset points", xytext=(0, 4),
                            ha="center", va="bottom", fontsize=6, color=color)
        ax.set_xticks(range(len(cands)))
        ax.set_xticklabels([ABBR.get(m, m) for m, *_ in cands], fontsize=7)
        ax.set_ylim(min(-0.02, 1.2 * ymin), 1.3 * ymax)
        ax.spines[["top", "right"]].set_visible(False)
        if col == 1:
            ax.set_ylabel("cumulative bits/spike (over baseline)", fontsize=8)
            ax.legend([Line2D([0], [0], color="#2ca02c", lw=2.4)],
                      ["accepted Δ (step)"], fontsize=6, loc="upper left",
                      framealpha=0.9)

    fig.suptitle(
        f"{probe_id}  cluster {cluster_id}  trial {tid}  —  forward selection "
        f"(condition-stratified 10-fold, one-sided Wilcoxon signed-rank, "
        f"α={config.selection_alpha}; History in baseline)",
        fontsize=10, fontweight="bold",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)


# --------------------------------------------------------------------------- #
# Fig 2 — population forward-selection summary + acid-stack (evolution of
# forward_selection_summary + acid_stack_sorted).
#   (1,1) histogram of #selected predictors per cluster (single colour).
#   (1,2) main effects: VIOLIN of per-cluster LOO-unique cv-bps (all points
#         overlaid), over the clusters that selected each term; linear axis.
#   (1,3) same for interaction terms (LOO-unique computed on the fly to match
#         variance_partition's full-additive / condition-stratified-5-fold defn).
#   (2,1) per-cluster stacked unique Δ cv-bps, all main effects except History &
#         Onset, sorted by stack total (symlog full / linear+zoom split).
#   (3,1) Onset LOO-unique Δ cv-bps per cluster, SAME sort as (2,1), single gray,
#         signed; mirrors the (2,1) total + zoom layout.
# --------------------------------------------------------------------------- #
# Single per-predictor colour map shared by (1,2), (1,3), and (2,1). History and
# Onset are not in those panels (baseline/nuisance terms); Onset gets its own
# row (3,1), History is NaN in the histbase run so it is not shown.
PREDICTOR_COLORS = {"Speed": "tab:green", "Acceleration": "tab:cyan", "TF": "tab:orange",
                    "SF": "tab:olive", "OR": "tab:red", "ME_face": "tab:purple"}
DISPLAY = {"Speed": "Speed", "Acceleration": "Accel", "TF": "TF", "SF": "SF",
           "OR": "OR", "ME_face": "ME"}
# (canonical label, variance_partition column) — main effects only, no History/Onset.
ME_UNIQUE = [("Speed", "unique_Speed"), ("Acceleration", "unique_Accel"),
             ("TF", "unique_TF"), ("SF", "unique_SF"), ("OR", "unique_OR"),
             ("ME_face", "unique_ME")]
# Acid-stack vars (display label, column, colour) — History AND Onset excluded.
ACID_VARS = [("Speed", "unique_Speed", "tab:green"), ("Accel", "unique_Accel", "tab:cyan"),
             ("TF", "unique_TF", "tab:orange"), ("SF", "unique_SF", "tab:olive"),
             ("OR", "unique_OR", "tab:red"), ("ME", "unique_ME", "tab:purple")]
# History stays in the full additive model when measuring uniques (it's a real
# nuisance regressor); it's only excluded from the DISPLAY.
FULL_ADD = ["Speed", "TF", "SF", "OR", "ME_face", "Acceleration", "History"]


def _interaction_colors(name):
    """(facecolor, hatch-colour) for an interaction bar = its two constituents."""
    parts = name.split("_x_")
    cols = [PREDICTOR_COLORS.get(p, "0.5") for p in parts]
    return cols[0], (cols[1] if len(cols) > 1 else cols[0])


def _sel_terms(s):
    s = str(s)
    return [] if s in ("", "nan") else s.split("+")


def _cv_with_folds(prep, vars_, folds, config, backend="irls"):
    X, _ = _design(prep, vars_, config)
    if X.shape[1] == 0 or X.shape[1] >= prep["y"].size:
        return float("nan")
    from rc2_glm.cross_validation import cross_validate_glm as _cv
    return float(_cv(X, prep["y"], prep["offset"], folds,
                     lambda_ridge=config.lambda_ridge, backend=backend).cv_bits_per_spike)


def compute_interaction_uniques(config, probes, mc, backend="irls"):
    """Per-cluster LOO-unique cv-bps per interaction term, over the clusters that
    selected it. unique_I = cv(full_additive ∪ I) − cv(full_additive),
    condition-stratified 5-fold (matches variance_partition_accel.partition).
    Returns ``{term: [per-cluster unique, ...]}`` (RAW Δ, negatives kept — the
    violin panel shows the full distribution)."""
    from rc2_glm.cross_validation import make_trial_folds
    per = {}
    for probe_id in probes:
        sub = mc[mc["probe_id"] == probe_id]
        need = {int(r.cluster_id): [v for v in _sel_terms(r.time_selected_vars) if "_x_" in v]
                for r in sub.itertuples() if any("_x_" in v for v in _sel_terms(r.time_selected_vars))}
        if not need:
            continue
        probe = load_probe_data(FORMATTED_DIR / f"{probe_id}.mat", config=config,
                                stimulus_lookup=_lookup(), cluster_set="selected")
        rf = load_rf_sf_or(config.rf_sf_or_parquet_dir, probe.probe_id,
                           min_concentration=getattr(config, "rf_min_concentration", 0.0))
        by = {c.cluster_id: c for c in probe.clusters}
        for cid, inters in need.items():
            if cid not in by:
                continue
            df = bin_cluster(probe, by[cid], rf_lookup=rf)
            if df.empty:
                continue
            prep = prepare_cluster_design(df, config)
            folds = make_trial_folds(
                df["trial_id"].to_numpy(np.int64), config.n_folds, config.cv_seed,
                condition_labels_per_bin=df["condition"].to_numpy(object),
                strategy="condition-stratified")
            cv_full = _cv_with_folds(prep, FULL_ADD, folds, config, backend)
            for I in inters:
                u = _cv_with_folds(prep, FULL_ADD + [I], folds, config, backend) - cv_full
                per.setdefault(I, []).append(u)
            log.info("interaction uniques: %s cluster %d done", probe_id, cid)
    return per


def _draw_acid_stack(ax, df_sub, present, *, ylim=None, order=None, signed=False):
    """Per-cluster stacked unique Δ cv-bps for a (sub)set of clusters.
    ``present`` = [(label, column, colour)]. By default negatives are clipped to 0
    and the columns are sorted ascending by stack total; pass ``signed=True`` to
    keep negative bars (used for the single-series Onset row, so every datapoint
    shows) and pass an explicit ``order`` (positional indices) to reuse another
    panel's cluster ordering so columns align. Returns ``(tot, order)``."""
    def _col(col):
        v = pd.to_numeric(df_sub[col], errors="coerce").fillna(0).to_numpy()
        return v if signed else np.clip(v, 0, None)
    heights = {lab: _col(col) for lab, col, _ in present}
    tot = np.sum(list(heights.values()), axis=0) if heights else np.zeros(len(df_sub))
    if order is None:
        order = np.argsort(tot)
    x = np.arange(len(df_sub)); bottom = np.zeros(len(df_sub))
    for lab, col, c in present:
        h = heights[lab][order]
        ax.bar(x, h, bottom=bottom, width=1.0, color=c, linewidth=0, label=lab)
        bottom = bottom + h
    ax.set_xlim(-0.5, len(df_sub) - 0.5)
    if ylim is not None:
        ax.set_ylim(*ylim)
    return tot, order


def plot_fig2_summary(mc, vp, inter_per, out_path, config, *, zoom_split=False):
    n = len(mc)
    fig = plt.figure(figsize=(16, 13) if zoom_split else (15, 13), constrained_layout=True)
    gs = fig.add_gridspec(3, 6, height_ratios=[1.0, 1.0, 1.0])
    ax11 = fig.add_subplot(gs[0, 0:2])
    ax12 = fig.add_subplot(gs[0, 2:4])
    ax13 = fig.add_subplot(gs[0, 4:6])
    if zoom_split:
        ax2 = fig.add_subplot(gs[1, 0:4]); ax2z = fig.add_subplot(gs[1, 4:6])
        ax3 = fig.add_subplot(gs[2, 0:4]); ax3z = fig.add_subplot(gs[2, 4:6])
    else:
        ax2 = fig.add_subplot(gs[1, :]); ax2z = None
        ax3 = fig.add_subplot(gs[2, :]); ax3z = None

    # (1,1) model-size histogram, single colour. History excluded from the count
    # (nuisance term, not shown anywhere); interactions counted.
    sizes = mc["time_selected_vars"].map(
        lambda s: len([t for t in _sel_terms(s) if t != "History"]))
    vc = sizes.value_counts().sort_index()
    ax11.bar(vc.index, vc.values, color="#4c72b0", edgecolor="0.3", width=0.8)
    for x, c in zip(vc.index, vc.values):
        ax11.text(x, c, str(int(c)), ha="center", va="bottom", fontsize=8)
    ax11.set_xlabel("# selected predictors (excl. History)"); ax11.set_ylabel("# clusters")
    ax11.set_title(f"Model complexity (n={n})", fontsize=10)
    ax11.set_xticks(vc.index)

    # (1,2) main effects: per-cluster LOO-unique cv-bps distribution (violin + all
    # points), over the clusters that SELECTED the term (same basis as (1,3)).
    sel_map = {(r.probe_id, int(r.cluster_id)): set(_sel_terms(r.time_selected_vars))
               for r in mc.itertuples()}
    me_labels, me_data = [], []
    for lab, col in ME_UNIQUE:
        if col not in vp.columns:
            continue
        mask = np.array([lab in sel_map.get((r.probe_id, int(r.cluster_id)), set())
                         for r in vp.itertuples()])
        vals = pd.to_numeric(vp[col], errors="coerce").to_numpy()
        me_labels.append(lab); me_data.append(vals[mask])
    _violin_panel(ax12, [DISPLAY[m] for m in me_labels], me_data, "main effects",
                  [PREDICTOR_COLORS[m] for m in me_labels])

    # (1,3) interactions: per-cluster LOO-unique cv-bps distribution (sorted by
    # Σ of positive uniques, so the strongest interaction is leftmost).
    it_labels = sorted(inter_per, key=lambda k: -float(np.nansum(np.clip(inter_per[k], 0, None))))
    _violin_panel(ax13, [_INT_ABBR(k) for k in it_labels],
                  [inter_per[k] for k in it_labels], "interaction terms",
                  [_interaction_colors(k)[0] for k in it_labels])

    # (2,1) acid stack, main effects (History & Onset excluded). Capture the
    # cluster ORDER so the Onset row (3,1) below aligns column-for-column.
    present = [(lab, col, c) for lab, col, c in ACID_VARS if col in vp.columns]
    tot, order_all = _draw_acid_stack(ax2, vp, present)
    ax2.set_xlabel("cluster (sorted; History & Onset excluded)")
    ax2.legend(ncol=len(present), fontsize=8, frameon=False, loc="upper left")

    # (3,1) Onset row: SAME cluster sort as row 2, single gray, signed (negatives
    # kept — Onset is a small ±0.025 nuisance term and we want every datapoint).
    onset_present = [("Onset", "unique_Onset", "0.45")]
    has_onset = "unique_Onset" in vp.columns

    if zoom_split:
        ax2.set_ylabel("stacked unique Δ cv-bps (linear)")
        ax2.set_title(f"Per-cluster unique contributions — all clusters (n={len(vp)})",
                      fontsize=10, fontweight="bold")
        mask = tot < 0.2
        vp_zoom = vp[mask]
        _, order_zoom = _draw_acid_stack(ax2z, vp_zoom, present, ylim=(0, 0.2))
        ax2z.set_xlabel(f"cluster (cumulative < 0.2, n={int(mask.sum())})")
        ax2z.set_ylabel("stacked unique Δ cv-bps")
        ax2z.set_title("magnified: cumulative < 0.2 bits/spike", fontsize=10)
        if has_onset:
            _draw_acid_stack(ax3, vp, onset_present, order=order_all, signed=True)
            _draw_acid_stack(ax3z, vp_zoom, onset_present, order=order_zoom, signed=True)
            ax3z.set_xlabel(f"cluster (same {int(mask.sum())} as above)")
            ax3z.set_ylabel("Onset unique Δ cv-bps (signed)")
            ax3z.set_title("magnified: same low-cumulative clusters", fontsize=10)
    else:
        ax2.set_yscale("symlog", linthresh=0.01)
        ax2.set_ylabel("stacked unique Δ cv-bps (symlog)")
        ax2.set_title(f"Per-cluster unique contributions — except History & Onset (n={len(vp)})",
                      fontsize=10, fontweight="bold")
        if has_onset:
            _draw_acid_stack(ax3, vp, onset_present, order=order_all, signed=True)

    if has_onset:
        ax3.axhline(0, color="0.6", lw=0.6)
        ax3.set_xlabel("cluster (same sort as row above; Onset only)")
        ax3.set_ylabel("Onset unique Δ cv-bps (signed)")
        ax3.set_title(f"Per-cluster Onset contribution — baseline nuisance (n={len(vp)})",
                      fontsize=10, fontweight="bold")
    else:
        ax3.text(0.5, 0.5, "unique_Onset not in variance_partition.csv",
                 ha="center", va="center", transform=ax3.transAxes, fontsize=9)
        if ax3z is not None:
            ax3z.axis("off")

    fig.suptitle("FENS fig 2 — forward-selection summary + acid stack (goggles, "
                 f"{RUN_ALL_DIR.name})", fontsize=11, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)


# --------------------------------------------------------------------------- #
# FORWARD-Δ variant of fig 2: the SAME acid layout, but every value is the
# per-step ACCEPTED forward-selection gain (glm_selection_history_full.csv,
# ``added_this_round``) instead of the LOO variance partition. No LOO anywhere —
# the Onset/LOO reference row is dropped entirely.
# --------------------------------------------------------------------------- #
_FWD_MAIN_COL = {"Speed": "unique_Speed", "Acceleration": "unique_Accel",
                 "TF": "unique_TF", "SF": "unique_SF", "OR": "unique_OR",
                 "ME_face": "unique_ME"}


def _load_forward_deltas(run_dir, probes):
    """Per-cluster ACCEPTED forward-selection Δ cv-bps from each probe's
    glm_selection_history_full.csv (the sequential gain a term added when it was
    admitted: ``added_this_round`` True). Returns ``(fvp, inter_fwd)``:
      ``fvp`` — one row per cluster, ``unique_<main>`` columns (NaN if the term was
                not selected) so the existing acid/violin helpers consume it as-is;
      ``inter_fwd`` — ``{interaction_candidate: array of accepted Δ}`` over the
                clusters that admitted it."""
    main_rows, inter = [], {}
    for probe in probes:
        csv = run_dir / "_runs" / probe / "glm_selection_history_full.csv"
        if not csv.exists():
            continue
        df = pd.read_csv(csv)
        adm = df[df["added_this_round"].astype(str).str.lower().isin(("true", "1"))]
        for cid in sorted(df["cluster_id"].unique()):
            csub = adm[adm["cluster_id"] == cid]
            row = {"probe_id": probe, "cluster_id": int(cid)}
            for cand, col in _FWD_MAIN_COL.items():
                m = csub[csub["candidate"] == cand]
                row[col] = float(m["delta_bps"].iloc[0]) if len(m) else np.nan
            main_rows.append(row)
        for r in adm.itertuples():
            if "_x_" in str(r.candidate):
                inter.setdefault(str(r.candidate), []).append(float(r.delta_bps))
    fvp = pd.DataFrame(main_rows)
    inter_fwd = {k: np.asarray(v, float) for k, v in inter.items()}
    return fvp, inter_fwd


def plot_fig2_forward_deltas(mc, fvp, inter_fwd, out_path, *, log_heat=False, unit_heat=False):
    """Fig 2, FORWARD-Δ version (NO LOO): model-complexity histogram + accepted
    forward-Δ violins for main effects / interactions (row 1), and the per-cluster
    stacked forward-Δ acid bars + low-cumulative zoom (row 2). Row 3 is the same
    per-regressor data as a heatmap: ``log_heat=False`` → linear colour + a rescaled
    low-cumulative zoom panel; ``log_heat=True`` → a single full-width LOG-colour
    heatmap (big + small clusters legible at once, so no zoom)."""
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    gs = fig.add_gridspec(3, 6, height_ratios=[1.0, 1.0, 0.7])
    ax11 = fig.add_subplot(gs[0, 0:2]); ax12 = fig.add_subplot(gs[0, 2:4])
    ax13 = fig.add_subplot(gs[0, 4:6])
    ax2 = fig.add_subplot(gs[1, 0:4]); ax2z = fig.add_subplot(gs[1, 4:6])
    if log_heat or unit_heat:
        ax3 = fig.add_subplot(gs[2, :]); ax3z = None      # single full-width heatmap
    else:
        ax3 = fig.add_subplot(gs[2, 0:4]); ax3z = fig.add_subplot(gs[2, 4:6])

    # (1,1) model-complexity histogram (from the real model_comparison).
    sizes = mc["time_selected_vars"].map(
        lambda s: len([t for t in _sel_terms(s) if t != "History"]))
    vc = sizes.value_counts().sort_index()
    ax11.bar(vc.index, vc.values, color="#4c72b0", edgecolor="0.3", width=0.8)
    for x, c in zip(vc.index, vc.values):
        ax11.text(x, c, str(int(c)), ha="center", va="bottom", fontsize=8)
    ax11.set_xlabel("# selected predictors (excl. History)"); ax11.set_ylabel("# clusters")
    ax11.set_title(f"Model complexity (n={len(mc)})", fontsize=10); ax11.set_xticks(vc.index)

    # (1,2) main-effect accepted forward Δ, over the clusters that selected each term.
    me_labels, me_data = [], []
    for lab, col in ME_UNIQUE:
        if col not in fvp.columns:
            continue
        vals = pd.to_numeric(fvp[col], errors="coerce").to_numpy()
        me_labels.append(lab); me_data.append(vals[np.isfinite(vals)])
    _violin_panel(ax12, [DISPLAY[m] for m in me_labels], me_data, "main effects",
                  [PREDICTOR_COLORS[m] for m in me_labels])
    ax12.set_ylabel("accepted forward Δ bits/spike")
    ax12.set_title("Per-cluster forward Δ — main effects", fontsize=10)

    # (1,3) interaction accepted forward Δ, sorted by Σ of positive gains.
    it_labels = sorted(inter_fwd, key=lambda k: -float(np.nansum(np.clip(inter_fwd[k], 0, None))))
    _violin_panel(ax13, [_INT_ABBR(k) for k in it_labels],
                  [inter_fwd[k] for k in it_labels], "interaction terms",
                  [_interaction_colors(k)[0] for k in it_labels])
    ax13.set_ylabel("accepted forward Δ bits/spike")
    ax13.set_title("Per-cluster forward Δ — interaction terms", fontsize=10)

    # (2,*) per-cluster stacked accepted forward Δ + low-cumulative zoom.
    present = [(lab, col, c) for lab, col, c in ACID_VARS if col in fvp.columns]
    tot, order_all = _draw_acid_stack(ax2, fvp, present)
    ax2.set_xlabel("cluster (sorted; History & Onset are baseline, not forward-selected)")
    ax2.set_ylabel("stacked forward Δ bits/spike (linear)")
    ax2.set_title(f"Per-cluster forward-Δ contributions — all clusters (n={len(fvp)})",
                  fontsize=10, fontweight="bold")
    ax2.legend(ncol=len(present), fontsize=8, frameon=False, loc="upper left")
    mask = tot < 0.2
    _, order_zoom = _draw_acid_stack(ax2z, fvp[mask], present, ylim=(0, 0.2))
    ax2z.set_xlabel(f"cluster (cumulative < 0.2, n={int(mask.sum())})")
    ax2z.set_ylabel("stacked forward Δ bits/spike")
    ax2z.set_title("magnified: cumulative < 0.2 bits/spike", fontsize=10)

    # (3,*) SAME per-regressor forward Δ as row 2, as a HEATMAP: y = main-effect
    # regressors (no interactions), x = clusters in the SAME sort as row 2, colour =
    # forward Δ bits/spike (white = the term was not selected for that cluster). The
    # zoom reuses row-2's low-cumulative subset with the colour scale RESCALED to it.
    reg_labels = [lab for lab, _col, _c in present]
    M = np.vstack([pd.to_numeric(fvp[col], errors="coerce").to_numpy()
                   for _lab, col, _c in present])              # (n_regressors, n_clusters)
    cmap = plt.get_cmap("magma").copy(); cmap.set_bad("white")

    def _heat(ax, mat, xlabel, title, *, vmax=None, norm=None):
        kw = {"norm": norm} if norm is not None else {"vmin": 0.0, "vmax": vmax}
        im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", cmap=cmap,
                       interpolation="nearest", **kw)
        ax.set_yticks(range(len(reg_labels))); ax.set_yticklabels(reg_labels, fontsize=8)
        ax.set_xlabel(xlabel); ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.01, label="forward Δ bits/spike")

    finite = M[np.isfinite(M)]
    if unit_heat:
        # Per-cluster min-max to [0, 1]: divide each cluster's forward Δ by its own
        # max, so the dominant regressor = 1 and the profile is comparable across
        # clusters regardless of absolute magnitude. Not-selected = 0 (low end of
        # the colour map, same as a near-zero gain). Single panel, no zoom.
        M0 = np.nan_to_num(M, nan=0.0)
        cmax = M0.max(axis=0)
        Mn = np.zeros_like(M0)
        nz = cmax > 0
        Mn[:, nz] = M0[:, nz] / cmax[nz]
        # Order clusters by their TOP contributor — Speed, Accel, TF, SF, OR, ME (the
        # row order); no-main-effect clusters last; within a group the strongest top
        # contributor first. White lines mark the group boundaries.
        top_reg = np.where(nz, M0.argmax(axis=0), len(reg_labels))
        heat_order = np.lexsort((-cmax, top_reg))
        im = ax3.imshow(Mn[:, heat_order], aspect="auto", cmap="magma",
                        vmin=0.0, vmax=1.0, interpolation="nearest")
        for b in np.where(np.diff(top_reg[heat_order]) != 0)[0]:
            ax3.axvline(b + 0.5, color="white", lw=0.8, alpha=0.7)
        ax3.set_yticks(range(len(reg_labels))); ax3.set_yticklabels(reg_labels, fontsize=8)
        ax3.set_xlabel("cluster (grouped by top contributor: " + "→".join(reg_labels) + ")")
        ax3.set_title("Per-cluster forward Δ — heatmap (per-cluster 0–1, grouped by top predictor)",
                      fontsize=10)
        fig.colorbar(im, ax=ax3, fraction=0.04, pad=0.01,
                     label="forward Δ (per-cluster max = 1)")
    elif log_heat:
        # LOG colour over the full range → big + small clusters legible in ONE
        # panel, so no separate zoom is needed.
        from matplotlib.colors import LogNorm
        pos = M[np.isfinite(M) & (M > 0)]
        vmin = max(float(pos.min()), 1e-3) if pos.size else 1e-3
        vmax = float(pos.max()) if pos.size else 1.0
        _heat(ax3, M[:, order_all], "cluster (same sort as row 2)",
              "Per-cluster forward Δ — heatmap (main effects, LOG colour, all clusters)",
              norm=LogNorm(vmin=vmin, vmax=vmax))
    else:
        vmax_all = float(np.percentile(finite, 99)) if finite.size else 1.0
        _heat(ax3, M[:, order_all], "cluster (same sort as row 2)",
              "Per-cluster forward Δ — heatmap (main effects)", vmax=vmax_all)
        zoom_idx = np.where(mask)[0]
        Mz = M[:, zoom_idx][:, order_zoom]
        finite_z = Mz[np.isfinite(Mz)]
        vmax_z = float(np.percentile(finite_z, 99)) if finite_z.size else vmax_all
        _heat(ax3z, Mz, f"cluster (cumulative < 0.2, n={int(mask.sum())})",
              "magnified: same clusters (colour scale rescaled)", vmax=vmax_z)

    fig.suptitle("FENS fig 2 — forward-selection summary, FORWARD-Δ version "
                 "(stacked accepted sequential gains, not LOO-unique)",
                 fontsize=12, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", out_path)


def _INT_ABBR(name):
    return name.replace("ME_face", "ME").replace("Acceleration", "A").replace("_x_", "×")


def _violin_panel(ax, labels, data_lists, title, colors):
    """Per-term distribution of per-cluster LOO-unique cv-bps as a violin with
    every datapoint overlaid (jittered strip + median bar). Linear axis (no log),
    negatives kept. ``data_lists[i]`` = the per-cluster uniques for ``labels[i]``;
    ``n=`` is the count of finite points."""
    xs = np.arange(len(labels))
    for xi, (data, c) in enumerate(zip(data_lists, colors)):
        data = np.asarray([d for d in np.asarray(data, float) if np.isfinite(d)])
        if data.size == 0:
            continue
        if data.size >= 2 and np.ptp(data) > 0:
            parts = ax.violinplot([data], positions=[xi], widths=0.8,
                                  showmeans=False, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(c); b.set_edgecolor(c); b.set_alpha(0.30)
        jit = (np.random.RandomState(xi).rand(data.size) - 0.5) * 0.28
        ax.scatter(xi + jit, data, s=10, color=c, edgecolor="0.2",
                   linewidth=0.3, alpha=0.85, zorder=3)
        ax.hlines(np.median(data), xi - 0.22, xi + 0.22, color="0.1", lw=1.4, zorder=4)
        ax.annotate(f"n={data.size}", (xi, data.max()), textcoords="offset points",
                    xytext=(0, 3), ha="center", va="bottom", fontsize=7)
    ax.axhline(0, color="0.85", lw=0.6)
    ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("LOO-unique cv-bps (per cluster)")
    ax.set_title(f"Per-cluster cv-bps — {title}", fontsize=10)


# --------------------------------------------------------------------------- #
# Oracle: buildup cv-bps for the stored Selected set == CSV value.
# --------------------------------------------------------------------------- #
def oracle_check(prep, config, probe_id, cluster_id):
    csv = RUN_ALL_DIR / "_runs" / probe_id / "glm_model_comparison.csv"
    if not csv.exists():
        log.warning("oracle: no %s — skipping", csv)
        return
    mc = pd.read_csv(csv)
    row = mc[mc["cluster_id"] == cluster_id]
    if row.empty:
        log.warning("oracle: cluster %d not in %s — skipping", cluster_id, csv)
        return
    sel_str = str(row["time_selected_vars"].iloc[0])
    sel_vars = [] if sel_str in ("", "nan") else sel_str.split("+")
    stored = float(row["time_Selected_cv_bps"].iloc[0])
    mine = cv_bps_for(prep, sel_vars, config)
    # 5e-3 tolerance: the stored run used multi-threaded BLAS; even thread-pinned
    # our recompute lands within ~few×1e-3 of it (IRLS float-summation order).
    ok = np.isfinite(mine) and abs(mine - stored) < 5e-3
    log.info(
        "ORACLE cluster %d: Selected=[%s] stored=%.5f mine=%.5f Δ=%.2e -> %s",
        cluster_id, sel_str, stored, mine, mine - stored, "PASS" if ok else "FAIL",
    )
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--probe", default="CAA-1124370_rec1_rec2_rec3")
    ap.add_argument("--clusters", type=int, nargs="+", default=[20, 29, 375])
    ap.add_argument("--n-trials", type=int, default=6,
                    help="VT trials to render per cluster (separate figures).")
    ap.add_argument("--baseline-min", type=int, default=1)
    ap.add_argument("--baseline-max", type=int, default=6,
                    help="Pick VT trials with baseline-min..baseline-max stationary "
                         "spikes ('a few, not too active').")
    ap.add_argument("--fig2c", action="store_true",
                    help="Hardcastle Fig-2C forward-selection panel for ONE trial "
                         "(needs --trial).")
    ap.add_argument("--trial", type=int, default=None,
                    help="Trial id for --fig2c.")
    ap.add_argument("--accepted", default="Acceleration,SF,Speed",
                    help="--fig2c: accepted terms (columns), in order.")
    ap.add_argument("--rejected", default="OR",
                    help="--fig2c: the rejected example term (last column).")
    ap.add_argument("--fig2-summary", action="store_true",
                    help="Fig 2: population forward-selection summary + acid stack "
                         "(all clusters, both probes). Ignores --clusters/--trial.")
    ap.add_argument("--fig2-forward-deltas", action="store_true",
                    help="Fig 2, FORWARD-Δ version: same acid layout driven by the "
                         "accepted forward-selection gains (glm_selection_history_full.csv) "
                         "instead of the LOO variance partition. No LOO; no Onset row.")
    ap.add_argument("--trial-structure", action="store_true",
                    help="Poster column 1: trial-structure schematic (VF, T, "
                         "TF, SF, OR, pink FR) for one trial per --clusters cluster. "
                         "Auto-picks a trial of --condition unless --trial is given.")
    ap.add_argument("--fr-stat", choices=("median", "mean"), default="median",
                    help="--trial-structure: across-trial FR statistic in col 2. "
                         "'median' (default, the typical trial; no band) or 'mean' "
                         "(with a ± SEM band). On sparse cells the median collapses "
                         "toward 0; 'mean' reads the PSTH-like bump.")
    ap.add_argument("--condition", choices=("V", "VT"), default="V",
                    help="--trial-structure: trial condition. V (ReplayOnly, "
                         "default) → VF = visual command, T ≡ 0. VT (StageOnly) → "
                         "VF = visual command, T = stage TRANSLATION speed.")
    ap.add_argument("--match-clouds", action="store_true",
                    help="--trial-structure: per cluster, render BOTH a V and a VT "
                         "figure from the SAME motion cloud (the cloud where the "
                         "cluster is most active across V+VT). Ignores --condition/"
                         "--trial.")
    ap.add_argument("--top-fr", type=int, default=None,
                    help="Tuning-only browsing grids (TF/SF/OR × V|VT, no "
                         "time-series) for the top-N most active clusters (mean "
                         "firing rate) → FENS_figures_poster/tuning_browse/. For "
                         "picking cells.")
    ap.add_argument("--all-clusters", action="store_true",
                    help="With the tuning-browse path: render EVERY selected "
                         "cluster (not just top-N).")
    ap.add_argument("--probes", nargs="+", default=None,
                    help="Probe list for the tuning-browse path (default: --probe). "
                         "Goggles cohort = CAA-1124370_rec1_rec2_rec3 "
                         "CAA-1124371_rec1_rec2_rec3.")
    ap.add_argument("--display", choices=("median_iqr", "mean_sd", "mean_sem"),
                    default="median_iqr",
                    help="Tuning error bars: per-bin median+IQR (default), mean±SD, "
                         "or mean±SEM. Each writes to its own tuning_browse[_*]/ folder.")
    ap.add_argument("--select", choices=("bic", "rsq_mean"), default="bic",
                    help="Model-selection criterion for TF/SF: lowest BIC (default) "
                         "or highest R²-on-the-mean-curve (the MATLAB "
                         "ModelSelectionTuning rule). rsq_mean → a *_selR2mean/ folder.")
    ap.add_argument("--reuse-fits", default=None,
                    help="Path to an existing tuning_fits.csv; reuse its cached "
                         "model fits + bootstrap-p instead of recomputing (e.g. to "
                         "re-render the SAME fits with a different --display).")
    ap.add_argument("--tuning-n-reps", type=int, default=1000,
                    help="Bootstrap reps for the tuning-significance p (default 1000; "
                         "drop to ~200 for fast browsing).")
    ap.add_argument("--tuning-significance-summary", action="store_true",
                    help="Cohort tuning-significance summaries over an existing "
                         "tuning_fits.csv (NO re-fit): a V-only 2×3 tuned/not-tuned "
                         "strip (all + good-fit halves) and the V→VT transition plot "
                         "→ FENS_figures_poster/. CSV from --reuse-fits, else the "
                         "tuning_browse_mean_sem_selR2mean/ one.")
    args = ap.parse_args()

    config = make_config_histbase_all()
    backend = "irls"

    if args.tuning_significance_summary:
        csv_path = (Path(args.reuse_fits) if args.reuse_fits else
                    OUT_DIR / "tuning_browse_mean_sem_selR2mean" / "tuning_fits.csv")
        if not csv_path.exists():
            ap.error(f"tuning_fits.csv not found: {csv_path} "
                     "(run the tuning-browse path first, or pass --reuse-fits)")
        fits_df = pd.read_csv(csv_path)
        log.info("tuning-significance summary: %d fit rows from %s",
                 len(fits_df), csv_path)
        plot_tuning_significance_summary(
            fits_df, OUT_DIR / "V_tuning_significance_summary", condition="V")
        plot_tuning_transitions(
            fits_df, OUT_DIR / "VtoVT_tuning_significance_transitions")
        return 0

    if args.fig2_summary:
        from run_glm_goggles_rf_sfor_20ms import PROBES
        mc = pd.read_csv(RUN_ALL_DIR / "glm_model_comparison.csv")
        vp = pd.read_csv(RUN_ALL_DIR / "diagnostics" / "variance_partition.csv")
        log.info("fig2: %d clusters (model_comparison), %d (variance_partition)",
                 len(mc), len(vp))
        inter_per = compute_interaction_uniques(config, PROBES, mc, backend)
        plot_fig2_summary(mc, vp, inter_per,
                          OUT_DIR / "fig2_forward_selection_summary", config)
        plot_fig2_summary(mc, vp, inter_per,
                          OUT_DIR / "fig2_forward_selection_summary_zoom", config,
                          zoom_split=True)

    if args.fig2_forward_deltas:
        from run_glm_goggles_rf_sfor_20ms import PROBES
        mc = pd.read_csv(RUN_ALL_DIR / "glm_model_comparison.csv")
        fvp, inter_fwd = _load_forward_deltas(RUN_ALL_DIR, PROBES)
        log.info("fig2 forward-Δ: %d clusters (model_comparison), %d (selection history), "
                 "%d interaction terms", len(mc), len(fvp), len(inter_fwd))
        plot_fig2_forward_deltas(
            mc, fvp, inter_fwd,
            OUT_DIR / "fig2_forward_selection_summary_zoom_forwarddeltas")
        plot_fig2_forward_deltas(
            mc, fvp, inter_fwd,
            OUT_DIR / "fig2_forward_selection_summary_forwarddeltas_logheat",
            log_heat=True)
        plot_fig2_forward_deltas(
            mc, fvp, inter_fwd,
            OUT_DIR / "fig2_forward_selection_summary_forwarddeltas_unitheat",
            unit_heat=True)
        return 0
        return 0

    # Tuning-only browsing grids over one or more probes (loads each itself).
    # mean_sd → its own folder so the median+IQR set isn't overwritten.
    if args.top_fr is not None or args.all_clusters:
        probes = args.probes if args.probes else [args.probe]
        sub = {"mean_sd": "tuning_browse_mean_sd",
               "mean_sem": "tuning_browse_mean_sem"}.get(args.display, "tuning_browse")
        if args.select == "rsq_mean":
            sub += "_selR2mean"  # alternative selection → its own folder
        out_dir = OUT_DIR / sub
        fits_lookup = None
        if args.reuse_fits:
            fits_lookup = _load_fits_lookup(args.reuse_fits)
            log.info("reusing %d cached fits from %s (no re-fit/bootstrap)",
                     len(fits_lookup), args.reuse_fits)
        all_rows = []
        for probe_id in probes:
            all_rows += run_tuning_browse(
                probe_id, config, out_dir,
                all_clusters=args.all_clusters, top_fr=args.top_fr,
                display=args.display, select_criterion=args.select,
                fits_lookup=fits_lookup, tuning_n_reps=args.tuning_n_reps,
            )
        # One combined CSV of the BIC-selected fits (cluster × value × condition)
        # across all probes rendered this run.
        if all_rows:
            csv_path = out_dir / "tuning_fits.csv"
            pd.DataFrame(all_rows).to_csv(csv_path, index=False)
            log.info("wrote %d tuning-fit rows → %s", len(all_rows), csv_path)
        return 0

    log.info("loading %s ...", args.probe)
    probe = load_probe_data(
        FORMATTED_DIR / f"{args.probe}.mat", config=config,
        stimulus_lookup=_lookup(), cluster_set="selected",
    )
    rf_lookup = load_rf_sf_or(
        config.rf_sf_or_parquet_dir, probe.probe_id,
        min_concentration=getattr(config, "rf_min_concentration", 0.0),
    )
    by_id = {c.cluster_id: c for c in probe.clusters}
    trials_by_id = {t.trial_id: t for t in probe.trials}
    # MATLAB precomputed TF/Speed tuning cache (per-condition 20 quantile bins +
    # per-trial FR); TF tuning is read from here so it matches MATLAB's PDFs.
    pc = load_precomputed_bin_edges(FORMATTED_DIR / f"{args.probe}.mat")

    if args.trial_structure:
        # --match-clouds: ONE combined figure per condition. All --clusters become
        # the tuning columns; the FIRST is the primary (drives the trial-structure
        # block + its FR + the matched-cloud trial pick). FR y-limits are harmonised
        # across the cluster columns AND between the V and VT figures.
        if args.match_clouds:
            tcl = []  # (cid, cluster_df, cluster)
            for cid in args.clusters:
                if cid not in by_id:
                    log.warning("cluster %d not on probe %s — skipping", cid, args.probe)
                    continue
                cdf = bin_cluster(probe, by_id[cid], rf_lookup=rf_lookup)
                if cdf.empty:
                    log.warning("cluster %d: empty binned df — skipping", cid)
                    continue
                tcl.append((cid, cdf, by_id[cid]))
            if not tcl:
                log.warning("no usable clusters for --match-clouds")
                return 0
            primary_cid, primary_df, primary_cluster = tcl[0]
            pair = select_matched_cloud_trials(primary_df, trials_by_id)
            if pair is None:
                log.warning("primary cl %d: no cloud shared by V and VT — abort",
                            primary_cid)
                return 0
            cloud, v_tid, vt_tid, v_spk, vt_spk = pair
            log.info("matched cloud %s (primary cl %d) → V trial %d (%d spk) / VT "
                     "trial %d (%d spk); tuning cols %s", cloud, primary_cid, v_tid,
                     v_spk, vt_tid, vt_spk, [c for c, _, _ in tcl])
            tuning_clusters = [(c, d) for c, d, _ in tcl]
            ids = "-".join(str(c) for c, _, _ in tcl)
            rendered = []  # (fig, tun_axes, out_path)
            for cond, ctid in (("V", v_tid), ("VT", vt_tid)):
                trial = trials_by_id.get(ctid)
                if trial is None or np.flatnonzero(trial.motion_mask).size == 0:
                    log.warning("%s trial %s missing/no motion — skipping", cond, ctid)
                    continue
                out = (OUT_DIR / f"trial_structure_{probe.probe_id}_cl{ids}"
                                 f"_{cond}_trial_{ctid}")
                fig, taxes = plot_trial_structure(
                    probe.probe_id, primary_cid, primary_df, config, ctid, trial,
                    primary_cluster.spike_times, out, trials_by_id,
                    tuning_n_reps=args.tuning_n_reps, pc=pc, condition=cond,
                    tuning_clusters=tuning_clusters, save=False,
                    fr_stat=args.fr_stat)
                rendered.append((fig, taxes, out))
            # Shared per-regressor FR y-limit across the V and VT figures.
            ymax: dict[str, float] = {}
            for _fig, taxes, _out in rendered:
                for key, axlist in taxes.items():
                    for ax in axlist:
                        ymax[key] = max(ymax.get(key, 0.0), float(ax.get_ylim()[1]))
            for fig, taxes, out in rendered:
                for key, axlist in taxes.items():
                    for ax in axlist:
                        ax.set_ylim(0.0, ymax[key])
                out.parent.mkdir(parents=True, exist_ok=True)
                for ext in ("pdf", "png"):
                    fig.savefig(out.with_suffix(f".{ext}"), dpi=150)
                plt.close(fig)
                log.info("wrote %s.{pdf,png}", out)
            return 0

        # Single-cluster path (per --condition / --trial).
        for cid in args.clusters:
            if cid not in by_id:
                log.warning("cluster %d not on probe %s — skipping", cid, args.probe)
                continue
            cluster = by_id[cid]
            df = bin_cluster(probe, cluster, rf_lookup=rf_lookup)
            if df.empty:
                log.warning("cluster %d: empty binned df — skipping", cid)
                continue
            if args.trial is not None:
                tid = args.trial
            else:
                pick = select_v_trial(df, condition=args.condition)
                if pick is None:
                    log.warning("cluster %d: no %s trial available — skipping",
                                cid, args.condition)
                    continue
                tid, motion_spk, base_spk = pick
                log.info("cluster %d: %s trial %d (%d motion / %d baseline spk)",
                         cid, args.condition, tid, motion_spk, base_spk)
            trial = trials_by_id.get(tid)
            if trial is None or np.flatnonzero(trial.motion_mask).size == 0:
                log.warning("cluster %d: trial %s missing/no motion — skipping", cid, tid)
                continue
            plot_trial_structure(
                probe.probe_id, cid, df, config, tid, trial, cluster.spike_times,
                OUT_DIR / f"trial_structure_{probe.probe_id}_cluster_{cid}"
                          f"_{args.condition}_trial_{tid}",
                trials_by_id, tuning_n_reps=args.tuning_n_reps, pc=pc,
                condition=args.condition, fr_stat=args.fr_stat)
        return 0

    if args.fig2c:
        if args.trial is None:
            ap.error("--fig2c requires --trial")
        accepted = [s for s in args.accepted.split(",") if s]
        for cid in args.clusters:
            cluster = by_id[cid]
            df = bin_cluster(probe, cluster, rf_lookup=rf_lookup)
            prep = prepare_cluster_design(df, config)
            oracle_check(prep, config, probe.probe_id, cid)
            selected = load_selection_order(probe.probe_id, cid) or []
            baseline_vars = ["History"] if "History" in selected else []
            trial = trials_by_id[args.trial]
            midx = np.flatnonzero(trial.motion_mask)
            t_motion_start = float(trial.probe_t[int(midx[0])])
            plot_forward_panel(
                probe.probe_id, cid, df, prep, config, args.trial,
                cluster.spike_times, t_motion_start,
                baseline_vars, accepted, args.rejected,
                OUT_DIR / f"fig2c_{probe.probe_id}_cluster_{cid}_trial_{args.trial}",
            )
        return 0

    for cid in args.clusters:
        if cid not in by_id:
            log.warning("cluster %d not on probe %s — skipping", cid, args.probe)
            continue
        cluster = by_id[cid]
        df = bin_cluster(probe, cluster, rf_lookup=rf_lookup)
        if df.empty:
            log.warning("cluster %d: empty binned df — skipping", cid)
            continue
        prep = prepare_cluster_design(df, config)
        oracle_check(prep, config, probe.probe_id, cid)

        selected = load_selection_order(probe.probe_id, cid)
        if not selected:
            log.warning("cluster %d: no forward-selection record — skipping", cid)
            continue
        steps = buildup_steps_for(selected)
        log.info("cluster %d: forward-selected = %s | %d rows",
                 cid, selected, len(steps))

        # Buildup cv-bps + in-sample rates (cluster-wide; same for every trial).
        step_results = []
        for label, vars_ in steps:
            cv = cv_bps_for(prep, vars_, config, backend)
            rate = insample_rate_for(prep, vars_, config, backend)
            step_results.append(dict(label=label, vars=vars_, cv=cv, rate=rate))
            log.info("cluster %d %-26s cv-bps=%.4f", cid, label, cv)

        trials = select_vt_trials(df, args.n_trials, args.baseline_min, args.baseline_max)
        if not trials:
            log.warning("cluster %d: no VT trials with %d-%d baseline spikes — "
                        "skipping (too sparse in baseline)", cid,
                        args.baseline_min, args.baseline_max)
            continue
        log.info("cluster %d: %d VT trials selected: %s", cid, len(trials),
                 [(t, b) for t, _, b in trials])
        for tid, _motion, baseline in trials:
            trial = trials_by_id.get(tid)
            if trial is None:
                continue
            midx = np.flatnonzero(trial.motion_mask)
            if midx.size == 0:
                continue
            t_motion_start = float(trial.probe_t[int(midx[0])])
            plot_buildup(
                probe.probe_id, cid, df, config, step_results, tid,
                cluster.spike_times, t_motion_start,
                OUT_DIR / f"buildup_{probe.probe_id}_cluster_{cid}_trial_{tid}",
                baseline_spk=baseline,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
