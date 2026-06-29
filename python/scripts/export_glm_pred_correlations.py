#!/usr/bin/env python
"""Per-trial, per-cluster GLM predicted FR + correlation vs a 20 ms-Gaussian
observed FR, for the committed goggles run ``current_rf_sfor_20ms_histbase_goggles_all``.

What it produces (into ``<run>/predictions/``):
  1. Predicted FR data matrices  — predicted[cluster, trial, bin] (Hz), .npy + .mat
  2. Observed FR data matrices   — 20 ms-Gaussian smoothed obs FR, same shape (.npy + .mat)
  3. Per-trial correlations      — Pearson r(predicted, observed) over the WHOLE trial
                                   (baseline + motion), correlations[cluster, trial], .npy + .mat
  4. A per-cluster correlation boxplot (clusters on x, per-trial r on y).

How (faithfulness, no re-selection, no re-fit):
  The Selected-model prediction is reconstructed from the COMMITTED coefficients
  (``glm_coefficients.csv``). Bases are rebuilt with the production basis
  functions and the production config (``make_config_histbase_all``); the design
  is assembled with the production ``assemble_design_matrix_selected``; the
  committed beta is matched to the rebuilt design BY COLUMN NAME (fail-loud on any
  mismatch). To prove the rebuilt design is byte-identical to the run's, a sample
  of clusters is RE-FIT and the recovered beta is asserted to match the committed
  beta (beta is deterministic given data + lambda). If that parity check fails the
  script raises — a basis/design drift cannot ship silently.

  Predicted FR = exp(clip(X @ beta)) in Hz (the log(bin_width) offset converts the
  Poisson expected-count to a rate, exactly as pipeline._fit_plot_models does).

Run (from the repo ``python/`` dir, with the project venv):
    python -m scripts.export_glm_pred_correlations
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import savemat
from scipy.ndimage import gaussian_filter1d

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Mirror the driver: put python/ (parent of scripts/) on the path so both
# ``rc2_glm`` and ``scripts.*`` import whether launched as a file or ``-m``.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster
from rc2_glm.rf_sf_or import load_rf_sf_or
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.fitting import fit_poisson_glm
from rc2_glm.basis import (
    value_basis,
    onset_kernel_basis,
    history_basis,
    convolve_history,
    raised_cosine_basis_linear,
    circular_basis,
)

# Reuse the run's exact config + paths — the single source of truth for what the
# committed run was fit with (motion mask, bases, ranges, rf-local SF/OR, etc.).
from scripts.run_glm_goggles_rf_sfor_20ms import (  # noqa: E402
    make_config_histbase_all,
    _lookup,
    FORMATTED_DIR,
    RF_PARQUET_DIR,
    PROBES,
    OUT_HISTBASE_ALL,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("glm_pred_corr")

RUN_DIR = OUT_HISTBASE_ALL
OUT_DIR = RUN_DIR / "predictions"
MODEL = "Selected"

# Observed FR smoothing: a sigma = 20 ms Gaussian (Laura's "pink" MATLAB
# FiringRate convention) — deliberately NOT the repo's 100 ms display boxcar
# (_smoothed_obs_per_trial) and NOT the GLM's raw-count target.
OBS_GAUSS_SIGMA_S = 0.02

# Per-trial plot colours: pink = MEASURED (the 20 ms-Gaussian observed FR, the
# MATLAB "pink" FiringRate convention), black = GLM predicted rate.
PINK = "#E75480"

# Reconstruction parity gate: re-fit this many clusters per probe and assert the
# recovered beta matches the committed beta (max-abs) below this tolerance.
PARITY_N = 5
PARITY_TOL = 1e-3


# --------------------------------------------------------------------------- #
# Basis construction — replicates pipeline._fit_one_cluster (lines 493-612)
# exactly, calling the SAME production basis functions. Validated by the refit
# parity gate in run_probe().
# --------------------------------------------------------------------------- #
def build_cluster_bases(df: pd.DataFrame, config) -> dict:
    speed = df["speed"].to_numpy(np.float64)
    tf = df["tf"].to_numpy(np.float64)
    onset = df["time_since_onset"].to_numpy(np.float64)
    sf_vals = df["sf"].to_numpy(np.float64)
    or_vals = df["orientation"].to_numpy(np.float64)
    y = df["spike_count"].to_numpy(np.float64)
    trial_ids = df["trial_id"].to_numpy(np.int64)
    motion_rows = (df["condition"] != "stationary").to_numpy()

    spacing = getattr(config, "speed_tf_basis_spacing", "log")
    B_speed = value_basis(speed, config.n_speed_bases, *config.speed_range, spacing=spacing)
    B_tf = value_basis(tf, config.n_tf_bases, *config.tf_range, spacing=spacing)
    B_onset = onset_kernel_basis(onset, config.n_onset_bases, config.onset_range[1])

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

    B_me_face = None
    if "me_face_raw" in df.columns and getattr(config, "include_me_face", True):
        me_raw = df["me_face_raw"].to_numpy(np.float64)
        fin = np.isfinite(me_raw) & motion_rows
        if int(fin.sum()) >= 10:
            m = float(me_raw[fin].mean())
            s = float(me_raw[fin].std(ddof=0)) or 1.0
            me_z = np.where(np.isfinite(me_raw), (me_raw - m) / s, 0.0)
            B_me_face = raised_cosine_basis_linear(
                me_z, config.n_me_face_bases,
                config.me_face_range[0], config.me_face_range[1],
            )

    B_accel = None
    if "acceleration" in df.columns and getattr(config, "include_acceleration", False):
        acc_raw = df["acceleration"].to_numpy(np.float64)
        fin = np.isfinite(acc_raw) & motion_rows
        if int(fin.sum()) >= 10:
            am = float(acc_raw[fin].mean())
            asd = float(acc_raw[fin].std(ddof=0)) or 1.0
            acc_z = np.where(np.isfinite(acc_raw), (acc_raw - am) / asd, 0.0)
            acc_z = np.clip(acc_z, config.accel_range[0], config.accel_range[1])
            B_accel = raised_cosine_basis_linear(
                acc_z, config.n_accel_bases,
                config.accel_range[0], config.accel_range[1],
            )

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
    else:
        B_sf = None
        B_or = None

    sf_valid = sf_vals[(sf_vals != 0.0) & ~np.isnan(sf_vals)]
    sf_ref_levels = np.sort(np.unique(sf_valid)).tolist() if sf_valid.size else []
    or_valid = or_vals[(or_vals != 0.0) & ~np.isnan(or_vals)]
    or_ref_levels = np.sort(np.unique(or_valid)).tolist() if or_valid.size else []

    return dict(
        B_speed=B_speed, B_tf=B_tf, B_onset=B_onset,
        sf_vals=sf_vals, or_vals=or_vals, y=y,
        B_history=B_history, B_me_face=B_me_face, B_accel=B_accel,
        B_sf=B_sf, B_or=B_or,
        sf_ref_levels=sf_ref_levels, or_ref_levels=or_ref_levels,
    )


def assemble_selected(bases: dict, selected_vars: list[str], config):
    return assemble_design_matrix_selected(
        bases["B_speed"], bases["B_tf"], bases["B_onset"],
        bases["sf_vals"], bases["or_vals"], selected_vars,
        sf_ref_levels=bases["sf_ref_levels"], or_ref_levels=bases["or_ref_levels"],
        B_history=bases["B_history"], B_me_face=bases["B_me_face"],
        B_accel=bases["B_accel"], B_sf=bases["B_sf"], B_or=bases["B_or"],
        include_onset_kernel=getattr(config, "include_onset_kernel", True),
        history_in_baseline=getattr(config, "history_in_baseline", False),
    )


def parse_selected_vars(s) -> list[str]:
    """``time_selected_vars`` is a '+'-joined string; interactions use '_x_'
    (no '+'), so splitting on '+' is unambiguous. NaN / '' → Null (no vars)."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return []
    s = str(s).strip()
    if not s or s.lower() == "nan":
        return []
    return [v for v in s.split("+") if v]


# --------------------------------------------------------------------------- #
# Per-probe processing
# --------------------------------------------------------------------------- #
def run_probe(probe: str, config) -> dict:
    run_probe_dir = RUN_DIR / "_runs" / probe
    coef = pd.read_csv(run_probe_dir / "glm_coefficients.csv")
    comp = pd.read_csv(run_probe_dir / "glm_model_comparison.csv")
    coef_sel = coef[coef["model"] == MODEL]

    # Cohort = clusters that were actually fit (have a Selected model). Spike-floor
    # exclusions never wrote coefficients, so they drop out here automatically.
    cohort = sorted(set(coef_sel["cluster_id"]) & set(comp["cluster_id"]))
    log.info("%s: %d clusters with a %s model", probe, len(cohort), MODEL)

    sel_by_cluster = dict(zip(comp["cluster_id"], comp["time_selected_vars"]))

    probe_data = load_probe_data(
        FORMATTED_DIR / f"{probe}.mat", config=config,
        stimulus_lookup=_lookup(), cluster_set="selected",
    )
    rf_lookup = load_rf_sf_or(
        RF_PARQUET_DIR, probe_data.probe_id,
        min_concentration=getattr(config, "rf_min_concentration", 0.0),
    )
    clusters_by_id = {c.cluster_id: c for c in probe_data.clusters}

    offset = float(np.log(config.time_bin_width))
    sigma_bins = OBS_GAUSS_SIGMA_S / config.time_bin_width

    # First pass: per-cluster predicted/observed/correlation in a tidy dict, and
    # capture the shared (trial -> time axis) grid from the first cluster.
    per_cluster: dict[int, dict] = {}
    model_labels: dict[int, str] = {}   # cluster_id -> chosen-model string (for plot titles)
    trial_order: list[int] | None = None
    trial_times: dict[int, np.ndarray] = {}
    parity_done = 0

    for cid in cohort:
        cluster = clusters_by_id.get(cid)
        if cluster is None:
            log.warning("%s cluster %d: not in loaded probe, skipping", probe, cid)
            continue
        df = bin_cluster(probe_data, cluster, rf_lookup=rf_lookup)
        if df.empty:
            continue
        bases = build_cluster_bases(df, config)
        selected_vars = parse_selected_vars(sel_by_cluster.get(cid))
        model_labels[int(cid)] = "+".join(selected_vars) if selected_vars else "Null"
        X, names = assemble_selected(bases, selected_vars, config)

        # Committed beta, matched to the rebuilt design BY NAME (fail-loud).
        cc = coef_sel[coef_sel["cluster_id"] == cid]
        committed = dict(zip(cc["coefficient"], cc["estimate"]))
        if set(names) != set(committed):
            missing = sorted(set(names) - set(committed))
            extra = sorted(set(committed) - set(names))
            raise RuntimeError(
                f"{probe} cluster {cid}: design/coeff name mismatch. "
                f"rebuilt-but-not-in-CSV={missing} ; CSV-but-not-rebuilt={extra}"
            )
        beta = np.array([committed[n] for n in names], dtype=np.float64)

        # Parity gate: re-fit the first PARITY_N clusters and assert recovered
        # beta == committed beta (proves the rebuilt design matches the run's).
        if parity_done < PARITY_N:
            fit = fit_poisson_glm(
                X, bases["y"], offset,
                lambda_ridge=config.lambda_ridge, backend="irls",
            )
            max_abs = float(np.max(np.abs(np.asarray(fit.beta) - beta)))
            status = "OK" if max_abs <= PARITY_TOL else "FAIL"
            log.info("%s cluster %d: refit-vs-committed beta max|Δ|=%.2e [%s]",
                     probe, cid, max_abs, status)
            if max_abs > PARITY_TOL:
                raise RuntimeError(
                    f"{probe} cluster {cid}: reconstruction parity FAILED "
                    f"(max|Δβ|={max_abs:.3e} > {PARITY_TOL}). Design rebuild "
                    f"does not match the committed run; aborting."
                )
            parity_done += 1

        pred = np.exp(np.clip(X @ beta, -20.0, 20.0))  # Hz
        obs_rate = df["spike_count"].to_numpy(np.float64) / config.time_bin_width
        tids = df["trial_id"].to_numpy(np.int64)
        tsince = df["time_since_onset"].to_numpy(np.float64)

        # Per-trial: order bins by within-trial time, smooth observed (20 ms
        # Gaussian), correlate whole-trial predicted vs smoothed observed.
        trials = sorted(np.unique(tids).tolist())
        if trial_order is None:
            trial_order = trials
        cl_pred: dict[int, np.ndarray] = {}
        cl_obs: dict[int, np.ndarray] = {}
        cl_r: dict[int, float] = {}
        for t in trials:
            rows = np.where(tids == t)[0]
            order = rows[np.argsort(tsince[rows], kind="stable")]
            p = pred[order]
            o = gaussian_filter1d(obs_rate[order], sigma=sigma_bins, mode="nearest")
            cl_pred[t] = p
            cl_obs[t] = o
            if p.std() > 0 and o.std() > 0 and p.size >= 3:
                cl_r[t] = float(np.corrcoef(p, o)[0, 1])
            else:
                cl_r[t] = np.nan
            if t not in trial_times:
                trial_times[t] = tsince[order]
        per_cluster[cid] = dict(pred=cl_pred, obs=cl_obs, r=cl_r)

    # Assemble dense matrices (cluster, trial, bin) NaN-padded.
    cluster_ids = np.array(sorted(per_cluster), dtype=np.int64)
    trial_ids = np.array(trial_order, dtype=np.int64)
    n_bins_max = max(t.size for t in trial_times.values())
    C, T, B = cluster_ids.size, trial_ids.size, n_bins_max

    predicted = np.full((C, T, B), np.nan, dtype=np.float32)
    observed = np.full((C, T, B), np.nan, dtype=np.float32)
    correlations = np.full((C, T), np.nan, dtype=np.float64)
    bin_time_s = np.full((T, B), np.nan, dtype=np.float32)
    for ti, t in enumerate(trial_ids):
        tt = trial_times[int(t)]
        bin_time_s[ti, : tt.size] = tt
    for ci, cid in enumerate(cluster_ids):
        d = per_cluster[int(cid)]
        for ti, t in enumerate(trial_ids):
            t = int(t)
            if t in d["pred"]:
                p = d["pred"][t]
                predicted[ci, ti, : p.size] = p
                observed[ci, ti, : d["obs"][t].size] = d["obs"][t]
                correlations[ci, ti] = d["r"][t]

    return dict(
        probe=probe,
        cluster_ids=cluster_ids, trial_ids=trial_ids, bin_time_s=bin_time_s,
        predicted=predicted, observed=observed, correlations=correlations,
        model_labels=[model_labels.get(int(c), "Null") for c in cluster_ids],
    )


def save_probe(res: dict) -> None:
    probe = res["probe"]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # .npy (one array each — predicted, observed, correlations)
    np.save(OUT_DIR / f"predicted_fr_{probe}.npy", res["predicted"])
    np.save(OUT_DIR / f"observed_fr_{probe}.npy", res["observed"])
    np.save(OUT_DIR / f"correlations_{probe}.npy", res["correlations"])
    # (the per-.npy index JSON is written by write_index())
    # .mat (bundle with index arrays)
    savemat(OUT_DIR / f"predictions_{probe}.mat", {
        "probe": probe,
        "model": MODEL,
        "obs_gauss_sigma_s": OBS_GAUSS_SIGMA_S,
        "cluster_ids": res["cluster_ids"],
        "trial_ids": res["trial_ids"],
        "bin_time_s": res["bin_time_s"],
        "predicted_fr": res["predicted"],
        "observed_fr": res["observed"],
        "correlations": res["correlations"],
    }, do_compression=True)
    log.info("%s: wrote .npy + .mat (predicted %s)", probe, res["predicted"].shape)


def write_index(res: dict, bin_width: float) -> None:
    probe = res["probe"]
    (OUT_DIR / f"index_{probe}.json").write_text(json.dumps({
        "probe": probe,
        "model": MODEL,
        "obs_gauss_sigma_s": OBS_GAUSS_SIGMA_S,
        "bin_width_s": bin_width,
        "n_clusters": int(res["cluster_ids"].size),
        "n_trials": int(res["trial_ids"].size),
        "n_bins_max": int(res["predicted"].shape[2]),
        "cluster_ids": res["cluster_ids"].tolist(),
        "trial_ids": res["trial_ids"].tolist(),
        "axes": {
            "predicted_fr": "[cluster, trial, bin] GLM predicted firing rate (Hz), NaN-padded",
            "observed_fr": "[cluster, trial, bin] 20ms-Gaussian smoothed observed FR (Hz)",
            "correlations": "[cluster, trial] Pearson r: predicted rate vs observed, whole trial",
        },
    }, indent=2))


def make_boxplot(results: list[dict]) -> None:
    # One box per cluster across both probes; per-trial r on y.
    labels, data, colours = [], [], []
    palette = {0: "#4C72B0", 1: "#C44E52"}
    for pi, res in enumerate(results):
        for ci, cid in enumerate(res["cluster_ids"]):
            r = res["correlations"][ci]
            r = r[np.isfinite(r)]
            if r.size == 0:
                continue
            labels.append(f"{res['probe'].split('_')[0].replace('CAA-', '')}·{cid}")
            data.append(r)
            colours.append(palette[pi % 2])
    fig_w = max(8.0, 0.18 * len(data))
    fig, ax = plt.subplots(figsize=(fig_w, 5.0))
    bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6)
    for patch, col in zip(bp["boxes"], colours):
        patch.set_facecolor(col)
        patch.set_alpha(0.6)
    for med in bp["medians"]:
        med.set_color("black")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=90, fontsize=5)
    ax.axhline(0.0, color="grey", lw=0.7, ls="--")
    ax.set_ylabel("Pearson r  (predicted vs 20 ms-Gaussian observed FR)")
    ax.set_xlabel("cluster  (probe·id)")
    ax.set_title("Per-trial GLM prediction–observation correlation, by cluster\n"
                 "goggles · current_rf_sfor_20ms_histbase_goggles_all · Selected model")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=palette[i], alpha=0.6) for i in (0, 1)]
    ax.legend(handles, [r["probe"].split("_")[0] for r in results],
              fontsize=7, loc="lower right")
    fig.tight_layout()
    out = OUT_DIR / "correlation_boxplot.pdf"
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s", out)


def save_combined(results: list[dict]) -> None:
    # Combined .mat/.npy of the per-trial correlations (clusters stacked across
    # probes), plus a tidy long-form CSV summary (median r per cluster).
    rows = []
    for res in results:
        for ci, cid in enumerate(res["cluster_ids"]):
            r = res["correlations"][ci]
            rf = r[np.isfinite(r)]
            rows.append({
                "probe": res["probe"], "cluster_id": int(cid),
                "n_trials": int(np.isfinite(r).sum()),
                "median_r": float(np.nanmedian(r)) if rf.size else np.nan,
                "mean_r": float(np.nanmean(r)) if rf.size else np.nan,
                "q1_r": float(np.nanquantile(r, 0.25)) if rf.size else np.nan,
                "q3_r": float(np.nanquantile(r, 0.75)) if rf.size else np.nan,
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "correlation_summary.csv", index=False)
    savemat(OUT_DIR / "correlations_all.mat", {
        "probe": [res["probe"] for res in results],
        "summary_probe": summary["probe"].to_numpy(dtype=object),
        "summary_cluster_id": summary["cluster_id"].to_numpy(),
        "summary_median_r": summary["median_r"].to_numpy(),
    }, do_compression=True)
    log.info("wrote correlation_summary.csv (%d clusters) + correlations_all.mat",
             len(summary))


def make_cluster_mean_swarm(results: list[dict]) -> None:
    # ONE violin over the per-cluster MEAN correlations (each cluster = 1 point),
    # points overlaid. Single colour (no probe split); aesthetic matches the FENS
    # _violin_panel convention (make_fens_poster_figures.py: violin alpha 0.30 with
    # facecolor=edgecolor=colour, jittered strip s=10 edgecolor "0.2", median bar,
    # n= annotation, faint zero line).
    means = []
    for res in results:
        for ci in range(res["cluster_ids"].size):
            r = res["correlations"][ci]
            r = r[np.isfinite(r)]
            if r.size == 0:
                continue
            means.append(float(np.mean(r)))
    means = np.asarray(means)
    c = "#4C72B0"
    fig, ax = plt.subplots(figsize=(3.2, 5.2))
    if means.size >= 2 and np.ptp(means) > 0:
        parts = ax.violinplot([means], positions=[0], widths=0.8,
                              showmeans=False, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(c)
            b.set_edgecolor(c)
            b.set_alpha(0.30)
    jit = (np.random.RandomState(0).rand(means.size) - 0.5) * 0.28
    ax.scatter(jit, means, s=10, color=c, edgecolor="0.2",
               linewidth=0.3, alpha=0.85, zorder=3)
    ax.hlines(np.median(means), -0.22, 0.22, color="0.1", lw=1.4, zorder=4)
    ax.annotate(f"n={means.size}", (0, means.max()), textcoords="offset points",
                xytext=(0, 3), ha="center", va="bottom", fontsize=7)
    ax.axhline(0, color="0.85", lw=0.6)
    ax.set_xticks([0])
    ax.set_xticklabels(["all clusters"], fontsize=8)
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylabel("mean per-trial Pearson r (per cluster)")
    ax.set_title("Per-cluster mean GLM\nprediction–observation correlation", fontsize=10)
    fig.tight_layout()
    out = OUT_DIR / "correlation_clustermean_swarm.pdf"
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s (median of cluster means = %.3f)", out, float(np.median(means)))


def make_violin_plot(results: list[dict]) -> None:
    # Same per-cluster data as the boxplot, as violins (per-trial r distribution
    # per cluster), one violin per cluster, coloured by probe, medians marked.
    labels, data, colours = [], [], []
    palette = {0: "#4C72B0", 1: "#C44E52"}
    for pi, res in enumerate(results):
        for ci, cid in enumerate(res["cluster_ids"]):
            r = res["correlations"][ci]
            r = r[np.isfinite(r)]
            if r.size < 2:   # a violin needs spread
                continue
            labels.append(f"{res['probe'].split('_')[0].replace('CAA-', '')}·{cid}")
            data.append(r)
            colours.append(palette[pi % 2])
    x = np.arange(1, len(labels) + 1)
    fig_w = max(8.0, 0.18 * len(labels))
    fig, ax = plt.subplots(figsize=(fig_w, 5.0))
    parts = ax.violinplot(data, positions=x, widths=0.8, showmedians=True,
                          showextrema=False)
    for body, col in zip(parts["bodies"], colours):
        body.set_facecolor(col)
        body.set_edgecolor("none")
        body.set_alpha(0.6)
    if "cmedians" in parts:
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=5)
    ax.axhline(0.0, color="grey", lw=0.7, ls="--")
    ax.set_ylabel("Pearson r  (predicted vs 20 ms-Gaussian observed FR)")
    ax.set_xlabel("cluster  (probe·id)")
    ax.set_title("Per-trial GLM prediction–observation correlation, by cluster (violin)\n"
                 "goggles · current_rf_sfor_20ms_histbase_goggles_all · Selected model")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=palette[i], alpha=0.6) for i in (0, 1)]
    ax.legend(handles, [r["probe"].split("_")[0] for r in results],
              fontsize=7, loc="lower right")
    fig.tight_layout()
    out = OUT_DIR / "correlation_violin.pdf"
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s", out)


def make_mean_sem_plot(results: list[dict]) -> None:
    # Same per-cluster data as the boxplot, shown as mean ± SEM of the per-trial
    # Pearson r (one point per cluster; SEM = std(ddof=1)/sqrt(n_trials)).
    labels, means, sems, colours = [], [], [], []
    palette = {0: "#4C72B0", 1: "#C44E52"}
    for pi, res in enumerate(results):
        for ci, cid in enumerate(res["cluster_ids"]):
            r = res["correlations"][ci]
            r = r[np.isfinite(r)]
            if r.size == 0:
                continue
            labels.append(f"{res['probe'].split('_')[0].replace('CAA-', '')}·{cid}")
            means.append(float(np.mean(r)))
            sems.append(float(np.std(r, ddof=1) / np.sqrt(r.size)) if r.size > 1 else 0.0)
            colours.append(palette[pi % 2])
    x = np.arange(1, len(labels) + 1)
    fig_w = max(8.0, 0.18 * len(labels))
    fig, ax = plt.subplots(figsize=(fig_w, 5.0))
    ax.errorbar(x, means, yerr=sems, fmt="none", ecolor="grey", elinewidth=0.8,
                capsize=2, zorder=1)
    ax.scatter(x, means, c=colours, s=14, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=5)
    ax.axhline(0.0, color="grey", lw=0.7, ls="--")
    ax.set_ylabel("mean ± SEM Pearson r  (predicted vs 20 ms-Gaussian observed FR)")
    ax.set_xlabel("cluster  (probe·id)")
    ax.set_title("Per-trial GLM prediction–observation correlation, by cluster (mean ± SEM)\n"
                 "goggles · current_rf_sfor_20ms_histbase_goggles_all · Selected model")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=palette[i], alpha=0.8) for i in (0, 1)]
    ax.legend(handles, [r["probe"].split("_")[0] for r in results],
              fontsize=7, loc="lower right")
    fig.tight_layout()
    out = OUT_DIR / "correlation_mean_sem.pdf"
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s", out)


def make_per_trial_pdfs(results: list[dict], *, y_shared: bool = False,
                        subdir: str = "plots") -> None:
    """One multipage PDF per cluster (<subdir>/<probe>/cluster_<id>.pdf), one
    trial per page: pink = measured (20 ms-Gaussian observed FR), black = GLM
    predicted rate; per-trial Pearson r annotated; the chosen model in the title.

    y_shared=False: each page autoscales independently. y_shared=True: every page
    of a cluster shares ylim [0, 1.05 * max over ALL that cluster's trials of both
    predicted and measured FR] — a per-cluster common scale (literal max; note a
    rare history-driven predicted spike can dominate it)."""
    plots_dir = OUT_DIR / subdir
    fig, ax = plt.subplots(figsize=(7.0, 3.2))  # reused across all pages
    for res in results:
        probe = res["probe"]
        pdir = plots_dir / probe
        pdir.mkdir(parents=True, exist_ok=True)
        cids = res["cluster_ids"]
        tids = res["trial_ids"]
        bt = res["bin_time_s"]
        pred = res["predicted"]
        obs = res["observed"]
        corr = res["correlations"]
        labels = res["model_labels"]
        for ci, cid in enumerate(cids):
            model = labels[ci]
            out = pdir / f"cluster_{int(cid)}.pdf"
            # Per-cluster shared y-max = literal max of both traces over all trials.
            ymax = None
            if y_shared:
                cluster_max = np.nanmax([
                    np.nanmax(pred[ci]) if np.isfinite(pred[ci]).any() else np.nan,
                    np.nanmax(obs[ci]) if np.isfinite(obs[ci]).any() else np.nan,
                ])
                if np.isfinite(cluster_max) and cluster_max > 0:
                    ymax = 1.05 * float(cluster_max)
            n_pages = 0
            with PdfPages(out) as pdf:
                for ti, tid in enumerate(tids):
                    x = bt[ti]
                    m = np.isfinite(x)
                    if int(m.sum()) < 2:
                        continue
                    xx = x[m]
                    o = obs[ci, ti][m]
                    p = pred[ci, ti][m]
                    if not np.isfinite(o).any() or not np.isfinite(p).any():
                        continue
                    ax.clear()
                    ax.plot(xx, o, color=PINK, lw=1.3, label="measured (20 ms Gaussian)")
                    ax.plot(xx, p, color="black", lw=1.0, label="GLM predicted")
                    ax.axvline(0.0, color="grey", lw=0.6, ls="--")  # motion onset
                    if ymax is not None:
                        ax.set_ylim(0.0, ymax)
                    r = corr[ci, ti]
                    rtxt = f"r = {r:.3f}" if np.isfinite(r) else "r = n/a"
                    ax.text(0.98, 0.95, rtxt, transform=ax.transAxes,
                            ha="right", va="top", fontsize=8,
                            bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.85))
                    ax.set_title(f"{probe} · cluster {int(cid)} · trial {int(tid)}\n"
                                 f"model: {model}", fontsize=8)
                    ax.set_xlabel("time from motion onset (s)", fontsize=8)
                    ax.set_ylabel("firing rate (Hz)", fontsize=8)
                    ax.legend(fontsize=6, loc="upper left", framealpha=0.85)
                    fig.tight_layout()
                    pdf.savefig(fig)
                    n_pages += 1
            log.info("%s cluster %d: wrote %s/%s (%d trial pages)",
                     probe, int(cid), subdir, out.name, n_pages)
    plt.close(fig)


def load_results_from_disk() -> list[dict]:
    """Rebuild the `results` structure from the saved predictions_<probe>.mat +
    the run's selected-vars, so --plots-only can (re)draw without reconstructing."""
    from scipy.io import loadmat
    results = []
    for probe in PROBES:
        matp = OUT_DIR / f"predictions_{probe}.mat"
        if not matp.exists():
            raise FileNotFoundError(
                f"{matp} not found — run the full export (no --plots-only) first")
        m = loadmat(matp)
        comp = pd.read_csv(RUN_DIR / "_runs" / probe / "glm_model_comparison.csv")
        sel = dict(zip(comp["cluster_id"], comp["time_selected_vars"]))
        cluster_ids = np.asarray(m["cluster_ids"]).ravel().astype(int)
        labels = ["+".join(parse_selected_vars(sel.get(int(c)))) or "Null"
                  for c in cluster_ids]
        results.append(dict(
            probe=probe,
            cluster_ids=cluster_ids,
            trial_ids=np.asarray(m["trial_ids"]).ravel().astype(int),
            bin_time_s=np.asarray(m["bin_time_s"]),
            predicted=np.asarray(m["predicted_fr"]),
            observed=np.asarray(m["observed_fr"]),
            correlations=np.asarray(m["correlations"]),
            model_labels=labels,
        ))
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--plots-only", action="store_true",
                    help="skip reconstruction; read the saved predictions_<probe>.mat "
                         "and (re)generate only the per-cluster trial PDFs")
    ap.add_argument("--which", choices=("auto", "yfixed", "both", "none"), default="both",
                    help="which per-trial PDF variant(s) to render: auto = "
                         "per-page autoscale (plots/), yfixed = per-cluster shared "
                         "y-max (plots_yfixed/), both (default), none = skip the "
                         "per-trial PDFs (summary figures only)")
    args = ap.parse_args()
    config = make_config_histbase_all()

    def render_plots(results: list[dict]) -> None:
        if args.which in ("auto", "both"):
            make_per_trial_pdfs(results, y_shared=False, subdir="plots")
        if args.which in ("yfixed", "both"):
            make_per_trial_pdfs(results, y_shared=True, subdir="plots_yfixed")

    if args.plots_only:
        log.info("plots-only (%s): loading saved matrices from %s", args.which, OUT_DIR)
        results = load_results_from_disk()
        make_boxplot(results)
        make_mean_sem_plot(results)
        make_violin_plot(results)
        make_cluster_mean_swarm(results)
        render_plots(results)
        log.info("done → %s", OUT_DIR)
        return 0
    log.info("run=%s | model=%s | obs=%g ms Gaussian | bin=%g ms",
             RUN_DIR.name, MODEL, OBS_GAUSS_SIGMA_S * 1e3,
             config.time_bin_width * 1e3)
    results = []
    for probe in PROBES:
        res = run_probe(probe, config)
        save_probe(res)
        write_index(res, config.time_bin_width)
        results.append(res)
    make_boxplot(results)
    make_mean_sem_plot(results)
    make_violin_plot(results)
    make_cluster_mean_swarm(results)
    save_combined(results)
    render_plots(results)
    log.info("done → %s", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
