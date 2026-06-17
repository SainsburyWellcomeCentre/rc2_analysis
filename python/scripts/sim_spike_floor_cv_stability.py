"""Spike-count floor re-check under 10-fold signed-rank CV.

Question (2026-06-16): the spike-count quality floor (≥50 total spikes,
≥50% trial occupancy) was calibrated when the forward-selection gate was a
fixed Δ-bps threshold on 5-fold (and 2-fold speed-profile) CV. The goggles
rf_sfor run now gates with the **Hardcastle one-sided Wilcoxon signed-rank
test on PER-FOLD bits/spike over 10-fold condition-stratified CV**. Is 50 still
the right floor?

Two things the floor protects, which fold count touches differently:
  - The POOLED cv-bps you report: denominator = total spikes, ~fold-independent
    (10 folds give slightly more training data per fold → marginally tighter).
  - The DECISION: the signed-rank test reads per-fold bps, whose denominator is
    ~total/n_folds. With 10 folds a 50-spike cell has ~5 spikes/fold → per-fold
    bps is Poisson-noise dominated. So the gate is MORE sensitive to low spikes.

This sim drives the REAL ``cross_validate_glm`` + the REAL signed-rank helper
(``forward_selection._signed_rank_greater``) on synthetic Poisson-GLM cells,
sweeping total spike count, for 5- vs 10-fold, measuring:
  1. SD of the pooled Δ cv-bps across reps (reported-number noise);
  2. mean within-cell per-fold bps SD (the signed-rank's input noise);
  3. false-admit rate — a NULL candidate (β=0): how often signed-rank wrongly
     admits at α=0.05 (Type I; should sit near α if calibrated);
  4. power — a candidate with a fixed modest true effect: admit rate.

Outputs a long CSV + a 3-panel figure to
``<figures>/glm/exploration/spike_floor_cv_stability.{csv,pdf,png}``. Reuses the
production CV machinery; runs nothing on the pipeline. Thread-pinned BLAS for
reproducible cv-bps (rc2_glm IRLS is BLAS-thread-nondeterministic).

Usage:
    python scripts/sim_spike_floor_cv_stability.py [--reps 200] [--quick]
"""

from __future__ import annotations

# BLAS thread pinning MUST precede numpy import (cv-bps reproducibility).
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.forward_selection import _signed_rank_greater

# --- Fixed design constants (kept close to the goggles rf_sfor run) ---------
BIN_WIDTH = 0.02                 # 20 ms bins
N_TRIALS = 36                    # 18 per condition → clean condition-stratified folds
N_BINS_PER_TRIAL = 30
N_BASES = 5                      # candidate covariate basis columns (Speed-like)
LAMBDA_RIDGE = 1.0               # production ridge on non-intercept columns
ALPHA = 0.05                     # signed-rank admission level
EFFECT_BETA = 0.55               # modest true effect magnitude (β on basis col 0)
SPIKE_TARGETS = (15, 30, 50, 75, 100, 200, 500)
FOLD_COUNTS = (5, 10)


def _covariate_basis(n: int, rng: np.random.Generator) -> np.ndarray:
    """A smooth Speed-like raised-cosine-ish basis on a random covariate."""
    x = rng.uniform(0.0, 1.0, n)
    cols = [np.cos(2.0 * np.pi * x * (k + 1) + 0.3 * k) for k in range(N_BASES)]
    return np.stack(cols, axis=1)


def _penalty(n_cols_with_intercept: int) -> np.ndarray:
    """Ridge on every column except the intercept (column 0)."""
    d = np.full(n_cols_with_intercept, LAMBDA_RIDGE)
    d[0] = 0.0
    return np.diag(d)


def _simulate_cell(target_spikes: int, effect: bool, rng: np.random.Generator):
    """One synthetic Poisson-GLM cell with ~target_spikes total spikes."""
    n = N_TRIALS * N_BINS_PER_TRIAL
    trial_ids = np.repeat(np.arange(N_TRIALS), N_BINS_PER_TRIAL)
    cond = np.array(["A" if t % 2 == 0 else "B" for t in range(N_TRIALS)])[trial_ids]
    Xb = _covariate_basis(n, rng)
    offset = float(np.log(BIN_WIDTH))

    beta = np.zeros(N_BASES)
    if effect:
        beta[0] = EFFECT_BETA
    eta_cov = Xb @ beta
    # Calibrate the intercept so E[total spikes] = target (mean over the cov term).
    intercept = float(np.log(target_spikes) - np.log(np.sum(np.exp(offset + eta_cov))))
    mu = np.exp(intercept + offset + eta_cov)
    y = rng.poisson(mu).astype(float)

    X_null = np.ones((n, 1))
    X_cand = np.column_stack([np.ones(n), Xb])
    return dict(y=y, offset=offset, X_null=X_null, X_cand=X_cand,
                trial_ids=trial_ids, cond=cond, total=int(y.sum()))


def _one_rep(target_spikes: int, n_folds: int, effect: bool,
             rng: np.random.Generator) -> dict | None:
    cell = _simulate_cell(target_spikes, effect, rng)
    folds = make_trial_folds(
        cell["trial_ids"], n_folds=n_folds, seed=int(rng.integers(1 << 30)),
        condition_labels_per_bin=cell["cond"],
    )
    cv_null = cross_validate_glm(
        cell["X_null"], cell["y"], cell["offset"], folds,
        lambda_ridge=0.0, backend="irls",
    )
    cv_cand = cross_validate_glm(
        cell["X_cand"], cell["y"], cell["offset"], folds,
        lambda_ridge=LAMBDA_RIDGE, backend="irls",
        penalty_matrix=_penalty(cell["X_cand"].shape[1]),
    )
    per_fold_delta = cv_cand.fold_bits_per_spike - cv_null.fold_bits_per_spike
    valid = per_fold_delta[np.isfinite(per_fold_delta)]
    if valid.size < 2:
        return None
    pval = _signed_rank_greater(per_fold_delta)
    return dict(
        total=cell["total"],
        pooled_delta=cv_cand.cv_bits_per_spike - cv_null.cv_bits_per_spike,
        cand_fold_bps_sd=float(np.nanstd(cv_cand.fold_bits_per_spike)),
        pval=pval,
        admit=bool(pval < ALPHA),
    )


def run(reps: int) -> pd.DataFrame:
    rows = []
    master = np.random.default_rng(20260616)
    for effect in (False, True):
        for n_folds in FOLD_COUNTS:
            for target in SPIKE_TARGETS:
                recs = []
                for _ in range(reps):
                    r = _one_rep(target, n_folds, effect,
                                 np.random.default_rng(master.integers(1 << 31)))
                    if r is not None:
                        recs.append(r)
                if not recs:
                    continue
                d = pd.DataFrame(recs)
                rows.append(dict(
                    effect="effect" if effect else "null",
                    n_folds=n_folds,
                    target_spikes=target,
                    mean_total_spikes=float(d["total"].mean()),
                    n_reps=len(d),
                    mean_pooled_delta=float(d["pooled_delta"].mean()),
                    sd_pooled_delta=float(d["pooled_delta"].std()),
                    mean_perfold_bps_sd=float(d["cand_fold_bps_sd"].mean()),
                    admit_rate=float(d["admit"].mean()),
                ))
                tag = "power " if effect else "FPR   "
                print(f"[{tag}] folds={n_folds:2d} target={target:4d} "
                      f"(~{rows[-1]['mean_total_spikes']:6.0f}) "
                      f"admit={rows[-1]['admit_rate']:.3f} "
                      f"perfoldSD={rows[-1]['mean_perfold_bps_sd']:.3f} "
                      f"poolΔSD={rows[-1]['sd_pooled_delta']:.3f}")
    return pd.DataFrame(rows)


def make_figure(df: pd.DataFrame, out_stem: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    colors = {5: "#1f77b4", 10: "#d62728"}

    # Panel 1: estimation noise (null cells).
    ax = axes[0]
    nd = df[df["effect"] == "null"]
    for nf in FOLD_COUNTS:
        s = nd[nd["n_folds"] == nf].sort_values("target_spikes")
        ax.plot(s["target_spikes"], s["mean_perfold_bps_sd"], "-o",
                color=colors[nf], label=f"{nf}-fold per-fold bps SD")
        ax.plot(s["target_spikes"], s["sd_pooled_delta"], "--s",
                color=colors[nf], alpha=0.6, label=f"{nf}-fold pooled Δ SD")
    ax.axvline(50, color="k", ls=":", lw=1, label="floor = 50")
    ax.set_xscale("log"); ax.set_xlabel("total spikes"); ax.set_ylabel("bits/spike SD")
    ax.set_title("Estimation noise vs spikes"); ax.legend(fontsize=6)

    # Panel 2: false-admit (Type I) rate, null cells.
    ax = axes[1]
    for nf in FOLD_COUNTS:
        s = nd[nd["n_folds"] == nf].sort_values("target_spikes")
        ax.plot(s["target_spikes"], s["admit_rate"], "-o", color=colors[nf],
                label=f"{nf}-fold")
    ax.axhline(ALPHA, color="k", ls="--", lw=1, label=f"α={ALPHA}")
    ax.axvline(50, color="k", ls=":", lw=1)
    ax.set_xscale("log"); ax.set_xlabel("total spikes")
    ax.set_ylabel("false-admit rate"); ax.set_ylim(0, max(0.12, ax.get_ylim()[1]))
    ax.set_title("Type I (null candidate)"); ax.legend(fontsize=7)

    # Panel 3: power, effect cells.
    ax = axes[2]
    ed = df[df["effect"] == "effect"]
    for nf in FOLD_COUNTS:
        s = ed[ed["n_folds"] == nf].sort_values("target_spikes")
        ax.plot(s["target_spikes"], s["admit_rate"], "-o", color=colors[nf],
                label=f"{nf}-fold")
    ax.axvline(50, color="k", ls=":", lw=1, label="floor = 50")
    ax.set_xscale("log"); ax.set_xlabel("total spikes"); ax.set_ylabel("power (admit rate)")
    ax.set_title(f"Power (true β={EFFECT_BETA})"); ax.legend(fontsize=7)

    fig.suptitle("Spike-count floor under signed-rank CV — 5 vs 10 folds", fontsize=12)
    for ext in ("pdf", "png"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--quick", action="store_true", help="reps=40 smoke run")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    reps = 40 if args.quick else args.reps

    out_dir = Path(args.out) if args.out else (
        Path.home() / "local_data" / "motion_clouds" / "figures" / "glm"
        / "exploration"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "spike_floor_cv_stability"

    print(f"running {reps} reps/config "
          f"({len(SPIKE_TARGETS)} spike levels × {len(FOLD_COUNTS)} fold-counts "
          f"× null/effect)…")
    df = run(reps)
    df.to_csv(stem.with_suffix(".csv"), index=False)
    make_figure(df, stem)
    print(f"\nwrote {stem}.csv / .pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
