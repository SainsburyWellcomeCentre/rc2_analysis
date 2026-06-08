"""Split-by-condition GLMs (2026-06-05): one GLM per cluster per condition.

Motivation
----------
The production ``glm/current/`` fit pools all three motion-clouds
conditions into a single design matrix — the Speed-vs-TF discrimination
is carried by *cross-condition* variance (Speed≡0 in V, TF≡0 in
T_Vstatic, TF=gain·Speed in VT). This driver does the complementary
thing: fit each cluster **separately within each condition** and compare
the resulting tuning curves across conditions. That reframes the two
live questions

  - does TF tuning change from V → VT?      (V-model TF  vs VT-model TF)
  - does Speed tuning change from T → VT?    (T-model Speed vs VT-model Speed)

as a direct cross-model comparison, instead of relying on the pooled
fit's within-design variance.

Configuration (per Laura, 2026-06-05)
-------------------------------------
  * No motion energy  → all 4 probes participate equally (244's camera is
    a placeholder; dropping ME removes the only probe asymmetry).
  * No spike history.
  * 100 ms bins.
  * Ridge λ = 1.0 (the 2026-05-07 retuned default).
  * Selection CV = condition-stratified 10-seed / 7-of-10 admission
    (the 2026-05-08 production rule: a candidate is admitted iff its
    Δ-cv_bps clears threshold in ≥7 of 10 random fold partitions). This
    is what makes the binary selection counts seed-robust. Speed-profile
    2-fold CV is added as a post-hoc per-cluster *diagnostic*
    (profile_cv_diagnostic=True): its folds are fixed to profile_id so
    they carry no seed variability and can't drive 10-seed admission, but
    the profile_cv_bps_* columns still report whether the selected model
    generalises train-on-profile-1 / test-on-profile-2.
    (The earlier 2026-06-05 variant used speed-profile folds AS the
    selection CV — a leak-test by construction, but single-seed because
    those folds are deterministic; kept at current_splitted_by_condition_
    speedprofile_singleseed/.)
  * Always forward-select, always from the intercept+onset baseline.
  * Per-condition candidate sets (the non-degenerate regressors only):
        V         : TF, SF, OR          (Speed ≡ 0)
        T_Vstatic : Speed               (TF ≡ 0; SF/OR undefined)
        VT        : Speed, TF, SF, OR
  * Cohort = the 88 GLM-cohort clusters, reused as-is across all three
    conditions (pulled per-probe from glm/current/_runs/<probe>/
    glm_model_comparison.csv; prefilter disabled so the set is fixed).

Outputs
-------
  ``figures/glm/current_splitted_by_condition/_runs/<probe>/<condition>/``
      per-probe-per-condition pipeline output (glm_model_comparison.csv,
      selection_history, coefficients, per-cluster overview/tuning/kernel
      pdfs, diagnostics/tuning_curves.csv, run.log).
  ``.../glm_model_comparison_<condition>.csv``
      4-probe-concatenated model comparison, one per condition.
  ``.../figs/forward_selection_summary_<condition>.pdf``
      per-condition forward-selection summary across all probes.
  ``.../figs/cross_condition/<probe>_cluster_<id>.pdf``
      the headline comparison: TF gain (V vs VT) and Speed gain
      (T_Vstatic vs VT) overlaid per cluster, from each condition's
      Additive model.

Usage::

    python -m scripts.run_glm_split_by_condition            # run everything
    python -m scripts.run_glm_split_by_condition --dry-run  # trial-count check only

Run with the rc2_analysis environment's interpreter. HOME-relative paths
throughout (works on both metonymy and synecdoche).
"""
from __future__ import annotations

import argparse
import logging
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import run_pipeline

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("glm_split_by_condition")

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds"
FORMATTED_DIR = ROOT / "formatted_data"
# Cohort source: the preserved pre-fix production run (Laura moved the old
# current/ etc. into no_speed_cv/ when the profile_id bug was found). The 88
# cohort cluster ids are stable, so we reuse them for the corrected re-runs.
CURRENT_RUNS = ROOT / "figures" / "glm" / "no_speed_cv" / "current" / "_runs"
OUT_ROOT = ROOT / "figures" / "glm" / "current_splitted_by_condition"
POOLED_OUT = ROOT / "figures" / "glm" / "current"   # all conditions pooled
MC_SEQUENCE = ROOT / "motion_cloud_sequence_250414.mat"
MC_FOLDERS = ROOT / "image_folders.mat"

PROBES = (
    "CAA-1123243_rec1",
    "CAA-1123244_rec1",
    "CAA-1123466_rec1",
    "CAA-1123467_rec1",
)
CONDITIONS = ("V", "T_Vstatic", "VT")

# Non-degenerate candidate regressors per condition (see module docstring).
CONDITION_CANDIDATES: dict[str, dict[str, tuple[str, ...]]] = {
    "V": {
        "main_effects": ("TF", "SF", "OR"),
        "interactions": ("TF_x_SF", "TF_x_OR", "SF_x_OR"),
    },
    "T_Vstatic": {
        "main_effects": ("Speed",),
        "interactions": (),
    },
    "VT": {
        "main_effects": ("Speed", "TF", "SF", "OR"),
        "interactions": (
            "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
            "TF_x_SF", "TF_x_OR", "SF_x_OR",
        ),
    },
}


# Pooled (all-conditions-together) candidate set — the full design. Used for
# the current/ run, identical methodology to the per-condition runs but with
# every condition's trials fit jointly (fit_condition=None).
POOLED_CANDIDATES: dict[str, tuple[str, ...]] = {
    "main_effects": ("Speed", "TF", "SF", "OR"),
    "interactions": (
        "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
        "TF_x_SF", "TF_x_OR", "SF_x_OR",
    ),
}


# Set by main() from --cv. Drives fold strategy + seed count for make_config.
#   "condition-stratified": random k-fold over (trial_id, condition); 10-seed
#       / 7-of-10 admission (seed-robust counts; both trajectories in every fold).
#   "speed-profile": 2-fold train-trajectory-1 / test-trajectory-2 leak-test
#       (profile_id velocity-derived); deterministic so single-seed.
# profile_cv_diagnostic is always on → each run renders the speed-profile vs
# condition-stratified cv-bps comparison figure (speed_profile_cv_comparison.pdf).
CV_STRATEGY = "condition-stratified"


def make_config(condition: str | None) -> GLMConfig:
    """GLMConfig for one condition (or pooled when None), under CV_STRATEGY."""
    cand = POOLED_CANDIDATES if condition is None else CONDITION_CANDIDATES[condition]
    speed_profile = CV_STRATEGY == "speed-profile"
    return replace(
        GLMConfig(),
        fit_condition=condition,
        main_effects=cand["main_effects"],
        interactions=cand["interactions"],
        include_me_face=False,       # no motion energy → all 4 probes equal
        include_history=False,
        include_onset_kernel=True,   # intercept+onset is the Null baseline
        lambda_ridge=1.0,
        cv_strategy=CV_STRATEGY,
        n_selection_seeds=1 if speed_profile else 10,
        selection_threshold_count=1 if speed_profile else 7,
        profile_cv_diagnostic=True,
        apply_prefilter=False,        # cohort is fixed to the reused 88
    )


def load_cohort(probe: str) -> set[int]:
    """The 88-cohort cluster ids for one probe, from glm/current/."""
    csv = CURRENT_RUNS / probe / "glm_model_comparison.csv"
    if not csv.exists():
        raise FileNotFoundError(
            f"cohort source missing: {csv}. Run glm/current/ first or fix the path."
        )
    df = pd.read_csv(csv)
    return {int(c) for c in df["cluster_id"].dropna().unique()}


# ----------------------------------------------------------------------
# Dry-run sanity check: per-condition × per-profile trial counts
# ----------------------------------------------------------------------
def dry_run() -> int:
    """Report trials per (condition, profile_id) for each probe — no fitting.

    Speed-profile 2-fold CV needs both profiles present within each
    condition subset. This loads only trial metadata, so it's cheap and
    safe to run before committing to the full 12-run fit.
    """
    rows: list[dict] = []
    for probe in PROBES:
        mat = FORMATTED_DIR / f"{probe}.mat"
        lookup = StimulusLookup(str(MC_SEQUENCE), str(MC_FOLDERS))
        data = load_probe_data(mat, config=GLMConfig(), stimulus_lookup=lookup,
                               visp_only=True)
        for t in data.trials:
            rows.append({
                "probe_id": probe,
                "condition": t.condition,
                "profile_id": int(t.profile_id),
            })
    df = pd.DataFrame(rows)
    counts = (
        df[df["condition"].isin(CONDITIONS)]
        .groupby(["probe_id", "condition", "profile_id"])
        .size()
        .unstack("profile_id", fill_value=0)
    )
    print("\n=== trials per (probe, condition) × profile_id ===")
    print(counts.to_string())
    print(
        "\nspeed-profile CV needs both profile columns > 0 in every row.\n"
        "If any condition has a single profile, its 2-fold CV-bps is undefined."
    )
    return 0


# ----------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------
def run_probe_condition(probe: str, condition: str, cohort: set[int]) -> Path:
    out_dir = OUT_ROOT / "_runs" / probe / condition
    cmp_path = out_dir / "glm_model_comparison.csv"
    if cmp_path.exists():
        log.info("skip %s/%s — already at %s", probe, condition, cmp_path)
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    mat = FORMATTED_DIR / f"{probe}.mat"
    lookup = StimulusLookup(str(MC_SEQUENCE), str(MC_FOLDERS))
    log.info("fitting %s / %s (%d clusters)", probe, condition, len(cohort))
    run_pipeline(
        mat_path=mat,
        config=make_config(condition),
        output_dir=out_dir,
        stimulus_lookup=lookup,
        backend="irls",
        visp_only=True,
        make_plots=True,
        plot_format="pdf",
        n_jobs=4,
        cluster_filter=cohort,
    )
    return out_dir


def run_probe_pooled(probe: str, cohort: set[int]) -> Path:
    """Pooled (all-conditions-together) fit for one probe → current/_runs/<probe>."""
    out_dir = POOLED_OUT / "_runs" / probe
    cmp_path = out_dir / "glm_model_comparison.csv"
    if cmp_path.exists():
        log.info("skip pooled %s — already at %s", probe, cmp_path)
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    lookup = StimulusLookup(str(MC_SEQUENCE), str(MC_FOLDERS))
    log.info("fitting pooled %s (%d clusters)", probe, len(cohort))
    run_pipeline(
        mat_path=FORMATTED_DIR / f"{probe}.mat",
        config=make_config(None),
        output_dir=out_dir,
        stimulus_lookup=lookup,
        backend="irls",
        visp_only=True,
        make_plots=True,
        plot_format="pdf",
        n_jobs=4,
        cluster_filter=cohort,
    )
    return out_dir


# ----------------------------------------------------------------------
# Aggregation
# ----------------------------------------------------------------------
def aggregate_pooled() -> pd.DataFrame | None:
    """Concatenate per-probe pooled model_comparison → current/glm_model_comparison.csv."""
    frames = []
    for probe in PROBES:
        csv = POOLED_OUT / "_runs" / probe / "glm_model_comparison.csv"
        if not csv.exists():
            log.warning("missing %s — skipping in pooled aggregation", csv)
            continue
        df = pd.read_csv(csv)
        df["probe_id"] = probe
        frames.append(df)
    if not frames:
        return None
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(POOLED_OUT / "glm_model_comparison.csv", index=False)
    log.info("wrote %s (%d clusters)", POOLED_OUT / "glm_model_comparison.csv", len(out))
    return out


def aggregate_condition(condition: str) -> pd.DataFrame | None:
    """Concatenate per-probe model_comparison for one condition."""
    frames = []
    for probe in PROBES:
        csv = OUT_ROOT / "_runs" / probe / condition / "glm_model_comparison.csv"
        if not csv.exists():
            log.warning("missing %s — skipping in %s aggregation", csv, condition)
            continue
        df = pd.read_csv(csv)
        df["probe_id"] = probe
        frames.append(df)
    if not frames:
        return None
    out = pd.concat(frames, ignore_index=True)
    out_csv = OUT_ROOT / f"glm_model_comparison_{condition}.csv"
    out.to_csv(out_csv, index=False)
    log.info("wrote %s (%d clusters)", out_csv, len(out))
    return out


def render_forward_selection_summary(df: pd.DataFrame, condition: str) -> None:
    from rc2_glm.plots import plot_forward_selection_summary, save_figure
    fig = plot_forward_selection_summary(df)
    fig.suptitle(
        f"Forward selection — {condition} only "
        f"(4 probes, n={df['cluster_id'].count()} clusters)",
        fontsize=12, fontweight="bold",
    )
    figdir = OUT_ROOT / "figs"
    figdir.mkdir(parents=True, exist_ok=True)
    for p in save_figure(fig, figdir / f"forward_selection_summary_{condition}",
                         fmt="pdf"):
        log.info("wrote %s", p)
    plt.close(fig)


# ----------------------------------------------------------------------
# Cross-condition tuning overlay (the headline comparison)
# ----------------------------------------------------------------------
def load_tuning(probe: str, condition: str) -> pd.DataFrame | None:
    csv = OUT_ROOT / "_runs" / probe / condition / "diagnostics" / "tuning_curves.csv"
    if not csv.exists():
        return None
    return pd.read_csv(csv)


def render_cross_condition_overlays(model: str = "Additive") -> None:
    """Per-cluster overlay: TF (V vs VT) and Speed (T_Vstatic vs VT).

    Two curves per condition, colour-matched (blue=cond_a, red=cond_b):
      * SOLID  — the condition's GLM marginal gain exp(B·β) (left axis,
        basis-rotation invariant model tuning).
      * DASHED — the OBSERVED median tuning curve from the MATLAB cache
        (median across trials of the per-trial-per-bin FR, right Hz axis)
        for the same variable+condition. This is the data, not a model.

    TF is compared V vs VT; Speed T_Vstatic vs VT.
    """
    from rc2_glm.pipeline import _observed_median_tuning
    from rc2_glm.precomputed_bins import load_precomputed_bin_edges

    tuning = {
        (probe, cond): load_tuning(probe, cond)
        for probe in PROBES for cond in CONDITIONS
    }
    figdir = OUT_ROOT / "figs" / "cross_condition"
    figdir.mkdir(parents=True, exist_ok=True)

    overlays = [
        ("TF", "V", "VT"),          # does TF tuning change V → VT?
        ("Speed", "T_Vstatic", "VT"),  # does Speed tuning change T → VT?
    ]
    n_made = 0
    for probe in PROBES:
        pre = load_precomputed_bin_edges(FORMATTED_DIR / f"{probe}.mat")
        cohort = sorted(load_cohort(probe))
        for cid in cohort:
            fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
            drew_any = False
            for ax, (var, cond_a, cond_b) in zip(axes, overlays):
                ax2 = ax.twinx()
                obs_drawn = False
                for cond, colour in ((cond_a, "tab:blue"), (cond_b, "tab:red")):
                    # SOLID — model marginal gain
                    tdf = tuning.get((probe, cond))
                    if tdf is not None:
                        sub = tdf[
                            (tdf["cluster_id"] == cid)
                            & (tdf["model"] == model)
                            & (tdf["variable"] == var)
                        ].sort_values("grid_x")
                        if not sub.empty:
                            ax.plot(sub["grid_x"], sub["gain"], color=colour,
                                    lw=2, label=f"{cond} model gain")
                            drew_any = True
                    # DASHED median + shaded Q1-Q3 — observed tuning (empirical, Hz)
                    obs = _observed_median_tuning(pre, cond, cid)
                    if obs and var in obs:
                        ox, omed, oq1, oq3 = obs[var]
                        band = np.isfinite(oq1) & np.isfinite(oq3)
                        ax2.fill_between(ox, oq1, oq3, where=band, color=colour,
                                         alpha=0.3, linewidth=0)
                        ax2.plot(ox, omed, color=colour, lw=1.4, ls="--",
                                 label=f"{cond} obs FR")
                        obs_drawn = True
                        drew_any = True
                ax.set_title(f"{var}: {cond_a} vs {cond_b}")
                ax.set_xlabel(var)
                ax.set_ylabel("model gain  exp(B·β)")
                ax.axhline(1.0, color="0.7", lw=0.8, ls=":")
                if obs_drawn:
                    ax2.set_ylabel("obs median FR (Hz)", fontsize=8)
                else:
                    ax2.set_yticks([])
                h1, l1 = ax.get_legend_handles_labels()
                h2, l2 = ax2.get_legend_handles_labels()
                if h1 or h2:
                    ax.legend(h1 + h2, l1 + l2, fontsize=6, loc="best")
            if not drew_any:
                plt.close(fig)
                continue
            fig.suptitle(
                f"{probe}  cluster {cid}  "
                f"(solid: {model} model gain · dashed: observed median FR)",
                fontsize=10, fontweight="bold",
            )
            fig.tight_layout()
            fig.savefig(figdir / f"{probe}_cluster_{cid}.pdf")
            plt.close(fig)
            n_made += 1
    log.info("wrote %d cross-condition overlay figures to %s", n_made, figdir)


def _mean_norm(y: np.ndarray) -> np.ndarray:
    """Divide a curve by its mean over bins (gain-cancelling normalisation).

    A pure multiplicative gain (curve2 = a*curve1) cancels exactly, so two
    mean-normalised curves coincide iff they differ only by gain; any
    residual divergence is a genuine shape (incl. modulation-depth) change.
    """
    m = np.nanmean(y)
    return y / m if (np.isfinite(m) and m != 0) else y


def render_cross_condition_normalized(model: str = "Additive") -> None:
    """Mean-normalised single-axis version of the cross-condition overlays.

    Every curve (model gain + observed median, and the observed Q1-Q3 band)
    is divided by its own mean, so all sit on one dimensionless axis around
    1.0. Curves that overlap differ only by a multiplicative gain; curves
    that diverge have a different SHAPE. Lets you read gain-vs-shape
    directly, which the dual-axis (Hz vs gain) version can't.
    """
    from rc2_glm.pipeline import _observed_median_tuning
    from rc2_glm.precomputed_bins import load_precomputed_bin_edges

    tuning = {
        (probe, cond): load_tuning(probe, cond)
        for probe in PROBES for cond in CONDITIONS
    }
    figdir = OUT_ROOT / "figs" / "cross_condition_normalized"
    figdir.mkdir(parents=True, exist_ok=True)
    overlays = [("TF", "V", "VT"), ("Speed", "T_Vstatic", "VT")]
    n_made = 0
    for probe in PROBES:
        pre = load_precomputed_bin_edges(FORMATTED_DIR / f"{probe}.mat")
        for cid in sorted(load_cohort(probe)):
            fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
            drew_any = False
            for ax, (var, cond_a, cond_b) in zip(axes, overlays):
                for cond, colour in ((cond_a, "tab:blue"), (cond_b, "tab:red")):
                    tdf = tuning.get((probe, cond))
                    if tdf is not None:
                        sub = tdf[
                            (tdf["cluster_id"] == cid)
                            & (tdf["model"] == model)
                            & (tdf["variable"] == var)
                        ].sort_values("grid_x")
                        if not sub.empty:
                            ax.plot(sub["grid_x"], _mean_norm(sub["gain"].to_numpy()),
                                    color=colour, lw=2, label=f"{cond} model")
                            drew_any = True
                    obs = _observed_median_tuning(pre, cond, cid)
                    if obs and var in obs:
                        ox, omed, oq1, oq3 = obs[var]
                        divisor = np.nanmean(omed)
                        if np.isfinite(divisor) and divisor != 0:
                            band = np.isfinite(oq1) & np.isfinite(oq3)
                            ax.fill_between(ox, oq1 / divisor, oq3 / divisor,
                                            where=band, color=colour, alpha=0.3,
                                            linewidth=0)
                            ax.plot(ox, omed / divisor, color=colour, lw=1.4,
                                    ls="--", label=f"{cond} obs")
                            drew_any = True
                ax.set_title(f"{var}: {cond_a} vs {cond_b}")
                ax.set_xlabel(var)
                ax.set_ylabel("normalised tuning (curve / its mean)")
                ax.axhline(1.0, color="0.7", lw=0.8, ls=":")
                h, l = ax.get_legend_handles_labels()
                if h:
                    ax.legend(fontsize=6, loc="best")
            if not drew_any:
                plt.close(fig)
                continue
            fig.suptitle(
                f"{probe}  cluster {cid}  — mean-normalised "
                f"(overlap = gain-only · divergence = shape change)",
                fontsize=10, fontweight="bold",
            )
            fig.tight_layout()
            fig.savefig(figdir / f"{probe}_cluster_{cid}.pdf")
            plt.close(fig)
            n_made += 1
    log.info("wrote %d normalised cross-condition figures to %s", n_made, figdir)


# ----------------------------------------------------------------------
def _run_pooled() -> None:
    """current/ — all conditions pooled, under CV_STRATEGY."""
    for probe in PROBES:
        run_probe_pooled(probe, load_cohort(probe))
    pooled = aggregate_pooled()
    if pooled is not None and not pooled.empty:
        from rc2_glm.plots import plot_forward_selection_summary, save_figure
        fig = plot_forward_selection_summary(pooled)
        fig.suptitle(
            f"Forward selection — all conditions pooled "
            f"(4 probes, n={pooled['cluster_id'].count()} clusters)",
            fontsize=12, fontweight="bold",
        )
        (POOLED_OUT / "figs").mkdir(parents=True, exist_ok=True)
        for p in save_figure(fig, POOLED_OUT / "figs" / "forward_selection_summary",
                             fmt="pdf"):
            log.info("wrote %s", p)
        plt.close(fig)
    log.info("done — pooled (cv=%s) under %s", CV_STRATEGY, POOLED_OUT)


def _run_split() -> None:
    """current_splitted_by_condition/ — per-condition, under CV_STRATEGY."""
    for probe in PROBES:
        cohort = load_cohort(probe)
        for cond in CONDITIONS:
            run_probe_condition(probe, cond, cohort)
    for cond in CONDITIONS:
        agg = aggregate_condition(cond)
        if agg is not None and not agg.empty:
            render_forward_selection_summary(agg, cond)
    render_cross_condition_overlays()
    render_cross_condition_normalized()
    log.info("done — split (cv=%s) under %s", CV_STRATEGY, OUT_ROOT)


def main() -> int:
    global CV_STRATEGY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scope", choices=("pooled", "split"), default="split",
        help="pooled = current/ (all conditions together); "
             "split = current_splitted_by_condition/ (per-condition).",
    )
    parser.add_argument(
        "--cv", choices=("condition-stratified", "speed-profile"),
        default="condition-stratified",
        help="Selection CV fold strategy. condition-stratified → 10-seed/7-of-10 "
             "random k-fold; speed-profile → 2-fold train-trajectory-1/test-2 "
             "leak-test (single-seed). profile_cv_diagnostic is on either way.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Only report per-condition × per-profile trial counts, then exit.",
    )
    parser.add_argument(
        "--overlays-only", action="store_true",
        help="Skip fitting; just re-render the split cross-condition overlays.",
    )
    args = parser.parse_args()
    CV_STRATEGY = args.cv

    if args.dry_run:
        return dry_run()
    if args.overlays_only:
        render_cross_condition_overlays()
        render_cross_condition_normalized()
        log.info("done — overlays under %s", OUT_ROOT / "figs")
        return 0

    log.info("scope=%s | cv=%s", args.scope, args.cv)
    if args.scope == "pooled":
        _run_pooled()
    else:
        _run_split()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
