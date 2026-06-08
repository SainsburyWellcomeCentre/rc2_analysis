"""current/-analog GLM fit for the goggles motion-clouds probes.

The goggles cohort (CAA-1124370/371, mouse VR goggles, 2026-03) is the new
experimental modality: each formatted file bundles rec1 (sparse noise) / rec2
(motion clouds) / rec3 (sparse noise); only rec2 is the motion-clouds block the
reader selects (``FormattedDataReader(..., session_of_interest='rec2')`` via
``load_probe_data``).

This driver reproduces the pooled ``glm/current/`` methodology on the goggles
data:

  * all three conditions (V / T_Vstatic / VT) fit jointly (``fit_condition=None``),
    full candidate design;
  * λ = 1.0, raised-cosine Speed/TF bases, intercept+onset baseline, no motion
    energy, no history;
  * forward selection condition-stratified 10-seed / 7-of-10 admission, with the
    speed-profile 2-fold CV kept as a per-cluster diagnostic
    (``profile_cv_diagnostic=True``; profile_id is velocity-derived in io, so it
    tracks the true trajectory split — see the 2026-06-05/08 CHANGELOG entries).

Differences from the screens driver (``run_glm_split_by_condition.py``):

  * goggles formatted dir + the two goggles probes;
  * the goggles stimulus parameter grid (``GOGGLES_STIMULUS``: SF 008/016/032,
    re-parameterized VX, gain ladder preserved — see ``trial_conditions.py``);
  * the goggles motion-cloud sequence / image-folders mats;
  * **cohort built by the prefilter** off the goggles stationary-vs-motion table
    (``apply_prefilter=True``), because there is no prior goggles run to reuse.
    Clusters absent from the precomputed SVM table fall out as non-significant,
    so the cohort universe is the SVM table's selected units.

Usage:
    python scripts/run_glm_goggles.py --dry-run      # trial/cohort sanity, no fit
    python scripts/run_glm_goggles.py                # fit both probes, aggregate
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

from rc2_formatted_data_reader import GOGGLES_STIMULUS, StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import run_pipeline

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("glm_goggles")

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds"
FORMATTED_DIR = ROOT / "formatted_data_goggles"
OUT_ROOT = ROOT / "figures" / "glm" / "current_goggles"
# Per-condition split run (one GLM per cluster per condition), reusing the
# current_goggles/ cohort. Mirrors screens current_splitted_by_condition/.
OUT_SPLIT = ROOT / "figures" / "glm" / "current_splitted_by_condition_goggles"
MC_SEQUENCE = ROOT / "motion_clouds_goggles_sequence_260420.mat"
MC_FOLDERS = ROOT / "image_folders_goggles.mat"

PROBES = (
    "CAA-1124370_rec1_rec2_rec3",
    "CAA-1124371_rec1_rec2_rec3",
)
CONDITIONS = ("V", "T_Vstatic", "VT")

# Non-degenerate candidate regressors per condition (stimulus-design driven, so
# modality-independent — identical to the screens split driver). Within a single
# condition the design is degenerate: V has Speed≡0, T_Vstatic has TF≡0 and
# SF/OR=NaN, so each condition gets a restricted candidate set.
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

# Pooled (all-conditions-together) candidate set — the full design, matching
# the screens current/ run (run_glm_split_by_condition.POOLED_CANDIDATES).
POOLED_CANDIDATES: dict[str, tuple[str, ...]] = {
    "main_effects": ("Speed", "TF", "SF", "OR"),
    "interactions": (
        "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
        "TF_x_SF", "TF_x_OR", "SF_x_OR",
    ),
}


def make_config() -> GLMConfig:
    """Pooled current/-analog config with the prefilter ON (builds the cohort)."""
    return replace(
        GLMConfig(),
        fit_condition=None,
        main_effects=POOLED_CANDIDATES["main_effects"],
        interactions=POOLED_CANDIDATES["interactions"],
        include_me_face=False,       # no motion energy → probes comparable
        include_history=False,
        include_onset_kernel=True,   # intercept+onset is the Null baseline
        lambda_ridge=1.0,
        cv_strategy="condition-stratified",
        n_selection_seeds=10,
        selection_threshold_count=7,
        profile_cv_diagnostic=True,
        apply_prefilter=True,        # cohort from the goggles SVM table
    )


def make_config_condition(condition: str) -> GLMConfig:
    """Per-condition config: restricted candidates, cohort reused (no prefilter)."""
    cand = CONDITION_CANDIDATES[condition]
    return replace(
        GLMConfig(),
        fit_condition=condition,
        main_effects=cand["main_effects"],
        interactions=cand["interactions"],
        include_me_face=False,
        include_history=False,
        include_onset_kernel=True,
        lambda_ridge=1.0,
        cv_strategy="condition-stratified",
        n_selection_seeds=10,
        selection_threshold_count=7,
        profile_cv_diagnostic=True,
        apply_prefilter=False,        # cohort fixed to the current_goggles/ run
    )


def load_cohort(probe: str) -> set[int]:
    """Cohort cluster ids for one probe, from the pooled current_goggles/ run."""
    csv = OUT_ROOT / "_runs" / probe / "glm_model_comparison.csv"
    if not csv.exists():
        raise FileNotFoundError(
            f"cohort source missing: {csv}. Run the pooled fit first "
            f"(python scripts/run_glm_goggles.py)."
        )
    df = pd.read_csv(csv)
    return {int(c) for c in df["cluster_id"].dropna().unique()}


def _lookup() -> StimulusLookup:
    return StimulusLookup(str(MC_SEQUENCE), str(MC_FOLDERS), GOGGLES_STIMULUS)


# ----------------------------------------------------------------------
# Dry-run sanity check: per-condition × per-profile trial counts + cohort size
# ----------------------------------------------------------------------
def dry_run() -> int:
    """Report trial counts and prefilter cohort size — no GLM fit."""
    from rc2_glm.prefilter import prefilter_probe

    rows: list[dict] = []
    for probe in PROBES:
        data = load_probe_data(
            FORMATTED_DIR / f"{probe}.mat",
            config=GLMConfig(),
            stimulus_lookup=_lookup(),
            visp_only=True,
        )
        for t in data.trials:
            rows.append({
                "probe_id": probe,
                "condition": t.condition,
                "profile_id": int(t.profile_id),
            })
        pf = prefilter_probe(data, config=GLMConfig())
        n_keep = int(pf["should_run_glm"].sum())
        log.info("%s: %d VISp clusters → cohort %d (prefilter)", probe, len(data.clusters), n_keep)

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
        "\ncondition-stratified CV is the gating selection; speed-profile 2-fold "
        "is the diagnostic and needs both profile columns > 0 per condition."
    )
    return 0


# ----------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------
def run_probe(probe: str) -> Path:
    out_dir = OUT_ROOT / "_runs" / probe
    cmp_path = out_dir / "glm_model_comparison.csv"
    if cmp_path.exists():
        log.info("skip %s — already at %s", probe, cmp_path)
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("fitting %s (prefilter cohort)", probe)
    run_pipeline(
        mat_path=FORMATTED_DIR / f"{probe}.mat",
        config=make_config(),
        output_dir=out_dir,
        stimulus_lookup=_lookup(),
        backend="irls",
        visp_only=True,
        make_plots=True,
        plot_format="pdf",
        n_jobs=4,
    )
    return out_dir


def aggregate() -> pd.DataFrame | None:
    """Concatenate per-probe model_comparison → current_goggles/glm_model_comparison.csv."""
    frames = []
    for probe in PROBES:
        csv = OUT_ROOT / "_runs" / probe / "glm_model_comparison.csv"
        if not csv.exists():
            log.warning("missing %s — skipping in aggregation", csv)
            continue
        df = pd.read_csv(csv)
        df["probe_id"] = probe
        frames.append(df)
    if not frames:
        return None
    out = pd.concat(frames, ignore_index=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_ROOT / "glm_model_comparison.csv", index=False)
    log.info("wrote %s (%d clusters)", OUT_ROOT / "glm_model_comparison.csv", len(out))
    render_forward_selection_summary(out)
    return out


def render_forward_selection_summary(pooled: pd.DataFrame) -> None:
    """Combined forward-selection summary across both goggles probes."""
    if pooled is None or pooled.empty:
        return
    from rc2_glm.plots import plot_forward_selection_summary, save_figure

    fig = plot_forward_selection_summary(pooled)
    n_probes = pooled["probe_id"].nunique()
    fig.suptitle(
        f"Forward selection — goggles, all conditions pooled "
        f"({n_probes} probes, n={pooled['cluster_id'].count()} clusters)",
        fontsize=12, fontweight="bold",
    )
    figdir = OUT_ROOT / "figs"
    figdir.mkdir(parents=True, exist_ok=True)
    for p in save_figure(fig, figdir / "forward_selection_summary", fmt="pdf"):
        log.info("wrote %s", p)


def run_probe_condition(probe: str, condition: str, cohort: set[int]) -> Path:
    out_dir = OUT_SPLIT / "_runs" / probe / condition
    cmp_path = out_dir / "glm_model_comparison.csv"
    if cmp_path.exists():
        log.info("skip %s/%s — already at %s", probe, condition, cmp_path)
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("fitting %s / %s (%d clusters)", probe, condition, len(cohort))
    run_pipeline(
        mat_path=FORMATTED_DIR / f"{probe}.mat",
        config=make_config_condition(condition),
        output_dir=out_dir,
        stimulus_lookup=_lookup(),
        backend="irls",
        visp_only=True,
        make_plots=True,
        plot_format="pdf",
        n_jobs=4,
        cluster_filter=cohort,
    )
    return out_dir


def aggregate_condition(condition: str) -> pd.DataFrame | None:
    """Concatenate per-probe model_comparison for one condition + summary plot."""
    frames = []
    for probe in PROBES:
        csv = OUT_SPLIT / "_runs" / probe / condition / "glm_model_comparison.csv"
        if not csv.exists():
            log.warning("missing %s — skipping in %s aggregation", csv, condition)
            continue
        df = pd.read_csv(csv)
        df["probe_id"] = probe
        frames.append(df)
    if not frames:
        return None
    out = pd.concat(frames, ignore_index=True)
    OUT_SPLIT.mkdir(parents=True, exist_ok=True)
    out_csv = OUT_SPLIT / f"glm_model_comparison_{condition}.csv"
    out.to_csv(out_csv, index=False)
    log.info("wrote %s (%d clusters)", out_csv, len(out))

    from rc2_glm.plots import plot_forward_selection_summary, save_figure

    fig = plot_forward_selection_summary(out)
    fig.suptitle(
        f"Forward selection — goggles, {condition} "
        f"({out['probe_id'].nunique()} probes, n={out['cluster_id'].count()} clusters)",
        fontsize=12, fontweight="bold",
    )
    figdir = OUT_SPLIT / "figs"
    figdir.mkdir(parents=True, exist_ok=True)
    for p in save_figure(fig, figdir / f"forward_selection_summary_{condition}", fmt="pdf"):
        log.info("wrote %s", p)
    return out


# ----------------------------------------------------------------------
# Cross-condition tuning overlays (the headline split-by-condition comparison)
# ----------------------------------------------------------------------
def load_tuning(probe: str, condition: str) -> pd.DataFrame | None:
    csv = OUT_SPLIT / "_runs" / probe / condition / "diagnostics" / "tuning_curves.csv"
    if not csv.exists():
        return None
    return pd.read_csv(csv)


def _mean_norm(y: np.ndarray) -> np.ndarray:
    """Divide a curve by its mean over bins (gain-cancelling normalisation)."""
    m = np.nanmean(y)
    return y / m if (np.isfinite(m) and m != 0) else y


# TF compared V vs VT; Speed compared T_Vstatic vs VT.
_OVERLAYS = (("TF", "V", "VT"), ("Speed", "T_Vstatic", "VT"))


def render_cross_condition_overlays(model: str = "Additive") -> None:
    """Per-cluster overlay: TF (V vs VT) and Speed (T_Vstatic vs VT).

    SOLID = the condition's GLM marginal gain exp(B·β) (left axis); DASHED =
    observed median FR from the goggles tuning cache (right Hz axis, with a
    shaded Q1-Q3 band). Mirrors the screens split driver.
    """
    from rc2_glm.pipeline import _observed_median_tuning
    from rc2_glm.precomputed_bins import load_precomputed_bin_edges

    tuning = {(p, c): load_tuning(p, c) for p in PROBES for c in CONDITIONS}
    figdir = OUT_SPLIT / "figs" / "cross_condition"
    figdir.mkdir(parents=True, exist_ok=True)
    n_made = 0
    for probe in PROBES:
        pre = load_precomputed_bin_edges(FORMATTED_DIR / f"{probe}.mat")
        for cid in sorted(load_cohort(probe)):
            fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
            drew_any = False
            for ax, (var, cond_a, cond_b) in zip(axes, _OVERLAYS):
                ax2 = ax.twinx()
                obs_drawn = False
                for cond, colour in ((cond_a, "tab:blue"), (cond_b, "tab:red")):
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


def render_cross_condition_normalized(model: str = "Additive") -> None:
    """Mean-normalised single-axis version (overlap = gain-only, divergence = shape)."""
    from rc2_glm.pipeline import _observed_median_tuning
    from rc2_glm.precomputed_bins import load_precomputed_bin_edges

    tuning = {(p, c): load_tuning(p, c) for p in PROBES for c in CONDITIONS}
    figdir = OUT_SPLIT / "figs" / "cross_condition_normalized"
    figdir.mkdir(parents=True, exist_ok=True)
    n_made = 0
    for probe in PROBES:
        pre = load_precomputed_bin_edges(FORMATTED_DIR / f"{probe}.mat")
        for cid in sorted(load_cohort(probe)):
            fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
            drew_any = False
            for ax, (var, cond_a, cond_b) in zip(axes, _OVERLAYS):
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


def run_split() -> None:
    """Per-condition fits (V / T_Vstatic / VT) for both probes + summaries."""
    for probe in PROBES:
        cohort = load_cohort(probe)
        for condition in CONDITIONS:
            run_probe_condition(probe, condition, cohort)
    for condition in CONDITIONS:
        aggregate_condition(condition)
    render_cross_condition_overlays()
    render_cross_condition_normalized()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report per-condition × profile trial counts + cohort size, then exit.",
    )
    parser.add_argument(
        "--probe", choices=PROBES, default=None,
        help="Fit a single probe (default: both).",
    )
    parser.add_argument(
        "--split", action="store_true",
        help="Per-condition split run (V / T_Vstatic / VT) → "
             "current_splitted_by_condition_goggles/, reusing the current_goggles/ "
             "cohort. Requires the pooled fit to have run first.",
    )
    parser.add_argument(
        "--overlays-only", action="store_true",
        help="Skip fitting; re-render the split cross-condition overlays only.",
    )
    args = parser.parse_args()

    if args.dry_run:
        return dry_run()

    if args.overlays_only:
        render_cross_condition_overlays()
        render_cross_condition_normalized()
        return 0

    if args.split:
        run_split()
        return 0

    probes = (args.probe,) if args.probe else PROBES
    for probe in probes:
        run_probe(probe)
    aggregate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
