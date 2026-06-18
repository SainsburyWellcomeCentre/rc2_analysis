"""Goggles GLM with RF-LOCAL SF/OR — 20 ms bins, Speed+TF+SF+OR+Acc.

The headline run for the receptive-field-based SF/OR work. Unlike
``run_glm_goggles.py`` (which feeds SF/OR as the categorical stimulus *tokens*,
100 ms bins), this driver sets ``sf_or_source="rf_local"``: SF and OR are the
per-cluster, per-bin *local* spatial-frequency (cpd) and orientation (deg) that
each cluster's receptive field sees on the motion cloud — the Gabor extraction
(``gabor_extract_gpu.py`` → cohort parquet) mapped onto the GLM bins via the
velocity-locked cloud-frame clock. SF becomes a continuous raised-cosine basis,
OR a circular basis (NOT dummies). Only clusters with a goggle RF are fit, so
the cohort is the prefilter cohort ∩ RF clusters (~61).

Config (per Laura's spec, 2026-06-12):
  * 20 ms bins; pooled all-conditions fit (V / T_Vstatic / VT together);
  * candidates Speed + TF + SF + OR + Acceleration (+ the standard stimulus
    pairwise interactions for the "as usual" forward-selection summary);
  * no face ME, no history;
  * λ = 1.0, intercept+onset Null, condition-stratified 10-seed / 7-of-10
    admission, speed-profile 2-fold kept as the per-cluster diagnostic.

Screens and the token goggles runs are untouched — this is the only driver that
flips ``sf_or_source``.

Usage:
    python scripts/run_glm_goggles_rf_sfor_20ms.py --dry-run   # cohort sizes
    python scripts/run_glm_goggles_rf_sfor_20ms.py --probe CAA-1124370_rec1_rec2_rec3 --max-clusters 1  # wiring smoke
    python scripts/run_glm_goggles_rf_sfor_20ms.py             # both probes
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import replace
from pathlib import Path

# Allow ``import scripts.<sibling>`` when this driver is run as a file path
# (``python scripts/run_...py``): then sys.path[0] is scripts/, not its parent,
# so the late in-pipeline ``import scripts.diagnostics_rf_sfor`` (main(), end of
# a full run) fails with "No module named 'scripts'" and the diagnostics step is
# silently skipped — variance_partition / acid figures never regenerate
# (histme_all run, 2026-06-16). Putting python/ (the parent of scripts/) on the
# path makes the package importable either way (file path or ``python -m``).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib

matplotlib.use("Agg")
import pandas as pd

from rc2_formatted_data_reader import GOGGLES_STIMULUS, StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import run_pipeline
from rc2_glm.rf_sf_or import load_rf_sf_or

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("glm_goggles_rf")

HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds"
FORMATTED_DIR = ROOT / "formatted_data_goggles"
OUT_ROOT = ROOT / "figures" / "glm" / "current_rf_sfor_20ms_goggles"
OUT_HISTME = ROOT / "figures" / "glm" / "current_rf_sfor_20ms_histme_goggles"
OUT_HISTME_ALL = ROOT / "figures" / "glm" / "current_rf_sfor_20ms_histme_goggles_all"
# Same all-clusters config as histme, but History flips to the BASELINE-nuisance
# mode (history_in_baseline=True): always-on in null + every fitted model, never a
# forward-selection candidate, never displayed, mean folded into the marginal's
# intercept. The "control for refractory structure without interpreting it" run.
OUT_HISTBASE_ALL = ROOT / "figures" / "glm" / "current_rf_sfor_20ms_histbase_goggles_all"
# Condition-split variants of the _all run (fit V or T_Vstatic trials alone) —
# isolates pure-visual (TF/SF/OR) vs pure-vestibular (Speed) tuning, away from
# the VT TF=gain·Speed entanglement. One folder per condition.
OUT_ALL_BY_COND = {
    "V": ROOT / "figures" / "glm" / "current_rf_sfor_20ms_histme_goggles_all_V",
    "T_Vstatic": ROOT / "figures" / "glm" / "current_rf_sfor_20ms_histme_goggles_all_Tvstatic",
}
# History-in-baseline variants of the condition splits (the histbase analogue of
# OUT_ALL_BY_COND): fit V or T_Vstatic trials alone, History as baseline nuisance.
OUT_HISTBASE_BY_COND = {
    "V": ROOT / "figures" / "glm" / "current_rf_sfor_20ms_histbase_goggles_all_V",
    "T_Vstatic": ROOT / "figures" / "glm" / "current_rf_sfor_20ms_histbase_goggles_all_Tvstatic",
}
MC_SEQUENCE = ROOT / "motion_clouds_goggles_sequence_260420.mat"
MC_FOLDERS = ROOT / "image_folders_goggles.mat"
RF_PARQUET_DIR = str(ROOT / "saved_goggles" / "_extract" / "cohort")

PROBES = (
    "CAA-1124370_rec1_rec2_rec3",
    "CAA-1124371_rec1_rec2_rec3",
)

# Speed+TF+SF+OR+Acceleration main effects; the stimulus pairwise interactions.
MAIN_EFFECTS = ("Speed", "TF", "SF", "OR", "Acceleration")
INTERACTIONS = (
    "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
    "TF_x_SF", "TF_x_OR", "SF_x_OR",
)
# Full pairwise interaction set on the histme value/visual predictors
# {Speed, TF, SF, OR, ME_face, Acceleration} = the 6 stimulus pairs above +
# ME_face_x_Speed + the Acceleration / ME_face full-pairwise set (2026-06-15):
# every condition's _all run carries all interactions on its candidates.
HISTME_INTERACTIONS = INTERACTIONS + (
    "ME_face_x_Speed",
    "Speed_x_Acceleration", "TF_x_ME_face", "TF_x_Acceleration",
    "SF_x_ME_face", "SF_x_Acceleration", "OR_x_ME_face",
    "OR_x_Acceleration", "ME_face_x_Acceleration",
)


def _lookup() -> StimulusLookup:
    return StimulusLookup(str(MC_SEQUENCE), str(MC_FOLDERS), GOGGLES_STIMULUS)


def make_config() -> GLMConfig:
    return replace(
        GLMConfig(),
        time_bin_width=0.02,
        sf_or_source="rf_local",
        rf_sf_or_parquet_dir=RF_PARQUET_DIR,
        # Goggles motion = the PHOTODIODE visual-stimulus window for V/VT (the
        # cloud displays ~0.7 s after the velocity command — that latency is the
        # real stationary↔motion gap, excluded; T_Vstatic falls back to
        # velocity). Screens keep the velocity default. See
        # reference_motion_clouds_goggles_trial_structure.
        motion_window_source="photodiode",
        fit_condition=None,
        main_effects=MAIN_EFFECTS,
        interactions=INTERACTIONS,
        include_acceleration=True,
        n_accel_bases=5,
        accel_range=(-3.0, 3.0),
        include_me_face=False,
        include_history=False,
        include_onset_kernel=True,
        lambda_ridge=1.0,
        # Selection gate = Hardcastle (2026-06-16). Forward selection admits a
        # variable iff a one-sided Wilcoxon signed-rank test on the PER-FOLD
        # paired Δ bits/spike (vs the running model) clears alpha, over a
        # 10-fold condition-stratified partition — the field-standard rule with
        # real per-fold power (10 folds; a fixed Δ-bps threshold has none). The
        # speed-profile trajectory split is demoted to the POST-HOC diagnostic
        # (profile_cv_diagnostic=True → speed_profile_cv_comparison in
        # diagnostics/), the "does Speed transfer across velocity trajectory"
        # supplementary, no longer the gate. signed_rank requires a single fold
        # partition (the folds are the test sample) → n_selection_seeds=1.
        cv_strategy="condition-stratified",
        n_folds=10,
        selection_rule="signed_rank",
        selection_alpha=0.05,
        n_selection_seeds=1,
        selection_threshold_count=1,
        profile_cv_diagnostic=True,
        apply_prefilter=False,  # whole selected cohort; prefilter is a diagnostic only
    )


def make_config_histme() -> GLMConfig:
    """make_config() + History (40 ms = 2 identity lag-bins) + face ME — the
    'complete-picture' run, and the test that the History intercept-fold and the
    ME / ME×Speed marginal reconstruction are all in place (both were latent/off
    in the main run, so a clean ZERO-FILLED log here proves the fixes)."""
    return replace(
        make_config(),
        include_history=True,
        history_window_s=0.04,          # 2 lag bins at 20 ms ("2 bins")
        history_basis_kind="identity",  # 2 dummy lags (lag-1, lag-2)
        n_history_bases=2,
        include_me_face=True,
        motion_energy_camera="camera0",  # face
        main_effects=MAIN_EFFECTS + ("ME_face",),
        interactions=HISTME_INTERACTIONS,
    )


def make_config_histme_all() -> GLMConfig:
    """make_config_histme() but INCLUDE clusters without a clean RF: give their
    SF/OR the per-cloud cohort NOMINAL (constant per trial) instead of dropping
    them (rf_sf_or_nominal_fallback=True). The full-cohort 'complete picture' —
    Speed/TF/Accel/ME/History are RF-independent so the no-RF clusters
    contribute to those; SF/OR is RF-local for RF clusters and the nominal
    stand-in for the rest (a HYBRID regressor — read SF/OR tuning on RF
    clusters only)."""
    return replace(make_config_histme(), rf_sf_or_nominal_fallback=True)


def make_config_histbase_all() -> GLMConfig:
    """make_config_histme_all() but History runs in BASELINE-nuisance mode
    (history_in_baseline=True): the spike-History term is always present in the
    null model and every fitted model (like the onset kernel), is NEVER a
    forward-selection candidate, is NEVER displayed, and its mean is folded into
    the intercept for the marginal. We control for refractory/autocorrelation
    structure without interpreting it. include_history is turned OFF (the legacy
    candidate mode) so History is baseline-only; the history window/basis/lags
    carry over from make_config_histme(). All other terms — rf_local SF/OR, ME,
    Speed/TF/SF/OR/Accel + the full interaction set — are unchanged."""
    return replace(
        make_config_histme_all(),
        include_history=False,       # not a candidate; baseline flag forces it on
        history_in_baseline=True,
    )


# Non-degenerate candidate set per condition (mirrors run_glm_split_by_condition
# + the histme additions): in V the platform is static (Speed≡Accel≡0) so only
# the visual terms + ME vary; in T_Vstatic the screen is grey (TF≡0, SF/OR≡NaN)
# so only Speed/Accel + ME vary. ME/History apply in both. Each condition
# carries ALL pairwise interactions on its own non-degenerate candidates
# (2026-06-15): a fair, complete candidate space, not a hand-picked subset.
#   V  : {TF, SF, OR, ME_face}        → C(4,2)=6 pairs
#   T  : {Speed, ME_face, Acceleration} → C(3,2)=3 pairs
# (Previously T set include_acceleration=True but omitted "Acceleration" from
# main_effects, so it was built and never offered as a candidate — the missing
# accel count in the T forward-selection summary. Fixed here.)
COND_CANDIDATES = {
    "V": {
        "main_effects": ("TF", "SF", "OR", "ME_face"),
        "interactions": (
            "TF_x_SF", "TF_x_OR", "SF_x_OR",
            "TF_x_ME_face", "SF_x_ME_face", "OR_x_ME_face",
        ),
        "include_acceleration": False,   # platform static → accel ≡ 0
    },
    "T_Vstatic": {
        "main_effects": ("Speed", "ME_face", "Acceleration"),
        "interactions": (
            "ME_face_x_Speed", "Speed_x_Acceleration", "ME_face_x_Acceleration",
        ),
        "include_acceleration": True,
    },
}


def make_config_histme_all_cond(condition: str) -> GLMConfig:
    """The _all (all-clusters, nominal-fallback) histme config, fit on ONE
    condition's trials with only that condition's non-degenerate candidates."""
    c = COND_CANDIDATES[condition]
    return replace(
        make_config_histme_all(),
        fit_condition=condition,
        main_effects=c["main_effects"],
        interactions=c["interactions"],
        include_acceleration=c["include_acceleration"],
    )


def make_config_histbase_all_cond(condition: str) -> GLMConfig:
    """make_config_histme_all_cond(condition) but History in BASELINE-nuisance
    mode — the condition-split analogue of make_config_histbase_all(). Reuses the
    same non-degenerate per-condition candidate set (History is never in it); only
    the history mode flips (include_history off, history_in_baseline on)."""
    return replace(
        make_config_histme_all_cond(condition),
        include_history=False,
        history_in_baseline=True,
    )


# Selected by --hist-me / --all-clusters in main(); run_probe/aggregate read these.
_CONFIG_FN = make_config
_RUN_LABEL = "Speed+TF+SF+OR+Acc"


def dry_run() -> int:
    """Report the selected cohort, RF clusters and the fit set per probe.

    The goggles cohort is the file's curated ``selected_clusters`` (not VISp
    anatomy). With ``apply_prefilter=False`` the rf_local fit set is
    selected ∩ RF; once every selected cluster has a measured-or-imputed RF
    that equals the full selected cohort.
    """
    for probe in PROBES:
        data = load_probe_data(
            FORMATTED_DIR / f"{probe}.mat", config=GLMConfig(),
            stimulus_lookup=_lookup(), cluster_set="selected",
        )
        selected = {c.cluster_id for c in data.clusters}
        rf = load_rf_sf_or(RF_PARQUET_DIR, probe).clusters
        log.info("%s: selected %d | RF %d | fit set (selected ∩ RF) = %d",
                 probe, len(selected), len(rf), len(selected & rf))
    return 0


def run_probe(probe: str, max_clusters: int | None = None,
              backend: str = "irls", device: str = "auto",
              n_jobs: int = 1) -> Path:
    out_dir = OUT_ROOT / "_runs" / probe
    out_dir.mkdir(parents=True, exist_ok=True)
    cluster_filter = None
    if max_clusters is not None:
        # Wiring smoke: restrict to the first few selected∩RF clusters.
        data = load_probe_data(
            FORMATTED_DIR / f"{probe}.mat", config=GLMConfig(),
            stimulus_lookup=_lookup(), cluster_set="selected",
        )
        cohort = {c.cluster_id for c in data.clusters}
        rf = load_rf_sf_or(RF_PARQUET_DIR, probe).clusters
        if getattr(_CONFIG_FN(), "rf_sf_or_nominal_fallback", False):
            # _all smoke: lead with no-RF clusters so the nominal fallback runs.
            ordered = sorted(cohort - rf) + sorted(cohort & rf)
            cluster_filter = set(ordered[:max_clusters])
        else:
            cluster_filter = set(sorted(cohort & rf)[:max_clusters])
        log.info("smoke: %s restricted to %s", probe, cluster_filter)
    # The fitting engine is a CLI choice: irls (CPU, the established lineage) or
    # nemos (JAX, GPU-capable on the cluster). device routes JAX to cpu/gpu/auto.
    cfg = replace(_CONFIG_FN(), device=device)
    run_pipeline(
        mat_path=FORMATTED_DIR / f"{probe}.mat",
        config=cfg,
        output_dir=out_dir,
        stimulus_lookup=_lookup(),
        backend=backend,
        cluster_set="selected",
        make_plots=True,
        plot_format="pdf",
        n_jobs=n_jobs,  # process-parallel fitting (loky). Plotting is main-process
                   # regardless. With the signed-rank 10-fold forward selection the
                   # serial path is slow, so n_jobs>1 is worth it — BUT pin BLAS to
                   # 1 thread/worker (OMP/OPENBLAS/MKL/VECLIB=1) to avoid thread
                   # oversubscription AND make cv-bps deterministic (single-thread
                   # BLAS has a fixed reduction order; the multi-thread wander is the
                   # ~±0.004 cv-bps nondeterminism). Default 1 keeps old runs as-is.
        cluster_filter=cluster_filter,
    )
    return out_dir


def aggregate() -> pd.DataFrame | None:
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

    from rc2_glm.plots import plot_forward_selection_summary, save_figure

    fig = plot_forward_selection_summary(out)
    fig.suptitle(
        f"Forward selection — goggles, RF-local SF/OR, 20 ms "
        f"({_RUN_LABEL}; n={out['cluster_id'].count()} clusters)",
        fontsize=12, fontweight="bold",
    )
    figdir = OUT_ROOT / "figs"
    figdir.mkdir(parents=True, exist_ok=True)
    for p in save_figure(fig, figdir / "forward_selection_summary", fmt="pdf"):
        log.info("wrote %s", p)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--probe", choices=PROBES, default=None)
    ap.add_argument("--max-clusters", type=int, default=None,
                    help="Restrict to the first N prefilter∩RF clusters (wiring smoke).")
    ap.add_argument("--hist-me", action="store_true",
                    help="Complete-picture variant: + History (2 lag bins) + face ME "
                         "→ current_rf_sfor_20ms_histme_goggles/.")
    ap.add_argument("--all-clusters", action="store_true",
                    help="hist-me config but INCLUDE no-RF clusters via the per-cloud "
                         "nominal SF/OR fallback → current_rf_sfor_20ms_histme_goggles_all/. "
                         "Implies --hist-me.")
    ap.add_argument("--condition", choices=("V", "T_Vstatic"), default=None,
                    help="Condition-SPLIT of the _all run: fit only V (VF) or T_Vstatic "
                         "trials, with that condition's non-degenerate candidates → "
                         "current_rf_sfor_20ms_histme_goggles_all_<cond>/.")
    ap.add_argument("--history-in-baseline", action="store_true",
                    help="all-clusters config but History runs as an always-on BASELINE "
                         "nuisance (never selected, never displayed) → "
                         "current_rf_sfor_20ms_histbase_goggles_all/. Implies --all-clusters.")
    ap.add_argument("--backend", choices=("irls", "nemos"), default="irls",
                    help="Fitting engine: irls (CPU, default, the established lineage) or "
                         "nemos (JAX, GPU-capable on the cluster).")
    ap.add_argument("--device", choices=("auto", "cpu", "gpu"), default="auto",
                    help="JAX device for --backend nemos (cpu/gpu/auto). Ignored for irls.")
    ap.add_argument("--n-jobs", type=int, default=1,
                    help="Process-parallel cluster fits (loky). Default 1 (serial). "
                         "For a deterministic parallel run pin BLAS to 1 thread/worker: "
                         "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
                         "VECLIB_MAXIMUM_THREADS=1 python scripts/run_glm_...py --n-jobs 8.")
    args = ap.parse_args()

    global OUT_ROOT, _CONFIG_FN, _RUN_LABEL
    if args.history_in_baseline and args.condition:
        cond = args.condition
        OUT_ROOT = OUT_HISTBASE_BY_COND[cond]
        _CONFIG_FN = lambda: make_config_histbase_all_cond(cond)
        _RUN_LABEL = f"all clusters, {cond} only; History=baseline"
    elif args.history_in_baseline:
        OUT_ROOT = OUT_HISTBASE_ALL
        _CONFIG_FN = make_config_histbase_all
        _RUN_LABEL = "Speed+TF+SF+OR+Acc+ME (all clusters; History=baseline)"
    elif args.condition:
        cond = args.condition
        OUT_ROOT = OUT_ALL_BY_COND[cond]
        _CONFIG_FN = lambda: make_config_histme_all_cond(cond)
        _RUN_LABEL = f"all clusters, {cond} only"
    elif args.all_clusters:
        OUT_ROOT = OUT_HISTME_ALL
        _CONFIG_FN = make_config_histme_all
        _RUN_LABEL = "Speed+TF+SF+OR+Acc+ME+History (all clusters)"
    elif args.hist_me:
        OUT_ROOT = OUT_HISTME
        _CONFIG_FN = make_config_histme
        _RUN_LABEL = "Speed+TF+SF+OR+Acc+ME+History"

    if args.dry_run:
        return dry_run()

    probes = (args.probe,) if args.probe else PROBES
    for probe in probes:
        run_probe(probe, max_clusters=args.max_clusters,
                  backend=args.backend, device=args.device, n_jobs=args.n_jobs)
    if args.max_clusters is None:
        aggregate()
        # Diagnostics are PART OF THE PIPELINE — generated at the end of every
        # full run (both probes), never a separate manual step (Laura 2026-06-15).
        # Late import breaks the circular dependency (diagnostics imports this
        # driver). A diagnostics failure is logged, not fatal — the fits + the
        # root figs/ are already on disk.
        run_key = ("histbase_V" if args.history_in_baseline and args.condition == "V"
                   else "histbase_Tvstatic" if args.history_in_baseline and args.condition == "T_Vstatic"
                   else "histbase" if args.history_in_baseline
                   else "all_V" if args.condition == "V"
                   else "all_Tvstatic" if args.condition == "T_Vstatic"
                   else "all" if args.all_clusters
                   else "histme" if args.hist_me else "main")
        try:
            from scripts.diagnostics_rf_sfor import generate as _generate_diagnostics
            log.info("pipeline: generating diagnostics for the %s run", run_key)
            _generate_diagnostics(run_key)
        except Exception as exc:  # noqa: BLE001 — keep the run's outputs even if diag fails
            log.warning("diagnostics step failed (run outputs intact): %s", exc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
