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
MC_SEQUENCE = ROOT / "motion_clouds_goggles_sequence_260420.mat"
MC_FOLDERS = ROOT / "image_folders_goggles.mat"

PROBES = (
    "CAA-1124370_rec1_rec2_rec3",
    "CAA-1124371_rec1_rec2_rec3",
)
CONDITIONS = ("V", "T_Vstatic", "VT")

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
    return out


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
    args = parser.parse_args()

    if args.dry_run:
        return dry_run()

    probes = (args.probe,) if args.probe else PROBES
    for probe in probes:
        run_probe(probe)
    aggregate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
