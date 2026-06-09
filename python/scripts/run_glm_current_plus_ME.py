"""current_plus_ME/ — pooled GLM + face motion energy, speed-profile CV.

Pooled model (all conditions together) with ME_face added as a candidate, on
the 3 probes with real face-camera data (243/244/466; 467 excluded, ~50%
coverage). Corrected velocity-derived profile_id + speed-profile leak-test CV
+ profile_cv_diagnostic — same identifiable setup as current/, just +ME and
3 probes. Output kept separate from the goggles run as current_plus_ME/.

Purpose: ME_face as an interpretable trial-to-trial *state* covariate
(Stringer 2019 / Musall 2019 / Lohuis) to mop up across-trial rate variability,
then a variance partition with ME as a group (the test: Speed's unique
contribution must survive adding ME).
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.pipeline import run_pipeline

import scripts.run_glm_split_by_condition as drv

PROBES_ME = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME"


def make_config_me() -> GLMConfig:
    return replace(
        GLMConfig(),
        fit_condition=None,                      # pooled: all conditions
        main_effects=("Speed", "TF", "SF", "OR", "ME_face"),
        interactions=("Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
                      "TF_x_SF", "TF_x_OR", "SF_x_OR", "ME_face_x_Speed"),
        include_me_face=True, motion_energy_camera="camera0",
        include_history=False, include_onset_kernel=True,
        lambda_ridge=1.0,
        cv_strategy="speed-profile", n_selection_seeds=1, selection_threshold_count=1,
        profile_cv_diagnostic=True,
        apply_prefilter=False,                   # cohort fixed to the reused 88 (3-probe subset)
    )


def main() -> int:
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    cfg = make_config_me()
    for probe in PROBES_ME:
        out_dir = OUT / "_runs" / probe
        if (out_dir / "glm_model_comparison.csv").exists():
            print(f"skip {probe} — already at {out_dir}")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"fitting {probe} (+ME, pooled, speed-profile)")
        run_pipeline(
            mat_path=drv.FORMATTED_DIR / f"{probe}.mat",
            config=cfg, output_dir=out_dir, stimulus_lookup=lookup,
            backend="irls", visp_only=True, make_plots=True, plot_format="pdf",
            n_jobs=4, cluster_filter=drv.load_cohort(probe),
        )
    frames = []
    for probe in PROBES_ME:
        f = OUT / "_runs" / probe / "glm_model_comparison.csv"
        if f.exists():
            d = pd.read_csv(f); d["probe_id"] = probe; frames.append(d)
    if frames:
        out = pd.concat(frames, ignore_index=True)
        out.to_csv(OUT / "glm_model_comparison.csv", index=False)
        sv = out["time_selected_vars"].fillna("Null")
        n = len(out)
        print(f"wrote {OUT}/glm_model_comparison.csv ({n} clusters)")
        print(f"ME_face selected: {int(sv.str.contains('ME_face', regex=False).sum())}/{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
