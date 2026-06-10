"""current_plus_ME_20ms/ — the pooled +ME run refit at 20 ms bins, with spike
history ON, as a parallel comparison to the 100 ms current_plus_ME/.

Why 20 ms: at 100 ms the only history lags available (>=100 ms) overlap the slow
stimulus envelope, so a "history" term becomes a catch-all that absorbs Speed
(diagnosed 2026-04-30). At 20 ms history finally has FAST lags (refractory ~one
20 ms bin, bursting 2-3 bins) that are timescale-separated from the slow
stimulus. We keep the window SHORT (100 ms = 5 lag bins) so it cannot reach back
into the slow regime. Everything else matches current_plus_ME/ (3 face-camera
probes, pooled, ME_face, speed-profile CV, reused cohort).

The honest test (run separately): does Speed's unique contribution survive
adding History (and ME)? If yes, history is a clean point-process nuisance.
"""
from __future__ import annotations

import sys
from dataclasses import replace

import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.pipeline import run_pipeline

import scripts.run_glm_split_by_condition as drv

PROBES_ME = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
# WITH history (4-lag/80 ms) -> current_plus_ME_20ms; --no-history -> _nohist;
# --refractory (2-lag/40 ms, refractory-only) -> _refractory. New folders each.
OUT_HIST = drv.ROOT / "figures" / "glm" / "current_plus_ME_20ms"
OUT_NOHIST = drv.ROOT / "figures" / "glm" / "current_plus_ME_20ms_nohist"
OUT_REFRAC = drv.ROOT / "figures" / "glm" / "current_plus_ME_20ms_refractory"
OUT = OUT_HIST                                   # default (diagnostic scripts import this)


def make_config_me_20ms(include_history: bool = True,
                        history_window_s: float = 0.08) -> GLMConfig:
    n_lags = max(1, round(history_window_s / 0.02))
    return replace(
        GLMConfig(),
        time_bin_width=0.02,                     # 20 ms bins
        fit_condition=None,                      # pooled
        main_effects=("Speed", "TF", "SF", "OR", "ME_face"),
        interactions=("Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
                      "TF_x_SF", "TF_x_OR", "SF_x_OR", "ME_face_x_Speed"),
        include_me_face=True, motion_energy_camera="camera0",
        include_history=include_history,         # ON here; control run = OFF
        history_window_s=history_window_s,       # 0.08 = 4 lags; 0.04 = 2-lag refractory
        history_basis_kind="identity",           # per-lag dummies, not raised cosine
        n_history_bases=n_lags,                   # ignored for identity
        allow_history_interactions=False,
        include_onset_kernel=True,
        lambda_ridge=1.0,
        cv_strategy="speed-profile", n_selection_seeds=1, selection_threshold_count=1,
        profile_cv_diagnostic=True,
        apply_prefilter=False,                   # cohort fixed to the reused 88 (3-probe subset)
    )


def main() -> int:
    no_history = "--no-history" in sys.argv
    refractory = "--refractory" in sys.argv
    if refractory:
        OUT, win, tag = OUT_REFRAC, 0.04, "+ME +history(refractory 2-lag/40 ms) 20 ms"
    elif no_history:
        OUT, win, tag = OUT_NOHIST, 0.08, "no-history 20 ms control"
    else:
        OUT, win, tag = OUT_HIST, 0.08, "+ME +history(4-lag/80 ms) 20 ms"
    print(f"=== current_plus_ME_20ms run: {tag} -> {OUT.name} ===")
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    cfg = make_config_me_20ms(include_history=not no_history, history_window_s=win)
    for probe in PROBES_ME:
        out_dir = OUT / "_runs" / probe
        if (out_dir / "glm_model_comparison.csv").exists():
            print(f"skip {probe} — already at {out_dir}")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"fitting {probe} (+ME, +history, 20 ms, pooled, speed-profile)")
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
        print(f"  ME_face  selected: {int(sv.str.contains('ME_face', regex=False).sum())}/{n}")
        print(f"  History  selected: {int(sv.str.contains('History', regex=False).sum())}/{n}")
        print(f"  Speed    selected: {int(sv.str.contains('Speed', regex=False).sum())}/{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
