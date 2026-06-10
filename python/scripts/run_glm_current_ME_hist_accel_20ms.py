"""current_ME_hist_accel_20ms/ — canonical + ME + History + Acceleration, the
forward-selectable model (NO interactions, NO population). Full pipeline run
with the standard per-probe figures + _runs.

20 ms bins; Speed/TF/SF/OR + onset + ME_face + refractory History (2-lag/40 ms,
identity dummies) + Acceleration (signed binned translation accel, value-axis
candidate). interactions=() — skip any mixed multiplicative term in forward
selection. Speed-profile CV. 3 face-camera probes (ME_face needs camera0).
"""
from __future__ import annotations

from dataclasses import replace

import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.pipeline import run_pipeline

import scripts.run_glm_split_by_condition as drv

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_ME_hist_accel_20ms"


def make_config() -> GLMConfig:
    return replace(
        GLMConfig(),
        time_bin_width=0.02,
        fit_condition=None,
        main_effects=("Speed", "TF", "SF", "OR", "ME_face", "Acceleration"),
        interactions=(),                          # skip mixed multiplicative terms
        include_me_face=True, motion_energy_camera="camera0",
        include_acceleration=True, n_accel_bases=5, accel_range=(-3.0, 3.0),
        include_history=True, history_window_s=0.04,
        history_basis_kind="identity", n_history_bases=2,
        include_onset_kernel=True,
        lambda_ridge=1.0,
        cv_strategy="speed-profile", n_selection_seeds=1, selection_threshold_count=1,
        profile_cv_diagnostic=True,
        apply_prefilter=False,
    )


def main() -> int:
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    cfg = make_config()
    for probe in PROBES:
        out_dir = OUT / "_runs" / probe
        if (out_dir / "glm_model_comparison.csv").exists():
            print(f"skip {probe} — already at {out_dir}")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"fitting {probe} (canonical+ME+History+Accel, 20 ms, no interactions)")
        run_pipeline(
            mat_path=drv.FORMATTED_DIR / f"{probe}.mat",
            config=cfg, output_dir=out_dir, stimulus_lookup=lookup,
            backend="irls", visp_only=True, make_plots=True, plot_format="pdf",
            n_jobs=4, cluster_filter=drv.load_cohort(probe),
        )
    frames = []
    for probe in PROBES:
        f = OUT / "_runs" / probe / "glm_model_comparison.csv"
        if f.exists():
            d = pd.read_csv(f); d["probe_id"] = probe; frames.append(d)
    if frames:
        out = pd.concat(frames, ignore_index=True)
        out.to_csv(OUT / "glm_model_comparison.csv", index=False)
        sv = out["time_selected_vars"].fillna("Null")
        n = len(out)
        print(f"wrote {OUT}/glm_model_comparison.csv ({n} clusters)")
        for v in ("Speed", "TF", "SF", "OR", "ME_face", "Acceleration", "History"):
            print(f"  {v:13s} {int(sv.str.contains(v, regex=False).sum())}/{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
