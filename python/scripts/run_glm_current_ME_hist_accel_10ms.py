"""current_ME_hist_accel_10ms/ — same model as the 20 ms run but at 10 ms bins,
with the history window = 20 ms total (= 2 lags of 10 ms, identity dummies).

Canonical Speed/TF/SF/OR + onset + ME_face + History(2-lag/20 ms) + Acceleration,
no interactions, no population. Full pipeline run with standard figures + _runs.
"""
from __future__ import annotations

from dataclasses import replace

import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.pipeline import run_pipeline

import scripts.run_glm_split_by_condition as drv
from scripts.run_glm_current_ME_hist_accel_20ms import make_config as _make_20, PROBES

OUT = drv.ROOT / "figures" / "glm" / "current_ME_hist_accel_10ms"


def make_config() -> GLMConfig:
    # 10 ms bins; history window 20 ms → round(0.02/0.01) = 2 identity lag dummies.
    return replace(_make_20(), time_bin_width=0.01, history_window_s=0.02)


def main() -> int:
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    cfg = make_config()
    for probe in PROBES:
        out_dir = OUT / "_runs" / probe
        if (out_dir / "glm_model_comparison.csv").exists():
            print(f"skip {probe} — already at {out_dir}")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"fitting {probe} (canonical+ME+History(20ms)+Accel, 10 ms bins)")
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
