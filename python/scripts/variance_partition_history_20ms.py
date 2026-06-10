"""History acid test at 20 ms — does Speed's unique survive adding History
(and ME)? Companion to variance_partition_ME.py, refit at 20 ms with the short
(100 ms / 5-lag) history kernel.

Per cohort cluster (3 face-camera probes, pooled, speed-profile folds) compute
cv-bps for nested models and report:
  - unique_<g> (leave-one-out in the full Stim+ME+History model) per block
  - unique_Speed under {no nuisance / +ME / +History / +both}  ← the acid test
  - unique_History, unique_ME within the full model

The honest expectation: a SHORT history kernel is timescale-separated from the
slow stimulus, so unique_Speed should be ~flat as History is added (unlike
ME-within). If Speed collapses under History, the 20 ms history is still
absorbing stimulus and needs an even shorter window.

Output: current_plus_ME_20ms/diagnostics/variance_partition_history.{csv,pdf,png}
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import (convolve_history, history_basis, onset_kernel_basis,
                           raised_cosine_basis_linear, value_basis)
from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

import sys

import scripts.run_glm_split_by_condition as drv
from scripts.run_glm_current_plus_ME_20ms import make_config_me_20ms

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
# --refractory -> 2-lag/40 ms history into the _refractory folder; else 4-lag/80 ms.
REFRACTORY = "--refractory" in sys.argv
OUT = drv.ROOT / "figures" / "glm" / (
    "current_plus_ME_20ms_refractory" if REFRACTORY else "current_plus_ME_20ms")
HISTORY_WINDOW = 0.04 if REFRACTORY else 0.08
STIM = ["Speed", "TF", "SF", "OR"]


def _build_me(df, cfg):
    if "me_face_raw" not in df.columns:
        return None
    me = df["me_face_raw"].to_numpy(float)
    motion = (df["condition"] != "stationary").to_numpy()
    fin = np.isfinite(me) & motion
    if int(fin.sum()) < 10:
        return None
    mu = float(me[fin].mean()); sd = float(me[fin].std(ddof=0)) or 1.0
    mez = np.where(np.isfinite(me), (me - mu) / sd, 0.0)
    return raised_cosine_basis_linear(mez, cfg.n_me_face_bases, *cfg.me_face_range)


def partition(df, cfg):
    if df["spike_count"].sum() == 0:
        return None
    B_me = _build_me(df, cfg)
    if B_me is None:
        return None
    speed = df["speed"].to_numpy(float); tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float); orient = df["orientation"].to_numpy(float)
    y = df["spike_count"].to_numpy(float)
    trial_ids = df["trial_id"].to_numpy(np.int64)
    profile_ids = df["profile_id"].to_numpy(np.int64)
    cond = df["condition"].to_numpy(object)
    if len(np.unique(profile_ids)) < 2:
        return None
    offset = float(np.log(cfg.time_bin_width))
    sp = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=sp)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=sp)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    hb = history_basis(cfg.n_history_bases, cfg.history_window_s, cfg.time_bin_width,
                       kind=cfg.history_basis_kind)
    B_hist = convolve_history(y, trial_ids, hb)
    folds = make_trial_folds(trial_ids, cfg.n_folds, cfg.cv_seed,
                             condition_labels_per_bin=cond,
                             strategy="speed-profile", profile_ids_per_bin=profile_ids)

    def cv(vars_, me_on, hist_on, onset_on=True):
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, list(vars_),
            include_onset_kernel=onset_on,
            B_me_face=B_me if me_on else None,
            B_history=B_hist if hist_on else None)
        return cross_validate_glm(X, y, offset, folds,
                                  lambda_ridge=cfg.lambda_ridge).cv_bits_per_spike

    noSpeed = [v for v in STIM if v != "Speed"]
    cv_int = cv([], False, False, onset_on=False)
    cv_stim = cv(STIM, False, False)
    cv_ME = cv(STIM + ["ME_face"], True, False)
    cv_H = cv(STIM + ["History"], False, True)
    cv_full = cv(STIM + ["ME_face", "History"], True, True)
    out = {
        "cv_intercept": cv_int, "cv_stim": cv_stim, "cv_ME": cv_ME,
        "cv_H": cv_H, "cv_full": cv_full,
        "total_full": cv_full - cv_int,
        "unique_History": cv_full - cv_ME,            # add H given ME present
        "unique_ME": cv_full - cv_H,                  # add ME given H present
        "unique_History_alone": cv_H - cv_stim,       # H beyond stim only
        # acid test: Speed's unique under each nuisance set
        "uSpeed_none": cv_stim - cv(noSpeed, False, False),
        "uSpeed_ME": cv_ME - cv(noSpeed + ["ME_face"], True, False),
        "uSpeed_H": cv_H - cv(noSpeed + ["History"], False, True),
        "uSpeed_both": cv_full - cv(noSpeed + ["ME_face", "History"], True, True),
    }
    out["unique_Onset"] = cv_full - cv(STIM + ["ME_face", "History"], True, True,
                                       onset_on=False)
    for g in STIM:
        out[f"unique_{g}"] = cv_full - cv(
            [v for v in STIM if v != g] + ["ME_face", "History"], True, True)
    return out


def main() -> int:
    cfg = make_config_me_20ms(history_window_s=HISTORY_WINDOW)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    rows = []
    for probe in PROBES:
        cohort = drv.load_cohort(probe)
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        for cl in pdata.clusters:
            if cl.cluster_id not in cohort:
                continue
            res = partition(bin_cluster(pdata, cl), cfg)
            if res is None:
                continue
            res.update(probe_id=probe, cluster_id=int(cl.cluster_id))
            rows.append(res)
        print(f"{probe}: {sum(r['probe_id']==probe for r in rows)} clusters")

    df = pd.DataFrame(rows)
    out_dir = OUT / "diagnostics"; out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "variance_partition_history.csv", index=False)
    well = df["total_full"] > 0.005

    def m(k):
        return df.loc[well, k].median()
    print(f"\nn={len(df)} (well-fit {int(well.sum())}); medians Δcv-bps:")
    print(f"  unique History (given ME) = {m('unique_History'):+.4f}  "
          f"[>0.005: {int((df.loc[well,'unique_History']>0.005).sum())}/{int(well.sum())}]")
    print(f"  unique ME (given History) = {m('unique_ME'):+.4f}")
    print(f"  Speed unique  none={m('uSpeed_none'):+.4f}  +ME={m('uSpeed_ME'):+.4f}  "
          f"+H={m('uSpeed_H'):+.4f}  +both={m('uSpeed_both'):+.4f}")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    groups = [("Onset", "unique_Onset", "tab:blue"), ("Speed", "unique_Speed", "tab:green"),
              ("TF", "unique_TF", "tab:orange"), ("SF", "unique_SF", "tab:olive"),
              ("OR", "unique_OR", "tab:red"), ("ME", "unique_ME", "tab:purple"),
              ("History", "unique_History", "tab:brown")]
    data = [df.loc[well, c].to_numpy() for _, c, _ in groups]
    bp = axL.boxplot(data, tick_labels=[g for g, _, _ in groups], showfliers=False,
                     patch_artist=True)
    for patch, (_, _, c) in zip(bp["boxes"], groups):
        patch.set_facecolor(c); patch.set_alpha(0.45)
    axL.axhline(0, color="0.6", lw=0.8); axL.axhline(0.005, color="0.6", lw=0.8, ls=":")
    axL.set_ylabel("unique Δ cv-bps (full 20 ms model)")
    axL.set_title(f"Unique contribution per group, 20 ms +ME +History (n={int(well.sum())})",
                  fontsize=10, fontweight="bold")

    conds = [("none", "uSpeed_none"), ("+ME", "uSpeed_ME"),
             ("+History", "uSpeed_H"), ("+both", "uSpeed_both")]
    M = np.column_stack([df.loc[well, c].to_numpy() for _, c in conds])
    x = np.arange(len(conds))
    for row in M:
        axR.plot(x, row, color="tab:green", alpha=0.12, lw=0.8)
    axR.plot(x, np.median(M, axis=0), color="black", lw=2.4, marker="o")
    for xi, med in zip(x, np.median(M, axis=0)):
        axR.annotate(f"{med:+.3f}", (xi, med), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=9, fontweight="bold")
    axR.axhline(0, color="0.6", lw=0.8)
    axR.set_xticks(x); axR.set_xticklabels([c for c, _ in conds])
    axR.set_ylabel("unique Speed (Δ cv-bps)")
    axR.set_title("Acid test at 20 ms: does Speed survive History?\n"
                  "(flat none→+History = history is a clean point-process nuisance)",
                  fontsize=10, fontweight="bold")

    fig.suptitle(f"History acid test ({OUT.name}, "
                 f"{int(round(HISTORY_WINDOW*1000))} ms / {int(round(HISTORY_WINDOW/0.02))}-lag kernel)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = out_dir / "variance_partition_history"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
