"""ME nuisance variance partition — does Speed's unique survive adding ME?

On the 3 face-camera probes (pooled, speed-profile folds), per cluster compute
cv-bps for nested models with and without ME_face, and:
  - unique_<g>  (leave-one-out within the +ME model) for Onset/Speed/TF/SF/OR/ME_face
  - unique_Speed WITHOUT ME vs WITH ME  ← the acid test
  - unique_ME = cv(full+ME) - cv(full no ME)

If unique_Speed is ~invariant to adding ME, ME is a clean interpretable nuisance
(state covariate) and Speed tuning stands. If unique_Speed collapses, ME is
absorbing stimulus structure.

Output: current_plus_ME/diagnostics/variance_partition_ME.{csv,pdf}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import (onset_kernel_basis, raised_cosine_basis_linear,
                           value_basis)
from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME"
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
    motion = df[df["condition"] != "stationary"]
    if motion.empty or df["spike_count"].sum() == 0:
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
    folds = make_trial_folds(trial_ids, cfg.n_folds, cfg.cv_seed,
                             condition_labels_per_bin=cond,
                             strategy="speed-profile", profile_ids_per_bin=profile_ids)

    def cv(vars_, me_on, onset_on=True):
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, list(vars_),
            include_onset_kernel=onset_on,
            B_me_face=B_me if me_on else None)
        return cross_validate_glm(X, y, offset, folds,
                                  lambda_ridge=cfg.lambda_ridge).cv_bits_per_spike

    cv_int = cv([], False, onset_on=False)
    cv_noME = cv(STIM, False)
    cv_noME_noSpeed = cv([v for v in STIM if v != "Speed"], False)
    cv_ME = cv(STIM + ["ME_face"], True)
    out = {
        "cv_intercept": cv_int, "cv_full_noME": cv_noME, "cv_full_ME": cv_ME,
        "total_noME": cv_noME - cv_int, "total_ME": cv_ME - cv_int,
        "unique_ME": cv_ME - cv_noME,
        "unique_Speed_noME": cv_noME - cv_noME_noSpeed,
        "unique_Speed_ME": cv_ME - cv(["TF", "SF", "OR", "ME_face"], True),
    }
    # leave-one-out uniques within the +ME model
    out["unique_Onset"] = cv_ME - cv(STIM + ["ME_face"], True, onset_on=False)
    for g in STIM:
        out[f"unique_{g}"] = cv_ME - cv([v for v in STIM if v != g] + ["ME_face"], True)
    return out


def main() -> int:
    cfg = drv.make_config(None)
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
    df.to_csv(out_dir / "variance_partition_ME.csv", index=False)
    n = len(df)
    well = df["total_ME"] > 0.005

    def m(k):
        return df.loc[well, k].median()
    print(f"\nn={n} (well-fit {int(well.sum())}); medians Δcv-bps:")
    print(f"  unique_ME             = {m('unique_ME'):+.4f}  [>0.005: {int((df.loc[well,'unique_ME']>0.005).sum())}/{int(well.sum())}]")
    print(f"  unique_Speed (no ME)  = {m('unique_Speed_noME'):+.4f}")
    print(f"  unique_Speed (+ME)    = {m('unique_Speed_ME'):+.4f}   <-- acid test")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6))
    groups = [("Onset", "unique_Onset", "tab:blue"), ("Speed", "unique_Speed_ME", "tab:green"),
              ("TF", "unique_TF", "tab:orange"), ("SF", "unique_SF", "tab:olive"),
              ("OR", "unique_OR", "tab:red"), ("ME_face", "unique_ME", "tab:purple")]
    data = [df.loc[well, c].to_numpy() for _, c, _ in groups]
    bp = axL.boxplot(data, tick_labels=[g for g, _, _ in groups], showfliers=False,
                     patch_artist=True)
    for patch, (_, _, c) in zip(bp["boxes"], groups):
        patch.set_facecolor(c); patch.set_alpha(0.4)
    axL.axhline(0, color="0.6", lw=0.8); axL.axhline(0.005, color="0.6", lw=0.8, ls=":")
    axL.set_ylabel("unique Δ cv-bps (+ME model)")
    axL.set_title(f"Unique contribution per group, +ME model (n={int(well.sum())})",
                  fontsize=10, fontweight="bold")

    # Acid test: Speed unique without vs with ME (paired).
    a = df.loc[well, "unique_Speed_noME"].to_numpy()
    b = df.loc[well, "unique_Speed_ME"].to_numpy()
    axR.scatter(a, b, s=18, color="tab:green", alpha=0.6, edgecolor="none")
    lim = [min(a.min(), b.min(), 0) * 1.05, max(a.max(), b.max()) * 1.05]
    axR.plot(lim, lim, ls="--", color="0.5", lw=1)
    axR.axhline(0, color="0.8", lw=0.6); axR.axvline(0, color="0.8", lw=0.6)
    axR.set_xlabel("unique Speed — NO ME"); axR.set_ylabel("unique Speed — +ME")
    axR.set_title(f"Acid test: does Speed survive ME?\n"
                  f"median {np.median(a):+.3f} → {np.median(b):+.3f} "
                  f"(on diagonal = survives)", fontsize=10, fontweight="bold")
    fig.suptitle("ME as nuisance — variance partition (current_plus_ME, speed-profile)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = out_dir / "variance_partition_ME"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
