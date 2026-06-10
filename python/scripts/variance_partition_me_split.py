"""ME split: across-trial vs within-trial — captured inter-trial variability
vs deviance shared with Speed.

ME_face per bin carries two physically different signals:
  - ME-across  = the trial mean (one value per trial), broadcast to its bins.
      The stimulus is identical across trials of a profile, so Speed/TF/onset
      CANNOT represent it; any deviance it explains is genuine inter-trial
      variability and is orthogonal to Speed by construction.
  - ME-within  = ME minus its trial mean (moment-to-moment fluctuation).
      Rides the velocity trajectory, so it overlaps Speed/onset within a trial
      — the only place ME can "steal" from Speed.

Per cohort cluster (3 face-camera probes, pooled, speed-profile folds) we
compute cv-bps for nested models and report:
  - unique(ME-across), unique(ME-within)  within the +both model
  - unique(Speed) under no-ME / +across / +within / +both   ← which ME part
    eats Speed's unique contribution.

Output: current_plus_ME/diagnostics/variance_partition_me_split.{csv,pdf,png}
"""
from __future__ import annotations

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


def _zfill(x):
    """z-score over finite entries; non-finite -> 0."""
    f = np.isfinite(x)
    if f.sum() < 2:
        return np.zeros_like(x)
    mu = float(x[f].mean()); sd = float(x[f].std(ddof=0)) or 1.0
    return np.where(np.isfinite(x), (x - mu) / sd, 0.0)


def _me_components(df, cfg):
    """Return (B_across, B_within) bases, or None if ME unusable."""
    if "me_face_raw" not in df.columns:
        return None
    me = df["me_face_raw"].to_numpy(float)
    tid = df["trial_id"].to_numpy()
    fin = np.isfinite(me)
    if int(fin.sum()) < 10:
        return None
    trial_mean = np.full_like(me, np.nan)
    for t in np.unique(tid):
        m = tid == t
        vals = me[m & fin]
        if vals.size:
            trial_mean[m] = vals.mean()
    across = trial_mean                       # trial-level state
    within = me - trial_mean                   # within-trial fluctuation
    n, lo_hi = cfg.n_me_face_bases, cfg.me_face_range
    B_across = raised_cosine_basis_linear(_zfill(across), n, *lo_hi)
    B_within = raised_cosine_basis_linear(_zfill(within), n, *lo_hi)
    return B_across, B_within


def partition(df, cfg):
    if df["spike_count"].sum() == 0:
        return None
    comp = _me_components(df, cfg)
    if comp is None:
        return None
    B_across, B_within = comp
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

    def cv(stim_vars, me_parts, onset_on=True):
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, list(stim_vars),
            include_onset_kernel=onset_on, B_me_face=None)
        extra = []
        if "across" in me_parts:
            extra.append(B_across)
        if "within" in me_parts:
            extra.append(B_within)
        if extra:
            X = np.hstack([X] + extra)
        return cross_validate_glm(X, y, offset, folds,
                                  lambda_ridge=cfg.lambda_ridge).cv_bits_per_spike

    noSpeed = [v for v in STIM if v != "Speed"]
    cv_int = cv([], [], onset_on=False)
    cv_stim = cv(STIM, [])
    cv_across = cv(STIM, ["across"])
    cv_within = cv(STIM, ["within"])
    cv_both = cv(STIM, ["across", "within"])
    return {
        "cv_intercept": cv_int, "cv_stim": cv_stim, "cv_both": cv_both,
        "total_both": cv_both - cv_int,
        # which ME component carries explained deviance (within the +both model)
        "unique_ME_across": cv_both - cv_within,
        "unique_ME_within": cv_both - cv_across,
        "unique_ME_total": cv_both - cv_stim,
        # which ME component eats Speed's unique
        "unique_Speed_noME": cv_stim - cv(noSpeed, []),
        "unique_Speed_across": cv_across - cv(noSpeed, ["across"]),
        "unique_Speed_within": cv_within - cv(noSpeed, ["within"]),
        "unique_Speed_both": cv_both - cv(noSpeed, ["across", "within"]),
    }


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
    df.to_csv(out_dir / "variance_partition_me_split.csv", index=False)
    well = df["total_both"] > 0.005

    def m(k):
        return df.loc[well, k].median()
    print(f"\nn={len(df)} (well-fit {int(well.sum())}); medians Δcv-bps:")
    print(f"  unique ME-across   = {m('unique_ME_across'):+.4f}  "
          f"[>0.005: {int((df.loc[well,'unique_ME_across']>0.005).sum())}/{int(well.sum())}]")
    print(f"  unique ME-within   = {m('unique_ME_within'):+.4f}  "
          f"[>0.005: {int((df.loc[well,'unique_ME_within']>0.005).sum())}/{int(well.sum())}]")
    print(f"  unique Speed  noME = {m('unique_Speed_noME'):+.4f}")
    print(f"  unique Speed +acrs = {m('unique_Speed_across'):+.4f}  (across-only)")
    print(f"  unique Speed +wthn = {m('unique_Speed_within'):+.4f}  (within-only)")
    print(f"  unique Speed +both = {m('unique_Speed_both'):+.4f}")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))

    # Panel L: which ME component carries the captured deviance.
    groups = [("ME-across\n(inter-trial)", "unique_ME_across", "tab:purple"),
              ("ME-within\n(movement)", "unique_ME_within", "tab:pink")]
    data = [df.loc[well, c].to_numpy() for _, c, _ in groups]
    bp = axL.boxplot(data, tick_labels=[g for g, _, _ in groups], showfliers=False,
                     patch_artist=True, widths=0.5)
    for patch, (_, _, c) in zip(bp["boxes"], groups):
        patch.set_facecolor(c); patch.set_alpha(0.45)
    for i, d in enumerate(data, 1):
        axL.scatter(np.full(len(d), i) + np.random.uniform(-.08, .08, len(d)), d,
                    s=10, alpha=0.4, color="0.3", edgecolor="none")
    axL.axhline(0, color="0.6", lw=0.8); axL.axhline(0.005, color="0.6", lw=0.8, ls=":")
    axL.set_ylabel("unique Δ cv-bps (+both ME model)")
    axL.set_title(f"Which ME part carries the captured deviance? (n={int(well.sum())})",
                  fontsize=10, fontweight="bold")

    # Panel R: Speed's unique under each ME condition (paired lines).
    conds = [("no ME", "unique_Speed_noME"), ("+across", "unique_Speed_across"),
             ("+within", "unique_Speed_within"), ("+both", "unique_Speed_both")]
    M = np.column_stack([df.loc[well, c].to_numpy() for _, c in conds])
    x = np.arange(len(conds))
    for row in M:
        axR.plot(x, row, color="tab:green", alpha=0.12, lw=0.8)
    axR.plot(x, np.median(M, axis=0), color="black", lw=2.4, marker="o",
             label="median")
    for xi, med in zip(x, np.median(M, axis=0)):
        axR.annotate(f"{med:+.3f}", (xi, med), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=9, fontweight="bold")
    axR.axhline(0, color="0.6", lw=0.8)
    axR.set_xticks(x); axR.set_xticklabels([c for c, _ in conds])
    axR.set_ylabel("unique Speed (Δ cv-bps)")
    axR.set_title("Which ME part eats Speed's unique?\n"
                  "(flat across→within = stealing is within-trial)",
                  fontsize=10, fontweight="bold")
    axR.legend(fontsize=9, frameon=False)

    fig.suptitle("ME across- vs within-trial — captured inter-trial variability "
                 "vs shared-with-Speed", fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = out_dir / "variance_partition_me_split"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
