"""Acceleration acid test — Speed vs Acceleration attribution.

Base model = the settled latest config: 20 ms bins, Speed/TF/SF/OR + ME_face +
refractory (2-lag/40 ms) spike history. We add an ACCELERATION block and ask:
  - does Acceleration carry unique cv-bps beyond the latest model?
  - does Speed's unique SURVIVE adding Acceleration (acid test)?
  - per-cell: Speed vs Acceleration unique (which dominates).

Acceleration = signed per-trial derivative of the binned translation speed
(d|v|/dt at the 20 ms model resolution; 0 in V where the platform is static,
like Speed). The two velocity profiles differ in acceleration time-course, so
the design may decouple Speed↔Accel (unlike Speed↔TF).

Output: current_plus_ME_20ms_accel/diagnostics/variance_partition_accel.{csv,pdf,png}
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
# --bin10 → the 10 ms run's config + folder; default = 20 ms.
BIN10 = "--bin10" in sys.argv
BIN_MS = 10 if BIN10 else 20
RUN_NAME = "current_ME_hist_accel_10ms" if BIN10 else "current_plus_ME_20ms_accel"
OUT = drv.ROOT / "figures" / "glm" / RUN_NAME
STIM = ["Speed", "TF", "SF", "OR"]
HISTORY_WINDOW = 0.04          # refractory base (20 ms run)


def _make_cfg():
    if BIN10:
        from scripts.run_glm_current_ME_hist_accel_10ms import make_config
        return make_config()
    return make_config_me_20ms(history_window_s=HISTORY_WINDOW)


def _zfill(x):
    f = np.isfinite(x)
    if f.sum() < 2:
        return np.zeros_like(x)
    mu = float(x[f].mean()); sd = float(x[f].std(ddof=0)) or 1.0
    return np.where(np.isfinite(x), (x - mu) / sd, 0.0)


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


def _build_accel(df, cfg):
    """Real binned signed translation acceleration (time_binning 'acceleration'
    column = MATLAB-style 100·diff(v), bin-mean, 0 in V), z-scored."""
    if "acceleration" not in df.columns:
        raise KeyError("no 'acceleration' column — rebuild needs updated time_binning")
    a = df["acceleration"].to_numpy(float)
    az = np.clip(_zfill(a), -3, 3)
    return raised_cosine_basis_linear(az, 5, -3.0, 3.0)


def partition(df, cfg):
    if df["spike_count"].sum() == 0:
        return None
    B_me = _build_me(df, cfg)
    if B_me is None:
        return None
    B_accel = _build_accel(df, cfg)
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

    def cv(vars_, me_on, hist_on, accel_on, onset_on=True):
        # assemble only activates B_me_face/B_history when the NAME is in the
        # selected list — must add them, not just pass the basis.
        sel = list(vars_)
        if me_on:
            sel = sel + ["ME_face"]
        if hist_on:
            sel = sel + ["History"]
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, sel,
            include_onset_kernel=onset_on,
            B_me_face=B_me if me_on else None,
            B_history=B_hist if hist_on else None)
        if accel_on:
            X = np.hstack([X, B_accel])
        return cross_validate_glm(X, y, offset, folds,
                                  lambda_ridge=cfg.lambda_ridge).cv_bits_per_spike

    noSpeed = [v for v in STIM if v != "Speed"]
    cv_int = cv([], False, False, False, onset_on=False)
    cv_full = cv(STIM, True, True, True)            # latest + Accel
    out = {
        "cv_intercept": cv_int, "cv_full": cv_full, "total_full": cv_full - cv_int,
        "unique_Accel": cv_full - cv(STIM, True, True, False),
        "unique_ME": cv_full - cv(STIM, False, True, True),
        "unique_History": cv_full - cv(STIM, True, False, True),
        # Speed acid sequence: does Accel further erode Speed beyond ME+History?
        "uSpeed_none": cv(STIM, False, False, False) - cv(noSpeed, False, False, False),
        "uSpeed_ME": cv(STIM, True, False, False) - cv(noSpeed, True, False, False),
        "uSpeed_MEH": cv(STIM, True, True, False) - cv(noSpeed, True, True, False),
        "uSpeed_all": cv_full - cv(noSpeed, True, True, True),
    }
    out["unique_Speed"] = out["uSpeed_all"]
    out["unique_Onset"] = cv_full - cv(STIM, True, True, True, onset_on=False)
    for g in STIM:
        if g == "Speed":
            continue
        out[f"unique_{g}"] = cv_full - cv([v for v in STIM if v != g], True, True, True)
    return out


def main() -> int:
    cfg = _make_cfg()
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
    df.to_csv(out_dir / "variance_partition_accel.csv", index=False)
    well = df["total_full"] > 0.005

    def m(k):
        return df.loc[well, k].median()
    print(f"\nn={len(df)} (well-fit {int(well.sum())}); medians Δcv-bps:")
    print(f"  unique Acceleration = {m('unique_Accel'):+.4f}  "
          f"[>0.005: {int((df.loc[well,'unique_Accel']>0.005).sum())}/{int(well.sum())}]")
    print(f"  Speed unique  none={m('uSpeed_none'):+.4f}  +ME={m('uSpeed_ME'):+.4f}  "
          f"+ME+H={m('uSpeed_MEH'):+.4f}  +all(+Accel)={m('uSpeed_all'):+.4f}")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    groups = [("Onset", "unique_Onset", "tab:blue"), ("Speed", "unique_Speed", "tab:green"),
              ("TF", "unique_TF", "tab:orange"), ("SF", "unique_SF", "tab:olive"),
              ("OR", "unique_OR", "tab:red"), ("ME", "unique_ME", "tab:purple"),
              ("History", "unique_History", "tab:brown"), ("Accel", "unique_Accel", "tab:cyan")]
    data = [df.loc[well, c].to_numpy() for _, c, _ in groups]
    bp = axL.boxplot(data, tick_labels=[g for g, _, _ in groups], showfliers=False,
                     patch_artist=True)
    for patch, (_, _, c) in zip(bp["boxes"], groups):
        patch.set_facecolor(c); patch.set_alpha(0.45)
    axL.axhline(0, color="0.6", lw=0.8); axL.axhline(0.005, color="0.6", lw=0.8, ls=":")
    axL.set_ylabel(f"unique Δ cv-bps (full {BIN_MS} ms model)")
    axL.tick_params(axis="x", labelrotation=20)
    axL.set_title(f"Unique per predictor, {BIN_MS} ms +ME +refrac-History +Accel (n={int(well.sum())})",
                  fontsize=10, fontweight="bold")

    conds = [("none", "uSpeed_none"), ("+ME", "uSpeed_ME"),
             ("+ME+H", "uSpeed_MEH"), ("+all(+Accel)", "uSpeed_all")]
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
    axR.set_title(f"Acid test: does Speed survive Acceleration?\n"
                  f"(unique Accel median {m('unique_Accel'):+.3f})",
                  fontsize=10, fontweight="bold")

    fig.suptitle(f"Acceleration acid test ({RUN_NAME}; "
                 "Speed vs Acceleration attribution)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = out_dir / "variance_partition_accel"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
