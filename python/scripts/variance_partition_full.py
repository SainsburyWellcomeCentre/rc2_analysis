"""Full consolidated partition — Speed/TF/SF/OR + ME + refractory History +
Acceleration + LOO Population + ME×Speed. Acid test + deviance per cluster.

Population = leave-one-out mean spike rate of all OTHER simultaneously-recorded
VISp cells on the probe, same 20 ms bin, z-scored, ONE additive column
(log-link → shared multiplicative gain; Okun/Stringer population-coupling).
Computed per probe from the binned spike-count matrix of every VISp cluster
(bins are velocity-determined, identical across clusters), LOO per target.

Per cohort cluster: unique cv-bps per predictor (leave-one-out in the full
model), Speed acid sequence for the NEW terms (base → +Accel → +Pop → +both),
and pseudo-R² with/without Population (the deviance-ceiling question).

Output: current_full_20ms/diagnostics/variance_partition_full.{csv,pdf,png}
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

import scripts.run_glm_split_by_condition as drv
from scripts.run_glm_current_plus_ME_20ms import make_config_me_20ms

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_full_20ms"
STIM = ["Speed", "TF", "SF", "OR"]
HISTORY_WINDOW = 0.04


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


def partition(df, cfg, pop_loo):
    if df["spike_count"].sum() == 0:
        return None
    B_me = _build_me(df, cfg)
    if B_me is None:
        return None
    B_accel = raised_cosine_basis_linear(np.clip(_zfill(df["acceleration"].to_numpy(float)), -3, 3),
                                         5, -3.0, 3.0)
    B_pop = _zfill(pop_loo).reshape(-1, 1)          # single additive column
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
    spc = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=spc)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=spc)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    hb = history_basis(cfg.n_history_bases, cfg.history_window_s, cfg.time_bin_width,
                       kind=cfg.history_basis_kind)
    B_hist = convolve_history(y, trial_ids, hb)
    folds = make_trial_folds(trial_ids, cfg.n_folds, cfg.cv_seed,
                             condition_labels_per_bin=cond,
                             strategy="speed-profile", profile_ids_per_bin=profile_ids)

    def cv(stim_vars, me, hist, mexs, accel, pop, onset_on=True):
        sel = list(stim_vars)
        if me:
            sel = sel + ["ME_face"]
        if hist:
            sel = sel + ["History"]
        if mexs:
            sel = sel + ["ME_face_x_Speed"]
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, sel, include_onset_kernel=onset_on,
            B_me_face=B_me if (me or mexs) else None,
            B_history=B_hist if hist else None)
        extra = []
        if accel:
            extra.append(B_accel)
        if pop:
            extra.append(B_pop)
        if extra:
            X = np.hstack([X] + extra)
        return cross_validate_glm(X, y, offset, folds,
                                  lambda_ridge=cfg.lambda_ridge).cv_bits_per_spike

    noSpeed = [v for v in STIM if v != "Speed"]
    cv_int = cv([], False, False, False, False, False, onset_on=False)
    full = cv(STIM, True, True, True, True, True)
    out = {
        "cv_intercept": cv_int, "cv_full": full, "total_full": full - cv_int,
        # full-model deviance with vs without Population (ceiling question)
        "cv_noPop": cv(STIM, True, True, True, True, False),
        # per-group unique (leave-one-out in the full model)
        "unique_Onset": full - cv(STIM, True, True, True, True, True, onset_on=False),
        "unique_Speed": full - cv(noSpeed, True, True, False, True, True),   # drop Speed + its interaction
        "unique_ME": full - cv(STIM, False, True, False, True, True),        # drop ME + interaction
        "unique_History": full - cv(STIM, True, False, True, True, True),
        "unique_MExS": full - cv(STIM, True, True, False, True, True),       # drop interaction only
        "unique_Accel": full - cv(STIM, True, True, True, False, True),
        "unique_Pop": full - cv(STIM, True, True, True, True, False),
    }
    for g in ("TF", "SF", "OR"):
        out[f"unique_{g}"] = full - cv([v for v in STIM if v != g], True, True, True, True, True)
    # Speed acid sequence for the NEW terms (base = ME+Hist+MExS, no Accel/Pop)
    def uspeed(accel, pop):
        return (cv(STIM, True, True, True, accel, pop)
                - cv(noSpeed, True, True, False, accel, pop))
    out["uSpeed_base"] = uspeed(False, False)
    out["uSpeed_Accel"] = uspeed(True, False)
    out["uSpeed_Pop"] = uspeed(False, True)
    out["uSpeed_all"] = uspeed(True, True)
    return out


def main() -> int:
    cfg = make_config_me_20ms(history_window_s=HISTORY_WINDOW)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    rows = []
    for probe in PROBES:
        cohort = drv.load_cohort(probe)
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        # Bin every VISp cluster once; bins are velocity-determined (identical
        # across clusters) so spike-count vectors align row-for-row.
        binned = {}
        for cl in pdata.clusters:
            d = bin_cluster(pdata, cl)
            binned[cl.cluster_id] = d
        lengths = {len(d) for d in binned.values()}
        if len(lengths) != 1:
            print(f"  [{probe}] WARNING: bin lengths differ {lengths} — skipping pop alignment check")
        ref_len = max(lengths)
        counts = np.zeros(ref_len, dtype=np.float64)
        usable = {cid: d for cid, d in binned.items() if len(d) == ref_len}
        for cid, d in usable.items():
            counts += d["spike_count"].to_numpy(float)
        print(f"  [{probe}] population from {len(usable)}/{len(binned)} VISp clusters, {ref_len} bins")
        for cl in pdata.clusters:
            if cl.cluster_id not in cohort or cl.cluster_id not in usable:
                continue
            df = usable[cl.cluster_id]
            pop_loo = counts - df["spike_count"].to_numpy(float)
            res = partition(df, cfg, pop_loo)
            if res is None:
                continue
            res.update(probe_id=probe, cluster_id=int(cl.cluster_id))
            rows.append(res)
        print(f"{probe}: {sum(r['probe_id']==probe for r in rows)} cohort clusters")

    df = pd.DataFrame(rows)
    out_dir = OUT / "diagnostics"; out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "variance_partition_full.csv", index=False)
    well = df["total_full"] > 0.005
    r2 = 1 - df["cv_full"] / df["cv_intercept"]
    r2_noPop = 1 - df["cv_noPop"] / df["cv_intercept"]

    def m(k):
        return df.loc[well, k].median()
    print(f"\nn={len(df)} (well-fit {int(well.sum())}); medians:")
    for g in ("Speed", "ME", "History", "MExS", "Accel", "Pop", "Onset", "TF", "SF", "OR"):
        col = f"unique_{g}"
        print(f"  unique {g:8s} = {m(col):+.4f}  [>0.005: {int((df.loc[well,col]>0.005).sum())}/{int(well.sum())}]")
    print(f"  pseudo-R2 with Pop = {r2[well].median():.4f} | without Pop = {r2_noPop[well].median():.4f}")
    print(f"  Speed unique: base={m('uSpeed_base'):+.4f} +Accel={m('uSpeed_Accel'):+.4f} "
          f"+Pop={m('uSpeed_Pop'):+.4f} +both={m('uSpeed_all'):+.4f}")

    # Figure: per-group uniques (L) + Speed acid for new terms (R)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 4.8))
    groups = [("Onset", "tab:blue"), ("Speed", "tab:green"), ("TF", "tab:orange"),
              ("SF", "tab:olive"), ("OR", "tab:red"), ("ME", "tab:purple"),
              ("History", "tab:brown"), ("MExS", "tab:gray"), ("Accel", "tab:cyan"),
              ("Pop", "tab:pink")]
    data = [df.loc[well, f"unique_{g}"].to_numpy() for g, _ in groups]
    bp = axL.boxplot(data, tick_labels=[g for g, _ in groups], showfliers=False, patch_artist=True)
    for patch, (_, c) in zip(bp["boxes"], groups):
        patch.set_facecolor(c); patch.set_alpha(0.45)
    axL.axhline(0, color="0.6", lw=0.8); axL.axhline(0.005, color="0.6", lw=0.8, ls=":")
    axL.tick_params(axis="x", labelrotation=25)
    axL.set_ylabel("unique Δ cv-bps (full model)")
    axL.set_title(f"Unique per predictor — full 20 ms model (n={int(well.sum())})",
                  fontsize=10, fontweight="bold")

    conds = [("base", "uSpeed_base"), ("+Accel", "uSpeed_Accel"),
             ("+Pop", "uSpeed_Pop"), ("+both", "uSpeed_all")]
    M = np.column_stack([df.loc[well, c].to_numpy() for _, c in conds])
    x = np.arange(len(conds))
    for row in M:
        axR.plot(x, row, color="tab:green", alpha=0.12, lw=0.8)
    axR.plot(x, np.median(M, axis=0), color="black", lw=2.4, marker="o")
    for xi, md in zip(x, np.median(M, axis=0)):
        axR.annotate(f"{md:+.3f}", (xi, md), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=9, fontweight="bold")
    axR.axhline(0, color="0.6", lw=0.8)
    axR.set_xticks(x); axR.set_xticklabels([c for c, _ in conds])
    axR.set_ylabel("unique Speed (Δ cv-bps)")
    axR.set_title(f"Speed acid test vs new terms\n"
                  f"pseudo-R²: {r2_noPop[well].median():.3f} → {r2[well].median():.3f} (+Pop)",
                  fontsize=10, fontweight="bold")

    fig.suptitle("Full model partition — +Accel +Population(LOO) +ME×Speed "
                 "(current_full_20ms)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = out_dir / "variance_partition_full"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
