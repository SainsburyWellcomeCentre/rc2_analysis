"""cv-bps variance partitioning for the motion-clouds GLM.

Decomposes each cluster's cross-validated explained performance (bits/spike,
under the speed-profile leak-test folds) into the UNIQUE contribution of each
predictor group and the SHARED (redundant) variance between Onset and the
stimulus block — to answer "is Speed/TF cosmetic, or does it add something
Onset can't?".

Per cluster, refit + CV these nested models (intercept always in):
  intercept            : baseline B0
  onset                : + onset kernel            (= "Null"; onset block alone)
  stim                 : + Speed/TF/SF/OR, NO onset (stimulus block alone)
  full                 : onset + stimulus
  full_no_<g>          : full minus one stimulus group g

Unique_g  = cv_full - cv_full_no_g          (g in stimulus groups)
Unique_onset_block = cv_full - cv_stim
Unique_stim_block  = cv_full - cv_onset
Shared(onset, stim) = (cv_onset - B0) + (cv_stim - B0) - (cv_full - B0)
                    =  cv_onset + cv_stim - cv_full - B0

The big SHARED block is expected (Onset and the fixed-trajectory stimulus are
collinear): partitioning *quantifies* the confound; the UNIQUE terms are the
trustworthy part.

Default: pooled model (all conditions VT+T_Vstatic+V together → current/).
Pass --condition {V,T_Vstatic,VT} for the per-condition variant later.

Output: <run>/diagnostics/variance_partition_<scope>.{csv,pdf}
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import onset_kernel_basis, value_basis
from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import _subset_to_condition
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

STIM_BY_SCOPE = {
    "pooled": ("Speed", "TF", "SF", "OR"),
    "VT": ("Speed", "TF", "SF", "OR"),
    "V": ("TF", "SF", "OR"),
    "T_Vstatic": ("Speed",),
}


def _cvbps(X, y, offset, fold_ids, lam):
    return cross_validate_glm(X, y, offset, fold_ids, lambda_ridge=lam).cv_bits_per_spike


def partition_cluster(df, cfg, stim_vars) -> dict | None:
    motion = df[df["condition"] != "stationary"]
    if motion.empty or df["spike_count"].sum() == 0:
        return None
    speed = df["speed"].to_numpy(float)
    tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float)
    orient = df["orientation"].to_numpy(float)
    y = df["spike_count"].to_numpy(float)
    trial_ids = df["trial_id"].to_numpy(np.int64)
    profile_ids = df["profile_id"].to_numpy(np.int64)
    cond_labels = df["condition"].to_numpy(object)
    if len(np.unique(profile_ids)) < 2:
        return None
    offset = float(np.log(cfg.time_bin_width))
    sp = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=sp)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=sp)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    folds = make_trial_folds(trial_ids, cfg.n_folds, cfg.cv_seed,
                             condition_labels_per_bin=cond_labels,
                             strategy="speed-profile", profile_ids_per_bin=profile_ids)

    def cv(vars_, onset_on):
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, list(vars_),
            include_onset_kernel=onset_on,
        )
        return _cvbps(X, y, offset, folds, cfg.lambda_ridge)

    stim = list(stim_vars)
    cv_intercept = cv([], False)
    cv_onset = cv([], True)
    cv_stim = cv(stim, False)
    cv_full = cv(stim, True)
    out = {
        "cv_intercept": cv_intercept, "cv_onset": cv_onset,
        "cv_stim": cv_stim, "cv_full": cv_full,
        "total_over_intercept": cv_full - cv_intercept,
        "unique_onset_block": cv_full - cv_stim,
        "unique_stim_block": cv_full - cv_onset,
        "shared_onset_stim": cv_onset + cv_stim - cv_full - cv_intercept,
    }
    for g in stim:
        out[f"unique_{g}"] = cv_full - cv([v for v in stim if v != g], True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--condition", choices=("pooled", "V", "T_Vstatic", "VT"),
                    default="pooled")
    scope = ap.parse_args().condition
    stim_vars = STIM_BY_SCOPE[scope]
    cond = None if scope == "pooled" else scope
    cfg = drv.make_config(cond)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    out_root = drv.POOLED_OUT if scope == "pooled" else drv.OUT_ROOT

    rows = []
    for probe in drv.PROBES:
        cohort = drv.load_cohort(probe)
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        for cl in pdata.clusters:
            if cl.cluster_id not in cohort:
                continue
            df = bin_cluster(pdata, cl)
            if cond is not None:
                df = _subset_to_condition(df, cond)
            if df.empty:
                continue
            res = partition_cluster(df, cfg, stim_vars)
            if res is None:
                continue
            res.update(probe_id=probe, cluster_id=int(cl.cluster_id))
            rows.append(res)
        print(f"{probe}: done ({sum(r['probe_id']==probe for r in rows)} clusters)")

    df = pd.DataFrame(rows)
    out_dir = out_root / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv = out_dir / f"variance_partition_{scope}.csv"
    df.to_csv(csv, index=False)
    print(f"wrote {csv} ({len(df)} clusters)")

    # Population summary.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    # Panel 1: Onset vs Stimulus block — unique + shared (relative to intercept).
    keys1 = ["unique_onset_block", "shared_onset_stim", "unique_stim_block"]
    labels1 = ["unique\nOnset", "shared\nOnset∩Stim", "unique\nStim"]
    data1 = [df[k].to_numpy() for k in keys1]
    axes[0].boxplot(data1, labels=labels1, showfliers=False)
    for i, d in enumerate(data1, 1):
        axes[0].scatter(np.full(len(d), i) + np.random.uniform(-.12, .12, len(d)),
                        d, s=8, alpha=0.4, color="0.3")
    axes[0].axhline(0, color="0.7", lw=0.8)
    axes[0].set_ylabel("Δ cv-bps (vs intercept)")
    axes[0].set_title(f"Onset vs Stimulus block — {scope}\n"
                      f"(n={len(df)}; medians: "
                      f"uOn={df['unique_onset_block'].median():.3f} "
                      f"sh={df['shared_onset_stim'].median():.3f} "
                      f"uSt={df['unique_stim_block'].median():.3f})")
    # Panel 2: per-group unique variance.
    gkeys = ["unique_onset_block"] + [f"unique_{g}" for g in stim_vars]
    glabels = ["Onset"] + list(stim_vars)
    gdata = [df[k].to_numpy() for k in gkeys]
    axes[1].boxplot(gdata, labels=glabels, showfliers=False)
    for i, d in enumerate(gdata, 1):
        axes[1].scatter(np.full(len(d), i) + np.random.uniform(-.12, .12, len(d)),
                        d, s=8, alpha=0.4, color="0.3")
    axes[1].axhline(0, color="0.7", lw=0.8)
    axes[1].set_ylabel("unique Δ cv-bps")
    axes[1].set_title("Unique variance per group (leave-one-out)")
    fig.suptitle(f"cv-bps variance partition ({scope}, speed-profile folds)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    pdf = out_dir / f"variance_partition_{scope}.pdf"
    fig.savefig(pdf)
    plt.close(fig)
    print(f"wrote {pdf}")

    # Headline composition figure: median partition as a stacked bar +
    # per-group unique medians — the communicative version.
    fig2, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))
    med = {k: float(df[k].median()) for k in
           ("unique_onset_block", "shared_onset_stim", "unique_stim_block")}
    parts = [("unique Onset", med["unique_onset_block"], "tab:blue"),
             ("shared Onset∩Stim", med["shared_onset_stim"], "0.6"),
             ("unique Stimulus", med["unique_stim_block"], "tab:red")]
    total = sum(v for _, v, _ in parts)
    left = 0.0
    for name, val, color in parts:
        axL.barh(0, val, left=left, height=0.5, color=color, edgecolor="white")
        if val > 0.002:
            pct = 100 * val / total if total else 0
            axL.text(left + val / 2, 0, f"{val:.3f}\n({pct:.0f}%)", ha="center",
                     va="center", fontsize=9,
                     color="white" if color != "0.6" else "black")
        left += val
    axL.set_yticks([])
    axL.set_xlabel("median Δ cv-bps (vs intercept)")
    axL.set_title(f"{scope}: explained variance partition\n(median across n={len(df)})")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in parts]
    axL.legend(handles, [n for n, _, _ in parts], loc="upper center",
               bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=8, frameon=False)
    groups = ["Onset"] + list(stim_vars)
    gmed = [float(df["unique_onset_block"].median())] + \
           [float(df[f"unique_{g}"].median()) for g in stim_vars]
    axR.bar(groups, gmed, color=["tab:blue"] + ["tab:red"] * len(stim_vars))
    axR.axhline(0, color="0.7", lw=0.8)
    axR.set_ylabel("median unique Δ cv-bps")
    axR.set_title("Unique variance per group (median, leave-one-out)")
    for i, v in enumerate(gmed):
        axR.text(i, v, f"{v:.3f}", ha="center",
                 va="bottom" if v >= 0 else "top", fontsize=8)
    fig2.suptitle(f"Variance partition — {scope} (speed-profile folds): "
                  f"Onset's UNIQUE share ≈ {med['unique_onset_block']:.3f}; "
                  f"Stimulus carries the unique variance",
                  fontsize=11, fontweight="bold")
    fig2.tight_layout()
    pdf2 = out_dir / f"variance_partition_{scope}_composition.pdf"
    fig2.savefig(pdf2)
    plt.close(fig2)
    print(f"wrote {pdf2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
