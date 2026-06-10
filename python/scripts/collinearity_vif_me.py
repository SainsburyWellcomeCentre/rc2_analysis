"""Collinearity (generalised VIF) for the POOLED design with ME split into
its across-trial and within-trial components.

Adds ME-across (trial mean) and ME-within (residual) as two extra blocks to the
pooled Speed/TF/Onset/SF/OR design and reports GVIF^(1/2df) per block. The
question: does ME-within show collinearity with Speed (confirming the within-
trial Speed↔ME overlap is a design-level issue), while ME-across stays clean?

3 face-camera probes; the design is spike-independent so one representative
cohort cluster per probe suffices, GVIF medianed across probes.

Output: current_plus_ME/diagnostics/collinearity_vif_me.{pdf,png}
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
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv
from scripts.collinearity_vif import _band, BAND_COLOR, GOOD, WORRY


def _gvif(R: np.ndarray, idx: list[int]) -> float:
    """Generalised VIF (Fox & Monette) via log-determinants — stable for the
    larger ME-augmented design where naive det(R) underflows."""
    others = [i for i in range(R.shape[0]) if i not in idx]
    s_full, ld_full = np.linalg.slogdet(R)
    if s_full <= 0:                       # genuinely rank-deficient
        return np.inf
    s_b, ld_b = np.linalg.slogdet(R[np.ix_(idx, idx)])
    if others:
        s_o, ld_o = np.linalg.slogdet(R[np.ix_(others, others)])
    else:
        s_o, ld_o = 1.0, 0.0
    if s_b <= 0 or s_o <= 0:
        return np.inf
    return float(np.exp(ld_b + ld_o - ld_full))

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME" / "diagnostics"
BLOCKS = ["Speed", "TF", "Onset", "SF", "OR", "MEacross", "MEwithin"]
PREFIX = {b: b + "_" for b in BLOCKS}


def _zfill(x):
    f = np.isfinite(x)
    if f.sum() < 2:
        return np.zeros_like(x)
    mu = float(x[f].mean()); sd = float(x[f].std(ddof=0)) or 1.0
    return np.where(np.isfinite(x), (x - mu) / sd, 0.0)


def _me_components(df, cfg):
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
    n, lo_hi = cfg.n_me_face_bases, cfg.me_face_range
    return (raised_cosine_basis_linear(_zfill(trial_mean), n, *lo_hi),
            raised_cosine_basis_linear(_zfill(me - trial_mean), n, *lo_hi))


def _design(df, cfg):
    speed = df["speed"].to_numpy(float); tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float); orient = df["orientation"].to_numpy(float)
    sp = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=sp)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=sp)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf, orient, ["Speed", "TF", "SF", "OR"],
        include_onset_kernel=True, B_me_face=None)
    keep = [i for i, nm in enumerate(names) if nm != "Intercept"]
    X = X[:, keep]; names = [names[i] for i in keep]
    comp = _me_components(df, cfg)
    if comp is None:
        return None
    B_across, B_within = comp
    X = np.hstack([X, B_across, B_within])
    names += [f"MEacross_{i+1}" for i in range(B_across.shape[1])]
    names += [f"MEwithin_{i+1}" for i in range(B_within.shape[1])]
    return X, names


def main() -> int:
    cfg = drv.make_config(None)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    rows = []
    for probe in PROBES:
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        cohort = drv.load_cohort(probe)
        cl = next(c for c in pdata.clusters if c.cluster_id in cohort)
        out = _design(bin_cluster(pdata, cl), cfg)
        if out is None:
            continue
        X, names = out
        Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
        R = np.corrcoef(Xs, rowvar=False)
        rec = {"probe": probe, "cond_number": float(np.linalg.cond(Xs))}
        for b in BLOCKS:
            idx = [i for i, nm in enumerate(names) if nm.startswith(PREFIX[b])]
            if not idx:
                continue
            g = _gvif(R, idx)
            rec[f"gvif_{b}"] = g ** (1.0 / (2 * len(idx))) if np.isfinite(g) else np.inf
        rows.append(rec)

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "collinearity_vif_me.csv", index=False)

    def med(b):
        return df[f"gvif_{b}"].replace(np.inf, np.nan).median()
    print("\n=== pooled +ME-split GVIF (SE scale, GVIF^(1/2df)); 2.24 / 3.16 ===")
    print(f"cond# median = {df['cond_number'].median():.1f}")
    for b in BLOCKS:
        print(f"  {b:10s} = {med(b):.2f}")

    vals = [med(b) for b in BLOCKS]
    colors = [BAND_COLOR[_band(v)] for v in vals]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    bars = ax.bar(BLOCKS, vals, color=colors, edgecolor="white")
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 0.02, f"{v:.2f}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.axhline(GOOD, color="0.4", ls="--", lw=1)
    ax.axhline(WORRY, color="0.4", ls=":", lw=1)
    ax.text(len(BLOCKS) - 0.4, GOOD + 0.02, "worrying", fontsize=7.5, color="0.4")
    ax.text(len(BLOCKS) - 0.4, WORRY + 0.02, "severe", fontsize=7.5, color="0.4")
    ax.set_ylabel("generalised VIF (SE scale, GVIF$^{1/2df}$)")
    ax.set_ylim(0, max(WORRY + 0.3, max(vals) * 1.15))
    ax.set_title("Pooled design + ME split — collinearity per block\n"
                 "(median over 3 face-camera probes)",
                 fontsize=10, fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=BAND_COLOR[b])
               for b in ("good", "worrying", "severe")]
    ax.legend(handles, ["good (<2.24)", "worrying (2.24–3.16)", "severe (>3.16)"],
              fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    out = OUT / "collinearity_vif_me"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
