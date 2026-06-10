"""Illustrate why the 2026-05-28 "FR-ME decouples in motion" and today's
"ME-within tracks Speed" are consistent — single-trial vs trial-averaged.

Reads the fr_me_corr_cache (built by fr_me_correlation_lohuis_style.py). Motion
condition = VT, pooled across the 3 face-camera probes.

A. Trial-AVERAGED FR & ME (z) vs time from onset — they track (mean-trace r high).
B. A few SINGLE trials (one cluster) — FR & ME are noisy, per-trial r low.
C. Summary of the r quantities that reconcile old & new: single-trial coupling
   drops baseline→motion (noise), trial-averaged coupling stays high (shared
   trajectory = today's ME-within); genuine across-trial state coupling is small.

Output: current_plus_ME/diagnostics/fr_me_reconciliation.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scripts.run_glm_split_by_condition as drv
from scripts.fr_me_correlation_lohuis_style import (
    BIN_WIDTH, CACHE_ROOT, PROBES, _cell_trial_indices, _load_cache,
    _mean_trace_r, _pearson_r_per_bin, _per_cluster_window_r,
    _stack_cell_arrays, _stack_multi_animal)

OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME" / "diagnostics"
COND = "VT"            # a platform-motion condition (stage moves the animal)
# today's GLM numbers (variance_partition_me_split.py), for annotation
GLM_ME_ACROSS, GLM_ME_WITHIN = 0.003, 0.013


def main() -> int:
    pbs = [_load_cache(CACHE_ROOT / f"{p}.npz") for p in PROBES]

    # ---- pooled stacks for VT (motion condition) ----
    fr_z, me_z, fr_raw, me_raw, off, groups = _stack_multi_animal(pbs, COND)
    t = off.astype(float) * BIN_WIDTH
    motion = off >= 0
    mean_fr = np.nanmean(fr_z, axis=0)
    mean_me = np.nanmean(me_z, axis=0)
    # mean-trace r: full window (onset step dominates) vs motion window only
    mtr_full = _mean_trace_r(fr_z, me_z)
    mtr_motion = _mean_trace_r(fr_z[:, motion], me_z[:, motion])
    # per-bin trial-residual r, averaged over motion bins
    rbin, _ = _pearson_r_per_bin(fr_z, me_z)
    resid_motion = float(np.nanmean(rbin[motion]))

    # ---- per-trial r (baseline vs motion), population over clusters ----
    pt_base, pt_motion = [], []
    for pb in pbs:
        for ci in range(pb.cluster_ids.size):
            v = _per_cluster_window_r(pb, ci, COND)
            b = v["per_trial_r_baseline"]; m = v["per_trial_r_motion"]
            if np.isfinite(b).sum() >= 3:
                pt_base.append(np.nanmedian(b))
            if np.isfinite(m).sum() >= 3:
                pt_motion.append(np.nanmedian(m))
    pt_base = np.array(pt_base); pt_motion = np.array(pt_motion)

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16, 4.6),
                                        gridspec_kw={"width_ratios": [1.1, 1.1, 1]})

    # Panel A: trial-averaged FR & ME (z) — they track.
    axA.plot(t, mean_fr, color="black", lw=1.8, label="mean FR (z)")
    axA.plot(t, mean_me, color="C3", lw=1.8, label="mean ME (z)")
    axA.axvline(0, color="0.5", ls="--", lw=0.8)
    axA.axvspan(0, t[motion][-1], color="0.92", zorder=0)
    axA.set_xlim(t[0], t[-1])
    axA.set_xlabel("time from motion onset (s)"); axA.set_ylabel("z-score")
    axA.legend(fontsize=8, frameon=False, loc="upper left")
    axA.set_title(f"Trial-AVERAGED ({COND}, pooled)\n"
                  f"mean-trace r = {mtr_full:+.2f} full (onset step) / "
                  f"{mtr_motion:+.2f} motion-only",
                  fontsize=10, fontweight="bold")

    # Panel B: single trials from one cluster — noisy.
    pb0 = pbs[0]
    ci = int(np.argmax(pb0.fr_mean))           # an active cluster
    idxs = _cell_trial_indices(pb0, COND, profile=1)
    frc, mec, offc = _stack_cell_arrays(pb0, ci, idxs)
    tc = offc.astype(float) * BIN_WIDTH
    mmask = offc >= 0
    shown, rs = 0, []
    for row in range(frc.shape[0]):
        f = frc[row]; m = mec[row]
        ok = np.isfinite(f) & np.isfinite(m) & mmask
        if ok.sum() < 8 or shown >= 3:
            continue
        rr = np.corrcoef(f[ok], m[ok])[0, 1]
        rs.append(rr)
        ax = axB if shown == 0 else axB  # overlay with vertical offset
        offset = shown * 3.0
        # z each within the shown motion window for visual overlay
        def _z(a, msk):
            a = a.copy(); s = a[msk].std() or 1.0
            return (a - np.nanmean(a[msk])) / s
        axB.plot(tc, _z(f, ok) + offset, color="black", lw=1.0)
        axB.plot(tc, _z(m, ok) + offset, color="C3", lw=1.0)
        axB.text(tc[-1], offset, f" r={rr:+.2f}", fontsize=8, va="center",
                 color="0.2")
        shown += 1
    axB.axvline(0, color="0.5", ls="--", lw=0.8)
    axB.set_xlim(tc[0], tc[-1]); axB.set_yticks([])
    axB.set_xlabel("time from motion onset (s)")
    axB.set_title(f"SINGLE trials (cluster {int(pb0.cluster_ids[ci])}, {COND})\n"
                  f"FR (black) & ME (red) noisy → per-trial r ≈ {np.median(rs):+.2f}",
                  fontsize=10, fontweight="bold")

    # Panel C: the reconciling quantities.
    labels = ["per-trial r\nBASELINE", "per-trial r\nMOTION",
              "mean-trace r\nFULL", "mean-trace r\nMOTION", "per-bin resid r\nMOTION"]
    vals = [np.median(pt_base), np.median(pt_motion), mtr_full, mtr_motion,
            resid_motion]
    colors = ["0.5", "0.5", "tab:green", "tab:olive", "tab:orange"]
    bars = axC.bar(range(len(vals)), vals, color=colors, edgecolor="white")
    for b, v in zip(bars, vals):
        axC.text(b.get_x() + b.get_width() / 2, v + (0.01 if v >= 0 else -0.03),
                 f"{v:+.2f}", ha="center", va="bottom" if v >= 0 else "top",
                 fontsize=9, fontweight="bold")
    axC.axhline(0, color="0.6", lw=0.8)
    axC.set_xticks(range(len(vals))); axC.set_xticklabels(labels, fontsize=7.5)
    axC.set_ylabel("Pearson r")
    axC.set_title("Single-trial coupling drops in motion (noise);\n"
                  "trial-averaged stays high (shared trajectory)",
                  fontsize=10, fontweight="bold")
    axC.text(0.5, -0.32,
             f"today's GLM (motion bins): unique ME-across +{GLM_ME_ACROSS:.3f}, "
             f"ME-within +{GLM_ME_WITHIN:.3f} cv-bps\n"
             "(ME-within ≈ the high mean-trace r; ME-across ≈ the small residual r)",
             transform=axC.transAxes, ha="center", va="top", fontsize=8,
             color="0.25", bbox=dict(facecolor="0.95", edgecolor="none", pad=4))

    fig.suptitle("FR–ME coupling: single-trial vs trial-averaged "
                 "(reconciling 2026-05-28 decoupling with today's ME-within↔Speed)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "fr_me_reconciliation"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"mean-trace r (motion) = {mtr_motion:+.3f}")
    print(f"per-bin resid r (motion) = {resid_motion:+.3f}")
    print(f"per-trial r: baseline {np.median(pt_base):+.3f} -> motion {np.median(pt_motion):+.3f}")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
