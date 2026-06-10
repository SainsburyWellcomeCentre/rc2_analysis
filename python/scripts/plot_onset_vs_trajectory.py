"""Decompose the FR-ME mean-trace coupling into the shared ONSET STEP vs the
within-motion TRAJECTORY — the part that most changes the interpretation.

The full-window mean-trace r (~0.84 in VT) is high mostly because FR and ME both
STEP UP at motion onset. Removing that step (centring within the motion window)
leaves the trajectory co-variation, ~0.21 — which is the ME-within ↔ Speed
signal. This figure shows the two pieces side by side, for VT and T_Vstatic.

Output: current_plus_ME/diagnostics/onset_vs_trajectory.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.fr_me_correlation_lohuis_style import (
    BIN_WIDTH, CACHE_ROOT, PROBES, _load_cache, _mean_trace_r,
    _stack_multi_animal)
import scripts.run_glm_split_by_condition as drv

OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME" / "diagnostics"
CONDS = ["VT", "T_Vstatic"]


def _r(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan
    x = x[ok] - x[ok].mean(); y = y[ok] - y[ok].mean()
    d = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / d) if d > 0 else np.nan


def main() -> int:
    pbs = [_load_cache(CACHE_ROOT / f"{p}.npz") for p in PROBES]
    fig, axes = plt.subplots(3, len(CONDS), figsize=(12, 10.2))
    for ci, cond in enumerate(CONDS):
        fr_z, me_z, _, _, off, _ = _stack_multi_animal(pbs, cond)
        t = off.astype(float) * BIN_WIDTH
        mfr = np.nanmean(fr_z, axis=0); mme = np.nanmean(me_z, axis=0)
        motion = off >= 0
        pre = off < 0
        r_full = _mean_trace_r(fr_z, me_z)
        # within-motion residual = motion trace minus its own post-onset mean
        fr_post = np.nanmean(mfr[motion]); me_post = np.nanmean(mme[motion])
        fr_res = mfr.copy(); me_res = mme.copy()
        fr_res[motion] -= fr_post; me_res[motion] -= me_post
        fr_res[pre] -= np.nanmean(mfr[pre]); me_res[pre] -= np.nanmean(mme[pre])
        r_motion = _r(mfr[motion], mme[motion])
        r_base = _r(mfr[pre], mme[pre])

        # Top: full mean traces with baseline level + onset step marked.
        axT = axes[0, ci]
        axT.plot(t, mfr, color="black", lw=1.8, label="mean FR (z)")
        axT.plot(t, mme, color="C3", lw=1.8, label="mean ME (z)")
        axT.axvline(0, color="0.5", ls="--", lw=0.9)
        axT.axvspan(t[motion][0], t[motion][-1], color="0.93", zorder=0)
        for m, c in [(mfr, "black"), (mme, "C3")]:
            axT.hlines(np.nanmean(m[pre]), t[0], 0, color=c, ls=":", lw=1.2)
            axT.hlines(np.nanmean(m[motion]), 0, t[motion][-1], color=c, ls=":", lw=1.2)
        axT.set_xlim(t[0], t[-1])
        axT.set_title(f"{cond}: full mean traces — onset STEP dominates\n"
                      f"full mean-trace r = {r_full:+.2f}",
                      fontsize=10, fontweight="bold")
        if ci == 0:
            axT.set_ylabel("z-score"); axT.legend(fontsize=8, frameon=False)

        # Bottom: within-motion residual (step removed) = the trajectory.
        axB = axes[1, ci]
        axB.plot(t[motion], fr_res[motion], color="black", lw=1.8)
        axB.plot(t[motion], me_res[motion], color="C3", lw=1.8)
        axB.axhline(0, color="0.6", lw=0.8)
        axB.set_xlim(t[motion][0], t[motion][-1])
        axB.set_xlabel("time from motion onset (s)")
        axB.set_title(f"{cond}: within-motion TRAJECTORY (step removed)\n"
                      f"trajectory mean-trace r = {r_motion:+.2f}  (= ME-within↔Speed)",
                      fontsize=10, fontweight="bold")
        if ci == 0:
            axB.set_ylabel("z (centred on post-onset mean)")

        # Third row: BASELINE window residual (step removed), same y-scale idea.
        axb = axes[2, ci]
        axb.plot(t[pre], fr_res[pre], color="black", lw=1.8)
        axb.plot(t[pre], me_res[pre], color="C3", lw=1.8)
        axb.axhline(0, color="0.6", lw=0.8)
        axb.set_xlim(t[pre][0], t[pre][-1])
        axb.set_xlabel("time from motion onset (s)  [baseline, pre-onset]")
        axb.set_title(f"{cond}: BASELINE window (pre-onset, centred)\n"
                      f"baseline mean-trace r = {r_base:+.2f}",
                      fontsize=10, fontweight="bold")
        if ci == 0:
            axb.set_ylabel("z (centred on pre-onset mean)")

    fig.suptitle("FR-ME coupling decomposed: shared ONSET STEP (big) vs "
                 "within-motion TRAJECTORY (modest = ME-within)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "onset_vs_trajectory"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
