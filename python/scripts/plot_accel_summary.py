"""One-glance summary of the Speed-vs-Acceleration finding.

From variance_partition_accel.csv (no refit):
  speed_only  = uSpeed_all           (Speed unique WITH acceleration present)
  shared      = uSpeed_MEH - uSpeed_all  (Speed's unique that Accel absorbs = Speed∩Accel)
  accel_only  = unique_Accel         (Acceleration unique WITH speed present)

Left:  population deviance split of the Speed+Accel joint contribution
       (Speed-only | shared | Accel-only) — Accel-only ≈ 0.
Right: per-cell unique Speed vs unique Accel (full model) — the cloud hugs the
       Speed axis: Speed carries unique deviance, Acceleration ≈ 0 regardless.

Output: current_plus_ME_20ms_accel/diagnostics/accel_summary.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys

_G = Path.home() / "local_data/motion_clouds/figures/glm"
ROOT = _G / ("current_ME_hist_accel_10ms" if "--bin10" in sys.argv
             else "current_plus_ME_20ms_accel")
CSV = ROOT / "diagnostics" / "variance_partition_accel.csv"


def main() -> int:
    df = pd.read_csv(CSV)
    df = df[df["total_full"] > 0.005].copy()
    df["speed_only"] = df["uSpeed_all"]
    df["shared"] = df["uSpeed_MEH"] - df["uSpeed_all"]
    df["accel_only"] = df["unique_Accel"]
    n = len(df)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8),
                                   gridspec_kw={"width_ratios": [1, 1]})

    # Left: stacked deviance split (medians, clipped at 0 for the bar).
    parts = [("Speed-only", df["speed_only"].median(), "tab:green"),
             ("shared\nSpeed∩Accel", df["shared"].median(), "0.6"),
             ("Accel-only", df["accel_only"].median(), "tab:cyan")]
    left = 0.0
    for name, val, c in parts:
        v = max(val, 0.0)
        axL.barh(0, v, left=left, height=0.5, color=c, edgecolor="white")
        if v > 0.0008:
            axL.text(left + v / 2, 0, f"{name}\n{val:+.4f}", ha="center", va="center",
                     fontsize=9, color="white" if c != "0.6" else "black",
                     fontweight="bold")
        else:
            axL.text(left, -0.42, f"{name} {val:+.4f}", ha="left", va="top",
                     fontsize=8, color="0.3")
        left += v
    axL.set_ylim(-0.6, 0.6); axL.set_yticks([])
    axL.set_xlabel("median Δ cv-bps (Speed + Acceleration joint contribution)")
    axL.set_title("Acceleration's contribution is ~entirely SHARED with Speed\n"
                  "— it has no own part; Speed has shared + a unique part",
                  fontsize=10, fontweight="bold")

    # Right: per-cell unique Speed vs unique Accel (full model).
    x = df["unique_Speed"].to_numpy(); y = df["unique_Accel"].to_numpy()
    axR.scatter(x, y, s=20, color="tab:cyan", alpha=0.6, edgecolor="none")
    lim = [min(x.min(), y.min(), 0) * 1.05, max(x.max(), y.max(), 0.02) * 1.05]
    axR.plot(lim, lim, ls="--", color="0.5", lw=1)
    axR.axhline(0, color="0.7", lw=0.8); axR.axvline(0, color="0.7", lw=0.8)
    axR.axhline(0.005, color="0.8", lw=0.7, ls=":")
    axR.set_xlabel("unique Speed (Δ cv-bps)")
    axR.set_ylabel("unique Acceleration (Δ cv-bps)")
    ns, na = int((x > 0.005).sum()), int((y > 0.005).sum())
    axR.set_title(f"Per-cell unique (full model, n={n})\n"
                  f"Speed>0.005: {ns}  ·  Accel>0.005: {na}  "
                  f"(cloud hugs the Speed axis)", fontsize=10, fontweight="bold")

    fig.suptitle("Speed vs Acceleration — acceleration is redundant; velocity is "
                 "the carrier (current_plus_ME_20ms_accel)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = ROOT / "diagnostics" / "accel_summary"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"n={n}  speed_only={df['speed_only'].median():+.4f}  "
          f"shared={df['shared'].median():+.4f}  accel_only={df['accel_only'].median():+.4f}")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
