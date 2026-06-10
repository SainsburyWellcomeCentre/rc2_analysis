"""Median deviance (McFadden pseudo-R²) and its scatter across clusters,
for the current models: 100 ms +ME vs 20 ms +ME+History.

pseudo-R² = 1 - cv_full / cv_intercept (speed-profile folds). Reads the
variance-partition csvs (no refitting).

Output: current_plus_ME_20ms/diagnostics/deviance_scatter.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "local_data/motion_clouds/figures/glm"
OUT = ROOT / "current_plus_ME_20ms" / "diagnostics"


def main() -> int:
    me = pd.read_csv(ROOT / "current_plus_ME/diagnostics/variance_partition_ME.csv")
    hi = pd.read_csv(ROOT / "current_plus_ME_20ms/diagnostics/variance_partition_history.csv")
    series = [
        ("100 ms\n+ME", 1 - me["cv_full_ME"] / me["cv_intercept"], "tab:purple"),
        ("20 ms\n+ME+History", 1 - hi["cv_full"] / hi["cv_intercept"], "tab:brown"),
    ]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.6),
                                   gridspec_kw={"width_ratios": [1, 1.2]})

    # Left: strip + box per model.
    for i, (lab, r2, c) in enumerate(series):
        x = np.full(len(r2), i) + np.random.uniform(-.12, .12, len(r2))
        axL.scatter(x, r2, s=14, color=c, alpha=0.5, edgecolor="none")
        bp = axL.boxplot([r2.to_numpy()], positions=[i], widths=0.45,
                         showfliers=False, patch_artist=True)
        bp["boxes"][0].set_facecolor(c); bp["boxes"][0].set_alpha(0.25)
        axL.text(i, r2.median() + 0.005, f"med {r2.median():.3f}", ha="center",
                 fontsize=9, fontweight="bold")
    axL.axhline(0, color="0.6", lw=0.8)
    axL.set_xticks(range(len(series)))
    axL.set_xticklabels([s[0] for s in series])
    axL.set_ylabel("cv deviance explained (McFadden pseudo-R²)")
    axL.set_title("Per-cluster deviance explained", fontsize=10, fontweight="bold")

    # Right: cumulative distribution (how it scatters across clusters).
    for lab, r2, c in series:
        xs = np.sort(r2.to_numpy())
        ys = np.arange(1, len(xs) + 1) / len(xs)
        axR.plot(xs, ys, color=c, lw=2,
                 label=f"{lab.replace(chr(10),' ')} (med {r2.median():.3f}, "
                       f"{100*(r2>0).mean():.0f}% >0)")
    axR.axvline(0, color="0.6", lw=0.8)
    axR.set_xlabel("pseudo-R²"); axR.set_ylabel("cumulative fraction of clusters")
    axR.set_title("Distribution across clusters", fontsize=10, fontweight="bold")
    axR.legend(fontsize=8, frameon=False, loc="lower right")

    fig.suptitle("Current models — deviance explained (low because single-trial "
                 "Poisson V1 is noisy; long right tail = the well-fit cells)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "deviance_scatter"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    for lab, r2, c in series:
        print(f"{lab.replace(chr(10),' '):20s} median {r2.median():.4f}  "
              f"IQR [{r2.quantile(.25):.4f}, {r2.quantile(.75):.4f}]  "
              f"{100*(r2>0).mean():.0f}% >0")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
