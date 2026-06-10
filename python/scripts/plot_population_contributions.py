"""Population view: per-cluster deviance + acid-test (unique contributions),
color-coded by predictor. From the 20 ms +ME+History variance partition.

Left:  per-cluster stacked unique cv-bps (one bar per cluster, sorted by total
       deviance), colored by predictor; black line = McFadden pseudo-R². Shows
       what carries the explainable structure as you move up the population.
Right: per-predictor unique-contribution distribution across clusters (the acid
       test, population summary).

Negative uniques are clipped to 0 in the stack (small overfit wiggles); the
boxplot keeps the raw values.

Output: current_plus_ME_20ms/diagnostics/population_contributions.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "local_data/motion_clouds/figures/glm/current_plus_ME_20ms"
CSV = ROOT / "diagnostics" / "variance_partition_history.csv"

# Consistent colours used across the session's figures.
VARS = [("Onset", "unique_Onset", "tab:blue"), ("Speed", "unique_Speed", "tab:green"),
        ("TF", "unique_TF", "tab:orange"), ("SF", "unique_SF", "tab:olive"),
        ("OR", "unique_OR", "tab:red"), ("ME", "unique_ME", "tab:purple"),
        ("History", "unique_History", "tab:brown")]


def main() -> int:
    df = pd.read_csv(CSV)
    df = df[df["total_full"] > 0.005].copy()        # well-fit cells
    df["pseudoR2"] = 1 - df["cv_full"] / df["cv_intercept"]
    df = df.sort_values("pseudoR2").reset_index(drop=True)
    n = len(df)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 5.2),
                                   gridspec_kw={"width_ratios": [1.7, 1]})

    # Left: stacked unique contributions per cluster (clip neg), + deviance line.
    x = np.arange(n)
    bottom = np.zeros(n)
    for lab, col, c in VARS:
        v = np.clip(df[col].to_numpy(), 0, None)
        axL.bar(x, v, bottom=bottom, width=1.0, color=c, label=lab, linewidth=0)
        bottom += v
    axL.set_xlim(-0.5, n - 0.5)
    axL.set_xlabel("cluster (sorted by deviance explained)")
    axL.set_ylabel("stacked unique Δ cv-bps")
    axL.legend(ncol=4, fontsize=8, frameon=False, loc="upper left")
    ax2 = axL.twinx()
    ax2.plot(x, df["pseudoR2"], color="black", lw=1.6)
    ax2.set_ylabel("McFadden pseudo-R² (black line)")
    ax2.set_ylim(0, max(0.05, df["pseudoR2"].max() * 1.05))
    axL.set_title(f"Per-cluster contribution composition (n={n}, 20 ms +ME+History)",
                  fontsize=11, fontweight="bold")

    # Right: per-predictor unique distribution (acid test, population).
    data = [df[col].to_numpy() for _, col, _ in VARS]
    bp = axR.boxplot(data, vert=True, tick_labels=[v[0] for v in VARS],
                     showfliers=False, patch_artist=True, widths=0.6)
    for patch, (_, _, c) in zip(bp["boxes"], VARS):
        patch.set_facecolor(c); patch.set_alpha(0.5)
    axR.axhline(0, color="0.6", lw=0.8); axR.axhline(0.005, color="0.6", lw=0.7, ls=":")
    axR.set_ylabel("unique Δ cv-bps")
    axR.tick_params(axis="x", labelrotation=30)
    axR.set_title("Unique contribution per predictor\n(acid test, population)",
                  fontsize=11, fontweight="bold")

    fig.suptitle("Population contributions — deviance + acid test per predictor "
                 "(full additive model, leave-one-out)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = ROOT / "diagnostics" / "population_contributions"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"n={n}; median pseudo-R²={df['pseudoR2'].median():.3f}")
    for lab, col, _ in VARS:
        print(f"  {lab:8s} unique median {df[col].median():+.4f} "
              f"(>0.005: {int((df[col]>0.005).sum())}/{n})")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
