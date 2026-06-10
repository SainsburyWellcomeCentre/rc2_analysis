"""Deviance progression across the model builds — shows the Population jump.

Median McFadden pseudo-R² (1 - cv_full/cv_intercept) at each stage, from the
partition CSVs. Stages 3→4 (full model without vs with Population) are the same
clusters, directly comparable; stages 1-2 are context (their own well-fit sets).

Output: current_full_20ms/diagnostics/model_progression.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

G = Path.home() / "local_data/motion_clouds/figures/glm"


def med_r2(csv, full_col, well_col):
    df = pd.read_csv(csv)
    w = df[well_col] > 0.005
    r2 = 1 - df[full_col] / df["cv_intercept"]
    return r2[w].to_numpy()


def main() -> int:
    stages = [
        ("100 ms\n+ME", med_r2(G / "current_plus_ME/diagnostics/variance_partition_ME.csv",
                               "cv_full_ME", "total_ME"), "tab:purple"),
        ("20 ms\n+ME+History", med_r2(G / "current_plus_ME_20ms/diagnostics/variance_partition_history.csv",
                                      "cv_full", "total_full"), "tab:brown"),
        ("20 ms full\n(no Pop)", med_r2(G / "current_full_20ms/diagnostics/variance_partition_full.csv",
                                        "cv_noPop", "total_full"), "0.6"),
        ("20 ms full\n+Population", med_r2(G / "current_full_20ms/diagnostics/variance_partition_full.csv",
                                           "cv_full", "total_full"), "tab:pink"),
    ]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    for i, (lab, r2, c) in enumerate(stages):
        med = np.median(r2)
        x = np.full(len(r2), i) + np.random.uniform(-.12, .12, len(r2))
        ax.scatter(x, r2, s=10, color=c, alpha=0.35, edgecolor="none")
        ax.bar(i, med, width=0.6, color=c, alpha=0.45, edgecolor=c, lw=1.5)
        ax.text(i, med + 0.004, f"{med:.3f}", ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(range(len(stages))); ax.set_xticklabels([s[0] for s in stages])
    ax.set_ylabel("cv deviance explained (median McFadden pseudo-R²)")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.annotate("", xy=(3, np.median(stages[3][1])), xytext=(2, np.median(stages[2][1])),
                arrowprops=dict(arrowstyle="->", color="tab:pink", lw=2))
    ax.text(2.5, (np.median(stages[2][1]) + np.median(stages[3][1])) / 2 + 0.01,
            "≈2× with\npopulation", color="tab:pink", fontsize=9, fontweight="bold", ha="center")
    ax.set_title("Deviance explained across model builds — the LOO population\n"
                 "term roughly doubles it (the shared-variability lever)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    out = G / "current_full_20ms/diagnostics/model_progression"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    for lab, r2, _ in stages:
        print(f"{lab.replace(chr(10),' '):22s} median pseudo-R² {np.median(r2):.4f} (n={len(r2)})")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
