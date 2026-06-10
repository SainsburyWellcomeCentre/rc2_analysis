"""Per-cluster acid test + deviance for the full model (current_full_20ms).

Left:  per-cluster stacked unique cv-bps (acid test), colored by predictor,
       sorted by deviance; black line = pseudo-R². Population dominates.
Right: per-cluster pseudo-R² WITHOUT vs WITH the population term — shows the
       deviance roughly doubling when the LOO population is included.

Reads variance_partition_full.csv (no refit). Negatives clipped to 0 in the stack.
Output: current_full_20ms/diagnostics/full_per_cluster.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "local_data/motion_clouds/figures/glm/current_full_20ms"
CSV = ROOT / "diagnostics" / "variance_partition_full.csv"

VARS = [("Onset", "unique_Onset", "tab:blue"), ("Speed", "unique_Speed", "tab:green"),
        ("TF", "unique_TF", "tab:orange"), ("SF", "unique_SF", "tab:olive"),
        ("OR", "unique_OR", "tab:red"), ("ME", "unique_ME", "tab:purple"),
        ("History", "unique_History", "tab:brown"), ("ME×Speed", "unique_MExS", "tab:gray"),
        ("Accel", "unique_Accel", "tab:cyan"), ("Population", "unique_Pop", "tab:pink")]


def main() -> int:
    df = pd.read_csv(CSV)
    df = df[df["total_full"] > 0.005].copy()
    df["r2"] = 1 - df["cv_full"] / df["cv_intercept"]
    df["r2_noPop"] = 1 - df["cv_noPop"] / df["cv_intercept"]
    df = df.sort_values("r2").reset_index(drop=True)
    n = len(df)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 5.2),
                                   gridspec_kw={"width_ratios": [1.7, 1]})

    x = np.arange(n); bottom = np.zeros(n)
    for lab, col, c in VARS:
        v = np.clip(df[col].to_numpy(), 0, None)
        axL.bar(x, v, bottom=bottom, width=1.0, color=c, label=lab, linewidth=0)
        bottom += v
    axL.set_xlim(-0.5, n - 0.5)
    axL.set_xlabel("cluster (sorted by deviance)")
    axL.set_ylabel("stacked unique Δ cv-bps")
    axL.legend(ncol=5, fontsize=7.5, frameon=False, loc="upper left")
    ax2 = axL.twinx(); ax2.plot(x, df["r2"], color="black", lw=1.6)
    ax2.set_ylabel("pseudo-R² (black line)"); ax2.set_ylim(0, max(0.05, df["r2"].max() * 1.05))
    axL.set_title(f"Per-cluster contribution (acid test), full model (n={n})",
                  fontsize=11, fontweight="bold")

    a = df["r2_noPop"].to_numpy(); b = df["r2"].to_numpy()
    axR.scatter(a, b, s=20, color="tab:pink", alpha=0.65, edgecolor="none")
    lim = [min(a.min(), b.min(), 0) * 1.05, max(a.max(), b.max()) * 1.05]
    axR.plot(lim, lim, ls="--", color="0.5", lw=1)
    axR.axhline(0, color="0.85", lw=0.6); axR.axvline(0, color="0.85", lw=0.6)
    axR.set_xlabel("pseudo-R²  WITHOUT population")
    axR.set_ylabel("pseudo-R²  WITH population")
    axR.set_title(f"Population doubles deviance explained\n"
                  f"median {np.median(a):.3f} → {np.median(b):.3f}",
                  fontsize=11, fontweight="bold")

    fig.suptitle("Full model — per-cluster acid test & deviance "
                 "(current_full_20ms; +Accel +Population +ME×Speed)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = ROOT / "diagnostics" / "full_per_cluster"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"n={n}; pseudo-R² noPop {np.median(a):.3f} -> withPop {np.median(b):.3f}")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
