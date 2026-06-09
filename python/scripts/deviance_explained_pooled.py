"""Population deviance-explained for the full pooled GLM (cross-validated).

Converts the cv-bps numbers (from variance_partition_pooled.csv) into the
GLM analogue of "variance explained":

  full-model cross-validated McFadden pseudo-R²  =  1 - LL_full / LL_null
     = 1 - cv_full / cv_intercept   (LL = bps · spikes · ln2; spikes/ln2 cancel)

  null = intercept-only (mean-rate) model; folds = speed-profile leak-test.

Also the deviance-EXPLAINED budget (fractions of the full model's gain over
intercept): unique-Onset / shared(Onset∩Stim) / unique-Stim, with Speed's share.

NOTE: this is explained DEVIANCE (pseudo-R²), NOT Gaussian variance; it's low
because single-trial Poisson spiking is noisy. Cross-validated, so it can be
negative for overfit cells.

Output: current/diagnostics/deviance_explained_pooled.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "local_data/motion_clouds/figures/glm/current"
CSV = ROOT / "diagnostics" / "variance_partition_pooled.csv"


def main() -> int:
    df = pd.read_csv(CSV)
    n = len(df)
    # Full-model cross-validated McFadden pseudo-R² (vs intercept-only null).
    r2 = 1.0 - df["cv_full"] / df["cv_intercept"]
    well = df["total_over_intercept"] > 0.005  # cells the full model actually helps

    print(f"FULL-model cv deviance-explained (McFadden pseudo-R²), pooled, n={n}")
    print(f"  median (all)        = {r2.median():.4f}")
    print(f"  median (well-fit)   = {r2[well].median():.4f}  (n={int(well.sum())})")
    print(f"  fraction R²>0       = {(r2>0).mean():.2%}")

    # Deviance-explained budget (fractions of total), well-fit cells.
    tot = df.loc[well, "total_over_intercept"]
    frac_on = (df.loc[well, "unique_onset_block"] / tot).median()
    frac_sh = (df.loc[well, "shared_onset_stim"] / tot).median()
    frac_st = (df.loc[well, "unique_stim_block"] / tot).median()
    frac_speed = (df.loc[well, "unique_Speed"] / tot).median()

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6),
                                   gridspec_kw={"width_ratios": [1, 1.1]})

    # Panel L: full-model pseudo-R² distribution.
    axL.boxplot([r2[well].to_numpy()], tick_labels=["full model"], showfliers=False,
                widths=0.5)
    axL.scatter(np.ones(int(well.sum())) + np.random.uniform(-.12, .12, int(well.sum())),
                r2[well], s=12, alpha=0.5, color="tab:blue", edgecolor="none")
    axL.axhline(0, color="0.6", lw=0.8)
    axL.set_ylabel("cv deviance explained (McFadden pseudo-R²)")
    axL.set_title(f"Full pooled GLM — cross-validated deviance explained\n"
                  f"(speed-profile folds; well-fit n={int(well.sum())}/{n}; "
                  f"median {r2[well].median():.3f})", fontsize=10, fontweight="bold")

    # Panel R: deviance-explained budget (median, % of full-model gain).
    parts = [("unique Onset", frac_on, "tab:blue"),
             ("shared Onset∩Stim", frac_sh, "0.6"),
             ("unique Stim", frac_st, "tab:red")]
    s = sum(p[1] for p in parts)
    left = 0.0
    for name, val, color in parts:
        axR.barh(0, val, left=left, height=0.5, color=color, edgecolor="white")
        if val > 0.03:
            axR.text(left + val / 2, 0, f"{100*val/s:.0f}%", ha="center", va="center",
                     color="white" if color != "0.6" else "black", fontsize=10)
        left += val
    axR.set_yticks([]); axR.set_xlim(0, s)
    axR.set_xlabel("median share of full-model explained deviance")
    axR.set_title(f"Deviance budget (well-fit cells)\n"
                  f"of which unique-Speed ≈ {100*frac_speed/s:.0f}% of total",
                  fontsize=10, fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in parts]
    axR.legend(handles, [p[0] for p in parts], loc="upper center",
               bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=8, frameon=False)

    fig.suptitle("Full GLM explained deviance (pseudo-R²) — population, pooled model",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = ROOT / "diagnostics" / "deviance_explained_pooled"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
