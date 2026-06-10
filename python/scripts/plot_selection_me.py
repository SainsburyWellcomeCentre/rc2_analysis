"""Forward-selection summary for the current_plus_ME run (3 face-camera probes).

Aggregate of which variables forward selection admitted, across the 78 cohort
clusters (243/244/466), pooled model + ME_face, speed-profile CV. Mirrors the
split-by-condition selection-fraction bars but for this single pooled+ME run.

Left:  fraction selecting each main effect, grouped by probe (consistency check).
Right: fraction with each interaction term (incl. ME_face×Speed).

Output: current_plus_ME/diagnostics/selection_me.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "local_data/motion_clouds/figures/glm/current_plus_ME"
CSV = ROOT / "glm_model_comparison.csv"

MAIN = [("Speed", "time_is_speed_tuned"), ("TF", "time_is_tf_tuned"),
        ("SF", "time_is_sf_tuned"), ("OR", "time_is_or_tuned"),
        ("ME_face", "time_is_me_face_tuned")]
INTER = [("Speed×TF", "time_has_speed_x_tf"), ("Speed×SF", "time_has_speed_x_sf"),
         ("Speed×OR", "time_has_speed_x_or"), ("TF×SF", "time_has_tf_x_sf"),
         ("TF×OR", "time_has_tf_x_or"), ("SF×OR", "time_has_sf_x_or"),
         ("ME×Speed", "time_has_me_face_x_speed")]


def main() -> int:
    df = pd.read_csv(CSV)
    probes = sorted(df["probe_id"].unique())
    n = len(df)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    # Panel L: main-effect selection fraction, grouped by probe + pooled.
    x = np.arange(len(MAIN))
    groups = probes + ["all"]
    w = 0.8 / len(groups)
    cmap = plt.get_cmap("tab10")
    for gi, g in enumerate(groups):
        sub = df if g == "all" else df[df["probe_id"] == g]
        frac = [sub[col].astype(bool).mean() for _, col in MAIN]
        color = "0.2" if g == "all" else cmap(gi)
        bars = axL.bar(x + (gi - (len(groups) - 1) / 2) * w, frac, w,
                       label=(f"all (n={n})" if g == "all"
                              else f"{g.split('_')[0]} (n={len(sub)})"),
                       color=color, edgecolor="white", linewidth=0.4)
        if g == "all":
            for b, f in zip(bars, frac):
                axL.text(b.get_x() + b.get_width() / 2, f + 0.01, f"{f:.0%}",
                         ha="center", va="bottom", fontsize=8, fontweight="bold")
    axL.set_xticks(x); axL.set_xticklabels([m for m, _ in MAIN])
    axL.set_ylim(0, 1.05); axL.set_ylabel("fraction of clusters selecting")
    axL.set_title(f"Forward-selection: main effects (n={n})",
                  fontsize=11, fontweight="bold")
    axL.legend(fontsize=8, frameon=False, ncol=2)

    # Panel R: interaction-term fraction (pooled).
    xi = np.arange(len(INTER))
    fraci = [df[col].astype(bool).mean() for _, col in INTER]
    colors = ["tab:purple" if "ME" in lab else "0.55" for lab, _ in INTER]
    bars = axR.bar(xi, fraci, 0.7, color=colors, edgecolor="white")
    for b, f in zip(bars, fraci):
        axR.text(b.get_x() + b.get_width() / 2, f + 0.005, f"{f:.0%}",
                 ha="center", va="bottom", fontsize=8)
    axR.set_xticks(xi); axR.set_xticklabels([m for m, _ in INTER], rotation=35,
                                            ha="right", fontsize=8)
    axR.set_ylim(0, max(0.2, max(fraci) * 1.25))
    axR.set_ylabel("fraction of clusters with interaction")
    axR.set_title("Selected interaction terms", fontsize=11, fontweight="bold")

    fig.suptitle("current_plus_ME — forward-selection summary "
                 "(pooled, speed-profile CV)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    out_dir = ROOT / "diagnostics"; out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "selection_me"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)

    print(f"n={n} across {len(probes)} probes")
    for lab, col in MAIN:
        print(f"  {lab:8s} selected {df[col].astype(bool).mean():.0%} "
              f"({int(df[col].astype(bool).sum())}/{n})")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
