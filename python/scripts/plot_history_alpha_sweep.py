"""Per-cluster history-α sweep — does the trial-level Pearson trade off
against the history-on cv_bps gain as we depotentiate history?

Reads:
    glm/current/diagnostics/history_alpha_sweep.csv
        — per-cluster Pearson r at α ∈ {0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5}
          where α scales the history β contribution multiplicatively in
          log-prediction space (η_α = η_no_history + α·(η_full − η_no_history))

Writes:
    glm/current/figs/history_alpha_sweep.{pdf,png}
        — top: median Pearson r vs α across clusters (with IQR band)
        — bottom: per-cluster faint lines, with the population median
                  highlighted; vertical at α=0 (no history) and α=1
                  (production)
        — annotation: argmax-α per cluster (median across clusters)
        — bottom-right: histogram of per-cluster argmax α

Useful for the "potentiate / depotentiate history" question Laura raised:
where does each cluster's predictive Pearson peak? At α=1 (production)?
At α<1 (over-fit)? At α>1 (history under-weighted)?
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CURRENT = Path("/Users/lauraporta/local_data/motion_clouds/figures/glm/current")
INPUT = CURRENT / "diagnostics" / "history_alpha_sweep.csv"
OUT_PDF = CURRENT / "figs" / "history_alpha_sweep.pdf"


def main() -> int:
    if not INPUT.is_file():
        print(f"ERROR: {INPUT} missing — run rc2-glm to generate", file=sys.stderr)
        return 1
    df = pd.read_csv(INPUT)
    if df.empty:
        print(f"ERROR: {INPUT} is empty", file=sys.stderr)
        return 1

    # Pivot to (cluster, α) → r
    pivot = df.pivot_table(
        index=["probe_id", "cluster_id"],
        columns="alpha",
        values="pearson_r_overall",
    )
    alphas = sorted(pivot.columns.tolist())
    pivot = pivot[alphas]

    fig = plt.figure(figsize=(13, 8), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    ax_main = fig.add_subplot(gs[0, :])
    ax_lines = fig.add_subplot(gs[1, :2])
    ax_argmax = fig.add_subplot(gs[1, 2])

    # row 0: median + IQR band across clusters
    median_r = pivot.median(axis=0).to_numpy()
    q25 = pivot.quantile(0.25, axis=0).to_numpy()
    q75 = pivot.quantile(0.75, axis=0).to_numpy()
    ax_main.fill_between(alphas, q25, q75, color="C0", alpha=0.2, label="IQR")
    ax_main.plot(alphas, median_r, "o-", color="C0", linewidth=2.5, markersize=8,
                 label=f"median across {len(pivot)} clusters")
    ax_main.axvline(1.0, color="grey", linestyle="--", linewidth=1, alpha=0.5)
    ax_main.text(1.02, ax_main.get_ylim()[0] if ax_main.get_ylim()[0] else 0,
                 "production (α=1)", fontsize=9, color="grey", verticalalignment="bottom")
    ax_main.axvline(0.0, color="red", linestyle=":", linewidth=1, alpha=0.5)
    ax_main.set_xlabel("α (scaling on history β contribution)")
    ax_main.set_ylabel("trial-level Pearson r (predicted vs observed counts)")
    ax_main.set_title(
        f"History-α sweep — per-cluster trial-level Pearson r vs α scaling "
        f"({len(pivot)} clusters)"
    )
    ax_main.legend()
    ax_main.grid(True, alpha=0.3)

    # row 1 left: per-cluster faint lines
    for _, row in pivot.iterrows():
        ax_lines.plot(alphas, row.to_numpy(), "-", color="grey",
                      alpha=0.15, linewidth=0.7)
    ax_lines.plot(alphas, median_r, "o-", color="C3", linewidth=2.5,
                  markersize=8, label="median")
    ax_lines.axvline(1.0, color="grey", linestyle="--", linewidth=1, alpha=0.5)
    ax_lines.axvline(0.0, color="red", linestyle=":", linewidth=1, alpha=0.5)
    ax_lines.set_xlabel("α")
    ax_lines.set_ylabel("Pearson r")
    ax_lines.set_title("per-cluster trajectories")
    ax_lines.legend()
    ax_lines.grid(True, alpha=0.3)

    # row 1 right: histogram of per-cluster argmax α
    argmax_alpha = np.array([alphas[i] for i in pivot.values.argmax(axis=1)])
    bins = np.array(alphas + [alphas[-1] + (alphas[-1] - alphas[-2])]) - \
           (alphas[1] - alphas[0]) / 2
    ax_argmax.hist(argmax_alpha, bins=bins, color="purple", alpha=0.7,
                   edgecolor="black")
    ax_argmax.axvline(1.0, color="grey", linestyle="--", linewidth=1, alpha=0.5)
    ax_argmax.set_xlabel("α at peak Pearson r")
    ax_argmax.set_ylabel("# clusters")
    ax_argmax.set_title(f"per-cluster optimal α (median={np.median(argmax_alpha):.2f})")
    ax_argmax.grid(True, alpha=0.3)

    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PDF.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"wrote {OUT_PDF}")
    print()
    print(f"=== summary ({len(pivot)} clusters) ===")
    for a, m in zip(alphas, median_r):
        print(f"  α={a:.2f}: median Pearson r = {m:+.4f}")
    print()
    print(f"per-cluster argmax α: median = {np.median(argmax_alpha):.2f}, "
          f"mean = {np.mean(argmax_alpha):.2f}")
    print(f"  α=0   peak: {(argmax_alpha == 0).sum()}/{len(argmax_alpha)} clusters")
    print(f"  α=1   peak: {(argmax_alpha == 1).sum()}/{len(argmax_alpha)} clusters")
    print(f"  α>1   peak: {(argmax_alpha > 1).sum()}/{len(argmax_alpha)} clusters")
    print(f"  α<1   peak: {(argmax_alpha < 1).sum()}/{len(argmax_alpha)} clusters")
    return 0


if __name__ == "__main__":
    sys.exit(main())
