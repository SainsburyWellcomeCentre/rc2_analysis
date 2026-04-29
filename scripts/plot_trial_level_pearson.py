"""Per-cluster trial-level Pearson r — current model + per-cluster distribution.

Reads:
    glm/current/diagnostics/trial_level_metrics.csv
        — per (cluster, model) overall + per-trial-mean Pearson r between
          observed counts and predicted λ·Δt

Writes:
    glm/current/figs/trial_level_pearson.{pdf,png}
        — top row: histograms of per-cluster overall and per-trial-mean r
                   for the Selected model
        — bottom row: per-model comparison (Null / Selected / Additive / FullInteraction)
                      + scatter of overall vs per-trial-mean r per cluster

This addresses the question "did trial-level reconstruction improve with
the new model?" — the production default is `Selected = History+Speed+...`,
and we should see the Pearson r being substantially > 0 for most clusters.
The Null model row is the floor (intercept-only) — Selected should be clearly
above it.
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
INPUT = CURRENT / "diagnostics" / "trial_level_metrics.csv"
OUT_PDF = CURRENT / "figs" / "trial_level_pearson.pdf"


def main() -> int:
    if not INPUT.is_file():
        print(f"ERROR: {INPUT} missing", file=sys.stderr)
        return 1
    df = pd.read_csv(INPUT)
    if df.empty:
        print(f"ERROR: {INPUT} is empty", file=sys.stderr)
        return 1

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)

    # row 0 (col 0): histogram of overall Pearson r — Selected model
    ax = axes[0, 0]
    sel = df[df["model"] == "Selected"]
    ax.hist(sel["pearson_r_overall"].dropna(), bins=30, color="C3", alpha=0.7,
            edgecolor="black")
    med = sel["pearson_r_overall"].median()
    ax.axvline(med, color="black", linewidth=1.5, label=f"median={med:.3f}")
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_xlabel("overall Pearson r (Selected)")
    ax.set_ylabel("# clusters")
    ax.set_title(f"Trial-level r (overall) — Selected, n={len(sel)}")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # row 0 (col 1): histogram of per-trial-mean Pearson r
    ax = axes[0, 1]
    ax.hist(sel["pearson_r_per_trial_mean"].dropna(), bins=30, color="C0",
            alpha=0.7, edgecolor="black")
    med = sel["pearson_r_per_trial_mean"].median()
    ax.axvline(med, color="black", linewidth=1.5, label=f"median={med:.3f}")
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_xlabel("per-trial-mean Pearson r (Selected)")
    ax.set_ylabel("# clusters")
    ax.set_title(f"Trial-level r (per-trial avg) — Selected")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # row 0 (col 2): scatter overall vs per-trial-mean r
    ax = axes[0, 2]
    ax.scatter(sel["pearson_r_overall"], sel["pearson_r_per_trial_mean"],
               color="C3", alpha=0.6, s=25)
    ax.plot([-0.2, 1.0], [-0.2, 1.0], "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("overall r")
    ax.set_ylabel("per-trial-mean r")
    ax.set_xlim(-0.2, 1.0); ax.set_ylim(-0.2, 1.0)
    ax.set_title("per-cluster: overall vs per-trial-mean (Selected)")
    ax.grid(True, alpha=0.3)

    # row 1: per-model comparison
    models = ["Null", "Selected", "Additive", "FullInteraction"]
    colors = ["grey", "C3", "C0", "C2"]

    # row 1 col 0: box plot of overall r per model
    ax = axes[1, 0]
    data = [df[df["model"] == m]["pearson_r_overall"].dropna().to_numpy()
            for m in models]
    bp = ax.boxplot(data, tick_labels=models, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color); patch.set_alpha(0.5)
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_ylabel("overall Pearson r")
    ax.set_title("per-model trial-level r (overall)")
    ax.grid(True, alpha=0.3)

    # row 1 col 1: box plot of per-trial-mean r per model
    ax = axes[1, 1]
    data = [df[df["model"] == m]["pearson_r_per_trial_mean"].dropna().to_numpy()
            for m in models]
    bp = ax.boxplot(data, tick_labels=models, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color); patch.set_alpha(0.5)
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_ylabel("per-trial-mean Pearson r")
    ax.set_title("per-model trial-level r (per-trial avg)")
    ax.grid(True, alpha=0.3)

    # row 1 col 2: lift over null per cluster (Selected r - Null r)
    ax = axes[1, 2]
    pivot = df.pivot_table(
        index=["probe_id", "cluster_id"],
        columns="model",
        values="pearson_r_overall",
    )
    if "Selected" in pivot.columns and "Null" in pivot.columns:
        delta = (pivot["Selected"] - pivot["Null"]).dropna()
        ax.hist(delta, bins=30, color="purple", alpha=0.7, edgecolor="black")
        med = delta.median()
        ax.axvline(med, color="black", linewidth=1.5, label=f"median Δ = {med:+.3f}")
        ax.axvline(0, color="grey", linewidth=0.5)
        ax.set_xlabel("Δ Pearson r (Selected − Null)")
        ax.set_ylabel("# clusters")
        ax.set_title("lift over Null (Selected)")
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Trial-level Pearson r — current production (no onset, history-on)",
        fontsize=13,
    )
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PDF.with_suffix(".png"), dpi=150)
    plt.close(fig)

    print(f"wrote {OUT_PDF}")
    print()
    print(f"=== Selected model trial-level Pearson r (n={len(sel)}) ===")
    print(f"  overall            median = {sel['pearson_r_overall'].median():+.3f}, mean = {sel['pearson_r_overall'].mean():+.3f}")
    print(f"  per-trial-mean     median = {sel['pearson_r_per_trial_mean'].median():+.3f}, mean = {sel['pearson_r_per_trial_mean'].mean():+.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
