"""Per-cluster acid test (unique cv-bps) with forward-selection overlay.

Like the per-cluster contribution figure, but each predictor's unique-bps
segment is drawn at full opacity if FORWARD SELECTION picked that predictor for
that cluster, and at alpha 0.5 if it did NOT. Marries attribution (full-model
leave-one-out unique) with selection (which predictors the Hardcastle step kept)
— makes the selection-vs-attribution gap visible per cluster:
  - bright & tall  = selected AND uniquely contributes
  - faded & tall   = uniquely contributes but selection dropped it (redundant)
  - faded & ~0     = correctly left out

Acid test from variance_partition_accel.csv (= the canonical+ME+History+Accel
model); selection from current_ME_hist_accel_20ms/glm_model_comparison.csv (same
model, same cohort/CV). Negatives clipped to 0 in the stack.

Output: current_ME_hist_accel_20ms/diagnostics/acid_vs_selection.{pdf,png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

G = Path.home() / "local_data/motion_clouds/figures/glm"
# --bin10 → 10 ms run (partition + selection both from current_ME_hist_accel_10ms);
# default = 20 ms (partition reused from current_plus_ME_20ms_accel).
if "--bin10" in sys.argv:
    RUN = G / "current_ME_hist_accel_10ms"
    PART = RUN / "diagnostics/variance_partition_accel.csv"
    SEL = RUN / "glm_model_comparison.csv"
    OUT = RUN / "diagnostics"
else:
    PART = G / "current_plus_ME_20ms_accel/diagnostics/variance_partition_accel.csv"
    SEL = G / "current_ME_hist_accel_20ms/glm_model_comparison.csv"
    OUT = G / "current_ME_hist_accel_20ms/diagnostics"

# (label, unique col, forward-selection name [None = always-on baseline], colour)
VARS = [
    ("Onset", "unique_Onset", None, "tab:blue"),
    ("Speed", "unique_Speed", "Speed", "tab:green"),
    ("TF", "unique_TF", "TF", "tab:orange"),
    ("SF", "unique_SF", "SF", "tab:olive"),
    ("OR", "unique_OR", "OR", "tab:red"),
    ("ME", "unique_ME", "ME_face", "tab:purple"),
    ("History", "unique_History", "History", "tab:brown"),
    ("Accel", "unique_Accel", "Acceleration", "tab:cyan"),
]


def plot_acid_figure(df, out, *, forward_only=False):
    """Per-cluster stacked unique cv-bps. Default: bright=forward-selected,
    faded(α0.5)=uniquely contributes but selection dropped it (the
    attribution-vs-selection gap). ``forward_only=True``: show ONLY the
    forward-selected predictors' contributions — drop the faded segments and
    sum only what selection kept. ``df`` must carry the ``unique_*`` columns,
    a ``sel_set`` column, and ``r2`` (sorted)."""
    n = len(df)
    fig, ax = plt.subplots(figsize=(15, 5.4))
    x = np.arange(n)
    bottom = np.zeros(n)
    for lab, col, selname, c in VARS:
        if col not in df.columns:
            continue
        h = np.clip(df[col].to_numpy(), 0, None)
        picked = (np.ones(n, dtype=bool) if selname is None
                  else df["sel_set"].apply(lambda s, sn=selname: sn in s).to_numpy())
        if forward_only:
            h = np.where(picked, h, 0.0)               # only selected predictors
            ax.bar(x, h, bottom=bottom, width=1.0, color=c, linewidth=0)
        else:
            for mask, alpha in ((picked, 1.0), (~picked, 0.5)):
                if mask.any():
                    ax.bar(x[mask], h[mask], bottom=bottom[mask], width=1.0,
                           color=c, alpha=alpha, linewidth=0)
        bottom += h
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_xlabel("cluster (sorted by deviance explained)")
    ax.set_ylabel("stacked unique Δ cv-bps"
                  + (" (forward-selected only)" if forward_only else " (acid test)"))
    ax2 = ax.twinx(); ax2.plot(x, df["r2"], color="black", lw=1.4)
    ax2.set_ylabel("pseudo-R² (black line)"); ax2.set_ylim(0, max(0.05, df["r2"].max() * 1.05))

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, label=lab) for lab, col, _, c in VARS if col in df.columns]
    if not forward_only:
        handles += [Patch(facecolor="0.4", alpha=1.0, label="forward-selected"),
                    Patch(facecolor="0.4", alpha=0.5, label="NOT selected (α0.5)")]
    ax.legend(handles=handles, ncol=5, fontsize=8, frameon=False, loc="upper left")
    title = (f"Per-cluster forward-selected contributions (n={n})" if forward_only else
             f"Per-cluster acid test × forward selection (n={n}) — "
             "faded = uniquely contributes but selection dropped it")
    ax.set_title(title, fontsize=11, fontweight="bold")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")


def plot_acid_stack_sorted(df, out):
    """Per-cluster stacked unique Δ cv-bps (acid test, ALL regressors), sorted by
    the stack total. Two subplots: all regressors (top) and all-except-History
    (bottom, re-sorted on its own total) — History dominates the stack, so
    dropping it exposes the Speed/Accel/ME/TF/SF/OR structure. Reuses VARS;
    negatives clipped to 0. ``df`` carries the per-cluster ``unique_*`` columns."""
    fig, (axT, axB) = plt.subplots(2, 1, figsize=(15, 8.5))
    n = len(df)
    for ax, vars_set, title in (
        (axT, VARS, "all regressors"),
        (axB, [v for v in VARS if v[0] != "History"], "all regressors except History"),
    ):
        present = [(lab, col, c) for lab, col, _, c in vars_set if col in df.columns]
        heights = {lab: np.clip(pd.to_numeric(df[col], errors="coerce").fillna(0).to_numpy(),
                                0, None) for lab, col, c in present}
        stack_tot = np.sum(list(heights.values()), axis=0)
        order = np.argsort(stack_tot)                  # ascending: small stacks left
        x = np.arange(n); bottom = np.zeros(n)
        for lab, col, c in present:
            h = heights[lab][order]
            ax.bar(x, h, bottom=bottom, width=1.0, color=c, linewidth=0, label=lab)
            bottom = bottom + h
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylabel("stacked unique Δ cv-bps")
        ax.set_title(f"Per-cluster unique contributions — {title} "
                     f"(n={n}, sorted by stack total)", fontsize=10, fontweight="bold")
        ax.legend(ncol=len(present), fontsize=7.5, frameon=False, loc="upper left")
    axB.set_xlabel("cluster (re-sorted by stacked unique Δ cv-bps within each subplot)")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")


def main() -> int:
    part = pd.read_csv(PART)
    sel = pd.read_csv(SEL)[["probe_id", "cluster_id", "time_selected_vars"]]
    df = part.merge(sel, on=["probe_id", "cluster_id"], how="inner")
    df = df[df["total_full"] > 0.005].copy()
    df["sel_set"] = df["time_selected_vars"].fillna("Null").apply(lambda s: set(s.split("+")))
    df["r2"] = 1 - df["cv_full"] / df["cv_intercept"]
    df = df.sort_values("r2").reset_index(drop=True)
    n = len(df)

    OUT.mkdir(parents=True, exist_ok=True)
    plot_acid_figure(df, OUT / "acid_vs_selection")

    # quick report: per predictor, how often unique>0.005 but NOT selected
    print(f"n={n}; selection-vs-attribution mismatches (unique>0.005 but NOT forward-selected):")
    for lab, col, selname, _ in VARS:
        if selname is None:
            continue
        u = df[col] > 0.005
        notsel = ~df["sel_set"].apply(lambda s, sn=selname: sn in s)
        print(f"  {lab:8s} unique>0.005 in {int(u.sum())}, of which NOT selected: {int((u & notsel).sum())}")
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
