"""Consolidated 'defensible statements' figure for the motion-clouds GLM.

Brings together what survives the identifiability analysis: the pooled-model
unique contribution per predictor group (cv-bps, speed-profile leak-test folds),
annotated with the VIF feasibility verdict and the honest caveats. This is the
panel we can stand behind, instead of the (non-identifiable) VT-split results.

Reads current/diagnostics/variance_partition_pooled.csv.
Writes current/diagnostics/honest_framing_summary.{pdf,png}.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "local_data" / "motion_clouds" / "figures" / "glm" / "current"
CSV = ROOT / "diagnostics" / "variance_partition_pooled.csv"

GROUPS = [("Onset", "unique_onset_block", "tab:blue"),
          ("Speed", "unique_Speed", "tab:green"),
          ("TF", "unique_TF", "tab:orange"),
          ("SF", "unique_SF", "tab:olive"),
          ("OR", "unique_OR", "tab:red")]


def main() -> int:
    df = pd.read_csv(CSV)
    n = len(df)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.2),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    # Panel L: pooled unique Δ cv-bps per group (speed-profile leak-test).
    data, labels, colors, meds = [], [], [], []
    for name, col, color in GROUPS:
        if col not in df:
            continue
        d = df[col].to_numpy()
        data.append(d); labels.append(name); colors.append(color)
        meds.append(np.median(d))
    bp = axL.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.35)
    for i, (d, c, m) in enumerate(zip(data, colors, meds), 1):
        axL.scatter(np.full(len(d), i) + np.random.uniform(-.12, .12, len(d)),
                    d, s=9, color=c, alpha=0.55, edgecolor="none")
        axL.text(i, max(d) , f"med {m:+.3f}", ha="center", va="bottom", fontsize=7)
    axL.axhline(0, color="0.6", lw=0.8)
    axL.axhline(0.005, color="0.6", lw=0.8, ls=":")
    axL.text(len(data) + 0.3, 0.005, "0.005\nthresh", fontsize=6, va="center")
    axL.set_ylabel("unique Δ cv-bps (pooled, speed-profile folds)")
    axL.set_title(f"What uniquely contributes — pooled model (n={n})\n"
                  "(leave-one-group-out; Speed carries it, Onset is redundant)",
                  fontsize=10, fontweight="bold")

    # Panel R: the honest framing, as text.
    axR.axis("off")
    txt = (
        "HONEST FRAMING — what we can defensibly say\n"
        "─────────────────────────────────────────\n\n"
        "FEASIBILITY (VIF / cond#):\n"
        "  • pooled  Speed 1.17 · TF 1.23 · Onset 1.17   → separable ✓\n"
        "  • T_Vstatic Speed 1.36 ✓   • V  TF 1.24 ✓\n"
        "  • VT  cond# 205, VIF→∞  → Speed/TF NOT separable ✗\n\n"
        "DEFENSIBLE (pooled + speed-profile leak-test):\n"
        f"  • Speed uniquely contributes (med "
        f"{df['unique_Speed'].median():+.3f} bps, "
        f"{int((df['unique_Speed']>0.005).sum())}/{n} cells),\n"
        "    and it GENERALISES across the two speed trajectories\n"
        "    → real speed encoding, not a time/onset artefact.\n"
        f"  • Onset unique ≈ 0 (med {df['unique_onset_block'].median():+.3f}) "
        "→ redundant, not a carrier.\n"
        "  • TF/SF/OR add little UNIQUE beyond Speed across conditions.\n\n"
        "NOT DEFENSIBLE (drop these):\n"
        "  • VT-split Speed-vs-TF counts (35/33) and the\n"
        "    'Speed wins integration' plot — tie-break artefacts of an\n"
        "    unidentifiable (TF=gain·Speed) design.\n"
        "  • Any 'how tuning CHANGES across conditions' claim —\n"
        "    pooling assumes shared tuning; VT can't isolate the change.\n\n"
        "WHY: 3 conditions de-correlate Speed/TF/Onset (Speed-no-TF in T,\n"
        "TF-no-Speed in V, Onset-no-Speed in V). The 2 speed profiles\n"
        "guard Speed↔time; only the 3 gains touch Speed↔TF (too weak).\n"
        "Cross-condition INTEGRATION needs the active paradigm."
    )
    axR.text(0.0, 1.0, txt, ha="left", va="top", family="monospace",
             fontsize=7.4, transform=axR.transAxes)
    fig.suptitle("Motion-clouds GLM — consolidated honest framing",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = ROOT / "diagnostics" / "honest_framing_summary"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
