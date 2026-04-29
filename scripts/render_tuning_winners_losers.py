"""Side-by-side tuning curves for the top 3 IMPROVED + top 3 REGRESSED
clusters from current (no-onset+history) vs legacy (with-onset).

For each cluster:
  - row 1: Speed tuning curve (legacy + current Selected lines on same axes)
  - row 2: TF tuning curve
  - cache observed (per-trial mean) overlaid as a third line

Reads tuning_pearson_by_cluster.csv to identify the winners/losers,
reads the per-probe tuning_curves.csv from both run dirs for the
predicted curves, and reads the MATLAB cache for observed.

Output: glm/current/figs/tuning_winners_losers.{pdf,png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python" / "src"))
from rc2_glm.precomputed_bins import load_precomputed_bin_edges

ROOT = Path("/Users/lauraporta/local_data/motion_clouds/figures")
CURRENT = ROOT / "glm" / "current"
LEGACY = ROOT / "glm" / "legacy_with_onset"
FORMATTED = Path("/Users/lauraporta/local_data/motion_clouds/formatted_data")
OUT = CURRENT / "figs" / "tuning_winners_losers.pdf"


def _cache(probe: str):
    return load_precomputed_bin_edges(FORMATTED / f"{probe}.mat")


def _predicted_curve(run_dir: Path, probe: str, cid: int, var: str) -> np.ndarray | None:
    p = run_dir / "_runs" / probe / "diagnostics" / "tuning_curves.csv"
    if not p.is_file():
        return None
    df = pd.read_csv(p)
    sub = df[(df.cluster_id == cid) & (df.variable == var) & (df.model == "Selected")]
    if sub.empty:
        return None
    grid = sub.grid_x.to_numpy(); gain = sub.gain.to_numpy()
    o = np.argsort(grid)
    return grid[o], gain[o]


def _observed_per_cond(probe: str, cid: int, var: str):
    """Returns dict {cond: (centres, mean_obs)}."""
    pc = _cache(probe)
    out = {}
    for cond in ("VT", "V", "T_Vstatic"):
        centres = pc.speed_centres(cond) if var == "Speed" else pc.tf_centres(cond)
        rec = (pc.speed_by_group if var == "Speed" else pc.tf_by_group).get(cond)
        if centres is None or rec is None:
            continue
        arr = rec.tuning.get(int(cid))
        if arr is None:
            continue
        obs = np.nanmean(arr, axis=0)
        valid = np.isfinite(obs)
        if valid.sum() < 2:
            continue
        out[cond] = (np.asarray(centres)[valid], obs[valid])
    return out


def render_cluster(ax, probe: str, cid: int, var: str, label_pearson: dict):
    obs = _observed_per_cond(probe, cid, var)
    cur = _predicted_curve(CURRENT, probe, cid, var)
    leg = _predicted_curve(LEGACY, probe, cid, var)
    cond_colors = {"VT": "C0", "V": "C2", "T_Vstatic": "C1"}
    for cond, (c, o) in obs.items():
        ax.plot(c, o, "o-", color=cond_colors[cond], alpha=0.7,
                markersize=5, linewidth=1.2, label=f"obs {cond}")
    if leg is not None:
        ax.plot(leg[0], leg[1], "--", color="grey", linewidth=1.8,
                label=f"legacy r={label_pearson.get('legacy', float('nan')):.2f}")
    if cur is not None:
        ax.plot(cur[0], cur[1], "-", color="black", linewidth=1.8,
                label=f"current r={label_pearson.get('current', float('nan')):.2f}")
    ax.set_xlabel(f"{var} ({'cm/s' if var == 'Speed' else 'Hz'})", fontsize=8)
    ax.set_ylabel("FR (Hz)", fontsize=8)
    ax.set_title(f"{probe} c{cid} — {var}", fontsize=9)
    ax.legend(fontsize=6, loc="best")
    ax.grid(True, alpha=0.3)


def main() -> int:
    pp = pd.read_csv(CURRENT / "diagnostics" / "tuning_pearson_by_cluster.csv")
    piv = pp[pp.model == "Selected"].pivot_table(
        index=["probe_id", "cluster_id", "variable"],
        columns="run", values="pearson_r", aggfunc="mean",
    ).reset_index()
    piv["delta"] = piv["current"] - piv["legacy"]
    piv = piv.dropna(subset=["delta"])
    top_imp = piv.nlargest(3, "delta")
    top_reg = piv.nsmallest(3, "delta")
    rows = pd.concat([top_imp, top_reg], ignore_index=True)
    rows["category"] = ["improved"] * 3 + ["regressed"] * 3

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for i, (_, r) in enumerate(rows.iterrows()):
        ax = axes[i // 3, i % 3]
        render_cluster(
            ax, r["probe_id"], int(r["cluster_id"]), r["variable"],
            label_pearson={"legacy": r["legacy"], "current": r["current"]},
        )
        ax.text(0.5, 1.10, f"Δr = {r['delta']:+.2f} ({r['category']})",
                transform=ax.transAxes, ha="center",
                fontsize=10, fontweight="bold",
                color="green" if r["category"] == "improved" else "red")

    fig.suptitle(
        "Tuning-curve examples — TOP 3 IMPROVED (Δ Pearson r > 0) on top, "
        "TOP 3 REGRESSED on bottom",
        fontsize=12,
    )
    fig.savefig(OUT)
    fig.savefig(OUT.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
