"""Boundary × activity + variable-inclusion analysis on the face_me run.

Variant of ``python/notebooks/boundary_venn_exploration.py`` (2026-04-24)
re-pointed at the prompt-06 face_me run. Same fragile-band logic
(Δ-bps ∈ [0.002, 0.008]); same UpSet / Venn structure; same per-probe
breakdown. Question: does adding ME_face change the fragile-row pattern
that bit us in legacy and current?

Inputs (read-only):
  - glm/exploration/face_me_onset_no_history/_runs/<probe>/glm_selection_history.csv
  - glm/exploration/face_me_onset_no_history/_runs/<probe>/glm_model_comparison.csv
  (Aggregated across the 4 probes.)

Outputs in glm/exploration/face_me_onset_no_history/boundary_venn/:
  - boundary_vs_activity.{pdf,png}
  - upset_5way.{pdf,png}              # now 5 main effects: +ME_face
  - venn_per_probe_stacked.{pdf,png}
  - boundary_clusters.csv
  - fragile_by_candidate.csv
  - venn_cells.csv, venn_cells_by_probe.csv

Comparison vs the legacy boundary_venn (in
glm/exploration/boundary_venn/) is the headline: did adding ME_face
move clusters in or out of the fragile band?

Usage::

    python -m scripts.boundary_venn_face_me
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, spearmanr

HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds" / "figures" / "glm" / "exploration" / "face_me_onset_no_history"
RUNS_DIR = ROOT / "_runs"
OUT_DIR = ROOT / "boundary_venn"
OUT_DIR.mkdir(parents=True, exist_ok=True)

THRESH = 0.005
FRAGILE_LOW = 0.002
FRAGILE_HIGH = 0.008
GLM_TYPE = "time"


def _load_aggregated() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Concatenate per-probe selection_history + model_comparison."""
    history_frames = []
    comparison_frames = []
    for probe_dir in sorted(RUNS_DIR.glob("CAA-*_rec1")):
        probe_id = probe_dir.name
        h = probe_dir / "glm_selection_history.csv"
        c = probe_dir / "glm_model_comparison.csv"
        if not (h.exists() and c.exists()):
            print(f"skip {probe_id}: missing CSVs")
            continue
        hdf = pd.read_csv(h)
        hdf["probe_id"] = probe_id
        history_frames.append(hdf)
        cdf = pd.read_csv(c)
        if "probe_id" not in cdf.columns:
            cdf["probe_id"] = probe_id
        comparison_frames.append(cdf)
    if not history_frames:
        raise RuntimeError("no per-probe CSVs found under " + str(RUNS_DIR))
    history = pd.concat(history_frames, ignore_index=True)
    comparison = pd.concat(comparison_frames, ignore_index=True)
    return history, comparison


def _per_round_long(history: pd.DataFrame) -> pd.DataFrame:
    """Long-format: one row per (probe, cluster, round, candidate, delta_bps)."""
    # Columns observed 2026-04-30: probe_id, cluster_id, round, phase,
    # best_candidate, delta_bps, added, cv_bps_after. One row per (cluster,
    # round) with the round's WINNING candidate's Δ-bps.
    keep = ["probe_id", "cluster_id", "round", "phase",
            "best_candidate", "delta_bps", "added"]
    return history[[c for c in keep if c in history.columns]].rename(
        columns={"best_candidate": "candidate"}
    )


def _boundary_clusters(per_round: pd.DataFrame) -> pd.DataFrame:
    """Each row's delta_bps is in the fragile band → could flip on noise."""
    mask = (
        (per_round["delta_bps"] >= FRAGILE_LOW)
        & (per_round["delta_bps"] <= FRAGILE_HIGH)
    )
    return per_round[mask].copy()


def main() -> int:
    history, comparison = _load_aggregated()
    print(f"loaded {len(history)} history rows, {len(comparison)} clusters")
    per_round = _per_round_long(history)
    boundary = _boundary_clusters(per_round)
    print(f"fragile rows (Δ-bps ∈ [{FRAGILE_LOW}, {FRAGILE_HIGH}]): {len(boundary)}")
    boundary.to_csv(OUT_DIR / "boundary_clusters.csv", index=False)

    cand_tbl = (
        boundary.groupby("candidate").size().rename("n_fragile_rows")
        .reset_index().sort_values("n_fragile_rows", ascending=False)
    )
    cand_tbl.to_csv(OUT_DIR / "fragile_by_candidate.csv", index=False)
    print("\nfragile by candidate:")
    print(cand_tbl.to_string(index=False))

    # --- Boundary × activity scatter (uses time_n_spikes from comparison) ---
    n_spikes_by_id = comparison.set_index(["probe_id", "cluster_id"])[
        "time_n_spikes"
    ].to_dict()
    boundary["n_spikes"] = boundary.apply(
        lambda r: n_spikes_by_id.get((r["probe_id"], r["cluster_id"]), np.nan),
        axis=1,
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    # Panel 1: per-round Δ-bps vs n_spikes (log)
    ax = axes[0]
    ax.scatter(
        boundary["n_spikes"], boundary["delta_bps"],
        c="C3", alpha=0.6, edgecolors="none", s=30,
    )
    ax.axhline(THRESH, color="k", linestyle="--", linewidth=0.8, label=f"threshold = {THRESH}")
    ax.axhspan(FRAGILE_LOW, FRAGILE_HIGH, color=(0.85, 0.85, 0.85), alpha=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("cluster n_spikes (log)")
    ax.set_ylabel("Δ-bps at fragile round")
    ax.set_title("Fragile rounds — activity × Δ-bps")
    ax.legend(fontsize=9)

    # Panel 2: quartile fragility rate
    valid = comparison[~comparison["time_n_spikes"].isna()]
    quartiles = pd.qcut(valid["time_n_spikes"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    fragile_set = set(boundary.set_index(["probe_id", "cluster_id"]).index)
    is_fragile = valid.apply(
        lambda r: (r["probe_id"], r["cluster_id"]) in fragile_set, axis=1,
    )
    rate_by_q = pd.DataFrame({"q": quartiles.values, "fragile": is_fragile.values})
    rate = rate_by_q.groupby("q", observed=True)["fragile"].mean().reset_index()
    ax2 = axes[1]
    ax2.bar(rate["q"], rate["fragile"], color="C3", edgecolor="black")
    ax2.set_ylabel("fraction of clusters with ≥1 fragile round")
    ax2.set_title("Fragility by activity quartile")

    # Panel 3: per-candidate fragile count (bar)
    ax3 = axes[2]
    ax3.barh(cand_tbl["candidate"], cand_tbl["n_fragile_rows"], color="C0")
    ax3.set_xlabel("# fragile rounds")
    ax3.set_title("Which candidates are fragile?")
    fig.suptitle(
        f"Boundary × activity (face_me run, n={len(comparison)} clusters, "
        f"{len(boundary)} fragile rounds)", fontsize=12, fontweight="bold",
    )
    fig.savefig(OUT_DIR / "boundary_vs_activity.pdf")
    fig.savefig(OUT_DIR / "boundary_vs_activity.png", dpi=150)
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'boundary_vs_activity.pdf'}")

    # --- Variable-inclusion UpSet across {Speed, TF, SF, OR, ME_face} ---
    candidates = ["Speed", "TF", "SF", "OR", "ME_face"]
    bool_cols = {}
    sv = comparison["time_selected_vars"].fillna("Null").astype(str)
    for c in candidates:
        bool_cols[c] = sv.apply(lambda s, c=c: c in s.split("+"))
    cell = pd.DataFrame(bool_cols).astype(int)
    cell["probe_id"] = comparison["probe_id"].values
    cell["cluster_id"] = comparison["cluster_id"].values
    cell["combo"] = cell[candidates].apply(
        lambda r: "+".join([c for c in candidates if r[c]]) or "(none)", axis=1,
    )
    cell_counts = cell["combo"].value_counts().reset_index()
    cell_counts.columns = ["combo", "n_clusters"]
    cell_counts.to_csv(OUT_DIR / "venn_cells.csv", index=False)
    per_probe = (
        cell.groupby(["probe_id", "combo"]).size()
        .rename("n_clusters").reset_index()
    )
    per_probe.to_csv(OUT_DIR / "venn_cells_by_probe.csv", index=False)
    print(f"\n{len(cell_counts)} non-empty inclusion cells across 5 candidates")
    print(cell_counts.head(15).to_string(index=False))

    # Simple horizontal bar plot of inclusion combos (top 15)
    top = cell_counts.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax.barh(top["combo"], top["n_clusters"], color="C0")
    ax.set_xlabel("# clusters")
    ax.set_title(
        f"Variable-inclusion combos (face_me run, top 15 of {len(cell_counts)})",
        fontsize=11,
    )
    fig.savefig(OUT_DIR / "upset_5way.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "upset_5way.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'upset_5way.pdf'}")

    # Per-probe stacked bar — same combos as the upset, broken out by probe.
    pivot = (
        per_probe.pivot(index="combo", columns="probe_id", values="n_clusters")
        .fillna(0).astype(int)
    )
    top_combos = cell_counts.head(12)["combo"].tolist()
    pivot = pivot.reindex(top_combos[::-1])
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    bottom = np.zeros(len(pivot))
    for col in pivot.columns:
        ax.barh(pivot.index, pivot[col], left=bottom, label=col)
        bottom += pivot[col].to_numpy()
    ax.set_xlabel("# clusters")
    ax.set_title("Per-probe inclusion combos (face_me run)", fontsize=11)
    ax.legend(loc="lower right", fontsize=9, ncol=2)
    fig.savefig(OUT_DIR / "venn_per_probe_stacked.pdf")
    fig.savefig(OUT_DIR / "venn_per_probe_stacked.png", dpi=150)
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'venn_per_probe_stacked.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
