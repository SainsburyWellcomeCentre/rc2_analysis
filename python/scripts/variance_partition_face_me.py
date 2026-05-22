"""Variance partition for the face_me run — v1 (Pearson-r based).

Reads the per-probe ``diagnostics/trial_level_metrics.csv`` files and
computes per-cluster:

  - ``r_selected``           : Selected model trial-level Pearson r (overall)
  - ``r_selected_no_me_face``: same but with ME_face β zeroed at prediction
  - ``unique_me_face_r``     : r_selected − r_selected_no_me_face (Δ Pearson
                                attributable to ME_face)
  - ``r_null``, ``r_additive``, ``r_full_interaction`` for context
  - ``me_face_selected``     : whether ME_face was in the cluster's selected_vars

Caveats:
  - This is **Pearson-r based**, not cv-bps based. A proper bps-level
    partition needs per-cluster CV refits with Selected_no_me / Selected_no_speed
    as fixed-model labels — that's queued for the next pipeline iteration
    alongside the ME_face tuning column. Pearson r drop tracks bps drop
    qualitatively but not exactly (Pearson is on motion bins; bps weights
    by per-bin Poisson likelihood).
  - ``Selected_no_speed`` is NOT currently emitted by the pipeline, so a
    Speed-side unique partition isn't computable here. Coming with the
    next pipeline edit.
  - ME_face is a session-level signal — clusters from probes without
    real camera data (CAA-1123244 placeholder; back half of CAA-1123467)
    have no ``Selected_no_me_face`` entry and are excluded from the
    partition.

Outputs (under
``~/local_data/motion_clouds/figures/glm/exploration/face_me_onset_no_history/diagnostics/``):
  - ``variance_partition_pearson.csv``: per-cluster table.
  - ``variance_partition_pearson.pdf``: 2-panel summary (per-cluster
    bars of Δ Pearson; histogram of unique_me_face_r).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HOME = Path.home()
ROOT = (
    HOME / "local_data" / "motion_clouds" / "figures" / "glm"
    / "exploration" / "face_me_onset_no_history"
)
RUNS_DIR = ROOT / "_runs"
OUT_DIR = ROOT / "diagnostics"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _load_per_probe() -> pd.DataFrame:
    frames = []
    for probe_dir in sorted(RUNS_DIR.glob("CAA-*_rec1")):
        f = probe_dir / "diagnostics" / "trial_level_metrics.csv"
        if not f.exists():
            print(f"skip {probe_dir.name}: {f} absent")
            continue
        df = pd.read_csv(f)
        if "probe_id" not in df.columns:
            df["probe_id"] = probe_dir.name
        frames.append(df)
    if not frames:
        raise RuntimeError("no per-probe trial_level_metrics.csv found")
    return pd.concat(frames, ignore_index=True)


def main() -> int:
    metrics = _load_per_probe()
    print(f"loaded {len(metrics)} (cluster, model) rows")
    print(f"models: {sorted(metrics['model'].unique())}")

    # Pivot so each cluster gets one row, with Pearson r for each model.
    pivot = metrics.pivot_table(
        index=["probe_id", "cluster_id"],
        columns="model",
        values="pearson_r_overall",
        aggfunc="first",
    ).reset_index()

    keep_models = [
        "Null", "Selected", "Additive", "FullInteraction",
        "Selected_no_me_face",
    ]
    for m in keep_models:
        if m not in pivot.columns:
            pivot[m] = np.nan

    pivot = pivot.rename(columns={
        "Null": "r_null",
        "Selected": "r_selected",
        "Additive": "r_additive",
        "FullInteraction": "r_full_interaction",
        "Selected_no_me_face": "r_selected_no_me_face",
    })

    # ME_face contribution to trial-level Pearson r. Defined only for
    # clusters where the pipeline emitted a Selected_no_me_face variant
    # (i.e. ME_face columns existed in the Selected fit).
    pivot["unique_me_face_r"] = (
        pivot["r_selected"] - pivot["r_selected_no_me_face"]
    )

    # Bring in the cluster's selected_vars to flag ME_face selection.
    sv_frames = []
    for probe_dir in sorted(RUNS_DIR.glob("CAA-*_rec1")):
        f = probe_dir / "glm_model_comparison.csv"
        if not f.exists():
            continue
        cmp = pd.read_csv(f)[["probe_id", "cluster_id", "time_selected_vars"]]
        sv_frames.append(cmp)
    sv = pd.concat(sv_frames, ignore_index=True)
    pivot = pivot.merge(sv, on=["probe_id", "cluster_id"], how="left")
    pivot["me_face_selected"] = (
        pivot["time_selected_vars"].fillna("Null").astype(str)
        .apply(lambda s: "ME_face" in s.split("+"))
    )

    out_csv = OUT_DIR / "variance_partition_pearson.csv"
    pivot.to_csv(out_csv, index=False)
    print(f"wrote {out_csv} ({len(pivot)} clusters)")

    # --- Summary stats ---
    has_me_var = pivot["unique_me_face_r"].notna()
    n_clusters_with_me = int(has_me_var.sum())
    print()
    print(f"clusters with Selected_no_me_face emission: {n_clusters_with_me}/{len(pivot)}")
    if n_clusters_with_me > 0:
        sub = pivot.loc[has_me_var, "unique_me_face_r"]
        print(f"  unique_me_face_r: median={sub.median():+.4f}, "
              f"mean={sub.mean():+.4f}, std={sub.std():+.4f}")
        print(f"  Δr > 0 (ME helps Pearson): {int((sub > 0).sum())}/{n_clusters_with_me}")
        print(f"  Δr ≤ 0 (ME doesn't help / hurts): {int((sub <= 0).sum())}/{n_clusters_with_me}")

    # Compare ME-selected vs ME-not-selected clusters
    me_sel_mask = pivot["me_face_selected"] & has_me_var
    me_not_mask = (~pivot["me_face_selected"]) & has_me_var
    if me_sel_mask.any():
        print(f"\n  ME-selected clusters (n={int(me_sel_mask.sum())}):")
        s = pivot.loc[me_sel_mask, "unique_me_face_r"]
        print(f"    median Δr = {s.median():+.4f}, mean = {s.mean():+.4f}")
    if me_not_mask.any():
        print(f"  ME-NOT-selected clusters (n={int(me_not_mask.sum())}):")
        s = pivot.loc[me_not_mask, "unique_me_face_r"]
        print(f"    median Δr = {s.median():+.4f}, mean = {s.mean():+.4f}")
        print("    (should be ~0 — these clusters have no ME_face columns to remove)")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    # Panel 1: histogram of unique_me_face_r
    ax = axes[0]
    finite = pivot["unique_me_face_r"].dropna().to_numpy()
    if finite.size:
        ax.hist(finite, bins=30, color="C0", edgecolor="black", linewidth=0.5)
        ax.axvline(0.0, color="k", linewidth=1, linestyle="--")
        ax.axvline(np.median(finite), color="C3", linewidth=1.5,
                   label=f"median = {np.median(finite):+.4f}")
        ax.set_xlabel("unique ME_face contribution to trial-level Pearson r")
        ax.set_ylabel("# clusters")
        ax.set_title(
            f"ME_face Pearson-r contribution\n"
            f"(Selected − Selected_no_me_face, n={len(finite)} clusters)",
            fontsize=10,
        )
        ax.legend(fontsize=9)

    # Panel 2: per-cluster bar — sorted by Δr
    ax = axes[1]
    ranked = pivot.dropna(subset=["unique_me_face_r"]).sort_values(
        "unique_me_face_r", ascending=False,
    ).reset_index(drop=True)
    bar_colors = ["C2" if s else "C7" for s in ranked["me_face_selected"]]
    ax.bar(np.arange(len(ranked)), ranked["unique_me_face_r"],
           color=bar_colors, edgecolor="none")
    ax.axhline(0.0, color="k", linewidth=0.6)
    ax.set_xlabel("cluster (sorted by Δr)")
    ax.set_ylabel("unique ME_face Δ Pearson r")
    ax.set_title(
        "Per-cluster ME_face contribution (green=ME selected, grey=not selected)",
        fontsize=10,
    )
    out_pdf = OUT_DIR / "variance_partition_pearson.pdf"
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
