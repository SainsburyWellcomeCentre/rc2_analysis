"""Prompt-03 §Ablation: spike history vs onset kernel — 4-variant comparison.

Runs four pipeline invocations on the 4-probe filtered set, each with a
different combination of {include_history, include_onset_kernel}:

    A — baseline                 : onset ✓, history ✗   (parity-protected)
    B — baseline + history       : onset ✓, history ✓   (prompt-03 main deliverable)
    C — history-only             : onset ✗, history ✓   (history replaces onset)
    D — minus-onset              : onset ✗, history ✗   (sanity ceiling)

Per cluster we extract `time_Selected_cv_bps` from each variant's
`glm_model_comparison.csv` and assemble:

    glm_ablation_history_onset.csv
      probe_id, cluster_id,
      cv_bps_A, cv_bps_B, cv_bps_C, cv_bps_D,
      selected_A, selected_B, selected_C, selected_D,
      delta_history (B - A), delta_history_replaces_onset (C - A)

A summary plot ``ablation_history_vs_onset.pdf`` is rendered alongside
(per-cluster paired lines, with a bold median line).

Usage::

    /path/to/rc2_analysis/python/.../python -m \
        scripts.run_ablation_history_onset

Outputs land in:
    ~/local_data/motion_clouds/figures/ablation_history_onset/
"""
from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("ablation")

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
ROOT = Path("/Users/lauraporta/local_data/motion_clouds")
FORMATTED_DIR = ROOT / "formatted_data"
CLUSTER_FILTER = ROOT / "figures" / "glm_single_cluster" / "prefilter_decision_tree.csv"
OUT_ROOT = ROOT / "figures" / "glm" / "exploration" / "history_ablation"
PROBES = (
    "CAA-1123243_rec1",
    "CAA-1123244_rec1",
    "CAA-1123466_rec1",
    "CAA-1123467_rec1",
)
PYTHON = "/Users/lauraporta/miniforge3/envs/rc2_analysis/bin/python"

VARIANTS: dict[str, dict[str, bool]] = {
    "A_baseline":               {"include_history": False, "include_onset_kernel": True},
    "B_baseline_plus_history":  {"include_history": True,  "include_onset_kernel": True},
    "C_history_only":           {"include_history": True,  "include_onset_kernel": False},
    "D_minus_onset":            {"include_history": False, "include_onset_kernel": False},
}


def run_variant(variant_name: str, flags: dict[str, bool]) -> Path:
    """Run rc2-glm with the variant's config flags. Returns the output dir."""
    out_dir = OUT_ROOT / variant_name
    if (out_dir / "glm_model_comparison.csv").exists():
        log.info("skipping %s — output already present at %s",
                 variant_name, out_dir)
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cli_args = []
    # Defaults flipped 2026-04-29: history is now ON, onset is now OFF.
    # Use the affirmative-opt-out / opt-in flags to override the new
    # defaults explicitly.
    if not flags["include_history"]:
        cli_args.append("--no-history")
    if flags["include_onset_kernel"]:
        cli_args.append("--with-onset-kernel")

    cmd = [
        PYTHON, "-m", "rc2_glm.pipeline",
        # No positionals — both mat (None → iterate FORMATTED_DATA_DIR) and
        # out_dir (None → use RC2_GLM_OUTPUT_DIR env) come from env.
        "--cluster-filter-csv", str(CLUSTER_FILTER),
        "--plot-clusters", "0",  # skip per-cluster plots; we only need CSVs
        *cli_args,
    ]
    import os
    env = os.environ.copy()
    env["RC2_FORMATTED_DATA_DIR"] = str(FORMATTED_DIR)
    env["RC2_GLM_OUTPUT_DIR"] = str(out_dir)
    log.info("running %s — flags=%s", variant_name, flags)
    log.info("  out_dir=%s", out_dir)
    log.info("  CLI: %s", " ".join(cmd))
    proc = subprocess.run(
        cmd, env=env, capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        log.error("variant %s FAILED (exit %d)", variant_name, proc.returncode)
        log.error("stderr tail:\n%s", proc.stderr[-2000:])
        raise RuntimeError(f"{variant_name} pipeline run failed")
    log.info("  done %s", variant_name)
    return out_dir


def aggregate(out_dirs: dict[str, Path]) -> pd.DataFrame:
    """Pull per-cluster Selected CV-bps + selected_vars from each variant."""
    frames = []
    for variant, out_dir in out_dirs.items():
        cmp_path = out_dir / "glm_model_comparison.csv"
        if not cmp_path.exists():
            log.warning("variant %s: missing %s — skipping in aggregation",
                        variant, cmp_path)
            continue
        df = pd.read_csv(cmp_path)[
            ["probe_id", "cluster_id", "time_Selected_cv_bps",
             "time_Null_cv_bps", "time_selected_vars"]
        ].copy()
        df = df.rename(columns={
            "time_Selected_cv_bps": f"cv_bps_{variant[0]}",
            "time_Null_cv_bps": f"null_bps_{variant[0]}",
            "time_selected_vars": f"selected_{variant[0]}",
        })
        frames.append(df)
    out = frames[0]
    for f in frames[1:]:
        out = out.merge(f, on=["probe_id", "cluster_id"], how="outer")
    # Deltas of interest
    out["delta_history"] = out["cv_bps_B"] - out["cv_bps_A"]
    out["delta_history_replaces_onset"] = out["cv_bps_C"] - out["cv_bps_A"]
    out["delta_minus_onset"] = out["cv_bps_D"] - out["cv_bps_A"]
    return out


def render_plot(df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    variants = ["A", "B", "C", "D"]
    labels = ["baseline\n(onset, no hist)",
              "baseline + history\n(onset + hist)",
              "history-only\n(no onset, hist)",
              "minus-onset\n(no onset, no hist)"]
    x = np.arange(4)
    # Per-cluster faint lines
    cv_cols = [f"cv_bps_{v}" for v in variants]
    for _, row in df.iterrows():
        ys = [row[c] for c in cv_cols]
        if any(np.isnan(ys)):
            continue
        ax.plot(x, ys, "-", color="grey", alpha=0.18, linewidth=0.6)
    # Median across clusters per variant
    medians = [np.nanmedian(df[c]) for c in cv_cols]
    ax.plot(x, medians, "o-", color="C3", linewidth=2.5, markersize=10,
            label="median across 88 clusters")
    # Boxplot for distribution
    box_data = [df[c].dropna().to_numpy() for c in cv_cols]
    ax.boxplot(box_data, positions=x, widths=0.4, showfliers=False,
               patch_artist=True,
               boxprops={"facecolor": "C0", "alpha": 0.25},
               medianprops={"color": "C0"})
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Selected-model CV bits/spike (per cluster)")
    ax.set_title(
        "Spike-history × onset-kernel ablation — "
        f"n={int(df[cv_cols].dropna().shape[0])} clusters across 4 probes"
    )
    # Headline summary on the figure
    delta_h = df["delta_history"].median()
    delta_c = df["delta_history_replaces_onset"].median()
    txt = (
        f"median Δ(B−A) = {delta_h:+.4f} bps  ← does history add value on top of onset?\n"
        f"median Δ(C−A) = {delta_c:+.4f} bps  ← can history substitute for onset?"
    )
    ax.text(
        0.02, 0.98, txt, transform=ax.transAxes, fontsize=9,
        verticalalignment="top",
        bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "grey"},
    )
    ax.legend(loc="lower left", fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    log.info("wrote %s", out_path)


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out_dirs: dict[str, Path] = {}
    for name, flags in VARIANTS.items():
        out_dirs[name] = run_variant(name, flags)
    df = aggregate(out_dirs)
    csv_path = OUT_ROOT / "glm_ablation_history_onset.csv"
    df.to_csv(csv_path, index=False)
    log.info("wrote %s (%d rows)", csv_path, len(df))
    render_plot(df, OUT_ROOT / "ablation_history_vs_onset.pdf")
    # Headline numbers in the log
    print()
    print("=== Median per-cluster CV-bps ===")
    for v in ("A", "B", "C", "D"):
        print(f"  {v}: {np.nanmedian(df[f'cv_bps_{v}']):.4f} bps")
    print()
    print(f"=== Per-cluster deltas ===")
    print(f"  median Δ(B - A) [history adds value]:     "
          f"{df['delta_history'].median():+.4f} bps")
    print(f"  median Δ(C - A) [history replaces onset]: "
          f"{df['delta_history_replaces_onset'].median():+.4f} bps")
    print(f"  median Δ(D - A) [no temporal at all]:     "
          f"{df['delta_minus_onset'].median():+.4f} bps")
    return 0


if __name__ == "__main__":
    sys.exit(main())
