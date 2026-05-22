"""Prompt-06 follow-up: 3-way selection-pattern comparison.

Compares the prompt-06 face_me_onset_no_history run against:
  - ``glm/current/`` — 2026-04-29 production (history ON, onset OFF, no ME)
  - ``glm/legacy_with_onset/`` — pre-2026-04-29 (history OFF, onset ON, no ME)

The question this answers — Laura's "how does the hierarchical
classification across all probes change when we add face motion energy?"
— lives in the resulting ``selection_pattern_shift.csv`` and the
sankey-style figure ``selection_shift.pdf``.

Per-cluster columns emitted:
  - ``probe_id``, ``cluster_id``
  - ``selected_vars_legacy`` / ``cv_bps_legacy``
  - ``selected_vars_current`` / ``cv_bps_current``
  - ``selected_vars_face_me`` / ``cv_bps_face_me``
  - ``me_face_selected``, ``me_face_x_speed_selected``
  - ``has_ME_face_main``, ``has_ME_face_x_Speed``

Plus aggregate summaries printed to stdout:
  - ME_face selection rate
  - prefilter category × selected_vars cross-tab
  - cv_bps medians per run

Usage::

    /Users/laura/mambaforge3/envs/rc2_analysis/bin/python -m \\
        scripts.compare_face_me_vs_current
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("compare_face_me")

HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds"
GLM_DIR = ROOT / "figures" / "glm"

LEGACY_CSV = GLM_DIR / "legacy_with_onset" / "glm_model_comparison.csv"
CURRENT_CSV = GLM_DIR / "current" / "glm_model_comparison.csv"
FACE_ME_CSV = (
    GLM_DIR / "exploration" / "face_me_onset_no_history"
    / "glm_model_comparison_aggregated.csv"
)
OUT_DIR = GLM_DIR / "exploration" / "face_me_onset_no_history"


def _load(csv_path: Path, label: str) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"{label} CSV not found at {csv_path}")
    df = pd.read_csv(csv_path)
    keep = ["probe_id", "cluster_id", "time_selected_vars", "time_Selected_cv_bps"]
    missing = [c for c in keep if c not in df.columns]
    if missing:
        raise ValueError(f"{label} CSV missing columns: {missing}")
    df = df[keep].rename(columns={
        "time_selected_vars": f"selected_vars_{label}",
        "time_Selected_cv_bps": f"cv_bps_{label}",
    })
    return df


def merge_three_way() -> pd.DataFrame:
    legacy = _load(LEGACY_CSV, "legacy")
    current = _load(CURRENT_CSV, "current")
    face_me = _load(FACE_ME_CSV, "face_me")
    out = legacy.merge(current, on=["probe_id", "cluster_id"], how="outer")
    out = out.merge(face_me, on=["probe_id", "cluster_id"], how="outer")
    sv = out["selected_vars_face_me"].fillna("Null")
    out["me_face_selected"] = sv.str.contains("ME_face", regex=False, na=False)
    out["me_face_x_speed_selected"] = sv.str.contains(
        "ME_face_x_Speed", regex=False, na=False,
    )
    return out


def _normalise_selection(s: str | float) -> str:
    """Sort selected-var tokens so 'Speed+TF' and 'TF+Speed' compare equal."""
    if pd.isna(s) or s in ("Null", ""):
        return "Null"
    parts = sorted(s.split("+"))
    return "+".join(parts)


def render_sankey(df: pd.DataFrame, out_path: Path) -> None:
    """Approximate sankey: rank top selection groupings; flow lines between runs."""
    df = df.copy()
    for col in ("selected_vars_legacy", "selected_vars_current",
                "selected_vars_face_me"):
        df[col] = df[col].apply(_normalise_selection)

    runs = ("legacy", "current", "face_me")
    fig, ax = plt.subplots(figsize=(11, 7))

    # Use top-K most common selections in each run for tidy labels.
    K = 8

    def _top_k(col: str) -> list[str]:
        return list(df[col].value_counts().head(K).index)

    layers: dict[str, list[str]] = {r: _top_k(f"selected_vars_{r}") for r in runs}
    # Map non-top-K to "other".
    for r in runs:
        col = f"selected_vars_{r}"
        keep = set(layers[r])
        df[col] = df[col].where(df[col].isin(keep), other="other")
        if "other" not in layers[r]:
            layers[r].append("other")

    # Layout: each layer is a vertical column; nodes spaced evenly.
    n_runs = len(runs)
    x_positions = {r: i / (n_runs - 1) for i, r in enumerate(runs)}
    y_pos: dict[tuple[str, str], float] = {}
    for r in runs:
        for i, lab in enumerate(layers[r]):
            y_pos[(r, lab)] = 1.0 - (i / max(1, len(layers[r]) - 1)) * 0.9 - 0.05

    # Draw flows: legacy → current and current → face_me
    pairs = [(runs[0], runs[1]), (runs[1], runs[2])]
    for src, dst in pairs:
        ct = pd.crosstab(df[f"selected_vars_{src}"], df[f"selected_vars_{dst}"])
        for src_lab in layers[src]:
            for dst_lab in layers[dst]:
                n = int(ct.loc[src_lab, dst_lab]) if (
                    src_lab in ct.index and dst_lab in ct.columns
                ) else 0
                if n == 0:
                    continue
                lw = 0.5 + 0.4 * np.sqrt(n)
                ax.plot(
                    [x_positions[src], x_positions[dst]],
                    [y_pos[(src, src_lab)], y_pos[(dst, dst_lab)]],
                    "-", color="grey", alpha=0.4, linewidth=lw,
                )

    # Draw nodes + labels
    for r in runs:
        for lab in layers[r]:
            count = int((df[f"selected_vars_{r}"] == lab).sum())
            ax.scatter(
                x_positions[r], y_pos[(r, lab)],
                s=80 + 4 * count, color="C0", edgecolor="black", zorder=3,
            )
            ax.text(
                x_positions[r] + 0.02, y_pos[(r, lab)],
                f"{lab} (n={count})", fontsize=8, va="center",
            )

    for r in runs:
        ax.text(
            x_positions[r], 1.04, r, fontsize=11, fontweight="bold",
            ha="center",
        )

    ax.set_xlim(-0.15, 1.4)
    ax.set_ylim(-0.05, 1.1)
    ax.set_axis_off()
    ax.set_title(
        "Selection-pattern shift across three configurations\n"
        "legacy (history off, onset on, no ME) → current (history on, onset off, no ME) "
        "→ face_me (history off, onset on, ME_face on)"
    )
    plt.tight_layout()
    fig.savefig(out_path)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    log.info("wrote %s", out_path)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = merge_three_way()
    csv_path = OUT_DIR / "selection_pattern_shift.csv"
    df.to_csv(csv_path, index=False)
    log.info("wrote %s (%d rows)", csv_path, len(df))

    print()
    print("=== Cluster counts (intersection across the 3 runs) ===")
    has_all_3 = (
        df["selected_vars_legacy"].notna()
        & df["selected_vars_current"].notna()
        & df["selected_vars_face_me"].notna()
    )
    print(f"clusters in all 3 runs: {int(has_all_3.sum())}")
    print(f"clusters in face_me only: "
          f"{int(df['selected_vars_face_me'].notna().sum() - has_all_3.sum())}")

    print()
    print("=== ME_face selection rate (face_me run) ===")
    n = int(df["selected_vars_face_me"].notna().sum())
    n_me = int(df["me_face_selected"].sum())
    n_int = int(df["me_face_x_speed_selected"].sum())
    print(f"  ME_face main:        {n_me}/{n} ({100*n_me/max(n,1):.1f}%)")
    print(f"  ME_face_x_Speed:     {n_int}/{n} ({100*n_int/max(n,1):.1f}%)")

    print()
    print("=== Median Selected cv_bps (per cluster, common subset) ===")
    sub = df[has_all_3]
    for run in ("legacy", "current", "face_me"):
        col = f"cv_bps_{run}"
        print(f"  {run:<10}: median = {sub[col].median():+.4f} bps  "
              f"(mean = {sub[col].mean():+.4f}, n = {sub[col].notna().sum()})")

    print()
    print("=== Top-10 most-common selected_vars per run ===")
    for run in ("legacy", "current", "face_me"):
        col = f"selected_vars_{run}"
        print(f"\n  {run}:")
        print(sub[col].apply(_normalise_selection).value_counts().head(10).to_string())

    render_sankey(sub, OUT_DIR / "selection_shift.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
