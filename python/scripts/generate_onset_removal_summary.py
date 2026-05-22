"""Per-cluster CV-bps comparison: legacy with-onset vs current no-onset+history.

Writes ``glm/current/onset_removal_summary.csv``:
    probe_id, cluster_id,
    cv_bps_legacy_with_onset, cv_bps_current_no_onset_history,
    delta_cv_bps, sign,
    selected_vars_legacy, selected_vars_current

The "expected median improvement matches ablation variant B (≈ +0.07 bps)"
prediction in the plan can be directly checked from this CSV.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

LEGACY_CSV = Path(
    "/Users/lauraporta/local_data/motion_clouds/figures/glm/legacy_with_onset"
    "/glm_model_comparison.csv"
)
CURRENT_CSV = Path(
    "/Users/lauraporta/local_data/motion_clouds/figures/glm/current"
    "/glm_model_comparison.csv"
)
OUT_CSV = CURRENT_CSV.parent / "onset_removal_summary.csv"


def _load(path: Path, suffix: str) -> pd.DataFrame:
    if not path.is_file():
        raise SystemExit(f"missing: {path}")
    df = pd.read_csv(path)
    return df[["probe_id", "cluster_id",
               "time_Selected_cv_bps", "time_selected_vars"]].rename(
        columns={"time_Selected_cv_bps": f"cv_bps_{suffix}",
                 "time_selected_vars": f"selected_vars_{suffix}"}
    )


def main() -> int:
    legacy = _load(LEGACY_CSV, "legacy_with_onset")
    current = _load(CURRENT_CSV, "current_no_onset_history")

    df = legacy.merge(current, on=["probe_id", "cluster_id"], how="outer")
    df["delta_cv_bps"] = (
        df["cv_bps_current_no_onset_history"] - df["cv_bps_legacy_with_onset"]
    )
    df["sign"] = np.where(
        df["delta_cv_bps"] > 0, "improved",
        np.where(df["delta_cv_bps"] < 0, "regressed", "unchanged"),
    )
    df = df[[
        "probe_id", "cluster_id",
        "cv_bps_legacy_with_onset", "cv_bps_current_no_onset_history",
        "delta_cv_bps", "sign",
        "selected_vars_legacy_with_onset", "selected_vars_current_no_onset_history",
    ]]
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(df)} rows)")
    print()
    n_paired = int(df["delta_cv_bps"].notna().sum())
    print(f"=== summary across {n_paired} paired clusters ===")
    deltas = df["delta_cv_bps"].dropna()
    print(f"  median Δ = {deltas.median():+.4f} bps")
    print(f"  mean Δ   = {deltas.mean():+.4f} bps")
    print(f"  improved: {(deltas > 0).sum()}/{n_paired}")
    print(f"  regressed: {(deltas < 0).sum()}/{n_paired}")
    print(f"  exactly equal: {(deltas == 0).sum()}/{n_paired}")
    print()
    print("=== top 5 improvements ===")
    print(df.nlargest(5, "delta_cv_bps")[
        ["probe_id", "cluster_id", "delta_cv_bps", "selected_vars_current_no_onset_history"]
    ].to_string(index=False))
    print()
    print("=== top 5 regressions ===")
    print(df.nsmallest(5, "delta_cv_bps")[
        ["probe_id", "cluster_id", "delta_cv_bps", "selected_vars_current_no_onset_history"]
    ].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
