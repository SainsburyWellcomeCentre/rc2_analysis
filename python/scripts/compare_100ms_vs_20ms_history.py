"""Compare 100 ms (production) vs 20 ms (Phase E sweep) history runs.

Reads:
  /Users/lauraporta/local_data/motion_clouds/figures/glm/current/
    glm_{model_comparison,coefficients}.csv                    ← 100 ms baseline
  /Users/lauraporta/local_data/motion_clouds/figures/glm/exploration/
    bin_width_20ms_history/glm_{model_comparison,coefficients}.csv  ← 20 ms

Writes:
  ${20ms_dir}/comparison_100ms_vs_20ms.csv     ← per-cluster table
  ${20ms_dir}/history_shape_comparison.pdf     ← representative-clusters plot

Per-cluster columns:
  probe_id, cluster_id,
  cv_bps_100ms, cv_bps_20ms, delta_cv_bps,
  selected_vars_100ms, selected_vars_20ms,
  history_in_selected_100ms, history_in_selected_20ms,
  hist_lag1_coef_20ms, hist_peak_lag_idx_20ms

Headline numbers printed to stdout.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/Users/lauraporta/local_data/motion_clouds/figures")
DIR_100MS = ROOT / "glm" / "current"
DIR_20MS = ROOT / "glm" / "exploration" / "bin_width_20ms_history"

OUT_CSV = DIR_20MS / "comparison_100ms_vs_20ms.csv"
OUT_PDF = DIR_20MS / "history_shape_comparison.pdf"


def _history_coefs(coefs: pd.DataFrame, probe: str, cluster_id: int) -> np.ndarray:
    """Return the History_1 ... History_N coefficient vector for one cluster.

    Returns an empty array if no History_* columns are present (i.e. history
    wasn't selected for this cluster's Selected model).
    """
    sub = coefs[
        (coefs["probe_id"] == probe)
        & (coefs["cluster_id"] == cluster_id)
        & (coefs["glm_type"] == "time")
        & (coefs["coefficient"].str.startswith("History_"))
    ].copy()
    if sub.empty:
        return np.array([])
    sub["lag"] = sub["coefficient"].str.split("_").str[-1].astype(int)
    sub = sub.sort_values("lag")
    return sub["estimate"].to_numpy()


def _history_in_selected(selected_vars: str) -> bool:
    if not isinstance(selected_vars, str):
        return False
    return "History" in selected_vars.split("+")


def _load_pair() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cmp_100 = pd.read_csv(DIR_100MS / "glm_model_comparison.csv")
    cmp_20 = pd.read_csv(DIR_20MS / "glm_model_comparison.csv")
    coef_100 = pd.read_csv(DIR_100MS / "glm_coefficients.csv")
    coef_20 = pd.read_csv(DIR_20MS / "glm_coefficients.csv")
    return cmp_100, cmp_20, coef_100, coef_20


def main() -> int:
    if not DIR_20MS.is_dir() or not (DIR_20MS / "glm_model_comparison.csv").is_file():
        print(f"ERROR: 20 ms run output missing at {DIR_20MS}", file=sys.stderr)
        return 1
    cmp_100, cmp_20, coef_100, coef_20 = _load_pair()

    # Per-cluster merge on (probe_id, cluster_id)
    df = cmp_100[[
        "probe_id", "cluster_id",
        "time_Selected_cv_bps", "time_selected_vars",
    ]].rename(columns={
        "time_Selected_cv_bps": "cv_bps_100ms",
        "time_selected_vars": "selected_vars_100ms",
    })
    df_20 = cmp_20[[
        "probe_id", "cluster_id",
        "time_Selected_cv_bps", "time_selected_vars",
    ]].rename(columns={
        "time_Selected_cv_bps": "cv_bps_20ms",
        "time_selected_vars": "selected_vars_20ms",
    })
    df = df.merge(df_20, on=["probe_id", "cluster_id"], how="outer")
    df["delta_cv_bps"] = df["cv_bps_20ms"] - df["cv_bps_100ms"]
    df["history_in_selected_100ms"] = df["selected_vars_100ms"].apply(_history_in_selected)
    df["history_in_selected_20ms"] = df["selected_vars_20ms"].apply(_history_in_selected)

    # Per-row (lag-1 + peak) summary of the 20ms history shape
    lag1: list[float] = []
    peak_lag: list[int | None] = []
    for _, row in df.iterrows():
        coefs = _history_coefs(coef_20, row["probe_id"], int(row["cluster_id"]))
        if coefs.size == 0:
            lag1.append(np.nan)
            peak_lag.append(None)
            continue
        lag1.append(float(coefs[0]))
        # Peak = strongest absolute deviation from 0
        peak_lag.append(int(np.argmax(np.abs(coefs))))
    df["hist_lag1_coef_20ms"] = lag1
    df["hist_peak_lag_idx_20ms"] = peak_lag

    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(df)} rows)")

    # ----- headlines -----
    print()
    print(f"=== 4-probe summary across {len(df)} clusters ===")
    n_paired = int(df["delta_cv_bps"].notna().sum())
    print(f"Paired (both runs returned a Selected fit): {n_paired}/{len(df)}")
    deltas = df["delta_cv_bps"].dropna()
    print(f"  median Δ cv_bps (20 ms − 100 ms): {deltas.median():+.4f}")
    print(f"  mean Δ cv_bps                  : {deltas.mean():+.4f}")
    print(f"  20ms-improved: {(deltas > 0).sum()}/{n_paired}")
    print(f"  20ms-regressed: {(deltas < 0).sum()}/{n_paired}")
    print()
    h100 = df["history_in_selected_100ms"].sum()
    h20 = df["history_in_selected_20ms"].sum()
    print(f"History selected: 100ms {h100}/{len(df)} ({100*h100/len(df):.1f}%)  "
          f"vs 20ms {h20}/{len(df)} ({100*h20/len(df):.1f}%)")
    print()
    print("=== top 5 changes (20 ms helped) ===")
    print(df.nlargest(5, "delta_cv_bps")[
        ["probe_id", "cluster_id", "delta_cv_bps",
         "selected_vars_100ms", "selected_vars_20ms"]
    ].to_string(index=False))
    print()
    print("=== top 5 changes (20 ms hurt) ===")
    print(df.nsmallest(5, "delta_cv_bps")[
        ["probe_id", "cluster_id", "delta_cv_bps",
         "selected_vars_100ms", "selected_vars_20ms"]
    ].to_string(index=False))

    # ----- history-shape comparison plot -----
    # Pick 6 representative clusters: 2 from each of the highest |Δ|, lowest
    # |Δ|, and middle, with history selected at 20 ms.
    have_hist = df[df["history_in_selected_20ms"]].copy()
    if len(have_hist) >= 6:
        have_hist["abs_delta"] = have_hist["delta_cv_bps"].abs()
        top = have_hist.nlargest(2, "abs_delta")
        bot = have_hist.nsmallest(2, "abs_delta")
        med = have_hist.iloc[[len(have_hist) // 3, 2 * len(have_hist) // 3]]
        sample = pd.concat([top, med, bot]).drop_duplicates(["probe_id", "cluster_id"])
    else:
        sample = have_hist.head(6)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=False)
    axes = axes.flatten()
    for ax, (_, row) in zip(axes, sample.iterrows()):
        coefs_20 = _history_coefs(coef_20, row["probe_id"], int(row["cluster_id"]))
        coefs_100 = _history_coefs(coef_100, row["probe_id"], int(row["cluster_id"]))
        if coefs_20.size:
            x_20 = np.arange(1, len(coefs_20) + 1)
            ax.plot(x_20, coefs_20, "o-", color="C3",
                    label=f"20ms (n={len(coefs_20)} bases)")
        if coefs_100.size:
            x_100 = np.arange(1, len(coefs_100) + 1)
            ax.plot(x_100, coefs_100, "s-", color="C0", alpha=0.6,
                    label=f"100ms (n={len(coefs_100)} bases)")
        ax.axhline(0, color="grey", linewidth=0.5, alpha=0.5)
        ax.set_title(
            f"{row['probe_id']} cluster {int(row['cluster_id'])}\n"
            f"Δcv_bps={row['delta_cv_bps']:+.3f}",
            fontsize=10,
        )
        ax.set_xlabel("History basis index (lag bin)")
        ax.set_ylabel("β coefficient")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    plt.suptitle("Spike-history filter shape — 100 ms vs 20 ms bin width",
                 fontsize=13)
    plt.tight_layout()
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PDF.with_suffix(".png"), dpi=150)
    print()
    print(f"wrote {OUT_PDF}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
