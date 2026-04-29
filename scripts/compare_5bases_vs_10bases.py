"""5-basis history vs 10-basis production — does halving the bases pay off?

Reads:
    glm/current/                                    ← 10-basis production (100 ms)
    glm/exploration/history_5bases_100ms/           ← 5-basis at 100 ms

Writes:
    glm/exploration/history_5bases_100ms/comparison_vs_production.csv
        per-cluster: cv_bps_10bases, cv_bps_5bases, delta,
                     selected_vars_{10,5}, history_in_selected_{10,5},
                     pearson_overall_{10,5}

    glm/exploration/history_5bases_100ms/figs/
        cv_bps_paired_comparison.pdf — per-cluster paired lines
        history_kernel_5_vs_10.pdf   — representative-cluster kernel shapes

A small Δ cv_bps (≤ ±0.05 bps) at half the basis count is the win
— same predictive performance with a more interpretable filter shape.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/Users/lauraporta/local_data/motion_clouds/figures")
DIR_10 = ROOT / "glm" / "current"
DIR_5 = ROOT / "glm" / "exploration" / "history_5bases_100ms"
OUT_CSV = DIR_5 / "comparison_vs_production.csv"
OUT_FIGS = DIR_5 / "figs"


def _load_run(d: Path, label: str) -> pd.DataFrame:
    cmp = pd.read_csv(d / "glm_model_comparison.csv")
    out = cmp[["probe_id", "cluster_id",
               "time_Selected_cv_bps", "time_selected_vars"]].copy()
    out = out.rename(columns={
        "time_Selected_cv_bps": f"cv_bps_{label}",
        "time_selected_vars": f"selected_vars_{label}",
    })
    out[f"history_in_{label}"] = out[f"selected_vars_{label}"].apply(
        lambda s: isinstance(s, str) and "History" in s.split("+")
    )
    # Add per-cluster trial-level Pearson r if available
    tl_path = d / "diagnostics" / "trial_level_metrics.csv"
    if tl_path.is_file():
        tl = pd.read_csv(tl_path)
        sel = tl[tl["model"] == "Selected"][
            ["probe_id", "cluster_id", "pearson_r_overall"]
        ].rename(columns={"pearson_r_overall": f"pearson_overall_{label}"})
        out = out.merge(sel, on=["probe_id", "cluster_id"], how="left")
    return out


def _history_coefs(coefs: pd.DataFrame, probe: str, cid: int) -> np.ndarray:
    sub = coefs[
        (coefs["probe_id"] == probe)
        & (coefs["cluster_id"] == cid)
        & (coefs["glm_type"] == "time")
        & (coefs["model"] == "Selected")
        & (coefs["coefficient"].str.startswith("History_"))
    ].copy()
    if sub.empty:
        return np.array([])
    sub["lag"] = sub["coefficient"].str.split("_").str[-1].astype(int)
    sub = sub.sort_values("lag")
    return sub["estimate"].to_numpy()


def render_cv_bps_paired(df: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)

    paired = df.dropna(subset=["cv_bps_10bases", "cv_bps_5bases"])
    n = len(paired)

    # Paired lines
    ax = axes[0]
    for _, row in paired.iterrows():
        ax.plot([0, 1], [row["cv_bps_10bases"], row["cv_bps_5bases"]],
                "-", color="grey", alpha=0.2, linewidth=0.6)
    ax.plot([0, 1],
            [paired["cv_bps_10bases"].median(), paired["cv_bps_5bases"].median()],
            "o-", color="C3", linewidth=3, markersize=10, label="median")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["10 bases\n(production)", "5 bases"])
    ax.set_ylabel("Selected cv_bps (per cluster)")
    ax.set_title(f"Per-cluster paired comparison (n={n})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Δ histogram
    ax = axes[1]
    delta = paired["cv_bps_5bases"] - paired["cv_bps_10bases"]
    ax.hist(delta, bins=30, color="purple", alpha=0.7, edgecolor="black")
    med = delta.median()
    ax.axvline(med, color="black", linewidth=1.5, label=f"median Δ = {med:+.4f}")
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_xlabel("Δ cv_bps (5 bases − 10 bases)")
    ax.set_ylabel("# clusters")
    ax.set_title(f"Δ distribution (improved: {(delta > 0).sum()}/{n}, "
                 f"regressed: {(delta < 0).sum()}/{n})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle("5-basis vs 10-basis history at 100 ms", fontsize=12)
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=150)
    plt.close(fig)


def render_kernel_comparison(out: Path, bin_width_s: float = 0.1,
                             history_window_s: float = 0.2) -> None:
    """6 representative clusters: history kernel SHAPE at 5 vs 10 bases.

    Plotted in lag space (kernel(lag) = sum_k β_k · basis_k(lag)) rather
    than raw β indexed by basis. The β values look different between 5
    and 10 bases (basis-rotation ambiguity — see findings-so-far §1b),
    but the lag-space kernel is the directly-comparable function.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python" / "src"))
    from rc2_glm.basis import history_basis

    coefs_10 = pd.read_csv(DIR_10 / "glm_coefficients.csv")
    coefs_5 = pd.read_csv(DIR_5 / "glm_coefficients.csv")
    ref_clusters = [
        ("CAA-1123243_rec1", 116), ("CAA-1123243_rec1", 339),
        ("CAA-1123244_rec1", 80),  ("CAA-1123466_rec1", 162),
        ("CAA-1123466_rec1", 70),  ("CAA-1123467_rec1", 564),
    ]
    n_lag_bins = max(1, int(round(history_window_s / bin_width_s)))
    lag_ms = np.arange(1, n_lag_bins + 1) * bin_width_s * 1000.0
    B_10 = history_basis(10, history_window_s, bin_width_s)
    B_5 = history_basis(5, history_window_s, bin_width_s)

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), constrained_layout=True)
    for ax, (probe, cid) in zip(axes.flatten(), ref_clusters):
        b10 = _history_coefs(coefs_10, probe, cid)
        b5 = _history_coefs(coefs_5, probe, cid)
        if b10.size:
            kernel_10 = B_10 @ b10
            ax.plot(lag_ms, kernel_10, "o-", color="C0",
                    linewidth=2, markersize=6, label=f"10 bases (n={len(b10)})")
        if b5.size:
            kernel_5 = B_5 @ b5
            ax.plot(lag_ms, kernel_5, "s-", color="C3",
                    linewidth=2, markersize=6, label=f"5 bases (n={len(b5)})")
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.set_title(f"{probe} c{cid}", fontsize=10)
        ax.set_xlabel("lag (ms post-spike)")
        ax.set_ylabel("kernel(lag) — log λ contribution")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle(
        f"History kernel in lag space — 5 vs 10 bases at "
        f"{int(bin_width_s * 1000)} ms (n_lag_bins resolved = {n_lag_bins})",
        fontsize=12,
    )
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=150)
    plt.close(fig)


def main() -> int:
    if not (DIR_5 / "glm_model_comparison.csv").is_file():
        print(f"ERROR: 5-basis run not yet at {DIR_5}", file=sys.stderr)
        return 1
    df_10 = _load_run(DIR_10, "10bases")
    df_5 = _load_run(DIR_5, "5bases")
    df = df_10.merge(df_5, on=["probe_id", "cluster_id"], how="outer")
    df["delta_cv_bps"] = df["cv_bps_5bases"] - df["cv_bps_10bases"]
    if "pearson_overall_10bases" in df.columns and "pearson_overall_5bases" in df.columns:
        df["delta_pearson"] = df["pearson_overall_5bases"] - df["pearson_overall_10bases"]

    OUT_FIGS.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(df)} rows)")

    print()
    print("=== 5 bases vs 10 bases (production) at 100 ms ===")
    paired = df.dropna(subset=["cv_bps_10bases", "cv_bps_5bases"])
    print(f"  paired clusters: {len(paired)}")
    print(f"  median cv_bps:")
    print(f"    10 bases: {paired['cv_bps_10bases'].median():+.4f}")
    print(f"    5  bases: {paired['cv_bps_5bases'].median():+.4f}")
    print(f"    Δ      : {(paired['cv_bps_5bases'] - paired['cv_bps_10bases']).median():+.4f}")
    print(f"  improved: {(paired['cv_bps_5bases'] > paired['cv_bps_10bases']).sum()}/{len(paired)}")
    print(f"  regressed: {(paired['cv_bps_5bases'] < paired['cv_bps_10bases']).sum()}/{len(paired)}")
    print()
    h10 = int(df["history_in_10bases"].sum())
    h5 = int(df["history_in_5bases"].sum())
    print(f"  History selected: 10bases {h10}/{len(df)} ({100*h10/len(df):.0f}%) "
          f"| 5bases {h5}/{len(df)} ({100*h5/len(df):.0f}%)")
    if "delta_pearson" in df.columns:
        dp = df["delta_pearson"].dropna()
        print(f"  trial-level Pearson r Δ (5 − 10 bases): median {dp.median():+.4f}")

    render_cv_bps_paired(df, OUT_FIGS / "cv_bps_paired_comparison.pdf")
    render_kernel_comparison(OUT_FIGS / "history_kernel_5_vs_10.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
