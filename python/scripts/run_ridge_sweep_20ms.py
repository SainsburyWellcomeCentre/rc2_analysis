"""20 ms ridge-sweep — does heavier regularisation tame the History-filter noise?

At 20 ms bin width the History filter resolves all 10 lag bins distinctly,
but with the production ridge λ=10⁻³ (tuned for 100 ms) the fitted
coefficients oscillate between ±5 — classic under-regularisation for the
higher-dim fit. This script sweeps λ ∈ {1e-3, 1e-2, 1e-1, 1.0} at
20 ms across all 4 probes (88 clusters) and emits a comparison table +
a representative-clusters plot.

Outputs:
    glm/exploration/ridge_sweep_20ms/
        lambda_1e-03/  lambda_1e-02/  lambda_1e-01/  lambda_1e+00/
            (full per-probe rc2-glm output for each λ)
        comparison.csv             — per-cluster CV-bps per λ
        history_shape_per_lambda.pdf  — 6 representative clusters

Each λ is its own ~25 min run; total wall-clock ~100 min. Runs sequentially
because each rc2-glm already saturates all cores via joblib.
"""
from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("ridge_sweep")

ROOT = Path("/Users/lauraporta/local_data/motion_clouds")
FORMATTED_DIR = ROOT / "formatted_data"
CLUSTER_FILTER = ROOT / "figures" / "matlab_reference" / "prefilter_decision_tree.csv"
OUT_ROOT = ROOT / "figures" / "glm" / "exploration" / "ridge_sweep_20ms"
PYTHON = "/Users/lauraporta/miniforge3/envs/rc2_analysis/bin/python"

LAMBDAS: tuple[float, ...] = (1e-3, 1e-2, 1e-1, 1.0)


def _lambda_dirname(lam: float) -> str:
    # Render lambdas like "lambda_1e-03" so the dirs sort lexically.
    return f"lambda_{lam:.0e}".replace("+0", "+").replace("-0", "-")


def run_one(lam: float) -> Path:
    out_dir = OUT_ROOT / _lambda_dirname(lam)
    if (out_dir / "glm_model_comparison.csv").exists():
        log.info("λ=%g already run — skipping", lam)
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        PYTHON, "-m", "rc2_glm.pipeline",
        "--bin-width", "0.02",
        "--lambda-ridge", str(lam),
        "--cluster-filter-csv", str(CLUSTER_FILTER),
        "--plot-clusters", "0",  # skip per-cluster plots
    ]
    env = os.environ.copy()
    env["RC2_FORMATTED_DATA_DIR"] = str(FORMATTED_DIR)
    env["RC2_GLM_OUTPUT_DIR"] = str(out_dir)
    log.info("λ=%g — running %s", lam, " ".join(cmd))
    proc = subprocess.run(
        cmd, env=env, capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        log.error("λ=%g FAILED (exit %d). stderr tail:\n%s",
                  lam, proc.returncode, proc.stderr[-1500:])
        raise RuntimeError(f"λ={lam} failed")
    log.info("  done λ=%g", lam)
    return out_dir


def aggregate(out_dirs: dict[float, Path]) -> pd.DataFrame:
    """Per-cluster table with one column per λ."""
    frames = []
    for lam, d in out_dirs.items():
        cmp_path = d / "glm_model_comparison.csv"
        if not cmp_path.is_file():
            log.warning("λ=%g: missing %s", lam, cmp_path)
            continue
        df = pd.read_csv(cmp_path)[
            ["probe_id", "cluster_id", "time_Selected_cv_bps",
             "time_selected_vars"]
        ].rename(columns={
            "time_Selected_cv_bps": f"cv_bps_lam_{lam:g}",
            "time_selected_vars": f"selected_vars_lam_{lam:g}",
        })
        frames.append(df)
    out = frames[0]
    for f in frames[1:]:
        out = out.merge(f, on=["probe_id", "cluster_id"], how="outer")
    return out


def _history_coefs(coefs: pd.DataFrame, probe: str, cluster_id: int) -> np.ndarray:
    sub = coefs[
        (coefs["probe_id"] == probe)
        & (coefs["cluster_id"] == cluster_id)
        & (coefs["glm_type"] == "time")
        & (coefs["model"] == "Selected")
        & (coefs["coefficient"].str.startswith("History_"))
    ].copy()
    if sub.empty:
        return np.array([])
    sub["lag"] = sub["coefficient"].str.split("_").str[-1].astype(int)
    sub = sub.sort_values("lag")
    return sub["estimate"].to_numpy()


def render_history_shape_per_lambda(out_dirs: dict[float, Path]) -> None:
    """Pick 6 reference clusters and overlay history filter shape per λ."""
    # Use 100 ms output as canonical 'reference' cluster ID list (spans
    # the cv_bps range).
    ref_clusters = [
        ("CAA-1123243_rec1", 116),
        ("CAA-1123243_rec1", 339),
        ("CAA-1123244_rec1", 80),
        ("CAA-1123466_rec1", 162),
        ("CAA-1123466_rec1", 70),
        ("CAA-1123467_rec1", 564),
    ]
    # Load coefs for each λ
    coefs_by_lam: dict[float, pd.DataFrame] = {}
    for lam, d in out_dirs.items():
        p = d / "glm_coefficients.csv"
        if p.is_file():
            coefs_by_lam[lam] = pd.read_csv(p)

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for ax, (probe, cid) in zip(axes.flatten(), ref_clusters):
        for i, (lam, coefs) in enumerate(sorted(coefs_by_lam.items())):
            beta = _history_coefs(coefs, probe, cid)
            if beta.size == 0:
                continue
            x = np.arange(1, len(beta) + 1)
            ax.plot(
                x, beta, "o-",
                color=cmap(i / max(len(coefs_by_lam) - 1, 1)),
                label=f"λ={lam:g}",
                linewidth=1.5, markersize=4,
            )
        ax.axhline(0, color="grey", linewidth=0.6, alpha=0.5)
        ax.set_title(f"{probe} c{cid}", fontsize=10)
        ax.set_xlabel("History basis index", fontsize=8)
        ax.set_ylabel("β coefficient", fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    fig.suptitle("History-filter shape vs ridge λ at 20 ms", fontsize=12)
    out = OUT_ROOT / "history_shape_per_lambda.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=150)
    plt.close(fig)
    log.info("wrote %s", out)


def render_cv_bps_summary(df: pd.DataFrame) -> None:
    """Box-plot per λ of per-cluster cv_bps."""
    cols = sorted([c for c in df.columns if c.startswith("cv_bps_lam_")],
                  key=lambda s: float(s.replace("cv_bps_lam_", "")))
    if not cols:
        return
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    data = [df[c].dropna().to_numpy() for c in cols]
    labels = [c.replace("cv_bps_lam_", "λ=") for c in cols]
    ax.boxplot(data, tick_labels=labels, showfliers=False,
               patch_artist=True,
               boxprops={"facecolor": "C0", "alpha": 0.3},
               medianprops={"color": "C0"})
    medians = [np.nanmedian(d) for d in data]
    ax.plot(range(1, len(data) + 1), medians, "o-", color="C3", linewidth=2,
            label="median")
    ax.set_ylabel("Selected CV-bps (per cluster)")
    ax.set_title(f"Ridge sweep at 20 ms — {len(df)} clusters across 4 probes")
    ax.set_xlabel("ridge λ")
    ax.legend()
    ax.grid(True, alpha=0.3)
    out = OUT_ROOT / "cv_bps_per_lambda.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=150)
    plt.close(fig)
    log.info("wrote %s", out)


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out_dirs: dict[float, Path] = {}
    for lam in LAMBDAS:
        out_dirs[lam] = run_one(lam)
    df = aggregate(out_dirs)
    df.to_csv(OUT_ROOT / "comparison.csv", index=False)
    log.info("wrote %s (%d rows)", OUT_ROOT / "comparison.csv", len(df))

    print()
    print("=== median CV-bps per λ ===")
    for lam in sorted(out_dirs):
        col = f"cv_bps_lam_{lam:g}"
        if col in df.columns:
            print(f"  λ={lam:g}: median = {df[col].median():+.4f} bps")

    render_cv_bps_summary(df)
    render_history_shape_per_lambda(out_dirs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
