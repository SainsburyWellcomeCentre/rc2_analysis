"""VF + T = VT additive-tuning test, per cluster — goggles rebuild.

Rebuild of the lost screens ``speed_two_models_per_cluster`` /
``tf_two_models_per_cluster`` figures, for the goggles cohort. For each of the
top-20 Speed+TF-selected clusters (by Δcv_bps in the pooled current_goggles
run), test whether the VT tuning curve is the additive sum of the two unimodal
contributions:

  * Speed axis  (Speed_T_to_VT): observed VT speed-tuning vs
        Model A (gain+offset) : VT ~ a·T(s) + b               [black]
        Model B (additive)     : VT ~ a·T(s) + b·V(ḡ·s) + c    [red]
    where T = T_Vstatic speed-tuning, and V(ḡ·s) is the V-condition TF-tuning
    sampled at the TF the visual flow carries at speed s (TF = ḡ·s).

  * TF axis     (TF_V_to_VT):    observed VT TF-tuning vs
        Model A : VT ~ a·V(f) + b
        Model B : VT ~ a·V(f) + b·T(f/ḡ) + c
    where V = V TF-tuning, and T(f/ḡ) is the T_Vstatic speed-tuning at the
    speed matching TF f.

Per-bin tuning = MEAN across trials from the MATLAB tuning cache
(``csvs/{tuning_curves,tf_tuning_curves}/<probe>.mat`` beside the formatted
file). ḡ = 2/30 Hz/(cm/s), the middle rung of the goggles gain ladder. Model B
beating Model A (and the red sum tracking the data) = additivity holds.

NB: this is a clean reconstruction (the original screens script was lost, only
its gain_offset_per_cluster.csv survived). Conventions chosen 2026-06-09 with
Laura: mean per-bin tuning, ḡ = middle gain, additive 2-regressor Model B.

Usage:
    python scripts/plot_two_models_per_cluster_goggles.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_glm.precomputed_bins import load_precomputed_bin_edges

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("two_models_goggles")

HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds"
FORMATTED_DIR = ROOT / "formatted_data_goggles"
POOLED_CMP = ROOT / "figures" / "glm" / "current_goggles" / "glm_model_comparison.csv"
OUT_DIR = ROOT / "figures" / "glm" / "exploration" / "subspace_population_goggles"

GBAR = 2.0 / 30.0   # Hz per cm/s — middle rung of the goggles gain ladder
N_TOP = 20


# ----------------------------------------------------------------------
def select_clusters() -> list[tuple[str, int]]:
    """Top-N (probe, cluster_id) selected for BOTH Speed and TF, by Δcv_bps."""
    df = pd.read_csv(POOLED_CMP)
    sv = df["time_selected_vars"].fillna("").astype(str)
    both = df[sv.str.contains(r"\bSpeed\b") & sv.str.contains(r"\bTF\b")].copy()
    both = both.sort_values("time_delta_selected_vs_null", ascending=False).head(N_TOP)
    return [(r.probe_id, int(r.cluster_id)) for r in both.itertuples()]


def _mean_tuning(pre, cond: str, cid: int, var: str):
    """(centres, mean-over-trials tuning) for one condition/cluster/variable."""
    tun = pre.speed_tuning(cond, cid) if var == "Speed" else pre.tf_tuning(cond, cid)
    cen = pre.speed_centres(cond) if var == "Speed" else pre.tf_centres(cond)
    if tun is None or cen is None or np.size(tun) == 0:
        return None, None
    return np.asarray(cen, dtype=float), np.nanmean(tun, axis=0)


def _r2(y: np.ndarray, pred: np.ndarray) -> float:
    ok = np.isfinite(y) & np.isfinite(pred)
    if ok.sum() < 2:
        return float("nan")
    yv = y[ok]
    ss_res = np.sum((yv - pred[ok]) ** 2)
    ss_tot = np.sum((yv - np.mean(yv)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def _fit(cols: list[np.ndarray], y: np.ndarray):
    """Least-squares y ~ [cols..., 1]; returns (coef, prediction-on-full-grid, r2)."""
    design = cols + [np.ones_like(y)]
    stack = np.vstack(design).T            # (n_bins, n_terms)
    ok = np.all(np.isfinite(stack), axis=1) & np.isfinite(y)
    if ok.sum() < stack.shape[1] + 1:
        return None
    coef, *_ = np.linalg.lstsq(stack[ok], y[ok], rcond=None)
    pred = stack @ coef                    # full grid (NaN where any col NaN)
    return coef, pred, _r2(y, pred)


def _models_for(pre, cid: int, axis: str):
    """Build observed VT + Model A/B predictions on the VT grid for one axis.

    axis='Speed' → Speed_T_to_VT; axis='TF' → TF_V_to_VT.
    Returns (grid, y_obs, predA, r2A, predB, r2B) or None.
    """
    if axis == "Speed":
        grid, y = _mean_tuning(pre, "VT", cid, "Speed")          # target
        _, uni = _mean_tuning(pre, "T_Vstatic", cid, "Speed")    # unimodal T(s)
        cV, V = _mean_tuning(pre, "V", cid, "TF")                # cross V(f)
        if grid is None or uni is None or V is None:
            return None
        matched = np.interp(GBAR * grid, cV, V, left=np.nan, right=np.nan)
    else:  # TF
        grid, y = _mean_tuning(pre, "VT", cid, "TF")             # target
        _, uni = _mean_tuning(pre, "V", cid, "TF")              # unimodal V(f)
        cT, T = _mean_tuning(pre, "T_Vstatic", cid, "Speed")    # cross T(s)
        if grid is None or uni is None or T is None:
            return None
        matched = np.interp(grid / GBAR, cT, T, left=np.nan, right=np.nan)

    a = _fit([uni], y)                    # Model A: gain+offset
    b = _fit([uni, matched], y)           # Model B: additive (2 regressors)
    if a is None or b is None:
        return None
    return grid, y, a[1], a[2], b[1], b[2]


# ----------------------------------------------------------------------
def render(axis: str, clusters: list[tuple[str, int]], pres: dict) -> pd.DataFrame:
    """Per-cluster grid for one axis; returns the per-cluster r2 table."""
    n = len(clusters)
    ncol = 5
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.6 * nrow))
    axes = np.atleast_1d(axes).ravel()
    rows = []
    xlabel = "speed (cm/s)" if axis == "Speed" else "TF (Hz)"
    for ax, (probe, cid) in zip(axes, clusters):
        res = _models_for(pres[probe], cid, axis)
        if res is None:
            ax.set_visible(False)
            continue
        grid, y, predA, r2A, predB, r2B = res
        ax.plot(grid, y, "o", ms=3.5, color="0.25", label="VT observed")
        ax.plot(grid, predA, "-", color="black", lw=1.6, label=f"gain+offset r²={r2A:.2f}")
        ax.plot(grid, predB, "-", color="tab:red", lw=1.6, label=f"additive r²={r2B:.2f}")
        ax.set_title(f"{probe[:11]} cl{cid}", fontsize=8)
        ax.set_xlabel(xlabel, fontsize=7)
        ax.set_ylabel("FR (Hz)", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=5.5, loc="best")
        rows.append({"probe": probe, "cluster": cid,
                     "analysis": "Speed_T_to_VT" if axis == "Speed" else "TF_V_to_VT",
                     "r2_gain_offset": r2A, "r2_additive": r2B})
    for ax in axes[n:]:
        ax.set_visible(False)
    fig.suptitle(
        f"VF+T = VT — {axis} tuning, top {n} Speed+TF clusters (goggles)\n"
        f"black = gain+offset (Model A) · red = additive T+V (Model B)",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = "speed_two_models_per_cluster" if axis == "Speed" else "tf_two_models_per_cluster"
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}", dpi=150)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", OUT_DIR / stem)
    return pd.DataFrame(rows)


def main() -> int:
    clusters = select_clusters()
    log.info("top %d Speed+TF clusters: %s", len(clusters),
             ", ".join(f"{p[:11]}/{c}" for p, c in clusters))
    pres = {probe: load_precomputed_bin_edges(FORMATTED_DIR / f"{probe}.mat")
            for probe in {p for p, _ in clusters}}

    tables = [render(axis, clusters, pres) for axis in ("Speed", "TF")]
    out = pd.concat(tables, ignore_index=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_DIR / "two_models_r2_per_cluster.csv", index=False)
    log.info("wrote %s", OUT_DIR / "two_models_r2_per_cluster.csv")

    print("\n=== median R² (B beats A count) ===")
    for analysis, g in out.groupby("analysis"):
        beats = int((g["r2_additive"] > g["r2_gain_offset"]).sum())
        print(f"  {analysis:14s}: gain+offset {g['r2_gain_offset'].median():.3f} → "
              f"additive {g['r2_additive'].median():.3f}  (B>A {beats}/{len(g)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
