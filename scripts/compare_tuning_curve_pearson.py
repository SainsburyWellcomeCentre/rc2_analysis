"""Per-cluster tuning-curve Pearson r — current (no-onset+history) vs legacy (with-onset).

For each (probe, cluster, variable, condition, model), computes:
    Pearson r between the model's predicted tuning curve (resampled onto
    the MATLAB-cache 5%-quantile bins) and the cache's per-trial-mean
    observed firing rate at those bins.

Aggregates per-cluster means and emits:
    glm/current/diagnostics/tuning_pearson_by_cluster.csv
        probe_id, cluster_id, variable, condition, model,
        pearson_r_current, pearson_r_legacy, delta_pearson, n_bins

    glm/current/figs/tuning_pearson_comparison.pdf
        - histogram of per-cluster mean Pearson r (current vs legacy)
        - scatter (legacy_r, current_r) coloured by variable
        - per-variable summary table panel

This addresses the concern that the new no-onset+history model improved
trial-level prediction (cv_bps + autocorrelation) but degraded tuning-
curve fits, by quantifying the degradation per (cluster, variable).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python" / "src"))
from rc2_glm.precomputed_bins import load_precomputed_bin_edges

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("tuning_pearson")

ROOT = Path("/Users/lauraporta/local_data/motion_clouds")
FORMATTED = ROOT / "formatted_data"
CURRENT = ROOT / "figures" / "glm" / "current"
LEGACY = ROOT / "figures" / "glm" / "legacy_with_onset"

# Cache .mat paths — the same precomputed per-trial-per-bin firing-rate
# cache the Observed row of cluster_<id>_tuning.pdf reads from.
CACHE_DIR = FORMATTED / "csvs"


_CACHE_BY_PROBE: dict[str, object] = {}


def _get_precomputed(probe_id: str):
    if probe_id not in _CACHE_BY_PROBE:
        # The MATLAB cache lives next to formatted_data .mat files.
        mat_path = FORMATTED / f"{probe_id}.mat"
        try:
            _CACHE_BY_PROBE[probe_id] = load_precomputed_bin_edges(mat_path)
        except Exception as exc:
            log.warning("cache load failed for %s: %s", probe_id, exc)
            _CACHE_BY_PROBE[probe_id] = None
    return _CACHE_BY_PROBE[probe_id]


def _load_cache(probe_id: str, var: str) -> dict:
    """Return {cond: {'centres': np.ndarray, 'tuning': {cid: 2D arr}}}."""
    pc = _get_precomputed(probe_id)
    if pc is None:
        return {}
    by_cond: dict[str, dict] = {}
    for cond in ("VT", "V", "T_Vstatic"):
        if var == "Speed":
            centres = pc.speed_centres(cond)
            cluster_ids = pc.speed_by_group.get(cond)
        else:
            centres = pc.tf_centres(cond)
            cluster_ids = pc.tf_by_group.get(cond)
        if centres is None or cluster_ids is None:
            continue
        # tuning lives at .tuning dict on the per-condition record
        per_cluster = dict(cluster_ids.tuning)
        by_cond[cond] = {"centres": np.asarray(centres), "tuning": per_cluster}
    return by_cond


def _load_predicted(run_dir: Path, probe_id: str) -> pd.DataFrame:
    p = run_dir / "_runs" / probe_id / "diagnostics" / "tuning_curves.csv"
    if not p.is_file():
        return pd.DataFrame()
    return pd.read_csv(p)


def _pearson_for_cluster(
    pred_df: pd.DataFrame,
    cache_by_cond: dict,
    probe: str, cid: int, var: str,
) -> list[dict]:
    """Compute Pearson r per (model, condition) for one cluster/variable.

    Predicted is the gain on a 100-point grid (`grid_x`, `gain`); we
    interpolate onto the cache's bin centres and correlate against the
    per-trial-mean observed firing rate at the same bins.
    """
    rows: list[dict] = []
    sub_pred = pred_df[
        (pred_df["cluster_id"] == cid) & (pred_df["variable"] == var)
    ]
    for cond, cache in cache_by_cond.items():
        per_cluster = cache["tuning"].get(int(cid))
        if per_cluster is None:
            continue
        centres = cache["centres"]
        observed = np.nanmean(per_cluster, axis=0)  # mean across trials
        # Skip bins where no trials contributed
        valid = ~np.isnan(observed) & np.isfinite(observed)
        if valid.sum() < 4:
            continue
        observed_v = observed[valid]
        centres_v = centres[valid]
        for model in sub_pred["model"].unique():
            sm = sub_pred[sub_pred["model"] == model]
            if len(sm) < 2:
                continue
            grid = sm["grid_x"].to_numpy()
            gain = sm["gain"].to_numpy()
            order = np.argsort(grid)
            grid = grid[order]; gain = gain[order]
            if not np.all(np.diff(grid) > 0):
                # collapse duplicates
                _, uidx = np.unique(grid, return_index=True)
                grid = grid[uidx]; gain = gain[uidx]
            pred_at_centres = np.interp(centres_v, grid, gain)
            try:
                r, _ = pearsonr(observed_v, pred_at_centres)
            except (ValueError, RuntimeWarning):
                continue
            rows.append({
                "probe_id": probe, "cluster_id": int(cid),
                "variable": var, "condition": cond, "model": model,
                "pearson_r": float(r), "n_bins": int(valid.sum()),
            })
    return rows


def collect_run(run_dir: Path, label: str) -> pd.DataFrame:
    rows: list[dict] = []
    for probe_dir in sorted((run_dir / "_runs").iterdir()):
        if not probe_dir.is_dir():
            continue
        probe = probe_dir.name
        pred = _load_predicted(run_dir, probe)
        if pred.empty:
            log.warning("%s: no predicted tuning_curves.csv at %s", label, probe)
            continue
        for var in ("Speed", "TF"):
            cache = _load_cache(probe, var)
            if not cache:
                log.warning("%s: %s cache missing for %s", label, var, probe)
                continue
            for cid in pred["cluster_id"].unique():
                rows.extend(_pearson_for_cluster(pred, cache, probe, int(cid), var))
        log.info("%s probe %s: %d (cluster, var, cond, model) rows so far",
                 label, probe, len(rows))
    df = pd.DataFrame(rows)
    df["run"] = label
    return df


def render_plot(merged: pd.DataFrame, out_pdf: Path) -> None:
    """Histogram + scatter + per-variable summary."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)

    # row 0: histograms of per-cluster mean Pearson r per variable
    for col_idx, var in enumerate(("Speed", "TF")):
        ax = axes[0, col_idx]
        for label, color in (("legacy", "C0"), ("current", "C3")):
            mask = (merged["run"] == label) & (merged["variable"] == var) & (merged["model"] == "Selected")
            data = merged.loc[mask].groupby(["probe_id", "cluster_id"])["pearson_r"].mean()
            data = data.dropna().to_numpy()
            if data.size == 0:
                continue
            ax.hist(data, bins=25, alpha=0.6, color=color, label=f"{label} (median r={np.nanmedian(data):.2f})")
        ax.set_title(f"{var} — Selected model")
        ax.set_xlabel("per-cluster mean Pearson r")
        ax.set_ylabel("# clusters")
        ax.set_xlim(-0.5, 1.0)
        ax.axvline(0, color="grey", linewidth=0.5)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # row 0 col 2: legacy vs current scatter (Speed + TF combined)
    ax = axes[0, 2]
    pivot = merged[merged["model"] == "Selected"].groupby(
        ["probe_id", "cluster_id", "variable", "run"]
    )["pearson_r"].mean().unstack("run").reset_index()
    if "current" in pivot.columns and "legacy" in pivot.columns:
        for var, color in (("Speed", "green"), ("TF", "orange")):
            sub = pivot[pivot["variable"] == var].dropna()
            if not sub.empty:
                ax.scatter(sub["legacy"], sub["current"], color=color, alpha=0.6, label=var, s=25)
    ax.plot([-0.5, 1.0], [-0.5, 1.0], "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("legacy with-onset Pearson r")
    ax.set_ylabel("current no-onset+history Pearson r")
    ax.set_xlim(-0.5, 1.0); ax.set_ylim(-0.5, 1.0)
    ax.set_title("per-cluster scatter (Selected)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # row 1: per-condition stratified means
    for col_idx, var in enumerate(("Speed", "TF")):
        ax = axes[1, col_idx]
        for cond_idx, cond in enumerate(("VT", "V", "T_Vstatic")):
            for label, color, hatch in (("legacy", "C0", ""), ("current", "C3", "//")):
                mask = (merged["run"] == label) & (merged["variable"] == var) & (
                    merged["condition"] == cond) & (merged["model"] == "Selected")
                data = merged.loc[mask, "pearson_r"].dropna()
                if len(data) == 0:
                    continue
                pos = cond_idx * 3 + (0 if label == "legacy" else 1)
                ax.bar(pos, data.median(), color=color, alpha=0.7, hatch=hatch,
                       width=0.8, edgecolor="black", linewidth=0.5)
                ax.errorbar(pos, data.median(),
                            yerr=[[data.median() - data.quantile(0.25)],
                                  [data.quantile(0.75) - data.median()]],
                            fmt="none", color="black", capsize=3, linewidth=1)
        ax.set_xticks([0.5, 3.5, 6.5])
        ax.set_xticklabels(["VT", "V", "T_Vstatic"])
        ax.set_ylabel("median Pearson r (IQR err)")
        ax.set_title(f"{var} — by condition (left=legacy, right=current)")
        ax.set_ylim(-0.5, 1.0)
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.grid(True, alpha=0.3)

    # row 1 col 2: per-cluster delta histogram (current - legacy)
    ax = axes[1, 2]
    if "current" in pivot.columns and "legacy" in pivot.columns:
        deltas = (pivot["current"] - pivot["legacy"]).dropna()
        ax.hist(deltas, bins=30, color="purple", alpha=0.7)
        ax.axvline(0, color="black", linewidth=0.8)
        med = deltas.median()
        ax.axvline(med, color="red", linewidth=1.2, label=f"median Δ = {med:+.3f}")
        ax.set_xlabel("Δ Pearson r (current − legacy)")
        ax.set_ylabel("# (cluster, variable)")
        ax.set_title(f"per-cluster Δ (Selected) — n={len(deltas)}")
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Tuning-curve Pearson r — current (no onset, history-on) vs legacy (with onset)",
        fontsize=12,
    )
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.with_suffix(".png"), dpi=150)
    plt.close(fig)
    log.info("wrote %s", out_pdf)


def main() -> int:
    log.info("collecting current run")
    cur = collect_run(CURRENT, "current")
    log.info("collecting legacy run")
    leg = collect_run(LEGACY, "legacy")
    merged = pd.concat([cur, leg], ignore_index=True)
    out_csv = CURRENT / "diagnostics" / "tuning_pearson_by_cluster.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_csv, index=False)
    log.info("wrote %s (%d rows)", out_csv, len(merged))

    print()
    print("=== Per-cluster mean Pearson r (Selected model) ===")
    for var in ("Speed", "TF"):
        for label in ("legacy", "current"):
            sub = merged[(merged.run == label) & (merged.variable == var)
                         & (merged.model == "Selected")]
            data = sub.groupby(["probe_id", "cluster_id"])["pearson_r"].mean().dropna()
            print(f"  {var:5s}  {label:8s}: median = {data.median():+.3f}, mean = {data.mean():+.3f}, n = {len(data)}")
    out_pdf = CURRENT / "figs" / "tuning_pearson_comparison.pdf"
    out_pdf.parent.mkdir(exist_ok=True)
    render_plot(merged, out_pdf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
