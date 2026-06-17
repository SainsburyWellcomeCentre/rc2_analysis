"""RF-local diagnostics → ``<run>/diagnostics/`` — REUSING the canonical figures.

This does NOT define any plots of its own. It computes the per-cluster variance
partition on the rf_local continuous SF/OR design and feeds it to the existing
figure code:

* ``variance_partition.{csv,pdf,png}`` — ``variance_partition_accel.partition``
  (acid test: unique cv-bps per predictor + the Speed-survival sequence) rendered
  by ``variance_partition_accel.plot_partition_figure``.
* ``acid_vs_selection.{pdf,png}`` — ``plot_acid_vs_selection.plot_acid_figure``
  with ``forward_only=True`` (per-cluster stack of ONLY the forward-selected
  predictors' contributions).

``generate(run)`` is the entry point the run driver calls at the end of a run, so
diagnostics are part of the pipeline and never skipped. Runs: main / histme / all.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from rc2_glm.basis import circular_basis, raised_cosine_basis_linear
from rc2_glm.io import load_probe_data
from rc2_glm.prefilter import passes_spike_floor
from rc2_glm.rf_sf_or import load_rf_sf_or
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_goggles_rf_sfor_20ms as drv
from scripts.plot_acid_vs_selection import (plot_acid_figure, plot_acid_stack_sorted,
                                            plot_spike_count_diagnostics)
from scripts.variance_partition_accel import partition, plot_partition_figure

RUNS = {
    "main":   (drv.make_config,            drv.OUT_ROOT),
    "histme": (drv.make_config_histme,     drv.OUT_HISTME),
    "all":    (drv.make_config_histme_all, drv.OUT_HISTME_ALL),
    "all_V":  (lambda: drv.make_config_histme_all_cond("V"), drv.OUT_ALL_BY_COND["V"]),
    "all_Tvstatic": (lambda: drv.make_config_histme_all_cond("T_Vstatic"),
                     drv.OUT_ALL_BY_COND["T_Vstatic"]),
    "histbase": (drv.make_config_histbase_all, drv.OUT_HISTBASE_ALL),
}
_STIM = ("Speed", "TF", "SF", "OR")


def _rf_bases(df, cfg):
    """Continuous RF-local SF/OR bases (same construction as the production fit)."""
    sf = df["sf"].to_numpy(float)
    orient = df["orientation"].to_numpy(float)
    n = sf.size
    sf_fin = np.isfinite(sf)
    B_sf = np.zeros((n, cfg.n_sf_bases))
    if sf_fin.any():
        B_sf[sf_fin] = raised_cosine_basis_linear(sf[sf_fin], cfg.n_sf_bases, *cfg.sf_cpd_range)
    or_fin = np.isfinite(orient)
    B_or = np.zeros((n, cfg.n_or_bases))
    if or_fin.any():
        B_or[or_fin] = circular_basis(orient[or_fin], cfg.n_or_bases)
    return B_sf, B_or


def compute_partition(run: str, limit: int | None = None) -> tuple[pd.DataFrame, object, object]:
    cfg_fn, folder = RUNS[run]
    cfg = cfg_fn()
    # The condition's non-degenerate stimulus set (V has no Speed; T_Vstatic no
    # TF/SF/OR) — so the partition drops degenerate blocks instead of fitting
    # constant columns. Derived from the run's own main_effects.
    stim = tuple(v for v in cfg.main_effects if v in _STIM)
    rows = []
    spike_rows = []
    for probe in drv.PROBES:
        mc = folder / "_runs" / probe / "glm_model_comparison.csv"
        if not mc.exists():
            continue
        cohort = list(pd.read_csv(mc)["cluster_id"])
        if limit is not None:
            cohort = cohort[:limit]
        cohort = set(cohort)
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=drv._lookup(), cluster_set="selected")
        rf = load_rf_sf_or(drv.RF_PARQUET_DIR, probe,
                           min_concentration=getattr(cfg, "rf_min_concentration", 0.0))
        for cl in pdata.clusters:
            if cl.cluster_id not in cohort:
                continue
            df = bin_cluster(pdata, cl, rf_lookup=rf)
            # Record spike stats for EVERY cohort cluster (incl. floor-excluded)
            # so the spike-count diagnostic shows the filter, not just survivors.
            per_trial = df.groupby("trial_id")["spike_count"].sum().to_numpy(float)
            q = (np.percentile(per_trial, [0, 25, 50, 75, 100])
                 if per_trial.size else np.zeros(5))
            passed = passes_spike_floor(df, cfg.min_spikes_floor, cfg.min_trial_occupancy)
            spike_rows.append(dict(
                probe_id=probe, cluster_id=int(cl.cluster_id),
                n_spikes=int(per_trial.sum()), n_trials=int(per_trial.size),
                trial_occupancy=float((per_trial > 0).mean()) if per_trial.size else 0.0,
                passed_floor=bool(passed),
                pt_min=q[0], pt_q1=q[1], pt_med=q[2], pt_q3=q[3], pt_max=q[4]))
            if not passed:
                continue
            B_sf, B_or = _rf_bases(df, cfg)
            # Pass the run's ACTUAL predictor set so the partition matches what
            # was fit (the main run has no ME/History; histme/all have both).
            res = partition(
                df, cfg, B_sf=B_sf, B_or=B_or, stim=stim,
                include_me=getattr(cfg, "include_me_face", False),
                include_hist=getattr(cfg, "include_history", False),
                include_accel=getattr(cfg, "include_acceleration", False),
            )
            if res is None:
                continue
            res.update(probe_id=probe, cluster_id=int(cl.cluster_id))
            rows.append(res)
        print(f"{run}/{probe}: {sum(r['probe_id'] == probe for r in rows)} clusters partitioned")
    return pd.DataFrame(rows), pd.DataFrame(spike_rows), cfg, folder


def _render(df: pd.DataFrame, spike_df, cfg, folder, run: str) -> None:
    """Render every diagnostic figure into ``<run>/diagnostics/`` — shared by
    ``generate`` (after computing) and ``--figures-only`` (after reading the
    saved CSVs), so the figures never drift."""
    out_dir = folder / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    # Cohort spike-count quality diagnostics (the floor gate) — independent of
    # the partition, so render even when no cluster passed to be partitioned.
    if spike_df is not None and len(spike_df):
        plot_spike_count_diagnostics(
            spike_df, out_dir / "spike_count_diagnostics",
            min_spikes=cfg.min_spikes_floor, min_trial_frac=cfg.min_trial_occupancy)
    if df is None or df.empty:
        return
    well = df["total_full"] > 0.005
    bin_ms = int(round(cfg.time_bin_width * 1000))
    plot_partition_figure(
        df, well, out_dir / "variance_partition", bin_ms=bin_ms,
        suptitle=f"Variance partition (acid test) — rf_local {run} "
                 "(Speed vs Accel/ME/History attribution)")

    sel = pd.read_csv(folder / "glm_model_comparison.csv")[
        ["probe_id", "cluster_id", "time_selected_vars"]]
    m = df.merge(sel, on=["probe_id", "cluster_id"], how="inner")
    m = m[m["total_full"] > 0.005].copy()
    m["sel_set"] = m["time_selected_vars"].fillna("Null").apply(lambda s: set(s.split("+")))
    m["r2"] = 1 - m["cv_full"] / m["cv_intercept"]
    m = m.sort_values("r2").reset_index(drop=True)
    plot_acid_figure(m, out_dir / "acid_vs_selection", forward_only=True)         # forward-selected only
    plot_acid_stack_sorted(m, out_dir / "acid_stack_sorted")                       # all regressors, sorted by stack ± History


def generate(run: str, limit: int | None = None) -> int:
    """Compute the rf_local variance partition for ``run``, save the CSV, and
    render every diagnostic figure into ``<run>/diagnostics/``. Pipeline entry."""
    df, spike_df, cfg, folder = compute_partition(run, limit=limit)
    (folder / "diagnostics").mkdir(parents=True, exist_ok=True)
    if spike_df is not None and len(spike_df):
        spike_df.to_csv(folder / "diagnostics" / "spike_stats.csv", index=False)
    if not df.empty:
        df.to_csv(folder / "diagnostics" / "variance_partition.csv", index=False)
    else:
        print(f"[diagnostics] {run}: no clusters partitioned (needs ME+History+Accel "
              "in the config) — rendering spike-count diagnostics only")
    _render(df, spike_df, cfg, folder, run)
    return 0 if not df.empty else 1


def figures_only(run: str) -> int:
    """Re-render the diagnostic figures from the SAVED variance_partition.csv —
    no partition recompute. For adding a new figure to runs already computed."""
    cfg_fn, folder = RUNS[run]
    csv = folder / "diagnostics" / "variance_partition.csv"
    if not csv.exists():
        print(f"[diagnostics] {run}: {csv} not found — run the full partition first")
        return 1
    scsv = folder / "diagnostics" / "spike_stats.csv"
    spike_df = pd.read_csv(scsv) if scsv.exists() else None
    _render(pd.read_csv(csv), spike_df, cfg_fn(), folder, run)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", choices=list(RUNS), required=True)
    ap.add_argument("--limit", type=int, default=None, help="first N clusters/probe (smoke)")
    ap.add_argument("--figures-only", action="store_true",
                    help="re-render figures from the saved variance_partition.csv (no recompute)")
    args = ap.parse_args()
    if args.figures_only:
        return figures_only(args.run)
    return generate(args.run, limit=args.limit)


if __name__ == "__main__":
    raise SystemExit(main())
