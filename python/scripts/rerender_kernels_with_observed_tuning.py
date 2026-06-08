"""Re-render per-cluster kernel figures with the empirical tuning overlay.

For the split-by-condition run, regenerate each
``_runs/<probe>/<cond>/figs/cluster_<id>_kernels.pdf`` with the OBSERVED
median Speed/TF tuning curve (from the MATLAB tuning cache) overlaid as a
black dashed line on a twin Hz axis — so the fitted log-gain kernel and
the empirical tuning can be read together.

This rebuilds the figures from saved data only:
  * kernels      ← glm_coefficients.csv (plot_cluster_kernels reads betas)
  * dashed tuning← MATLAB cache (csvs/{tuning_curves,tf_tuning_curves}/<probe>.mat)
                   median across trials of per-trial-per-bin firing rate.
NO GLM re-fit. The dashed curve is the data, not a model prediction.

Usage::  python -m scripts.rerender_kernels_with_observed_tuning
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from rc2_glm.config import GLMConfig
from rc2_glm.pipeline import _observed_median_tuning
from rc2_glm.plots import plot_cluster_kernels, save_figure
from rc2_glm.precomputed_bins import load_precomputed_bin_edges

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("rerender_kernels")

HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds"
RUN = ROOT / "figures" / "glm" / "current_splitted_by_condition"
FORMATTED_DIR = ROOT / "formatted_data"
PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1", "CAA-1123467_rec1")
CONDITIONS = ("V", "T_Vstatic", "VT")
# Reuse the production split-by-condition knobs so the kernel grids
# (basis counts, ranges) match the fits that produced glm_coefficients.csv.
CONFIG = GLMConfig(include_me_face=False, include_history=False)


def main() -> int:
    n_fig = 0
    for probe in PROBES:
        pre = load_precomputed_bin_edges(FORMATTED_DIR / f"{probe}.mat")
        for cond in CONDITIONS:
            run_dir = RUN / "_runs" / probe / cond
            coef_csv = run_dir / "glm_coefficients.csv"
            if not coef_csv.exists():
                log.warning("missing %s — skipping", coef_csv)
                continue
            figs_dir = run_dir / "figs"
            coef = pd.read_csv(coef_csv)
            for cid, cdf in coef.groupby("cluster_id"):
                cid = int(cid)
                observed = _observed_median_tuning(pre, cond, cid)
                fig = plot_cluster_kernels(
                    probe_id=probe, cluster_id=cid,
                    coef_df_cluster=cdf, config=CONFIG,
                    observed_tuning=observed,
                )
                for p in save_figure(fig, figs_dir / f"cluster_{cid}_kernels", fmt="pdf"):
                    n_fig += 1
            log.info("re-rendered %s / %s", probe, cond)
    log.info("done — %d kernel figures rewritten under %s/_runs", n_fig, RUN)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
