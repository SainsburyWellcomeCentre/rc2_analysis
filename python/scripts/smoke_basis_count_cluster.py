"""Basis-count × ridge smoke test for one cluster / one condition.

Re-fits a single cluster's condition GLM over a grid of Speed/TF basis
counts × ridge λ, and overlays each fitted Additive-model kernel against
the observed median tuning (+Q1-Q3 band), mean-normalised so the SHAPE
match is readable. Used to see whether raising the basis count and/or
tuning regularisation lets the kernel capture a mid-range tuning peak the
5-basis λ=1 production blurs (and whether high-basis ringing is a
regularisation artefact).

Default target: CAA-1123243_rec1 cluster 376, VT,
n_bases ∈ {8,10,15,20} × λ ∈ {0.01,0.1,1,10,100}.
Output: /tmp/basis_ridge_smoke_<probe>_cl<cluster>_<cond>_<var>.pdf (one per var)
"""
from __future__ import annotations

import argparse
import tempfile
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.pipeline import _observed_median_tuning, run_pipeline
from rc2_glm.plots import _kernel_for_var
from rc2_glm.precomputed_bins import load_precomputed_bin_edges

import scripts.run_glm_split_by_condition as drv

PROBE = "CAA-1123243_rec1"
CLUSTER = 376
COND = "VT"
BASIS_COUNTS = (8, 10, 15, 20)
LAMBDAS = (0.01, 0.1, 1.0, 10.0, 100.0)
VARS = ("Speed", "TF")


def _mean_norm(y: np.ndarray) -> np.ndarray:
    m = np.nanmean(y)
    return y / m if (np.isfinite(m) and m != 0) else y


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spacing", choices=("log", "linear"), default="log",
                    help="Speed/TF basis spacing (default log/Weber).")
    spacing = ap.parse_args().spacing

    mat = drv.FORMATTED_DIR / f"{PROBE}.mat"
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    pre = load_precomputed_bin_edges(mat)
    observed = _observed_median_tuning(pre, COND, CLUSTER) or {}

    # Fit the cluster's VT GLM over the basis-count × ridge grid (single
    # seed — we read the Additive model, which is seed-independent).
    kernels: dict[tuple[int, float], dict[str, tuple]] = {}
    for n in BASIS_COUNTS:
        for lam in LAMBDAS:
            cfg = replace(
                drv.make_config(COND),
                n_speed_bases=n, n_tf_bases=n,
                lambda_ridge=lam, full_interaction_lambda=lam,
                speed_tf_basis_spacing=spacing,
                n_selection_seeds=1, selection_threshold_count=1,
            )
            with tempfile.TemporaryDirectory() as tmp:
                run_pipeline(
                    mat_path=mat, config=cfg, output_dir=Path(tmp),
                    stimulus_lookup=lookup, backend="irls", visp_only=True,
                    make_plots=False, n_jobs=1, cluster_filter={CLUSTER},
                )
                coef = pd.read_csv(Path(tmp) / "glm_coefficients.csv")
            add = coef[(coef["cluster_id"] == CLUSTER) & (coef["model"] == "Additive")]
            kernels[(n, lam)] = {
                var: _kernel_for_var(add, var, cfg) for var in VARS
            }
            print(f"n={n:2d} λ={lam:<6g}: fitted")

    # One grid figure per variable: rows = basis count, cols = ridge λ.
    # Solid = mean-normalised Additive model gain; dashed + band = observed.
    for var in VARS:
        obs = observed.get(var)
        fig, axes = plt.subplots(
            len(BASIS_COUNTS), len(LAMBDAS),
            figsize=(3.2 * len(LAMBDAS), 2.6 * len(BASIS_COUNTS)),
            squeeze=False, sharex=True,
        )
        for r, n in enumerate(BASIS_COUNTS):
            for c, lam in enumerate(LAMBDAS):
                ax = axes[r][c]
                kd = kernels[(n, lam)][var]
                if kd is not None:
                    x, kern, xlabel = kd
                    ax.plot(x, _mean_norm(np.exp(kern)), color="tab:red", lw=2,
                            label="model gain")
                    if r == len(BASIS_COUNTS) - 1:
                        ax.set_xlabel(xlabel, fontsize=8)
                if obs is not None:
                    ox, omed, oq1, oq3 = obs
                    div = np.nanmean(omed)
                    if np.isfinite(div) and div != 0:
                        band = np.isfinite(oq1) & np.isfinite(oq3)
                        ax.fill_between(ox, oq1 / div, oq3 / div, where=band,
                                        color="0.4", alpha=0.3, linewidth=0)
                        ax.plot(ox, omed / div, color="black", lw=1.4, ls="--",
                                label="obs median")
                ax.axhline(1.0, color="0.7", lw=0.8, ls=":")
                ax.tick_params(labelsize=7)
                if r == 0:
                    ax.set_title(f"λ = {lam:g}", fontsize=10, fontweight="bold")
                if c == 0:
                    ax.set_ylabel(f"{n} bases\n(÷ mean)", fontsize=9)
                if r == 0 and c == 0:
                    ax.legend(fontsize=6)
        fig.suptitle(
            f"{PROBE} cluster {CLUSTER} · {COND} · {var} kernel — "
            f"basis-count × ridge grid · {spacing} spacing "
            f"(solid: Additive model gain · dashed: observed median)",
            fontsize=11, fontweight="bold",
        )
        fig.tight_layout()
        out = Path(f"/tmp/basis_ridge_smoke_{PROBE}_cl{CLUSTER}_{COND}_{spacing}_{var}.pdf")
        fig.savefig(out)
        plt.close(fig)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
