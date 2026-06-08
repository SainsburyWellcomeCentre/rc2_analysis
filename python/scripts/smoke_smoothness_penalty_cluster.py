"""Hardcastle-style smoothness-penalty smoke test for one cluster / condition.

Instead of ridge (L2 on coefficient MAGNITUDE), penalise the CURVATURE of the
reconstructed Speed/TF kernel:

    penalty(beta) = lambda_s * || D2 @ (B_grid @ beta) ||^2
                  = beta.T (lambda_s * B_grid.T D2.T D2 B_grid) beta

(second-difference roughness — the Park & Pillow 2011 / Hardcastle 2017
smoothness prior, already used for the History block in rc2_glm.penalty).
This shrinks wiggliness, not amplitude, so a tall LOCALISED peak survives if
the data support it — the thing ridge could only flatten.

Fits the Additive model for one cluster/condition over a basis-count ×
penalty grid (first column = plain ridge λ=1 for contrast, rest = smoothness
λ_s) and overlays the reconstructed Speed kernel on the observed median
tuning (+Q1-Q3). TF is penalised the same way but Speed is plotted.

Target: CAA-1123243_rec1 cluster 376 VT.
Output: /tmp/smooth_penalty_smoke_<probe>_cl<cluster>_<cond>_<spacing>.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import onset_kernel_basis, value_basis
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.fitting import fit_poisson_glm
from rc2_glm.io import load_probe_data
from rc2_glm.penalty import second_difference_matrix
from rc2_glm.pipeline import _observed_median_tuning, _subset_to_condition
from rc2_glm.precomputed_bins import load_precomputed_bin_edges
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

PROBE = "CAA-1123243_rec1"
CLUSTER = 376
COND = "VT"
BASIS_COUNTS = (5,)
SMOOTH_LAMBDAS = (1.0, 10.0, 100.0, 1000.0)   # curvature-penalty strengths
GRID_N = 200
FLOOR = 1e-3                                   # tiny ridge floor (keep PD)


def _mean_norm(y: np.ndarray) -> np.ndarray:
    m = np.nanmean(y)
    return y / m if (np.isfinite(m) and m != 0) else y


def _curvature_block(B_grid: np.ndarray) -> np.ndarray:
    """Normalised curvature matrix S = Bᵀ D2ᵀ D2 B (mean-diagonal scaled to 1)."""
    D2 = second_difference_matrix(B_grid.shape[0])
    S = B_grid.T @ D2.T @ D2 @ B_grid
    scale = np.trace(S) / S.shape[0]
    return S / scale if scale > 0 else S


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spacing", choices=("log", "linear"), default="log")
    spacing = ap.parse_args().spacing

    cfg = drv.make_config(COND)
    mat = drv.FORMATTED_DIR / f"{PROBE}.mat"
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    pre = load_precomputed_bin_edges(mat)
    observed = (_observed_median_tuning(pre, COND, CLUSTER) or {}).get("Speed")

    # Bin cluster 376, restrict to VT trials (+ their stationary prelude).
    probe = load_probe_data(mat, config=cfg, stimulus_lookup=lookup, visp_only=True)
    cluster = next(c for c in probe.clusters if c.cluster_id == CLUSTER)
    df = _subset_to_condition(bin_cluster(probe, cluster), COND)
    speed = df["speed"].to_numpy(float)
    tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float)
    orient = df["orientation"].to_numpy(float)
    y = df["spike_count"].to_numpy(float)
    offset = float(np.log(cfg.time_bin_width))

    s_grid = np.linspace(*cfg.speed_range, GRID_N)
    t_grid = np.linspace(*cfg.tf_range, GRID_N)

    # columns: ridge λ=1 (reference) + smoothness λ_s sweep
    col_specs = [("ridge", 1.0)] + [("smooth", ls) for ls in SMOOTH_LAMBDAS]
    kernels: dict[tuple[int, int], np.ndarray] = {}  # (n, col) -> gain on s_grid

    for n in BASIS_COUNTS:
        B_speed = value_basis(speed, n, *cfg.speed_range, spacing=spacing)
        B_tf = value_basis(tf, n, *cfg.tf_range, spacing=spacing)
        B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
        X, names = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient,
            ["Speed", "TF", "SF", "OR"], include_onset_kernel=True,
        )
        # Some bases can be dropped as zero-variance over this cluster's
        # support, so map kept columns back to their basis index.
        sp = [(i, int(nm.split("_")[1]) - 1) for i, nm in enumerate(names)
              if nm.startswith("Speed_")]
        tf_ = [(i, int(nm.split("_")[1]) - 1) for i, nm in enumerate(names)
               if nm.startswith("TF_")]
        sp_idx, sp_basis = [i for i, _ in sp], [b for _, b in sp]
        tf_idx, tf_basis = [i for i, _ in tf_], [b for _, b in tf_]
        # Reconstruction grid restricted to the kept bases.
        Bsg = value_basis(s_grid, n, *cfg.speed_range, spacing=spacing)[:, sp_basis]
        Btg = value_basis(t_grid, n, *cfg.tf_range, spacing=spacing)[:, tf_basis]
        Ss, St = _curvature_block(Bsg), _curvature_block(Btg)

        for c, (kind, lam) in enumerate(col_specs):
            if kind == "ridge":
                fit = fit_poisson_glm(X, y, offset=offset, lambda_ridge=lam,
                                      backend="irls")
            else:
                P = FLOOR * np.eye(len(names))
                P[0, 0] = 0.0  # intercept unpenalised
                P[np.ix_(sp_idx, sp_idx)] = lam * Ss + FLOOR * np.eye(len(sp_idx))
                P[np.ix_(tf_idx, tf_idx)] = lam * St + FLOOR * np.eye(len(tf_idx))
                fit = fit_poisson_glm(X, y, offset=offset, backend="irls",
                                      penalty_matrix=P)
            beta_speed = fit.beta[sp_idx]
            kernels[(n, c)] = np.exp(Bsg @ beta_speed)
            print(f"n={n:2d} {kind} λ={lam:<7g}: fitted")

    # Plot grid: rows = basis count, cols = penalty spec.
    fig, axes = plt.subplots(len(BASIS_COUNTS), len(col_specs),
                             figsize=(3.1 * len(col_specs), 2.7 * len(BASIS_COUNTS)),
                             squeeze=False, sharex=True)
    for r, n in enumerate(BASIS_COUNTS):
        for c, (kind, lam) in enumerate(col_specs):
            ax = axes[r][c]
            ax.plot(s_grid, _mean_norm(kernels[(n, c)]), color="tab:red", lw=2,
                    label="model gain")
            if observed is not None:
                ox, omed, oq1, oq3 = observed
                div = np.nanmean(omed)
                if np.isfinite(div) and div != 0:
                    band = np.isfinite(oq1) & np.isfinite(oq3)
                    ax.fill_between(ox, oq1 / div, oq3 / div, where=band,
                                    color="0.4", alpha=0.3, linewidth=0)
                    ax.plot(ox, omed / div, color="black", lw=1.4, ls="--",
                            label="obs median")
            ax.axhline(1.0, color="0.7", lw=0.8, ls=":")
            ax.tick_params(labelsize=7)
            if r == len(BASIS_COUNTS) - 1:
                ax.set_xlabel("speed (cm/s)", fontsize=8)
            title = "ridge λ=1" if kind == "ridge" else f"smooth λ={lam:g}"
            if r == 0:
                ax.set_title(title, fontsize=10, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{n} bases\n(÷ mean)", fontsize=9)
            if r == 0 and c == 0:
                ax.legend(fontsize=6)
    fig.suptitle(
        f"{PROBE} cluster {CLUSTER} · {COND} · Speed kernel — ridge vs "
        f"curvature-smoothness penalty · {spacing} spacing "
        f"(solid: Additive model gain · dashed: observed median)",
        fontsize=10, fontweight="bold",
    )
    fig.tight_layout()
    out = Path(f"/tmp/smooth_penalty_smoke_{PROBE}_cl{CLUSTER}_{COND}_{spacing}.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
