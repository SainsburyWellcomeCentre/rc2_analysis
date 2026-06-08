"""Park-2014-style tuning build-up by REFITTING nested model subsets (VT).

For each chosen VT cluster, fits a sequence of nested GLMs and shows each
model's marginal predicted tuning against the observed median, so you can see
how progressively richer models reconstruct the tuning curve:

    Null (intercept+onset) → +Speed → +<other selected vars, in order>

Speed is ALWAYS the first stimulus step (a deliberate "is this cell explained
by Speed alone?" probe — included even when forward selection dropped it),
followed by the cell's actually-selected variables. Each curve is a SEPARATE
refit (ridge λ=1) — not a decomposition of one full model — so the {Speed}-only
curve is the genuine speed-only prediction. The final curve is the full
selected (∪Speed) model and should track the observed median (dashed).

Marginal predicted tuning = per-bin predicted rate (Hz) of the refit model,
binned over the cluster's motion bins into the MATLAB cache's Speed/TF bins
(median per bin) — same construction as the observed median.

Output: <run>/figs/component_buildup_vt/<probe>_cl<cluster>_VT.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import onset_kernel_basis, value_basis
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.fitting import fit_poisson_glm
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import _observed_median_tuning, _subset_to_condition
from rc2_glm.precomputed_bins import load_precomputed_bin_edges
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

COND = "VT"
ORDER = ("Speed", "TF", "SF", "OR")
CLUSTERS = [
    ("CAA-1123243_rec1", c) for c in (148, 248, 362, 376)
] + [
    ("CAA-1123244_rec1", c) for c in (34, 70, 80)
]
OUTDIR = drv.OUT_ROOT / "figs" / "component_buildup_vt"
PREFIX = {"Speed": "Speed_", "TF": "TF_", "SF": "SF_", "OR": "OR_"}


def _selected(probe: str, cluster: int) -> list[str]:
    df = pd.read_csv(drv.OUT_ROOT / "_runs" / probe / COND / "glm_model_comparison.csv")
    row = df.loc[df["cluster_id"] == cluster, "time_selected_vars"]
    sv = row.iloc[0] if not row.empty else None
    return [] if (not isinstance(sv, str) or sv == "Null") else sv.split("+")


def _marginalise(rate, var_vals, edges, n_centres):
    idx = np.digitize(var_vals, edges) - 1
    out = np.full(n_centres, np.nan)
    for b in range(n_centres):
        m = idx == b
        if m.any():
            out[b] = np.nanmedian(rate[m])
    return out


def buildup(probe_data, pre, cfg, probe: str, cluster: int) -> None:
    observed = _observed_median_tuning(pre, COND, cluster) or {}
    selected = _selected(probe, cluster)
    # Speed always first (explained-by-speed probe), then selected vars.
    build_vars = list(dict.fromkeys(["Speed"] + [v for v in ORDER if v in selected]))

    cl = next(c for c in probe_data.clusters if c.cluster_id == cluster)
    df = _subset_to_condition(bin_cluster(probe_data, cl), COND)
    speed = df["speed"].to_numpy(float)
    tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float)
    orient = df["orientation"].to_numpy(float)
    y = df["spike_count"].to_numpy(float)
    motion = (df["condition"] != "stationary").to_numpy()
    offset = float(np.log(cfg.time_bin_width))
    sp = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=sp)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=sp)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])

    # Refit each nested subset; store predicted rate (Hz) per bin.
    subsets = [("Null", [])]
    for i, v in enumerate(build_vars):
        subsets.append((f"+{v}", build_vars[: i + 1]))
    preds = []  # (label, rate_per_bin, is_selected_marker)
    for label, vars_ in subsets:
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, vars_, include_onset_kernel=True,
        )
        fit = fit_poisson_glm(X, y, offset=offset, lambda_ridge=cfg.lambda_ridge,
                              backend="irls")
        preds.append((label, np.exp(X @ fit.beta)))   # offset cancels → Hz

    panels = [("Speed", speed, "speed (cm/s)", pre.speed_edges, pre.speed_centres),
              ("TF", tf, "TF (Hz)", pre.tf_edges, pre.tf_centres)]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    cmap = plt.get_cmap("viridis")
    for ax, (var, vals, xlabel, edges_fn, centres_fn) in zip(axes, panels):
        centres, edges = centres_fn(COND), edges_fn(COND)
        if centres is None or edges is None:
            ax.set_axis_off(); continue
        n_c = len(centres)
        for k, (label, rate) in enumerate(preds):
            final = k == len(preds) - 1
            ax.plot(centres, _marginalise(rate[motion], vals[motion], edges, n_c),
                    color="tab:red" if final else cmap(k / max(len(preds) - 1, 1)),
                    lw=2.4 if final else 1.4,
                    label=label + (" (full)" if final else ""))
        obs = observed.get(var)
        if obs is not None:
            ox, omed, oq1, oq3 = obs
            band = np.isfinite(oq1) & np.isfinite(oq3)
            ax.fill_between(ox, oq1, oq3, where=band, color="0.4", alpha=0.25, linewidth=0)
            ax.plot(ox, omed, color="black", lw=1.6, ls="--", label="observed median")
        ax.set_xlabel(xlabel); ax.set_ylabel("firing rate (Hz)")
        ax.set_title(f"{var} tuning — build-up")
        ax.legend(fontsize=7)
    fig.suptitle(
        f"{probe} cluster {cluster} · VT · nested-refit build-up "
        f"(forward-selected: {'+'.join(selected) or 'Null'})",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"{probe}_cl{cluster}_VT.pdf"
    fig.savefig(out); plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    cfg = drv.make_config(COND)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    by_probe: dict[str, list[int]] = {}
    for probe, cl in CLUSTERS:
        by_probe.setdefault(probe, []).append(cl)
    for probe, clusters in by_probe.items():
        mat = drv.FORMATTED_DIR / f"{probe}.mat"
        pre = load_precomputed_bin_edges(mat)
        probe_data = load_probe_data(mat, config=cfg, stimulus_lookup=lookup,
                                     visp_only=True)
        for cl in clusters:
            buildup(probe_data, pre, cfg, probe, cl)
    print(f"done — {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
