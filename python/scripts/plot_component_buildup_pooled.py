"""Component build-up on the POOLED model (all conditions together).

Like plot_component_buildup_vt.py but fit on ALL conditions jointly (the
identifiable model), for the curated clusters. Nested refits
  Null(intercept+onset) → +Speed → +<other selected vars>
each refit on the full pooled data; shows the cumulative marginal predicted
tuning vs the observed, along the Speed and TF axes. Speed is always the first
stimulus step (the "does Speed add beyond onset?" probe).

Observed and predicted are marginalised the SAME way (per-bin rate binned by
speed/TF over the pooled motion bins, median per bin), so they're directly
comparable. Selected vars are read from the pooled run (current/).

Output: current/diagnostics/component_buildup_pooled/<probe>_cl<cluster>.pdf
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
from rc2_glm.pipeline import _observed_median_tuning
from rc2_glm.precomputed_bins import load_precomputed_bin_edges
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

ORDER = ("Speed", "TF", "SF", "OR")
CLUSTERS = [("CAA-1123243_rec1", c) for c in (148, 248, 362, 376)] + \
           [("CAA-1123244_rec1", c) for c in (34, 70, 80)]
OUTDIR = drv.POOLED_OUT / "diagnostics" / "component_buildup_pooled"
PREFIX = {"Speed": "Speed_", "TF": "TF_", "SF": "SF_", "OR": "OR_"}
DISPLAY_COND = "VT"   # display tuning in VT (has both Speed and TF, full range)


def _selected(probe: str, cluster: int) -> list[str]:
    df = pd.read_csv(drv.POOLED_OUT / "_runs" / probe / "glm_model_comparison.csv")
    row = df.loc[df["cluster_id"] == cluster, "time_selected_vars"]
    sv = row.iloc[0] if not row.empty else None
    return [] if (not isinstance(sv, str) or sv == "Null") else sv.split("+")


def _marg(rate, vals, edges):
    """Median rate per cache bin — bins predicted into the SAME 20 5%-quantile
    bins the cache uses for the observed tuning."""
    idx = np.digitize(vals, edges) - 1
    n = len(edges) - 1
    out = np.full(n, np.nan)
    for b in range(n):
        m = idx == b
        if m.any():
            out[b] = np.nanmedian(rate[m])
    return 0.5 * (edges[:-1] + edges[1:]), out


def buildup(pdata, pre, cfg, probe, cluster):
    cl = next(c for c in pdata.clusters if c.cluster_id == cluster)
    df = bin_cluster(pdata, cl)
    if df.empty:
        return
    selected = _selected(probe, cluster)
    build_vars = list(dict.fromkeys(["Speed"] + [v for v in ORDER if v in selected]))
    speed = df["speed"].to_numpy(float); tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float); orient = df["orientation"].to_numpy(float)
    y = df["spike_count"].to_numpy(float)
    # Fit on the FULL pooled data; DISPLAY the tuning in VT (its 20 cache bins).
    disp = (df["condition"] == DISPLAY_COND).to_numpy()
    offset = float(np.log(cfg.time_bin_width))
    sp = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=sp)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=sp)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])

    subsets = [("Null", [])] + [(f"+{v}", build_vars[:i + 1])
                                for i, v in enumerate(build_vars)]
    preds = []
    for label, vars_ in subsets:
        X, _ = assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf, orient, vars_, include_onset_kernel=True)
        fit = fit_poisson_glm(X, y, offset=offset, lambda_ridge=cfg.lambda_ridge,
                              backend="irls")
        preds.append((label, np.exp(X @ fit.beta)))  # Hz (offset cancels)

    observed = _observed_median_tuning(pre, DISPLAY_COND, cluster) or {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    cmap = plt.get_cmap("viridis")
    panels = [("Speed", speed, pre.speed_edges(DISPLAY_COND), "speed (cm/s)"),
              ("TF", tf, pre.tf_edges(DISPLAY_COND), "TF (Hz)")]
    for ax, (var, vals, edges, xlabel) in zip(axes, panels):
        if edges is None:
            ax.set_axis_off(); continue
        for k, (label, rate) in enumerate(preds):
            c, m = _marg(rate[disp], vals[disp], edges)
            final = k == len(preds) - 1
            ax.plot(c, m, color="tab:red" if final else cmap(k / max(len(preds) - 1, 1)),
                    lw=2.4 if final else 1.4, label=label + (" (full)" if final else ""))
        obs = observed.get(var)
        if obs is not None:
            ox, omed, oq1, oq3 = obs
            band = np.isfinite(oq1) & np.isfinite(oq3)
            ax.fill_between(ox, oq1, oq3, where=band, color="0.4", alpha=0.25, linewidth=0)
            ax.plot(ox, omed, color="black", lw=1.6, ls="--", label=f"observed ({DISPLAY_COND})")
        ax.set_xlabel(xlabel); ax.set_ylabel("firing rate (Hz)")
        ax.set_title(f"{var} tuning — pooled build-up")
        ax.legend(fontsize=7)
    fig.suptitle(f"{probe} cluster {cluster} · POOLED fit, shown in {DISPLAY_COND} "
                 f"(20 cache bins) · selected: {'+'.join(selected) or 'Null'}",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"{probe}_cl{cluster}.pdf"
    fig.savefig(out); fig.savefig(out.with_suffix(".png"), dpi=110)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    cfg = drv.make_config(None)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    by_probe: dict[str, list[int]] = {}
    for probe, cl in CLUSTERS:
        by_probe.setdefault(probe, []).append(cl)
    for probe, clusters in by_probe.items():
        mat = drv.FORMATTED_DIR / f"{probe}.mat"
        pre = load_precomputed_bin_edges(mat)
        pdata = load_probe_data(mat, config=cfg, stimulus_lookup=lookup, visp_only=True)
        for cl in clusters:
            buildup(pdata, pre, cfg, probe, cl)
    print(f"done — {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
