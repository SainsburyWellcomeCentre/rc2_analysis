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
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

ORDER = ("Speed", "TF", "SF", "OR")
CLUSTERS = [("CAA-1123243_rec1", c) for c in (148, 248, 362, 376)] + \
           [("CAA-1123244_rec1", c) for c in (34, 70, 80)]
OUTDIR = drv.POOLED_OUT / "diagnostics" / "component_buildup_pooled"
PREFIX = {"Speed": "Speed_", "TF": "TF_", "SF": "SF_", "OR": "OR_"}
NBINS = 18


def _selected(probe: str, cluster: int) -> list[str]:
    df = pd.read_csv(drv.POOLED_OUT / "_runs" / probe / "glm_model_comparison.csv")
    row = df.loc[df["cluster_id"] == cluster, "time_selected_vars"]
    sv = row.iloc[0] if not row.empty else None
    return [] if (not isinstance(sv, str) or sv == "Null") else sv.split("+")


def _marg(rate, vals, edges):
    idx = np.digitize(vals, edges) - 1
    n = len(edges) - 1
    out = np.full(n, np.nan)
    for b in range(n):
        m = idx == b
        if m.any():
            out[b] = np.nanmedian(rate[m])
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, out


def buildup(pdata, cfg, probe, cluster):
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
    motion = (df["condition"] != "stationary").to_numpy()
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

    s_edges = np.linspace(*cfg.speed_range, NBINS + 1)
    t_edges = np.linspace(*cfg.tf_range, NBINS + 1)
    obs = y[motion] / cfg.time_bin_width
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    cmap = plt.get_cmap("viridis")
    for ax, (var, vals, edges, xlabel) in zip(
            axes, [("Speed", speed, s_edges, "speed (cm/s)"),
                   ("TF", tf, t_edges, "TF (Hz)")]):
        for k, (label, rate) in enumerate(preds):
            c, m = _marg(rate[motion], vals[motion], edges)
            final = k == len(preds) - 1
            ax.plot(c, m, color="tab:red" if final else cmap(k / max(len(preds) - 1, 1)),
                    lw=2.4 if final else 1.4, label=label + (" (full)" if final else ""))
        c, mo = _marg(obs, vals[motion], edges)
        ax.plot(c, mo, color="black", lw=1.6, ls="--", label="observed (pooled)")
        ax.set_xlabel(xlabel); ax.set_ylabel("firing rate (Hz)")
        ax.set_title(f"{var} tuning — pooled build-up")
        ax.legend(fontsize=7)
    fig.suptitle(f"{probe} cluster {cluster} · POOLED · nested-refit build-up "
                 f"(selected: {'+'.join(selected) or 'Null'})",
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
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        for cl in clusters:
            buildup(pdata, cfg, probe, cl)
    print(f"done — {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
