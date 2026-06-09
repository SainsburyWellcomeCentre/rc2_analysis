"""Collinearity diagnostics for the motion-clouds GLM design.

Quantifies how separable each predictor block (Speed / TF / Onset / SF / OR)
is in each condition, BEFORE any fitting — a feasibility check for "can we
estimate a unique Speed contribution here?".

Reports, per (probe, condition) and pooled:
  - generalised VIF (Fox & Monette) per block, on the SE scale
    GVIF^(1/2df): compare to sqrt(5)=2.24 (worrying) / sqrt(10)=3.16 (severe).
    Variance-scale equivalent = the square.
  - design condition number (standardised columns): >30 moderate, >100 severe.

The design (Speed/TF/SF/OR/Onset bases) depends only on the stimulus+behaviour,
not the cluster's spikes, so one representative cluster per probe suffices.
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
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import _subset_to_condition
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

BLOCK_PREFIX = {"Speed": "Speed_", "TF": "TF_", "Onset": "Onset_",
                "SF": "SF_", "OR": "OR_"}


def _gvif(R: np.ndarray, idx: list[int]) -> float:
    others = [i for i in range(R.shape[0]) if i not in idx]
    detR = np.linalg.det(R)
    if detR <= 1e-12:
        return np.inf
    det_b = np.linalg.det(R[np.ix_(idx, idx)])
    det_o = np.linalg.det(R[np.ix_(others, others)]) if others else 1.0
    return det_b * det_o / detR


def _design(df, cfg):
    speed = df["speed"].to_numpy(float); tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float); orient = df["orientation"].to_numpy(float)
    sp = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=sp)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=sp)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf, orient, ["Speed", "TF", "SF", "OR"],
        include_onset_kernel=True)
    keep = [i for i, nm in enumerate(names) if nm != "Intercept"]
    return X[:, keep], [names[i] for i in keep]


def main() -> int:
    cfg = drv.make_config(None)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    scopes = [("pooled", None), ("T_Vstatic", "T_Vstatic"), ("V", "V"), ("VT", "VT")]
    rows = []
    for probe in drv.PROBES:
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        cohort = drv.load_cohort(probe)
        cl = next(c for c in pdata.clusters if c.cluster_id in cohort)
        df_full = bin_cluster(pdata, cl)
        for scope, cond in scopes:
            df = df_full if cond is None else _subset_to_condition(df_full, cond)
            if df.empty:
                continue
            X, names = _design(df, cfg)
            Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
            R = np.corrcoef(Xs, rowvar=False)
            blocks = {b: [i for i, nm in enumerate(names) if nm.startswith(p)]
                      for b, p in BLOCK_PREFIX.items()}
            rec = {"probe": probe, "scope": scope,
                   "cond_number": float(np.linalg.cond(Xs))}
            for b, idx in blocks.items():
                if not idx:
                    continue
                g = _gvif(R, idx)
                rec[f"gvif_{b}"] = g ** (1.0 / (2 * len(idx))) if np.isfinite(g) else np.inf
            rows.append(rec)

    df = pd.DataFrame(rows)
    print("\n=== generalised VIF (SE scale, GVIF^(1/2df)); compare 2.24 / 3.16; "
          "cond_number >30 moderate / >100 severe ===\n")
    for scope, _ in scopes:
        sub = df[df["scope"] == scope]
        if sub.empty:
            continue
        def med(col):
            return sub[col].replace(np.inf, np.nan).median() if col in sub else np.nan
        parts = [f"Speed={med('gvif_Speed'):.2f}" if "gvif_Speed" in sub else "",
                 f"TF={med('gvif_TF'):.2f}" if "gvif_TF" in sub else "",
                 f"Onset={med('gvif_Onset'):.2f}" if "gvif_Onset" in sub else "",
                 f"SF={med('gvif_SF'):.2f}" if "gvif_SF" in sub else "",
                 f"OR={med('gvif_OR'):.2f}" if "gvif_OR" in sub else ""]
        parts = [p for p in parts if p]
        print(f"{scope:10s} cond#={med('cond_number'):7.1f}  | " + "  ".join(parts))

    _plot_vif(df)
    return 0


# Bands on the SE scale (GVIF^(1/2df)); 2.24/3.16 = variance-scale 5/10.
GOOD, WORRY = 2.24, 3.16
BAND_COLOR = {"good": "#2ca02c", "worrying": "#ff8c00", "severe": "#d62728",
              "rank-deficient": "#6a1b9a", "n/a": "0.85"}


def _band(v: float) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "n/a"
    if not np.isfinite(v):
        return "rank-deficient"
    if v < GOOD:
        return "good"
    if v < WORRY:
        return "worrying"
    return "severe"


def _plot_vif(df: pd.DataFrame) -> None:
    scopes = ["pooled", "T_Vstatic", "V", "VT"]
    cols = ["Speed", "TF", "Onset", "SF", "OR"]

    def med(scope, col):
        sub = df[df["scope"] == scope]
        if col not in sub or sub.empty:
            return np.nan
        return sub[col].median()  # inf medians stay inf

    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    for r, scope in enumerate(scopes):
        for c, block in enumerate(cols):
            v = med(scope, f"gvif_{block}")
            band = _band(v)
            label = ("n/a" if band == "n/a"
                     else "∞" if band == "rank-deficient" else f"{v:.2f}")
            ax.add_patch(plt.Rectangle((c, len(scopes) - 1 - r), 1, 1,
                                       facecolor=BAND_COLOR[band], edgecolor="white",
                                       lw=2))
            tc = "white" if band in ("severe", "rank-deficient") else "black"
            ax.text(c + 0.5, len(scopes) - 1 - r + 0.5, label, ha="center",
                    va="center", fontsize=9, color=tc, fontweight="bold")
    ax.set_xlim(0, len(cols)); ax.set_ylim(0, len(scopes))
    ax.set_xticks(np.arange(len(cols)) + 0.5); ax.set_xticklabels(cols)
    ax.set_yticks(np.arange(len(scopes)) + 0.5)
    ax.set_yticklabels(list(reversed(scopes)))
    ax.set_xticks(np.arange(len(cols) + 1), minor=True)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BAND_COLOR[b]) for b in
               ("good", "worrying", "severe", "rank-deficient", "n/a")]
    labels = ["good (<2.24)", "worrying (2.24–3.16)", "severe (>3.16)",
              "rank-deficient (∞)", "n/a (degenerate)"]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=3, fontsize=7.5, frameon=False)
    ax.set_title("Collinearity per condition — generalised VIF (SE scale, GVIF^(1/2df))\n"
                 "Speed/TF separable in pooled & T_Vstatic; NOT in VT (rank-deficient)",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    out = Path.home() / "local_data/motion_clouds/figures/glm/current/diagnostics/collinearity_vif"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")


if __name__ == "__main__":
    raise SystemExit(main())
