"""Collinearity diagnostics for the GOGGLES motion-clouds GLM design, WITH ME_face.

Goggles port of collinearity_vif.py. Adds the face motion-energy block
(ME_face, camera0 = face video, de-interleaved 30 Hz) to the design and asks:
in the POOLED fit, how collinear is ME_face with Speed (and the rest)?

We report only the pooled scope as the headline (the per-condition split is
rank-deficient in VT — TF = gain·Speed — so its VIF is uninformative; kept in
the heatmap for parity but read pooled). Same generalised-VIF (Fox & Monette,
SE scale GVIF^(1/2df); 2.24 worrying / 3.16 severe) + design condition number.

The design depends only on stimulus+behaviour+camera, not the cluster's spikes,
so one representative cohort cluster per probe suffices.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis_linear, value_basis
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import _subset_to_condition
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_goggles as drv

BLOCK_PREFIX = {"Speed": "Speed_", "TF": "TF_", "Onset": "Onset_",
                "SF": "SF_", "OR": "OR_", "ME_face": "ME_face_"}


def _gvif(R: np.ndarray, idx: list[int]) -> float:
    others = [i for i in range(R.shape[0]) if i not in idx]
    detR = np.linalg.det(R)
    if detR <= 1e-12:
        return np.inf
    det_b = np.linalg.det(R[np.ix_(idx, idx)])
    det_o = np.linalg.det(R[np.ix_(others, others)]) if others else 1.0
    return det_b * det_o / detR


def _me_basis(df, cfg):
    """Replicate the pipeline's ME_face basis: z-score over motion rows, raised cosines."""
    if "me_face_raw" not in df.columns:
        return None
    me_raw = df["me_face_raw"].to_numpy(dtype=np.float64)
    motion = (df["condition"] != "stationary").to_numpy()
    fin = np.isfinite(me_raw) & motion
    if int(fin.sum()) < 10:
        return None
    mu = float(me_raw[fin].mean())
    sd = float(me_raw[fin].std(ddof=0)) or 1.0
    me_z = np.where(np.isfinite(me_raw), (me_raw - mu) / sd, 0.0)
    return raised_cosine_basis_linear(me_z, cfg.n_me_face_bases, *cfg.me_face_range)


def _design(df, cfg):
    speed = df["speed"].to_numpy(float); tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float); orient = df["orientation"].to_numpy(float)
    sp = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=sp)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=sp)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    B_me = _me_basis(df, cfg)
    selected = ["Speed", "TF", "SF", "OR"] + (["ME_face"] if B_me is not None else [])
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf, orient, selected,
        include_onset_kernel=True, B_me_face=B_me)
    keep = [i for i, nm in enumerate(names) if nm != "Intercept"]
    return X[:, keep], [names[i] for i in keep]


def _block_collinearity(Xs, names, a_prefix, b_prefix):
    """Max |pairwise corr| and first canonical correlation between two blocks."""
    ia = [i for i, nm in enumerate(names) if nm.startswith(a_prefix)]
    ib = [i for i, nm in enumerate(names) if nm.startswith(b_prefix)]
    if not ia or not ib:
        return np.nan, np.nan
    R = np.corrcoef(Xs, rowvar=False)
    max_pair = float(np.nanmax(np.abs(R[np.ix_(ia, ib)])))
    # first canonical correlation via SVD of cross-correlation of orthonormal bases
    Qa, _ = np.linalg.qr(Xs[:, ia]); Qb, _ = np.linalg.qr(Xs[:, ib])
    s = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
    return max_pair, float(s[0])


def main() -> int:
    cfg = drv.make_config_me()
    lookup = drv._lookup()
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
            mp, cc = _block_collinearity(Xs, names, "ME_face_", "Speed_")
            rec["me_speed_maxcorr"] = mp
            rec["me_speed_canoncorr"] = cc
            rows.append(rec)

    df = pd.DataFrame(rows)
    print("\n=== generalised VIF (SE scale GVIF^(1/2df); 2.24 worry / 3.16 severe); "
          "cond# >30 mod / >100 severe ===\n")
    for scope, _ in scopes:
        sub = df[df["scope"] == scope]
        if sub.empty:
            continue
        med = lambda c: (sub[c].replace(np.inf, np.nan).median() if c in sub else np.nan)
        parts = [f"{b}={med('gvif_'+b):.2f}" for b in BLOCK_PREFIX if f"gvif_{b}" in sub
                 and np.isfinite(med('gvif_'+b))]
        print(f"{scope:10s} cond#={med('cond_number'):8.1f} | " + "  ".join(parts))

    pooled = df[df["scope"] == "pooled"]
    print("\n=== POOLED ME_face ↔ Speed collinearity (the question) ===")
    print(f"  max |pairwise corr| ME_face×Speed : {pooled['me_speed_maxcorr'].median():.3f}")
    print(f"  first canonical corr (block-block): {pooled['me_speed_canoncorr'].median():.3f}")
    print("  (compare: ~0 independent, →1 collinear; GVIF Speed/ME_face above are the omnibus view)")

    _plot_vif(df)
    return 0


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
    cols = ["Speed", "TF", "Onset", "SF", "OR", "ME_face"]

    def med(scope, col):
        sub = df[df["scope"] == scope]
        if col not in sub or sub.empty:
            return np.nan
        return sub[col].median()

    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    for r, scope in enumerate(scopes):
        for c, block in enumerate(cols):
            v = med(scope, f"gvif_{block}")
            band = _band(v)
            label = ("n/a" if band == "n/a"
                     else "∞" if band == "rank-deficient" else f"{v:.2f}")
            ax.add_patch(plt.Rectangle((c, len(scopes) - 1 - r), 1, 1,
                                       facecolor=BAND_COLOR[band], edgecolor="white", lw=2))
            tc = "white" if band in ("severe", "rank-deficient") else "black"
            ax.text(c + 0.5, len(scopes) - 1 - r + 0.5, label, ha="center",
                    va="center", fontsize=9, color=tc, fontweight="bold")
    ax.set_xlim(0, len(cols)); ax.set_ylim(0, len(scopes))
    ax.set_xticks(np.arange(len(cols)) + 0.5); ax.set_xticklabels(cols)
    ax.set_yticks(np.arange(len(scopes)) + 0.5); ax.set_yticklabels(list(reversed(scopes)))
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BAND_COLOR[b]) for b in
               ("good", "worrying", "severe", "rank-deficient", "n/a")]
    labels = ["good (<2.24)", "worrying (2.24–3.16)", "severe (>3.16)",
              "rank-deficient (∞)", "n/a"]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=3, fontsize=7.5, frameon=False)
    ax.set_title("Goggles collinearity + ME_face — generalised VIF (SE scale)\n"
                 "read the POOLED row (split VT is rank-deficient: TF=gain·Speed)",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    out = (Path.home() /
           "local_data/motion_clouds/figures/glm/current_plus_ME_goggles/diagnostics/collinearity_vif")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out}.pdf / .png")


if __name__ == "__main__":
    raise SystemExit(main())
