"""Generalised VIF for the full current_full_20ms design — ALL predictors:
Speed / TF / Onset / SF / OR / ME_face / History / Acceleration / Population.

History, Acceleration and Population are per-cluster (History & Pop depend on
spikes), so computed per cohort cluster and medianed. Population = same-bin LOO
mean rate of all other VISp cells (one column), as in variance_partition_full.

Output: current_full_20ms/diagnostics/collinearity_vif_full.{pdf,png}
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import (convolve_history, history_basis, onset_kernel_basis,
                           raised_cosine_basis_linear, value_basis)
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv
from scripts.run_glm_current_plus_ME_20ms import make_config_me_20ms
from scripts.collinearity_vif import _band, BAND_COLOR, GOOD, WORRY
from scripts.collinearity_vif_me import _gvif
from scripts.variance_partition_full import _build_me, _zfill

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_full_20ms" / "diagnostics"
BLOCKS = ["Speed", "TF", "Onset", "SF", "OR", "ME_face", "History", "Accel", "Population"]
PREFIX = {b: b + "_" for b in BLOCKS}


def _design(df, cfg, pop_loo):
    speed = df["speed"].to_numpy(float); tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float); orient = df["orientation"].to_numpy(float)
    y = df["spike_count"].to_numpy(float); trial_ids = df["trial_id"].to_numpy(np.int64)
    spc = cfg.speed_tf_basis_spacing
    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=spc)
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=spc)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    X, names = assemble_design_matrix_selected(
        B_speed, B_tf, B_onset, sf, orient, ["Speed", "TF", "SF", "OR"],
        include_onset_kernel=True, B_me_face=None)
    keep = [i for i, nm in enumerate(names) if nm != "Intercept"]
    X = X[:, keep]; names = [names[i] for i in keep]
    B_me = _build_me(df, cfg)
    if B_me is None:
        return None
    hb = history_basis(cfg.n_history_bases, cfg.history_window_s, cfg.time_bin_width,
                       kind=cfg.history_basis_kind)
    B_hist = convolve_history(y, trial_ids, hb)
    B_accel = raised_cosine_basis_linear(np.clip(_zfill(df["acceleration"].to_numpy(float)), -3, 3),
                                         5, -3.0, 3.0)
    B_pop = _zfill(pop_loo).reshape(-1, 1)
    X = np.hstack([X, B_me, B_hist, B_accel, B_pop])
    names += [f"ME_face_{i+1}" for i in range(B_me.shape[1])]
    names += [f"History_{i+1}" for i in range(B_hist.shape[1])]
    names += [f"Accel_{i+1}" for i in range(B_accel.shape[1])]
    names += ["Population_1"]
    return X, names


def main() -> int:
    cfg = make_config_me_20ms(history_window_s=0.04)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    rows = []
    for probe in PROBES:
        cohort = drv.load_cohort(probe)
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        binned = {cl.cluster_id: bin_cluster(pdata, cl) for cl in pdata.clusters}
        ref_len = max(len(d) for d in binned.values())
        usable = {cid: d for cid, d in binned.items() if len(d) == ref_len}
        counts = np.zeros(ref_len)
        for d in usable.values():
            counts += d["spike_count"].to_numpy(float)
        for cl in pdata.clusters:
            if cl.cluster_id not in cohort or cl.cluster_id not in usable:
                continue
            d = usable[cl.cluster_id]
            out = _design(d, cfg, counts - d["spike_count"].to_numpy(float))
            if out is None:
                continue
            X, names = out
            Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
            R = np.corrcoef(Xs, rowvar=False)
            rec = {"cluster_id": int(cl.cluster_id), "probe": probe}
            for b in BLOCKS:
                idx = [i for i, nm in enumerate(names) if nm.startswith(PREFIX[b])]
                if idx:
                    g = _gvif(R, idx)
                    rec[f"gvif_{b}"] = g ** (1.0 / (2 * len(idx))) if np.isfinite(g) else np.inf
            rows.append(rec)
        print(f"{probe}: {sum(r['probe']==probe for r in rows)} clusters")

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "collinearity_vif_full.csv", index=False)
    vals = [df[f"gvif_{b}"].replace(np.inf, np.nan).median() for b in BLOCKS]
    print("\n=== full current_full_20ms GVIF (median over cohort) ===")
    for b, v in zip(BLOCKS, vals):
        print(f"  {b:11s} = {v:.2f}")

    colors = [BAND_COLOR[_band(v)] for v in vals]
    fig, ax = plt.subplots(figsize=(9.5, 4.3))
    bars = ax.bar(BLOCKS, vals, color=colors, edgecolor="white")
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 0.02, f"{v:.2f}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.axhline(GOOD, color="0.4", ls="--", lw=1); ax.axhline(WORRY, color="0.4", ls=":", lw=1)
    ax.text(len(BLOCKS) - 0.6, GOOD + 0.02, "worrying", fontsize=7.5, color="0.4")
    ax.text(len(BLOCKS) - 0.6, WORRY + 0.02, "severe", fontsize=7.5, color="0.4")
    ax.set_ylabel("generalised VIF (SE scale, GVIF$^{1/2df}$)")
    ax.set_ylim(0, max(WORRY + 0.3, max(vals) * 1.15))
    ax.tick_params(axis="x", labelrotation=20)
    ax.set_title(f"Full design VIF — all predictors (median over {len(df)} cohort clusters)",
                 fontsize=10, fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=BAND_COLOR[b]) for b in ("good", "worrying", "severe")]
    ax.legend(handles, ["good (<2.24)", "worrying (2.24–3.16)", "severe (>3.16)"],
              fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    out = OUT / "collinearity_vif_full"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
