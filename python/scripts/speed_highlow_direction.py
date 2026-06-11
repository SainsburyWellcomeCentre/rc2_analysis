"""Per-cell HIGH-vs-LOW speed tuning direction (marginal), to test whether the
population "motion-counts" hides heterogeneous single-cell speed tuning.

Companion to speed_shuffle_null.py. That script asks "does each cell carry
within-motion speed-value tuning above the platform-on step" (a magnitude +
significance). This one asks "in which DIRECTION" — does the cell fire more at
HIGH or LOW translation speed — read from the MARGINAL prediction (not the bare
partial kernel; see the partial-vs-marginal lesson):

  Fit the full model (Speed/TF/SF/OR + ME_face + refractory History + Accel) once
  per cell on all rows. Then override the SPEED column on motion rows to a low
  value (20th pct of observed motion speed) and to a high value (80th pct),
  rebuild B_speed, predict, and average the predicted rate over motion rows.
      direction = log(mu_high / mu_low)
  Everything else stays at its observed value, so this is the marginal speed
  effect holding the other regressors at the data.

Population read-out: the DISTRIBUTION of per-cell direction (not its mean). If
cells split into high- and low-preferring (wide, two-sided), the flat population
average is heterogeneous speed tuning averaging out — not "no speed tuning". Merge
with speed_shuffle_null.csv to colour cells that individually beat the trajectory
(within-motion) null.

Output: <run>/diagnostics/speed_highlow_direction.{csv,pdf,png}  (+ --bin10)
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import (convolve_history, history_basis, onset_kernel_basis,
                           value_basis)
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.fitting import fit_poisson_glm
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

import sys

import scripts.run_glm_split_by_condition as drv
from scripts.variance_partition_accel import (PROBES, STIM, OUT, _build_accel,
                                              _build_me, _make_cfg)

LO_PCT, HI_PCT = 20.0, 80.0
# --profile-time: replace the single shared onset/time kernel with a SEPARATE
# onset kernel per velocity profile (profile×time interaction). This fully
# saturates time-in-trajectory, common AND profile-specific, so any residual
# motion-onset transient can't masquerade as speed tuning. Speed then only
# explains cross-profile speed differences at MATCHED trajectory time -> the
# onset-confound control on the low-speed-preference result.
PROFILE_TIME = "--profile-time" in sys.argv
SUFFIX = "_profiletime" if PROFILE_TIME else ""


def cell_direction(df: pd.DataFrame, cfg):
    if df["spike_count"].sum() == 0:
        return None
    B_me = _build_me(df, cfg)
    if B_me is None:
        return None
    B_accel = _build_accel(df, cfg)

    speed = df["speed"].to_numpy(float)
    tf = df["tf"].to_numpy(float)
    onset = df["time_since_onset"].to_numpy(float)
    sf = df["sf"].to_numpy(float)
    orient = df["orientation"].to_numpy(float)
    y = df["spike_count"].to_numpy(float)
    trial_ids = df["trial_id"].to_numpy(np.int64)
    profile_ids = df["profile_id"].to_numpy(np.int64)
    cond = df["condition"].to_numpy(object)
    if PROFILE_TIME and len(np.unique(profile_ids)) < 2:
        return None

    motion = cond != "stationary"
    motion_idx = np.flatnonzero(motion)
    pos = motion & (speed > 0)
    if int(pos.sum()) < 10:
        return None
    s_lo = float(np.percentile(speed[pos], LO_PCT))
    s_hi = float(np.percentile(speed[pos], HI_PCT))
    if not (s_hi > s_lo):
        return None

    offset = float(np.log(cfg.time_bin_width))
    spacing = cfg.speed_tf_basis_spacing
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=spacing)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    hb = history_basis(cfg.n_history_bases, cfg.history_window_s,
                       cfg.time_bin_width, kind=cfg.history_basis_kind)
    B_hist = convolve_history(y, trial_ids, hb)
    sel = list(STIM) + ["ME_face", "History"]

    # Prediction mode (pass ref-levels) disables zero-variance column dropping,
    # so the design width is fixed when we later override speed to a constant
    # and reuse the fitted beta (otherwise a column gets re-pruned -> matmul
    # shape mismatch).
    sf_valid = sf[(sf != 0) & ~np.isnan(sf)]
    sf_ref = np.sort(np.unique(sf_valid)) if sf_valid.size else np.array([])
    or_valid = orient[(orient != 0) & ~np.isnan(orient)]
    or_ref = np.sort(np.unique(or_valid)) if or_valid.size else np.array([])

    # Profile-specific onset blocks (B_onset gated by each profile indicator).
    # Their sum is the shared onset, so we drop the shared kernel to avoid the
    # exact collinearity and let these saturate time per profile (ridge handles
    # any residual conditioning).
    onset_extra = None
    if PROFILE_TIME:
        onset_extra = np.hstack([B_onset * (profile_ids == p)[:, None]
                                 for p in np.unique(profile_ids)])

    def design(B_speed_use):
        X, _ = assemble_design_matrix_selected(
            B_speed_use, B_tf, B_onset, sf, orient, sel,
            sf_ref_levels=sf_ref, or_ref_levels=or_ref,
            B_me_face=B_me, B_history=B_hist,
            include_onset_kernel=not PROFILE_TIME)
        blocks = [X, onset_extra, B_accel] if PROFILE_TIME else [X, B_accel]
        return np.hstack(blocks)

    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=spacing)
    X_real = design(B_speed)
    if X_real.shape[1] >= y.size:
        return None
    fit = fit_poisson_glm(X_real, y, offset, lambda_ridge=cfg.lambda_ridge)
    beta = np.asarray(fit.beta, float)

    def marginal_rate(s_val):
        sp = speed.copy()
        sp[motion_idx] = s_val           # override only motion rows
        X = design(value_basis(sp, cfg.n_speed_bases, *cfg.speed_range, spacing=spacing))
        mu = np.exp(np.clip(X @ beta + offset, -20, 20))
        return float(mu[motion_idx].mean())

    mu_lo = marginal_rate(s_lo)
    mu_hi = marginal_rate(s_hi)
    return {
        "speed_lo": s_lo, "speed_hi": s_hi,
        "mu_lo": mu_lo, "mu_hi": mu_hi,
        "direction_logratio": float(np.log(max(mu_hi, 1e-12) / max(mu_lo, 1e-12))),
        "n_spikes": int(y.sum()),
    }


def _figure(df: pd.DataFrame, out_dir, run_name: str) -> None:
    nlc = "speed_shuffle_null.csv"
    null_path = out_dir / nlc
    if null_path.exists():
        nd = pd.read_csv(null_path)[["probe_id", "cluster_id", "p_traj", "total_full"]]
        df = df.merge(nd, on=["probe_id", "cluster_id"], how="left")
    else:
        df["p_traj"] = np.nan
        df["total_full"] = np.nan
    well = df[(df["total_full"] > 0.005) | df["total_full"].isna()].copy()
    d = well["direction_logratio"].to_numpy()
    sig = (well["p_traj"] < 0.05).to_numpy()
    n = len(well)
    n_hi = int((d > 0).sum())
    n_lo = int((d < 0).sum())
    n_sig = int(np.nansum(sig))
    sig_hi = int(((d > 0) & sig).sum())
    sig_lo = int(((d < 0) & sig).sum())

    fig, ax = plt.subplots(figsize=(8.5, 5))
    bins = np.linspace(np.floor(d.min() * 10) / 10, np.ceil(d.max() * 10) / 10, 31)
    ax.hist(d[~sig], bins=bins, color="0.7", alpha=0.8, label="n.s. vs within-motion null")
    ax.hist(d[sig], bins=bins, color="tab:green", alpha=0.8,
            label="beats within-motion null (p<0.05)")
    ax.axvline(0, color="black", lw=1.2)
    ax.axvline(np.median(d), color="tab:red", lw=1.4, ls="--",
               label=f"median {np.median(d):+.3f}")
    from scipy.stats import wilcoxon
    try:
        w_p = wilcoxon(d).pvalue           # is the population median ≠ 0?
    except ValueError:
        w_p = float("nan")
    shift = ("LOW-speed-preferring" if np.median(d) < 0 else "HIGH-speed-preferring")
    tag = " · per-profile time kernel (onset-confound control)" if PROFILE_TIME else ""
    ax.set_xlabel("per-cell speed tuning direction  log(rate@high / rate@low)")
    ax.set_ylabel("n cells")
    ax.set_title(
        f"Single-cell speed tuning direction ({run_name}, n={n}){tag}\n"
        f"low-pref {n_lo} · high-pref {n_hi} (Wilcoxon median≠0 p={w_p:.1e}) — "
        f"net {shift}; sig vs within-motion null {n_sig} (low {sig_lo} / high {sig_hi})",
        fontsize=9.5, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    out = out_dir / f"speed_highlow_direction{SUFFIX}"
    fig.savefig(f"{out}.pdf")
    fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"\nn={n}: high-pref {n_hi}, low-pref {n_lo}; "
          f"median direction {np.median(d):+.4f}")
    print(f"  beats within-motion null: {n_sig} (high {sig_hi} / low {sig_lo})")
    print(f"wrote {out}.pdf / .png")


def main() -> int:
    cfg = _make_cfg()
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    rows = []
    for probe in PROBES:
        cohort = drv.load_cohort(probe)
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        for cl in pdata.clusters:
            if cl.cluster_id not in cohort:
                continue
            res = cell_direction(bin_cluster(pdata, cl), cfg)
            if res is None:
                continue
            res.update(probe_id=probe, cluster_id=int(cl.cluster_id))
            rows.append(res)
        print(f"{probe}: {sum(r['probe_id'] == probe for r in rows)} cells")

    df = pd.DataFrame(rows)
    out_dir = OUT / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"speed_highlow_direction{SUFFIX}.csv", index=False)
    print(f"wrote {out_dir / f'speed_highlow_direction{SUFFIX}.csv'} ({len(df)} cells)")
    _figure(df, out_dir, OUT.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
