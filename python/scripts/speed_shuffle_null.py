"""Per-cell circular-shift NULL for Speed's unique cv-bps (the noise floor).

Answers the statistician's objection: the surviving unique-Speed (~+0.007 at
20 ms, ~+0.013 at 10 ms; positive in ~35/77 cells) has no reference scale, so it
can't be told apart from CV-estimation noise. This builds a per-cell null:

  Keep the full design fixed. Circularly shift the per-bin SPEED values WITHIN
  each trial's motion rows (independent random roll per trial). That preserves
  Speed's value distribution AND its within-trial autocorrelation, but destroys
  the alignment of the velocity TRAJECTORY to the spike train. Rebuild B_speed
  from the shifted speed, recompute

      unique_Speed = cv_full(shifted Speed) - cv_full(no Speed)

  N times -> per-cell null distribution. cv_full(no Speed) is identical to the
  observed case (Speed absent either way), so the test reduces to: does the REAL
  speed-spike alignment beat a time-scrambled one. Real if observed unique
  exceeds the 95th percentile of its own null.

Everything else (TF, SF, OR, ME_face, refractory History, Acceleration) stays at
its real, unshuffled value — so this is specifically Speed's UNIQUE temporal
information beyond the other regressors, matching variance_partition_accel's
uSpeed_all definition exactly (same cohort, folds, ridge, cv_full).

Population read-out: how many cells beat their own 95th-pct null vs the 5%
false-positive expectation (one-sided binomial).

Output: <run>/diagnostics/speed_shuffle_null.{csv,pdf,png}
  --bin10  : 10 ms run (current_ME_hist_accel_10ms) instead of 20 ms.
  --n N    : shuffles per cell (default 200).
  --quick  : N=20, first probe only, <=8 cells (smoke test).
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.basis import (convolve_history, history_basis, onset_kernel_basis,
                           value_basis)
from rc2_glm.cross_validation import cross_validate_glm, make_trial_folds
from rc2_glm.design_matrix import assemble_design_matrix_selected
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

# Reuse the EXACT config + ME/Accel basis builders the defended figure used, so
# observed_unique here == uSpeed_all in variance_partition_accel (no drift).
import scripts.run_glm_split_by_condition as drv
from scripts.variance_partition_accel import (PROBES, STIM, OUT, _build_accel,
                                              _build_me, _make_cfg)

N_SHUFFLES = 200
for i, a in enumerate(sys.argv):
    if a == "--n" and i + 1 < len(sys.argv):
        N_SHUFFLES = int(sys.argv[i + 1])
QUICK = "--quick" in sys.argv
if QUICK:
    N_SHUFFLES = 20
BASE_SEED = 20260611


def cell_null(df: pd.DataFrame, cfg, n_shuffles: int, seed: int):
    """Observed unique-Speed + its circular-shift null for one cluster."""
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
    if len(np.unique(profile_ids)) < 2:
        return None

    offset = float(np.log(cfg.time_bin_width))
    spacing = cfg.speed_tf_basis_spacing
    B_tf = value_basis(tf, cfg.n_tf_bases, *cfg.tf_range, spacing=spacing)
    B_onset = onset_kernel_basis(onset, cfg.n_onset_bases, cfg.onset_range[1])
    hb = history_basis(cfg.n_history_bases, cfg.history_window_s,
                       cfg.time_bin_width, kind=cfg.history_basis_kind)
    B_hist = convolve_history(y, trial_ids, hb)
    folds = make_trial_folds(trial_ids, cfg.n_folds, cfg.cv_seed,
                             condition_labels_per_bin=cond,
                             strategy="speed-profile", profile_ids_per_bin=profile_ids)

    def cv(B_speed_use, sel, *, accel_on=True, onset_on=True,
           me_on=True, hist_on=True):
        s = list(sel)
        if me_on:
            s = s + ["ME_face"]
        if hist_on:
            s = s + ["History"]
        X, _ = assemble_design_matrix_selected(
            B_speed_use, B_tf, B_onset, sf, orient, s,
            include_onset_kernel=onset_on,
            B_me_face=B_me if me_on else None,
            B_history=B_hist if hist_on else None)
        if accel_on:
            X = np.hstack([X, B_accel])
        return cross_validate_glm(X, y, offset, folds,
                                  lambda_ridge=cfg.lambda_ridge).cv_bits_per_spike

    B_speed = value_basis(speed, cfg.n_speed_bases, *cfg.speed_range, spacing=spacing)
    noSpeed = [v for v in STIM if v != "Speed"]
    cv_int = cv(B_speed, [], accel_on=False, onset_on=False, me_on=False, hist_on=False)
    cv_noSpeed = cv(B_speed, noSpeed)
    cv_full = cv(B_speed, STIM)
    observed = cv_full - cv_noSpeed

    motion = cond != "stationary"
    motion_idx = np.flatnonzero(motion)
    trial_motion_idx = [np.flatnonzero((trial_ids == t) & motion)
                        for t in np.unique(trial_ids)]
    trial_motion_idx = [ix for ix in trial_motion_idx if ix.size >= 2]
    rng = np.random.default_rng(seed)

    # TWO nulls, both rebuilding B_speed from a perturbed speed vector and
    # recomputing unique = cv_full(perturbed) - cv_noSpeed (cv_noSpeed is
    # shared; stationary rows stay speed=0 in both):
    #  - null_traj: circular roll WITHIN each trial's motion rows. Preserves
    #    the platform-on/off gating (a T/VT trial stays all-nonzero, a V trial
    #    all-zero) and the within-trial value distribution; scrambles only the
    #    velocity TRAJECTORY timing. => "does the time-course matter beyond
    #    'platform is moving'?"
    #  - null_floor: permute speed across ALL motion rows of the cell. Breaks
    #    the gating too => pure same-marginal estimation floor. => "is unique
    #    Speed above what a randomly-aligned same-marginal regressor scores?"
    null_traj = np.empty(n_shuffles)
    null_floor = np.empty(n_shuffles)
    for i in range(n_shuffles):
        sp = speed.copy()
        for ix in trial_motion_idx:
            k = int(rng.integers(1, ix.size))
            sp[ix] = np.roll(speed[ix], k)
        null_traj[i] = cv(value_basis(sp, cfg.n_speed_bases, *cfg.speed_range,
                                      spacing=spacing), STIM) - cv_noSpeed
        spf = speed.copy()
        spf[motion_idx] = rng.permutation(speed[motion_idx])
        null_floor[i] = cv(value_basis(spf, cfg.n_speed_bases, *cfg.speed_range,
                                       spacing=spacing), STIM) - cv_noSpeed

    p_traj = (1 + int((null_traj >= observed).sum())) / (n_shuffles + 1)
    p_floor = (1 + int((null_floor >= observed).sum())) / (n_shuffles + 1)
    return {
        "cv_intercept": cv_int, "cv_noSpeed": cv_noSpeed, "cv_full": cv_full,
        "total_full": cv_full - cv_int, "observed_unique": observed,
        "null_traj_mean": float(null_traj.mean()), "null_traj_std": float(null_traj.std()),
        "null_traj_p95": float(np.quantile(null_traj, 0.95)), "p_traj": p_traj,
        "null_floor_mean": float(null_floor.mean()), "null_floor_std": float(null_floor.std()),
        "null_floor_p95": float(np.quantile(null_floor, 0.95)), "p_floor": p_floor,
        "n_shuffles": n_shuffles, "n_spikes": int(y.sum()),
    }


def _figure(df: pd.DataFrame, out_dir, run_name: str) -> None:
    well = df[df["total_full"] > 0.005].copy().reset_index(drop=True)
    well = well.sort_values("observed_unique").reset_index(drop=True)
    n = len(well)
    sig_floor = (well["p_floor"] < 0.05).to_numpy()   # beats estimation floor
    sig_traj = (well["p_traj"] < 0.05).to_numpy()     # beats trajectory-timing null
    n_floor, n_traj = int(sig_floor.sum()), int(sig_traj.sum())
    bt_floor = binomtest(n_floor, n, 0.05, alternative="greater")
    bt_traj = binomtest(n_traj, n, 0.05, alternative="greater")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 4.8),
                                   gridspec_kw={"width_ratios": [1.7, 1]})

    x = np.arange(n)
    obs = well["observed_unique"].to_numpy()
    # two null bands per cell: floor (mean→p95) and trajectory (mean→p95)
    axL.vlines(x, well["null_floor_mean"], well["null_floor_p95"],
               color="0.8", lw=2.4, zorder=1, label="estimation-floor null (mean→p95)")
    axL.vlines(x, well["null_traj_mean"], well["null_traj_p95"],
               color="tab:orange", alpha=0.45, lw=2.4, zorder=2,
               label="trajectory-shuffle null (mean→p95)")
    axL.scatter(x[~sig_traj], obs[~sig_traj], s=20, color="0.45", zorder=4, label="obs (n.s. vs traj)")
    axL.scatter(x[sig_traj], obs[sig_traj], s=26, color="tab:green", zorder=4,
                label="obs > trajectory p95")
    axL.axhline(0, color="0.6", lw=0.8)
    axL.set_xlabel("cell (sorted by observed unique-Speed)")
    axL.set_ylabel("unique-Speed Δ cv-bps")
    axL.set_title("Observed unique-Speed vs two nulls\n"
                  "floor = same-marginal scramble · trajectory = within-trial circular shift",
                  fontsize=9.5, fontweight="bold")
    axL.legend(fontsize=7.5, frameon=False, loc="upper left")

    # Right: median decomposition floor + platform-on + trajectory = observed.
    floor = np.median(well["null_floor_mean"])
    platform_on = np.median(well["null_traj_mean"]) - floor
    trajectory = np.median(obs) - np.median(well["null_traj_mean"])
    parts = [("est. floor", floor, "0.7"),
             ("platform-on\ncontrast", platform_on, "tab:blue"),
             ("trajectory\ntuning", trajectory, "tab:green")]
    left = 0.0
    for name, val, c in parts:
        axR.barh(0, val, left=left, height=0.5, color=c, edgecolor="white")
        axR.text(left + val / 2, 0, f"{name}\n{val:+.4f}", ha="center", va="center",
                 fontsize=8.5, color="white" if c != "0.7" else "black", fontweight="bold")
        left += val
    axR.axvline(np.median(obs), color="black", lw=1.4, ls="--")
    axR.text(np.median(obs), 0.42, f"median observed\n{np.median(obs):+.4f}",
             ha="center", va="bottom", fontsize=8)
    axR.set_ylim(-0.6, 0.8); axR.set_yticks([])
    axR.set_xlabel("median Δ cv-bps")
    axR.set_title(f"unique-Speed decomposition (n={n})\n"
                  f"beats floor: {n_floor}/{n} (p={bt_floor.pvalue:.1e}) · "
                  f"beats traj: {n_traj}/{n} (p={bt_traj.pvalue:.1e})",
                  fontsize=9.5, fontweight="bold")

    fig.suptitle(f"Speed noise floor ({run_name}, "
                 f"N={int(well['n_shuffles'].iloc[0])} shuffles/cell, n={n})",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = out_dir / "speed_shuffle_null"
    fig.savefig(f"{out}.pdf")
    fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"\nn={n} well-fit")
    print(f"  beats estimation floor: {n_floor}/{n} (binom-greater p={bt_floor.pvalue:.2e})")
    print(f"  beats trajectory null : {n_traj}/{n} (binom-greater p={bt_traj.pvalue:.2e})")
    print(f"  median obs={np.median(obs):+.4f}  floor={floor:+.4f}  "
          f"platform-on={platform_on:+.4f}  trajectory={trajectory:+.4f}")
    print(f"wrote {out}.pdf / .png")


def main() -> int:
    cfg = _make_cfg()
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    probes = PROBES[:1] if QUICK else PROBES

    tasks = []
    for probe in probes:
        cohort = drv.load_cohort(probe)
        pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                                stimulus_lookup=lookup, visp_only=True)
        for cl in pdata.clusters:
            if cl.cluster_id not in cohort:
                continue
            tasks.append((probe, int(cl.cluster_id), bin_cluster(pdata, cl)))
        print(f"{probe}: {sum(t[0] == probe for t in tasks)} clusters binned")
    if QUICK:
        tasks = tasks[:8]

    def _run(probe, cid, df, seed):
        res = cell_null(df, cfg, N_SHUFFLES, seed)
        if res is None:
            return None
        res.update(probe_id=probe, cluster_id=cid)
        return res

    if QUICK:
        rows = [_run(p, c, d, BASE_SEED + i) for i, (p, c, d) in enumerate(tasks)]
    else:
        from joblib import Parallel, delayed
        rows = Parallel(n_jobs=-1)(
            delayed(_run)(p, c, d, BASE_SEED + i)
            for i, (p, c, d) in enumerate(tasks))
    rows = [r for r in rows if r is not None]

    df = pd.DataFrame(rows)
    out_dir = OUT / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "speed_shuffle_null.csv", index=False)
    print(f"wrote {out_dir / 'speed_shuffle_null.csv'} ({len(df)} cells)")
    _figure(df, out_dir, OUT.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
