"""Speed-resolved multisensory linearity test.

Question: does the rate-additive linearity (VF+T_vstatic ≈ a·T_vstatic + b·VF + c)
hold when trials are split into lower vs higher speed subsets, or does the
integration rule shift with stimulus strength?

Method:
  1. Load each probe via `rc2_glm.io.load_probe_data` to get per-trial mean motion-period speed.
     For T_Vstatic / VT the velocity channel is `stage` (or `filtered_teensy`); for V it is
     `multiplexer_output` (visual-flow replay speed). Units match by experimental design
     (V trials replay a prior `stage` trajectory).
  2. Median-split trials within each (probe, condition) into low / high speed subsets.
  3. Per (probe, cluster, condition, split): mean (motion_fr − stationary_fr) from the
     MATLAB-cached per-trial firing rate CSVs.
  4. For each split, run OLS: VF+T_vstatic = a·T_vstatic + b·VF + c on the 88 cluster rows.
  5. Plot 1×2: low-speed panel | high-speed panel, same axes.

Cache: per-trial mean motion speed table at
    OUT_DIR / "per_trial_motion_speed.csv"
(rebuilt automatically if missing or older than any .mat file).

Outputs:
    multisensory_linearity_speed_split.{pdf,png}
    multisensory_linearity_per_cluster_split.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/Users/lauraporta/local_data/motion_clouds")
FR_DIR = DATA / "formatted_data" / "csvs" / "stationary_vs_motion_fr"
PREFILTER = DATA / "figures" / "matlab_reference" / "prefilter_decision_tree.csv"
OUT_DIR = DATA / "figures" / "glm" / "exploration" / "multisensory_linearity"
MAT_DIR = DATA / "formatted_data"

PROBES = ["CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1", "CAA-1123467_rec1"]


def build_per_trial_speed():
    """Load .mat per probe and emit (probe_id, trial_id, condition, mean_motion_speed)."""
    from rc2_glm.io import load_probe_data
    rows = []
    for p in PROBES:
        probe = load_probe_data(MAT_DIR / f"{p}.mat")
        for tr in probe.trials:
            if tr.excluded:
                continue
            mask = tr.motion_mask.astype(bool)
            if mask.sum() == 0:
                rows.append((p, int(tr.trial_id), tr.condition, np.nan))
                continue
            mean_speed = float(np.mean(tr.velocity[mask]))
            rows.append((p, int(tr.trial_id), tr.condition, mean_speed))
    df = pd.DataFrame(rows, columns=["probe_id", "trial_id", "condition", "mean_motion_speed"])
    return df


def get_per_trial_speed():
    cache = OUT_DIR / "per_trial_motion_speed.csv"
    if cache.exists():
        cache_mtime = cache.stat().st_mtime
        mat_mtime = max((MAT_DIR / f"{p}.mat").stat().st_mtime for p in PROBES)
        if cache_mtime > mat_mtime:
            print(f"Reading cached per-trial speeds: {cache}")
            return pd.read_csv(cache)
    print("Computing per-trial speeds from .mat (this takes ~30 s)...")
    df = build_per_trial_speed()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    print(f"Cached per-trial speeds → {cache}")
    return df


def median_split(speed_df):
    """Add a 'speed_split' column ∈ {low, high} per (probe, condition) median."""
    speed_df = speed_df.dropna(subset=["mean_motion_speed"]).copy()
    medians = (speed_df.groupby(["probe_id", "condition"])["mean_motion_speed"]
               .median().rename("median_speed").reset_index())
    speed_df = speed_df.merge(medians, on=["probe_id", "condition"])
    speed_df["speed_split"] = np.where(
        speed_df["mean_motion_speed"] >= speed_df["median_speed"], "high", "low"
    )
    print("\nPer (probe, condition) median speed (cm/s; for V it's replay flow speed):")
    print(medians.to_string(index=False))
    print("\nTrial counts per (probe, condition, split):")
    print(speed_df.groupby(["probe_id", "condition", "speed_split"]).size().unstack(fill_value=0))
    return speed_df


def load_fr_long():
    frames = []
    for p in PROBES:
        df = pd.read_csv(FR_DIR / f"{p}.csv")
        df["probe_id"] = p
        df["evoked"] = df["motion_fr"] - df["stationary_fr"]
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def filter_to_glm_clusters(long):
    pf = pd.read_csv(PREFILTER)
    keep = pf.loc[pf["should_run_glm"] == 1, ["probe_id", "cluster_id"]]
    return long.merge(keep, on=["probe_id", "cluster_id"], how="inner"), keep


def build_wide(long):
    agg = (long.groupby(["probe_id", "cluster_id", "trial_group_label"])
                .agg(evoked_mean=("evoked", "mean"),
                     n_trials=("trial_id", "nunique"))
                .reset_index())
    wide = agg.pivot_table(
        index=["probe_id", "cluster_id"],
        columns="trial_group_label",
        values=["evoked_mean", "n_trials"],
    )
    wide.columns = [f"{a}__{b}" for a, b in wide.columns]
    return wide.reset_index()


def regress_one(sub, label):
    needed = ["evoked_mean__T_Vstatic", "evoked_mean__V", "evoked_mean__VT"]
    sub = sub.dropna(subset=needed).copy()
    dT = sub["evoked_mean__T_Vstatic"].to_numpy()
    dV = sub["evoked_mean__V"].to_numpy()
    dVT = sub["evoked_mean__VT"].to_numpy()
    X = np.column_stack([dT, dV, np.ones_like(dT)])
    coef, *_ = np.linalg.lstsq(X, dVT, rcond=None)
    a, b, c = coef
    fit = X @ coef
    sst = float(((dVT - dVT.mean()) ** 2).sum())
    sse = float(((dVT - fit) ** 2).sum())
    r2 = 1.0 - sse / sst
    n = len(dVT)
    sigma2 = sse / (n - 3)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    pa = 2 * (1 - stats.t.cdf(abs((a - 1.0) / se[0]), df=n - 3))
    pb = 2 * (1 - stats.t.cdf(abs((b - 1.0) / se[1]), df=n - 3))
    print(f"\n[{label}] n={n}  a={a:+.3f}±{se[0]:.3f} (p_vs_1={pa:.3e})  "
          f"b={b:+.3f}±{se[1]:.3f} (p_vs_1={pb:.3e})  c={c:+.3f}±{se[2]:.3f}  R²={r2:.3f}")
    out = sub.copy()
    out["dT"] = dT; out["dV"] = dV; out["dVT"] = dVT
    out["fit"] = fit; out["resid_fit"] = dVT - fit
    out["split"] = label
    return out, dict(label=label, a=a, b=b, c=c, r2=r2, n=n,
                     se_a=se[0], se_b=se[1], se_c=se[2], p_a=pa, p_b=pb)


def plot_two(low_sub, low_fp, high_sub, high_fp, out_pdf, out_png):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True, sharey=True)

    all_x = np.concatenate([low_sub["fit"].to_numpy(), high_sub["fit"].to_numpy(),
                             low_sub["dVT"].to_numpy(), high_sub["dVT"].to_numpy()])
    lo, hi = float(all_x.min()), float(all_x.max())
    pad = 0.05 * (hi - lo)

    for ax, sub, fp, color in [(axes[0], low_sub, low_fp, "C0"),
                                (axes[1], high_sub, high_fp, "C3")]:
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=1, alpha=0.5, label="y = x")
        ax.scatter(sub["fit"], sub["dVT"], s=22, alpha=0.7, edgecolor="k", linewidth=0.3, color=color)
        ax.set_xlabel(
            rf"${fp['a']:+.2f}\,\mathrm{{T_{{vstatic}}}} {fp['b']:+.2f}\,\mathrm{{VF}} {fp['c']:+.2f}$ "
            "(spk/s, evoked)"
        )
        ax.set_title(
            f"{fp['label']}-speed trials\n"
            rf"$\mathrm{{VF+T_{{vstatic}}}} = {fp['a']:+.2f}\,\mathrm{{T_{{vstatic}}}}"
            rf" {fp['b']:+.2f}\,\mathrm{{VF}} {fp['c']:+.2f}$"
            f",  R²={fp['r2']:.3f},  n={fp['n']}"
        )
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        ax.legend(loc="upper left", fontsize=9); ax.grid(alpha=0.3)

    axes[0].set_ylabel(r"$\mathrm{VF+T_{vstatic}}$ (spk/s, evoked)")
    fig.suptitle("Multisensory linearity, speed-resolved (per-probe-per-condition median split)", y=1.02)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    print(f"\nWrote: {out_pdf}\n       {out_png}")


def main():
    speed_df = get_per_trial_speed()
    speed_df = median_split(speed_df)

    long = load_fr_long()
    long, _keep = filter_to_glm_clusters(long)

    # condition column in CSV is `trial_group_label`; in speed_df it's `condition`
    long = long.merge(
        speed_df[["probe_id", "trial_id", "condition", "speed_split"]],
        left_on=["probe_id", "trial_id", "trial_group_label"],
        right_on=["probe_id", "trial_id", "condition"],
        how="inner",
    )
    print(f"\nFR rows after speed-split join: {len(long)}")

    low_long = long[long["speed_split"] == "low"]
    high_long = long[long["speed_split"] == "high"]

    low_wide = build_wide(low_long)
    high_wide = build_wide(high_long)

    low_sub, low_fp = regress_one(low_wide, "low")
    high_sub, high_fp = regress_one(high_wide, "high")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    combined = pd.concat([low_sub.assign(split="low"), high_sub.assign(split="high")],
                         ignore_index=True)
    csv_path = OUT_DIR / "multisensory_linearity_per_cluster_split.csv"
    combined.to_csv(csv_path, index=False)
    print(f"\nWrote per-cluster table: {csv_path}")

    pdf_path = OUT_DIR / "multisensory_linearity_speed_split.pdf"
    png_path = OUT_DIR / "multisensory_linearity_speed_split.png"
    plot_two(low_sub, low_fp, high_sub, high_fp, pdf_path, png_path)


if __name__ == "__main__":
    main()
