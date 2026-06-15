"""Multisensory linearity test on rc2 motion-clouds firing rates.

Question: per cluster, does ΔR_VT ≈ ΔR_T + ΔR_V (rate-additive integration of
vestibular and visual drive)? No GLM machinery — just per-trial mean rates from
the MATLAB-cached `stationary_vs_motion_fr/<probe>.csv` tables, restricted to the
88-cluster `should_run_glm=1` set. Per condition X ∈ {T_Vstatic, V, VT}, ΔR_X is
the whole firing-rate change from stationary to motion (not a model-decomposed
visual or vestibular component) — modality isolation is by experimental
paradigm: V holds the animal stationary, T_Vstatic drops the visual stimulus.

Per-cluster aggregation:
    ΔR_X = mean over trials in condition X ∈ {T_Vstatic, V, VT} of (motion_fr − stationary_fr)

Outputs (under ~/local_data/motion_clouds/figures/glm/exploration/):
    multisensory_linearity_scatter.pdf
    multisensory_linearity_scatter.png
    multisensory_linearity_per_cluster.csv

Run with the rc2_analysis conda env:
    /Users/lauraporta/miniforge3/envs/rc2_analysis/bin/python \
        scripts/multisensory_linearity_test.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/Users/lauraporta/local_data/motion_clouds")
FR_DIR = DATA / "formatted_data" / "csvs" / "stationary_vs_motion_fr"
PREFILTER = DATA / "figures" / "matlab_reference" / "prefilter_decision_tree.csv"
OUT_DIR = DATA / "figures" / "glm" / "exploration" / "multisensory_linearity"

PROBES = ["CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1", "CAA-1123467_rec1"]


def load_long():
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
    print(f"Prefilter clusters with should_run_glm=1: {len(keep)}")
    return long.merge(keep, on=["probe_id", "cluster_id"], how="inner")


def build_wide(long):
    agg = (
        long.groupby(["probe_id", "cluster_id", "trial_group_label"])
            .agg(evoked_mean=("evoked", "mean"),
                 motion_mean=("motion_fr", "mean"),
                 stationary_mean=("stationary_fr", "mean"),
                 n_trials=("trial_id", "nunique"))
            .reset_index()
    )
    wide = agg.pivot_table(
        index=["probe_id", "cluster_id"],
        columns="trial_group_label",
        values=["evoked_mean", "motion_mean", "stationary_mean", "n_trials"],
    )
    wide.columns = [f"{a}__{b}" for a, b in wide.columns]
    return wide.reset_index()


def regress(wide):
    needed = ["evoked_mean__T_Vstatic", "evoked_mean__V", "evoked_mean__VT"]
    sub = wide.dropna(subset=needed).copy()
    print(f"Clusters with all three condition means: {len(sub)} / {len(wide)}")
    dT = sub["evoked_mean__T_Vstatic"].to_numpy()
    dV = sub["evoked_mean__V"].to_numpy()
    dVT = sub["evoked_mean__VT"].to_numpy()

    print("\n=== Test 1: strict additivity dVT ≟ dT + dV (no fit) ===")
    pred = dT + dV
    resid = dVT - pred
    r, p = stats.pearsonr(dVT, pred)
    sst = float(((dVT - dVT.mean()) ** 2).sum())
    sse_id = float(((dVT - pred) ** 2).sum())
    r2_id = 1.0 - sse_id / sst
    print(f"  Pearson r(dT+dV, dVT): {r:.3f} (p={p:.2e})")
    print(f"  R^2 against y=x: {r2_id:.3f}")
    print(f"  Residual mean {resid.mean():+.3f}, median {np.median(resid):+.3f}, std {resid.std():.3f} spk/s")

    print("\n=== Test 2: best-fit dVT = a·dT + b·dV + c (least squares) ===")
    X = np.column_stack([dT, dV, np.ones_like(dT)])
    coef, *_ = np.linalg.lstsq(X, dVT, rcond=None)
    a, b, c = coef
    fit = X @ coef
    sse = float(((dVT - fit) ** 2).sum())
    r2 = 1.0 - sse / sst
    n = len(dVT)
    sigma2 = sse / (n - 3)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    print(f"  a (dT slope) = {a:+.3f} ± {se[0]:.3f}")
    print(f"  b (dV slope) = {b:+.3f} ± {se[1]:.3f}")
    print(f"  c (intercept) = {c:+.3f} ± {se[2]:.3f} spk/s")
    print(f"  R^2 = {r2:.3f}")
    for name, val, s in [("a=1", a - 1.0, se[0]), ("b=1", b - 1.0, se[1]), ("c=0", c, se[2])]:
        t = val / s
        pp = 2 * (1 - stats.t.cdf(abs(t), df=n - 3))
        print(f"  H0 {name}: t={t:+.2f}, p={pp:.3e}")

    print("\n=== Test 3: per-probe regression (sanity) ===")
    print(f"  {'probe':<22} {'n':>4}  {'a':>7}  {'b':>7}  {'c':>7}  {'R2':>5}")
    for probe, g in sub.groupby("probe_id"):
        x = g["evoked_mean__T_Vstatic"].to_numpy()
        y = g["evoked_mean__V"].to_numpy()
        z = g["evoked_mean__VT"].to_numpy()
        Xp = np.column_stack([x, y, np.ones_like(x)])
        cp, *_ = np.linalg.lstsq(Xp, z, rcond=None)
        fp = Xp @ cp
        sst_p = ((z - z.mean()) ** 2).sum()
        sse_p = ((z - fp) ** 2).sum()
        r2p = 1.0 - sse_p / sst_p
        print(f"  {probe:<22} {len(g):>4}  {cp[0]:+7.3f}  {cp[1]:+7.3f}  {cp[2]:+7.3f}  {r2p:5.3f}")

    sub["pred_strict"] = pred
    sub["resid_strict"] = resid
    sub["pred_fit"] = fit
    sub["resid_fit"] = dVT - fit
    return sub, dict(a=a, b=b, c=c, r2=r2, r2_strict=r2_id, n=n,
                     se_a=se[0], se_b=se[1], se_c=se[2])


def plot_scatter(sub, fp, out_pdf, out_png):
    import matplotlib.pyplot as plt

    dT = sub["evoked_mean__T_Vstatic"].to_numpy()
    dV = sub["evoked_mean__V"].to_numpy()
    dVT = sub["evoked_mean__VT"].to_numpy()
    fit = fp["a"] * dT + fp["b"] * dV + fp["c"]

    fig, ax = plt.subplots(1, 1, figsize=(6, 5.5))

    lo, hi = float(min(fit.min(), dVT.min())), float(max(fit.max(), dVT.max()))
    pad = 0.05 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=1, alpha=0.5, label="y = x")
    ax.scatter(fit, dVT, s=22, alpha=0.7, edgecolor="k", linewidth=0.3, color="C1")
    ax.set_xlabel(
        rf"${fp['a']:+.2f}\,\mathrm{{T_{{vstatic}}}} {fp['b']:+.2f}\,\mathrm{{VF}} {fp['c']:+.2f}$ "
        "(spk/s, evoked)"
    )
    ax.set_ylabel(r"$\mathrm{VF+T_{vstatic}}$ (spk/s, evoked)")
    ax.set_title(
        "Multisensory linearity: best-fit additive\n"
        rf"$\mathrm{{VF+T_{{vstatic}}}} = {fp['a']:+.2f}\,\mathrm{{T_{{vstatic}}}}"
        rf" {fp['b']:+.2f}\,\mathrm{{VF}} {fp['c']:+.2f}$"
        f",  R²={fp['r2']:.3f},  n={fp['n']}"
    )
    ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
    ax.legend(loc="upper left", fontsize=9); ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    print(f"\nWrote: {out_pdf}\n       {out_png}")


def main():
    long = load_long()
    print(f"Loaded {len(long)} per-trial rows across {long['probe_id'].nunique()} probes")
    long = filter_to_glm_clusters(long)
    wide = build_wide(long)
    sub, fp = regress(wide)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "multisensory_linearity_per_cluster.csv"
    pdf_path = OUT_DIR / "multisensory_linearity_scatter.pdf"
    png_path = OUT_DIR / "multisensory_linearity_scatter.png"

    sub.to_csv(csv_path, index=False)
    print(f"\nWrote per-cluster table: {csv_path}")
    plot_scatter(sub, fp, pdf_path, png_path)


if __name__ == "__main__":
    main()
