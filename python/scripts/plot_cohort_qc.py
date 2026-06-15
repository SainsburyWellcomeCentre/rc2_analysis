"""Cohort QC for the GPU Gabor extraction (SF/OR per RF across 36 clouds).

Cheap, no recompute — reads the per-cloud parquets and makes:
  (1) cohort_qc_overview.png — outlier spotter: RF positions/coverage on a frame,
      per-RF SF recovery vs the 3 tokens, per-RF OR offset-from-token, concentration.
  (2) cohort_qc_per_rf.png — one tiny SF-recovery panel per RF (obs vs token at the
      3 levels), red-bordered if anything looks off (SF, OR offset, or low concentration),
      so a weird single RF the pooled stats would hide is visible at a glance.
"""
from __future__ import annotations
import glob, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpi
from matplotlib.patches import Circle

GOG = os.path.expanduser("~/local_data/motion_clouds/saved_goggles")
COH = os.path.join(GOG, "_extract", "cohort")
FIG = os.path.join(GOG, "_figs"); os.makedirs(FIG, exist_ok=True)
DPP = 111.6 / 400.0
SF_WIN, OR_WIN = 25, 10


def cmean(d): return np.degrees(np.angle(np.mean(np.exp(1j * 2 * np.radians(d)))) / 2) % 180
def cstd(d):
    z = np.mean(np.exp(1j * 2 * np.radians(d))); return np.degrees(np.sqrt(max(-2 * np.log(abs(z)), 0))) / 2


def load():
    df = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(os.path.join(COH, "*.parquet")))],
                   ignore_index=True)
    df["sf_tok"] = (df["cloud"].str.extract(r"_sf(\d+p\d+)_")[0].str.replace("p", ".").astype(float) / DPP).round(3)
    df["th_tok"] = (np.degrees(df["cloud"].str.extract(r"theta(-?\d+p\d+)")[0].str.replace("p", ".").astype(float)) % 180).round()
    df["or_off"] = ((df["or_deg"] - df["th_tok"] + 90) % 180) - 90
    df["rf"] = df["probe"].str[-3:] + ":" + df["cluster"] + df["rf_type"].str[0].str.upper()
    return df


def per_rf(df):
    recs = []
    for rf, g in df.groupby("rf"):
        rec = dict(rf=rf, cx=g.cx.iloc[0], cy=g.cy.iloc[0], conc=g.concentration.median(),
                   edge=g.edge.max(), or_off_mean=cmean(g.or_deg.values - 0) if False else
                   (((cmean(((g.or_deg - g.th_tok) % 180).values) + 90) % 180) - 90),
                   or_sigma=cstd(g.or_off.values))
        for lv, gg in g.groupby("sf_tok"):
            rec[f"sf_{lv:.3f}"] = gg.sf_cpd.mean()
        recs.append(rec)
    return pd.DataFrame(recs)


def overview(df, R, sf_levels):
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
    # (A) RF coverage on a frame
    f0 = sorted(glob.glob(os.path.join(GOG, df.cloud.iloc[0], "*.png")))[650]
    im = mpi.imread(f0); im = im[..., :3].mean(-1) if im.ndim == 3 else im
    ax[0].imshow(im, cmap="gray")
    for _, r in R.iterrows():
        ax[0].add_patch(Circle((r.cx, r.cy), SF_WIN / DPP / 2, fill=False,
                               ec="red" if r.edge > 0 else "lime", lw=0.6, alpha=0.6))
    ax[0].invert_xaxis(); ax[0].set(xticks=[], yticks=[],
                                    title=f"RF coverage (n={len(R)}, SF windows)\nred=edge-clipped")
    # (B) SF recovery
    for lv in sf_levels:
        ax[1].scatter(np.full(len(R), lv), R[f"sf_{lv:.3f}"], s=10, alpha=0.5)
        ax[1].axhline(lv, color="0.6", ls=":", lw=0.8)
    lim = [0, max(sf_levels) * 1.3]; ax[1].plot(lim, lim, "k--", lw=0.8)
    ax[1].set(xlabel="SF token (cpd)", ylabel="per-RF obs SF (cpd)", xlim=lim, ylim=lim,
              title="(B) SF recovery per RF\n(should sit on y=x)")
    # (C) OR offset per RF
    ax[2].axvline(0, color="r", ls="--", lw=1)
    ax[2].errorbar(R.or_off_mean, np.arange(len(R)), xerr=R.or_sigma, fmt="o", ms=3, lw=0.5, alpha=0.5)
    ax[2].set(xlabel="OR − token (°)  (mean ± σ)", ylabel="RF index", xlim=(-90, 90),
              title="(C) OR offset per RF\n(centred at token=0)")
    # (D) concentration
    ax[3].hist(R.conc, bins=30, color="0.6"); ax[3].axvline(2, color="r", ls=":", lw=1)
    ax[3].set(xlabel="energy concentration (reliability)", ylabel="# RFs",
              title=f"(D) reliability\n{(R.conc < 2).sum()} RFs < 2 (red)")
    fig.suptitle(f"Cohort QC overview — {len(R)} RFs × 36 clouds (GPU extraction)", y=1.02, fontsize=12)
    fig.tight_layout(); out = os.path.join(FIG, "cohort_qc_overview.png")
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig); print(f"[saved] {out}")


def per_rf_grid(R, sf_levels):
    n = len(R); ncol = 10; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 1.6, nrow * 1.6))
    lim = [0, max(sf_levels) * 1.3]
    for i, (_, r) in enumerate(R.sort_values("rf").reset_index(drop=True).iterrows()):
        ax = axes.flat[i]
        obs = [r[f"sf_{lv:.3f}"] for lv in sf_levels]
        bad = (any(abs(o - lv) > 0.4 * lv for o, lv in zip(obs, sf_levels))
               or abs(r.or_off_mean) > 20 or r.conc < 2 or r.edge > 0)
        ax.plot(lim, lim, "k--", lw=0.5)
        ax.scatter(sf_levels, obs, s=14, c="C3" if bad else "C0")
        ax.set(xlim=lim, ylim=lim, xticks=[], yticks=[])
        ax.set_title(f"{r.rf}\nΔor{r.or_off_mean:+.0f} c{r.conc:.0f}", fontsize=5.5,
                     color="red" if bad else "black")
        for s in ax.spines.values():
            s.set_color("red" if bad else "0.7"); s.set_linewidth(1.4 if bad else 0.5)
    for j in range(n, nrow * ncol):
        axes.flat[j].axis("off")
    nbad = int(((R[[f"sf_{lv:.3f}" for lv in sf_levels]].sub(sf_levels).abs()
                 .gt([0.4 * lv for lv in sf_levels])).any(axis=1) | (R.or_off_mean.abs() > 20)
                | (R.conc < 2) | (R.edge > 0)).sum())
    fig.suptitle(f"Per-RF QC — SF obs vs token (3 levels) · red = flagged "
                 f"(SF off / |OR−tok|>20° / conc<2 / edge): {nbad}/{len(R)}", y=1.0, fontsize=11)
    fig.tight_layout(); out = os.path.join(FIG, "cohort_qc_per_rf.png")
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig); print(f"[saved] {out}")


def main():
    df = load()
    sf_levels = sorted(df.sf_tok.unique())
    R = per_rf(df)
    print(f"{len(R)} RFs; {(R.conc < 2).sum()} low-conc, {(R.edge > 0).sum()} edge, "
          f"{(R.or_off_mean.abs() > 20).sum()} |OR offset|>20°")
    overview(df, R, sf_levels)
    per_rf_grid(R, sf_levels)
    R.to_csv(os.path.join(COH, "per_rf_summary.csv"), index=False)
    print(f"[saved] {os.path.join(COH, 'per_rf_summary.csv')}")


if __name__ == "__main__":
    main()
