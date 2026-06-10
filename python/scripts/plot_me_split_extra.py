"""Illustrative plots for the ME across-/within-trial decomposition.

(1) me_decomposition_illustration: what the split IS, on real data —
    a single example trial's ME(t) split into trial-mean (across) + residual
    (within), with the velocity trace overlaid to show within-ME tracks the
    trajectory; plus the per-trial mean ME across trials (the across signal =
    actual inter-trial variability).

(2) speed_survival_across_within: the acid test made per-cell — Speed's unique
    cv-bps with no ME vs +ME-across (stays on the diagonal = survives) and vs
    +ME-within (falls below = eaten), plus a population deviance waterfall
    (intercept -> +stim -> +ME-across -> +ME-within).

Outputs in current_plus_ME/diagnostics/.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME" / "diagnostics"
CSV = OUT / "variance_partition_me_split.csv"
EX_PROBE = "CAA-1123243_rec1"


def _example_trial_df():
    """Binned df for one well-behaved cohort cluster on the example probe."""
    cfg = drv.make_config(None)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    cohort = drv.load_cohort(EX_PROBE)
    pdata = load_probe_data(drv.FORMATTED_DIR / f"{EX_PROBE}.mat", config=cfg,
                            stimulus_lookup=lookup, visp_only=True)
    for cl in pdata.clusters:
        if cl.cluster_id in cohort:
            return bin_cluster(pdata, cl)
    raise RuntimeError("no cohort cluster found")


def plot_decomposition_illustration() -> None:
    df = _example_trial_df()
    motion = df[(df["condition"] != "stationary") & np.isfinite(df["me_face_raw"])]
    # per-trial mean ME (the across component) over motion bins
    tmean = motion.groupby("trial_id")["me_face_raw"].mean()
    # example trial must come from a vestibular condition (stage moves -> speed
    # varies) so the |velocity| overlay is meaningful; V trials have speed==0.
    vest = motion[motion["condition"].isin(["T_Vstatic", "VT"])]
    srange = vest.groupby("trial_id")["speed"].agg(lambda s: s.max() - s.min())
    counts = vest.groupby("trial_id").size()
    cand = srange[(counts.reindex(srange.index) >= 8) & (srange > 0)]
    ex_tid = cand.idxmax()                      # clearest velocity trajectory
    ex = vest[vest["trial_id"] == ex_tid].sort_values("time_since_onset")
    t = ex["time_since_onset"].to_numpy()
    me = ex["me_face_raw"].to_numpy()
    spd = ex["speed"].to_numpy()
    mu = float(me.mean())

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.6),
                                   gridspec_kw={"width_ratios": [1.05, 1]})

    # Panel L: one trial — ME(t), trial-mean (across), within residual, speed.
    axL.plot(t, me, color="tab:purple", lw=1.8, marker="o", ms=3, label="ME(t) raw")
    axL.axhline(mu, color="tab:purple", ls="--", lw=1.6,
                label=f"trial mean = ME-across ({mu:.2f})")
    axL.fill_between(t, mu, me, color="tab:pink", alpha=0.35,
                     label="deviation = ME-within")
    axL.set_xlabel("time since motion onset (s)")
    axL.set_ylabel("face motion energy", color="tab:purple")
    axL.tick_params(axis="y", labelcolor="tab:purple")
    ax2 = axL.twinx()
    ax2.plot(t, spd, color="0.35", lw=1.6, ls="-", label="|velocity|")
    ax2.set_ylabel("stimulus speed", color="0.35")
    ax2.tick_params(axis="y", labelcolor="0.35")
    axL.set_title(f"One trial, {EX_PROBE.split('_')[0]} cl {int(ex_tid)} "
                  f"profile {int(ex['profile_id'].iloc[0])}\n"
                  "ME = (trial mean) + (within-trial deviation)",
                  fontsize=10, fontweight="bold")
    h1, l1 = axL.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    axL.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="upper right")

    # Panel R: per-trial mean ME across all trials = the inter-trial signal.
    prof = motion.groupby("trial_id")["profile_id"].first()
    order = np.arange(len(tmean))
    colors = ["tab:blue" if prof.loc[tid] == 1 else "tab:orange" for tid in tmean.index]
    axR.scatter(order, tmean.to_numpy(), c=colors, s=22, alpha=0.8, edgecolor="none")
    axR.axhline(tmean.mean(), color="0.5", ls="--", lw=1,
                label=f"grand mean ({tmean.mean():.2f})")
    axR.set_xlabel("trial (in order)")
    axR.set_ylabel("trial-mean ME (ME-across)")
    axR.set_title("ME-across = each trial's mean movement\n"
                  "(this trial-to-trial spread IS the inter-trial variability)",
                  fontsize=10, fontweight="bold")
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", ls="", color="tab:blue", label="profile 1"),
               Line2D([0], [0], marker="o", ls="", color="tab:orange", label="profile 2"),
               Line2D([0], [0], ls="--", color="0.5", label="grand mean")]
    axR.legend(handles=handles, fontsize=8, frameon=False)

    fig.suptitle("What the ME across-/within-trial split is (example data)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = OUT / "me_decomposition_illustration"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png  (example trial {int(ex_tid)})")


def plot_speed_survival() -> None:
    df = pd.read_csv(CSV)
    well = df["total_both"] > 0.005
    d = df[well]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6),
                             gridspec_kw={"width_ratios": [1, 1, 1.05]})

    def scatter(ax, ycol, title, color):
        x = d["unique_Speed_noME"].to_numpy(); y = d[ycol].to_numpy()
        ax.scatter(x, y, s=20, color=color, alpha=0.6, edgecolor="none")
        lim = [min(x.min(), y.min(), 0) * 1.05, max(x.max(), y.max()) * 1.05]
        ax.plot(lim, lim, ls="--", color="0.5", lw=1)
        ax.axhline(0, color="0.85", lw=0.6); ax.axvline(0, color="0.85", lw=0.6)
        ax.set_xlabel("unique Speed — no ME")
        ax.set_ylabel(f"unique Speed — {title}")
        n_below = int((y < x - 0.005).sum())
        ax.set_title(f"Speed vs {title}\nmedian {np.median(x):+.3f} → "
                     f"{np.median(y):+.3f}  ({n_below}/{len(d)} cells fall)",
                     fontsize=10, fontweight="bold")

    scatter(axes[0], "unique_Speed_across", "+ME-across", "tab:purple")
    scatter(axes[1], "unique_Speed_within", "+ME-within", "tab:pink")

    # Panel 3: cumulative deviance build-up over intercept, using PER-CELL
    # increment medians (difference-of-medians would overstate; pairing is
    # the correct measure). cv_across reconstructed = cv_both - unique_ME_within.
    cv_across = d["cv_both"] - d["unique_ME_within"]
    inc_stim = (d["cv_stim"] - d["cv_intercept"]).median()
    inc_across = (cv_across - d["cv_stim"]).median()      # add across to stim
    inc_within = (d["cv_both"] - cv_across).median()       # add within after across
    cum = [0.0, inc_stim, inc_stim + inc_across, inc_stim + inc_across + inc_within]
    labels = ["intercept", "+stim", "+ME-across", "+ME-within"]
    colors = ["0.7", "tab:green", "tab:purple", "tab:pink"]
    incs = [None, inc_stim, inc_across, inc_within]
    ax = axes[2]
    ax.bar(range(4), cum, color=colors, edgecolor="white")
    for i, v in enumerate(cum):
        ax.text(i, v + 0.002, f"{v:.3f}", ha="center", va="bottom", fontsize=9)
        if incs[i] is not None:
            ax.annotate(f"{incs[i]:+.3f}", (i, v), textcoords="offset points",
                        xytext=(0, 16), ha="center", fontsize=8, color="0.3")
    ax.set_xticks(range(4)); ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("cumulative cv-bps gain over intercept")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_title("Cumulative deviance build-up\n(median per-cell sequential gain)",
                 fontsize=10, fontweight="bold")

    fig.suptitle("Acid test by ME component — across leaves Speed on the diagonal, "
                 "within pulls it down", fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = OUT / "speed_survival_across_within"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    plot_speed_survival()
    plot_decomposition_illustration()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
