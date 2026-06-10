"""Per-trial within-trial STD of ME vs trial order, by condition.

Companion to plot_me_across_by_condition.py: instead of each trial's MEAN
(the across-trial signal), plot each trial's STD = how much ME fluctuates
WITHIN the trial. This is the magnitude of the ME-within component — the part
that overlaps Speed. Tests whether the within-trial fluctuation also splits by
condition (platform-motion vs visual) and whether it habituates.

Output: current_plus_ME/diagnostics/me_std_by_condition.{pdf,png}
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.io import load_probe_data
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME" / "diagnostics"
COLORS = {"V": "tab:blue", "T_Vstatic": "tab:red", "VT": "tab:green"}


def _trial_stds(probe, cfg, lookup):
    cohort = drv.load_cohort(probe)
    pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                            stimulus_lookup=lookup, visp_only=True)
    cl = next(c for c in pdata.clusters if c.cluster_id in cohort)
    df = bin_cluster(pdata, cl)
    m = df[(df["condition"] != "stationary") & np.isfinite(df["me_face_raw"])]
    g = (m.groupby("trial_id")
         .agg(stdME=("me_face_raw", lambda s: s.std(ddof=0)),
              n=("me_face_raw", "size"), cond=("condition", "first"))
         .reset_index().sort_values("trial_id").reset_index(drop=True))
    g = g[g["n"] >= 4]
    g["order"] = np.arange(len(g))
    return g


def main() -> int:
    cfg = drv.make_config(None)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, probe in zip(axes, PROBES):
        g = _trial_stds(probe, cfg, lookup)
        for cond, c in COLORS.items():
            sub = g[g["cond"] == cond]
            ax.scatter(sub["order"], sub["stdME"], s=16, color=c, alpha=0.6,
                       edgecolor="none", label=cond)
            if len(sub) > 6:
                r = sub.sort_values("order")
                roll = r["stdME"].rolling(9, center=True, min_periods=3).median()
                ax.plot(r["order"], roll, color=c, lw=2, alpha=0.9)
        ax.set_xlabel("trial (in order)")
        ax.set_title(probe.split("_")[0], fontsize=10, fontweight="bold")
        ax.axhline(g["stdME"].median(), color="0.7", ls="--", lw=0.8)
    axes[0].set_ylabel("within-trial STD of ME")
    axes[0].legend(fontsize=8, frameon=False, title="condition")
    fig.suptitle("Within-trial ME fluctuation (STD) by condition — the magnitude "
                 "of the ME-within signal that overlaps Speed",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "me_std_by_condition"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
