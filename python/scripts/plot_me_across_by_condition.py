"""Per-trial mean ME (the ME-across signal) vs trial order, by condition.

Tests Laura's observation: the across-trial ME splits into two bands that track
CONDITION (platform-motion T_Vstatic/VT = high; visual-only V = low) and the
high band HABITUATES away over the session. One panel per face-camera probe.

Output: current_plus_ME/diagnostics/me_across_by_condition.{pdf,png}
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

PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1")
OUT = drv.ROOT / "figures" / "glm" / "current_plus_ME" / "diagnostics"
COLORS = {"V": "tab:blue", "T_Vstatic": "tab:red", "VT": "tab:green"}


def _trial_means(probe, cfg, lookup):
    cohort = drv.load_cohort(probe)
    pdata = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat", config=cfg,
                            stimulus_lookup=lookup, visp_only=True)
    cl = next(c for c in pdata.clusters if c.cluster_id in cohort)
    df = bin_cluster(pdata, cl)
    m = df[(df["condition"] != "stationary") & np.isfinite(df["me_face_raw"])]
    g = (m.groupby("trial_id")
         .agg(meanME=("me_face_raw", "mean"), cond=("condition", "first"))
         .reset_index().sort_values("trial_id").reset_index(drop=True))
    g["order"] = np.arange(len(g))
    return g


def main() -> int:
    cfg = drv.make_config(None)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), sharey=False)
    for ax, probe in zip(axes, PROBES):
        g = _trial_means(probe, cfg, lookup)
        for cond, c in COLORS.items():
            sub = g[g["cond"] == cond]
            ax.scatter(sub["order"], sub["meanME"], s=16, color=c, alpha=0.6,
                       edgecolor="none", label=cond)
            # rolling median (in order) per condition to show habituation
            if len(sub) > 6:
                r = sub.sort_values("order")
                roll = r["meanME"].rolling(9, center=True, min_periods=3).median()
                ax.plot(r["order"], roll, color=c, lw=2, alpha=0.9)
        ax.set_xlabel("trial (in order)")
        ax.set_title(probe.split("_")[0], fontsize=10, fontweight="bold")
        ax.axhline(g["meanME"].median(), color="0.7", ls="--", lw=0.8)
    axes[0].set_ylabel("trial-mean ME (ME-across)")
    axes[0].legend(fontsize=8, frameon=False, title="condition")
    fig.suptitle("ME-across by condition: platform-motion (T_Vstatic/VT) evokes "
                 "high face-ME that habituates; V stays low",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "me_across_by_condition"
    fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out}.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
