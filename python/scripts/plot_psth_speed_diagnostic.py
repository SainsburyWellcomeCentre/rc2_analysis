"""FR-vs-time vs FR-vs-speed diagnostic, split by speed profile (VT).

Tests whether a cell's apparent Speed tuning is real Speed tuning or just the
onset/adaptation time course re-plotted against the (ramping) speed. Three
panels per cluster, each with the TWO speed profiles as separate lines:

  1. FR vs time-since-onset (PSTH)   — the response time course / adaptation
  2. speed vs time-since-onset        — shows the profiles differ in speed(t)
  3. FR vs speed (the tuning curve)   — THE ADJUDICATOR

Read: if FR-vs-speed overlaps for both profiles (despite different speed(t)),
it's genuine speed tuning; if instead FR-vs-time overlaps but FR-vs-speed
splits by profile, the "speed tuning" is really the time course.

All observed data (FR = spike_count / bin width); no model.
Output: <run>/figs/psth_speed_diagnostics/<probe>_cl<cluster>_VT.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.io import load_probe_data
from rc2_glm.pipeline import _subset_to_condition
from rc2_glm.time_binning import bin_cluster

import scripts.run_glm_split_by_condition as drv

COND = "VT"
CLUSTERS = [("CAA-1123243_rec1", 376), ("CAA-1123244_rec1", 34)]
DT = 0.1
PROFILE_COLORS = {1: "tab:blue", 2: "tab:orange"}


def _binned_mean(x, vals, edges):
    idx = np.digitize(x, edges) - 1
    centres = 0.5 * (edges[:-1] + edges[1:])
    mean = np.full(len(centres), np.nan)
    sem = np.full(len(centres), np.nan)
    for b in range(len(centres)):
        m = idx == b
        if m.sum() > 0:
            mean[b] = np.nanmean(vals[m])
            sem[b] = np.nanstd(vals[m]) / max(np.sqrt(m.sum()), 1)
    return centres, mean, sem


def plot_cluster(probe_data, cfg, probe, cluster):
    cl = next(c for c in probe_data.clusters if c.cluster_id == cluster)
    df = _subset_to_condition(bin_cluster(probe_data, cl), COND)
    mo = (df["condition"] != "stationary").to_numpy()
    t = df["time_since_onset"].to_numpy()[mo]
    sp = df["speed"].to_numpy()[mo]
    fr = df["spike_count"].to_numpy()[mo] / DT
    prof = df["profile_id"].to_numpy()[mo]

    t_edges = np.linspace(0, np.ceil(t.max()), 17)
    sp_edges = np.linspace(0, 50, 16)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for pid in sorted(np.unique(prof)):
        col = PROFILE_COLORS.get(int(pid), "0.4")
        s = prof == pid
        # 1. FR vs time
        c, m, e = _binned_mean(t[s], fr[s], t_edges)
        axes[0].plot(c, m, color=col, lw=2, label=f"profile {int(pid)}")
        axes[0].fill_between(c, m - e, m + e, color=col, alpha=0.2, linewidth=0)
        # 2. speed vs time
        c2, m2, _ = _binned_mean(t[s], sp[s], t_edges)
        axes[1].plot(c2, m2, color=col, lw=2, label=f"profile {int(pid)}")
        # 3. FR vs speed
        c3, m3, e3 = _binned_mean(sp[s], fr[s], sp_edges)
        axes[2].plot(c3, m3, color=col, lw=2, label=f"profile {int(pid)}")
        axes[2].fill_between(c3, m3 - e3, m3 + e3, color=col, alpha=0.2, linewidth=0)

    axes[0].set_xlabel("time since onset (s)"); axes[0].set_ylabel("FR (Hz)")
    axes[0].set_title("PSTH  (FR vs time)")
    axes[1].set_xlabel("time since onset (s)"); axes[1].set_ylabel("speed (cm/s)")
    axes[1].set_title("speed(t)  — profiles differ?")
    axes[2].set_xlabel("speed (cm/s)"); axes[2].set_ylabel("FR (Hz)")
    axes[2].set_title("FR vs speed  (tuning — adjudicator)")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.suptitle(
        f"{probe} cluster {cluster} · VT · is the speed tuning real or the PSTH? "
        f"(lines = the two speed profiles)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    outdir = drv.OUT_ROOT / "figs" / "psth_speed_diagnostics"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{probe}_cl{cluster}_VT.pdf"
    fig.savefig(out); plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    cfg = drv.make_config(COND)
    lookup = StimulusLookup(str(drv.MC_SEQUENCE), str(drv.MC_FOLDERS))
    by_probe: dict[str, list[int]] = {}
    for probe, cl in CLUSTERS:
        by_probe.setdefault(probe, []).append(cl)
    for probe, clusters in by_probe.items():
        probe_data = load_probe_data(drv.FORMATTED_DIR / f"{probe}.mat",
                                     config=cfg, stimulus_lookup=lookup,
                                     visp_only=True)
        for cl in clusters:
            plot_cluster(probe_data, cfg, probe, cl)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
