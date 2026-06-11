"""Overview of all extracted SF(t)/OR(t) regressors for one RF across the clouds.

Reads `_extract/<label>_allclouds_regressors.csv` (cloud, frame, sf_cpd, or_deg) and
draws, in one figure:
  row 1 — overlaid time-series traces (SF coloured by SF token; OR as offset-from-token,
          coloured by θ token) so the per-frame fluctuation is visible;
  row 2 — heatmaps, one row per cloud, sorted by token (SF: viridis; OR: twilight cyclic)
          so the extraction's separation into SF levels / orientations reads at a glance.
"""
from __future__ import annotations
import argparse, csv, os, re
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GOG = os.path.expanduser("~/local_data/goggle_clouds")
DPP = 111.6 / 400.0


def tokens(cloud):
    sf = float(re.search(r"_sf(\d+p\d+)_", cloud).group(1).replace("p", ".")) / DPP
    th = np.degrees(float(re.search(r"theta(-?\d+p\d+)", cloud).group(1).replace("p", "."))) % 180
    return round(sf, 3), round(th)


def off(a, tok):
    return ((np.asarray(a, float) - tok + 90) % 180) - 90


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="cl90")
    ap.add_argument("--out", default=os.path.join(GOG, "_figs", "cl90_allclouds_regressors_overview.png"))
    a = ap.parse_args()
    src = os.path.join(GOG, "_extract", f"{a.label}_allclouds_regressors.csv")
    sf, orr = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(src)):
        sf[r["cloud"]].append(float(r["sf_cpd"])); orr[r["cloud"]].append(float(r["or_deg"]))
    clouds = list(sf)
    SF = {c: np.array(sf[c]) for c in clouds}; OR = {c: np.array(orr[c]) for c in clouds}
    tok = {c: tokens(c) for c in clouds}
    sf_levels = sorted({t[0] for t in tok.values()}); th_levels = sorted({t[1] for t in tok.values()})
    sfcol = {lv: plt.cm.viridis(i / (len(sf_levels) - 1)) for i, lv in enumerate(sf_levels)}
    thcol = {lv: plt.cm.hsv(i / len(th_levels)) for i, lv in enumerate(th_levels)}

    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 2, wspace=0.2, hspace=0.28)
    # row 1 left: SF traces
    aS = fig.add_subplot(gs[0, 0])
    for c in clouds:
        aS.plot(SF[c], color=sfcol[tok[c][0]], lw=0.4, alpha=0.5)
    for lv in sf_levels:
        aS.axhline(lv, color=sfcol[lv], ls="--", lw=1.2, label=f"token {lv:.3f}")
    aS.set(xlabel="frame", ylabel="SF (cpd)", title="SF(t) — all clouds (colour = SF token)")
    aS.legend(fontsize=7, loc="upper right")
    # row 1 right: OR offset, split into the four θ groups; each group = a pair sharing the
    # OR (offset-from-token) y-axis — [marginal: counts] | [time series], centred at 0 = token
    inner = gs[0, 1].subgridspec(2, 2, hspace=0.5, wspace=0.42)
    nfr = max(len(OR[c]) for c in clouds)
    for gi, th in enumerate(th_levels):
        cell = inner[gi // 2, gi % 2].subgridspec(1, 2, width_ratios=[1, 3], wspace=0.05)
        aH = fig.add_subplot(cell[0, 0]); aT = fig.add_subplot(cell[0, 1], sharey=aH)
        grp = [c for c in clouds if tok[c][1] == th]
        for c in grp:
            aT.plot(off(OR[c], th), color="0.35", lw=0.4, alpha=0.45)
        pool = np.concatenate([off(OR[c], th) for c in grp])
        aH.hist(pool, bins=30, range=(-90, 90), orientation="horizontal", color="0.55")
        cmu = np.degrees(np.angle(np.mean(np.exp(1j * 2 * np.radians(pool)))) / 2)  # circular mean offset
        for ax in (aH, aT):
            ax.axhline(0, color="r", ls="--", lw=0.9)        # token = supposed mean
            ax.axhline(cmu, color="C0", ls="-", lw=0.9)      # observed CIRCULAR mean
            # ±90 are the SAME orientation (⊥ token): mark the wrap seam
            ax.axhline(90, color="k", lw=0.6, alpha=0.5); ax.axhline(-90, color="k", lw=0.6, alpha=0.5)
            ax.set_ylim(-90, 90); ax.set_yticks([-90, -45, 0, 45, 90]); ax.tick_params(labelsize=6)
        aH.invert_xaxis()                                    # counts grow left; OR axis abuts the time plot
        aH.set_ylabel("OR − token (°)", fontsize=8, labelpad=1)
        aH.set_yticklabels([-90, -45, 0, 45, 90], fontsize=6.5)
        aT.set_xlim(0, nfr); aT.set_yticklabels([])
        aT.set_title(f"θ = {th}°  (n={len(grp)})   circ-mean {cmu:+.1f}°", fontsize=7.5)
        if gi // 2 == 1:
            aH.set_xlabel("count", fontsize=6); aT.set_xlabel("frame", fontsize=7)
        else:
            aH.set_xticklabels([]); aT.set_xticklabels([])
    fig.text(0.74, 0.965, "OR offset from token — by θ group  (left: marginal counts · right: time; "
             "±90 wrap = ⊥token; red=token, blue=circ-mean)", ha="center", fontsize=7.5)
    # row 2: heatmaps (rows = clouds, sorted by token)
    hS = fig.add_subplot(gs[1, 0]); hO = fig.add_subplot(gs[1, 1])
    ord_sf = sorted(clouds, key=lambda c: (tok[c][0], tok[c][1]))
    ord_or = sorted(clouds, key=lambda c: (tok[c][1], tok[c][0]))
    MS = np.array([SF[c] for c in ord_sf]); MO = np.array([OR[c] for c in ord_or])
    im1 = hS.imshow(MS, aspect="auto", cmap="viridis", interpolation="nearest")
    hS.set(xlabel="frame", ylabel="cloud (sorted by SF then θ)", title="SF(t) heatmap")
    plt.colorbar(im1, ax=hS, label="SF (cpd)")
    # mark SF-level group boundaries
    for i in range(1, len(ord_sf)):
        if tok[ord_sf[i]][0] != tok[ord_sf[i - 1]][0]:
            hS.axhline(i - 0.5, color="w", lw=1.2)
    im2 = hO.imshow(MO, aspect="auto", cmap="twilight", vmin=0, vmax=180, interpolation="nearest")
    hO.set(xlabel="frame", ylabel="cloud (sorted by θ then SF)", title="OR(t) heatmap (cyclic)")
    plt.colorbar(im2, ax=hO, label="OR (deg)")
    for i in range(1, len(ord_or)):
        if tok[ord_or[i]][1] != tok[ord_or[i - 1]][1]:
            hO.axhline(i - 0.5, color="w", lw=1.2)
    fig.suptitle(f"{a.label} — extracted SF/OR regressors across {len(clouds)} clouds "
                 f"(adaptive 25°/10° windows, 0.279°/px)", fontsize=12, y=1.01)
    fig.tight_layout(); fig.savefig(a.out, dpi=130, bbox_inches="tight")
    print(f"[saved] {a.out}  ({len(clouds)} clouds)")


if __name__ == "__main__":
    main()
