"""Extraction-quality assessment across several real goggle RFs.

Runs the Gabor extractor at a spread of real RFs (different clusters, positions,
sizes; ON & OFF) on the smallest-TF sf032 cloud, and shows how quality varies.
Window per RF = its MIN-axis size (robust to the CSV elongation artefacts, e.g.
cl90's bogus 25°). Quality proxies: recovered SF vs token, Gabor energy-peak
concentration (sharp peak = reliable), and cycles-in-window = SF·size.

deg/px = 0.279 (wisecoco). OR uses the locked flip stim=(180−mod).
"""
from __future__ import annotations
import csv, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import windowed_cloud_stat_recovery as W
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

GOG = os.path.expanduser("~/local_data/goggle_clouds")
RFDIR = os.path.join(GOG, "_rfs")
DPP = 111.6 / 400.0
CLOUD = "theta0p000_Btheta0p785_sf00p032_Bsf0p005_VX0p095_BV0p100"  # token SF 0.115, OR 0
SF0 = 0.032 / DPP; BSF = 0.005 / DPP
SF_LO, SF_HI = max(0.01, SF0 - 3 * BSF), SF0 + 3 * BSF


def load_rfs():
    out = []
    for fp in glob.glob(os.path.join(RFDIR, "*_rf_metrics.csv")):
        probe = os.path.basename(fp).split("_")[0][-3:]
        for r in csv.DictReader(open(fp)):
            try:
                cx = float(r["centroid_azimuth_pixels"]); cy = float(r["centroid_elevation_pixels"])
                az = float(r["size_azimuth_deg"]); el = float(r["size_elevation_deg"])
            except (KeyError, ValueError):
                continue
            ms = min(az, el)
            if 50 < cx < 350 and 50 < cy < 350 and 5 <= ms <= 20:
                out.append(dict(probe=probe, cl=r["cluster_id"], pol=r["rf_type"],
                                cx=cx, cy=cy, ms=ms))
    return out


def select(rfs, n=6):
    rfs = sorted(rfs, key=lambda d: d["ms"])
    idx = np.linspace(0, len(rfs) - 1, n).astype(int)
    return [rfs[i] for i in idx]


def run_rf(fs, rf, bank_cache):
    cx, cy = int(round(rf["cx"])), int(round(rf["cy"]))
    h = max(8, int(round(rf["ms"] / DPP / 2)))            # half = min-axis radius (px)
    key = h
    if key not in bank_cache:
        bk, ors, sfs = W.build_gabor_bank(2 * h, 2 * h, DPP, sf_cpd=(SF_LO, SF_HI), n_sf=24)
        bank_cache[key] = (np.stack([b[0] for b in bk]), np.stack([b[1] for b in bk]), ors, sfs)
    GE, GO, ors, sfs = bank_cache[key]
    acc = None; sf_l, or_l = [], []
    patch_mean = None
    for f in fs:
        import matplotlib.image as mpi
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        p = im[cy - h:cy + h, cx - h:cx + h]
        patch_mean = p.astype(float) if patch_mean is None else patch_mean + p
        pdc = p - p.mean()
        re = np.einsum("kij,ij->k", GE, pdc); iv = np.einsum("kij,ij->k", GO, pdc)
        en = (re * re + iv * iv).reshape(len(ors), len(sfs))
        acc = en.astype(float) if acc is None else acc + en
        i, j = np.unravel_index(en.argmax(), en.shape)
        sf_l.append(sfs[j]); or_l.append((180 - np.degrees(ors[i])) % 180)
    acc /= len(fs); patch_mean /= len(fs)
    conc = float(acc.max() / acc.mean())                  # peak concentration (quality)
    sf_med = float(np.median(sf_l))
    or_med = float(np.degrees(np.angle(np.mean(np.exp(1j * 2 * np.radians(or_l)))) / 2) % 180)
    return dict(h=h, acc=acc, ors=ors, sfs=sfs, patch=patch_mean, conc=conc,
                sf=sf_med, orr=or_med, cyc=SF0 * 2 * h * DPP)


def main():
    fs = sorted(glob.glob(os.path.join(GOG, CLOUD, "*.png")))[::32]   # ~40 frames
    sel = select(load_rfs(), 6)
    print(f"cloud token SF={SF0:.3f} cpd, OR=0°; {len(fs)} frames")
    bank_cache = {}
    res = []
    for rf in sel:
        r = run_rf(fs, rf, bank_cache); r.update(rf); res.append(r)
        print(f"cl{rf['cl']:>3} {rf['pol']:>5} minsz={rf['ms']:>4.0f}° "
              f"({2*r['cyc']:.1f} cyc): SF={r['sf']:.4f}  OR={r['orr']:5.1f}°  conc={r['conc']:.1f}")

    fig, ax = plt.subplots(3, 6, figsize=(19, 9),
                           gridspec_kw=dict(height_ratios=[1, 1, 1.1]))
    for c, r in enumerate(res):
        a = ax[0, c]
        a.imshow(r["patch"], cmap="gray")
        a.add_patch(Circle((r["h"], r["h"]), r["h"], fill=False, ec="red", lw=1.5))
        a.set(xticks=[], yticks=[],
              title=f"{r['probe']} cl{r['cl']} {r['pol']}\n{r['ms']:.0f}° ({2*r['cyc']:.1f} cyc)")
        a.title.set_fontsize(8)
        b = ax[1, c]
        b.imshow(r["acc"].T, origin="lower", aspect="auto", cmap="magma")
        pk = np.unravel_index(r["acc"].argmax(), r["acc"].shape); b.plot(pk[0], pk[1], "co", ms=7, mfc="none", mew=2)
        b.set(xticks=[], yticks=[], title=f"SF={r['sf']:.3f} OR={r['orr']:.0f}°\nconc={r['conc']:.1f}")
        b.title.set_fontsize(8)
    # bottom: summary scatters spanning
    for k in range(6):
        ax[2, k].remove()
    gs = ax[2, 0].get_gridspec()
    axA = fig.add_subplot(gs[2, 0:3]); axB = fig.add_subplot(gs[2, 3:6])
    cyc = np.array([2 * r["cyc"] for r in res]); sfm = np.array([r["sf"] for r in res])
    conc = np.array([r["conc"] for r in res]); szs = np.array([r["ms"] for r in res])
    axA.scatter(szs, sfm, c=cyc, cmap="viridis", s=90, edgecolor="k")
    axA.axhline(SF0, ls="--", color="k", label=f"token {SF0:.3f}")
    axA.set(xlabel="RF min-axis size (deg)", ylabel="recovered SF (cpd)",
            title="recovered SF vs RF size (colour=cycles)"); axA.legend(fontsize=8)
    axB.scatter(cyc, conc, c=szs, cmap="plasma", s=90, edgecolor="k")
    axB.axvline(1.0, color="r", lw=1, ls=":"); axB.set(xlabel="cycles in window (SF·size)",
            ylabel="energy-peak concentration", title="extraction sharpness vs cycles-in-window")
    fig.suptitle(f"Extraction quality across real RFs — cloud {CLOUD} (token SF {SF0:.3f}, OR 0°), "
                 f"window=RF min-axis, {DPP:.3f}°/px", y=1.0, fontsize=10)
    fig.tight_layout()
    out = os.path.join(GOG, "_figs", "rf_quality_diag.png")
    fig.savefig(out, dpi=130); print(f"[saved] {out}")


if __name__ == "__main__":
    main()
