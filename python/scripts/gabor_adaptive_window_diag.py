"""Adaptive (separate) windows for SF and OR — do the marginals match expected?

SF wants a LARGE window (accurate mean, std → real B_sf); OR wants a SMALL window
(preserve the B_theta local-orientation spread). We extract per-RF temporal SF(t)
[large win] and OR(t) [small win, + Viterbi smoothing], and compare their marginals
(mean & std) to the stimulus's own SPATIAL marginals (sliding the SAME windows
across the full frame = the expected local distribution, by ergodicity).

cl90, sf032 cloud, deg/px = 0.279. OR uses locked flip stim=(180−mod).
"""
from __future__ import annotations
import glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import windowed_cloud_stat_recovery as W
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpi

GOG = os.path.expanduser("~/local_data/motion_clouds/saved_goggles")
DPP = 111.6 / 400.0
CX, CY = 267, 237
CLOUD = "theta0p000_Btheta0p785_sf00p032_Bsf0p005_VX0p095_BV0p100"
SF_WIN, OR_WIN = 25, 10                                  # deg (diameters): big for SF, small for OR
SF0 = 0.032 / DPP; BSF = 0.005 / DPP
SF_LO, SF_HI = max(0.01, SF0 - 4 * BSF), SF0 + 4 * BSF


def bank(win_deg):
    h = max(8, int(round(win_deg / DPP / 2)))
    bk, ors, sfs = W.build_gabor_bank(2 * h, 2 * h, DPP, sf_cpd=(SF_LO, SF_HI), n_sf=40)
    return h, np.stack([b[0] for b in bk]), np.stack([b[1] for b in bk]), ors, sfs


def est(patch, GE, GO, ors, sfs):
    pdc = patch - patch.mean()
    re = np.einsum("kij,ij->k", GE, pdc); iv = np.einsum("kij,ij->k", GO, pdc)
    en = (re * re + iv * iv).reshape(len(ors), len(sfs))
    i, j = np.unravel_index(en.argmax(), en.shape)
    return sfs[j], (180 - np.degrees(ors[i])) % 180, en


def circ_std(deg):
    z = np.mean(np.exp(1j * 2 * np.radians(deg)))
    return np.degrees(np.sqrt(max(-2 * np.log(abs(z)), 0))) / 2


def circ_mean(deg):
    return np.degrees(np.angle(np.mean(np.exp(1j * 2 * np.radians(deg)))) / 2) % 180


def main():
    fs = sorted(glob.glob(os.path.join(GOG, CLOUD, "*.png")))
    hS, GES, GOS, orsS, sfsS = bank(SF_WIN)
    hO, GEO, GOO, orsO, sfsO = bank(OR_WIN)

    # --- per-RF TEMPORAL marginals at cl90 ---
    sf_t, or_t = [], []
    for f in fs[::20]:                                   # ~65 frames
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        sfv, _, _ = est(im[CY - hS:CY + hS, CX - hS:CX + hS], GES, GOS, orsS, sfsS)
        _, orv, _ = est(im[CY - hO:CY + hO, CX - hO:CX + hO], GEO, GOO, orsO, sfsO)
        sf_t.append(sfv); or_t.append(orv)
    sf_t = np.array(sf_t); or_t = np.array(or_t)
    # OR Viterbi smoothing on the temporal sequence (re-extract energies for path)
    en_seq = []
    for f in fs[::20]:
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        _, _, en = est(im[CY - hO:CY + hO, CX - hO:CX + hO], GEO, GOO, orsO, sfsO)
        en_seq.append(en.sum(axis=1))
    from gabor_goggle_animation import viterbi_or
    path = viterbi_or(np.array(en_seq), np.degrees(orsO), 0.003)
    or_t_sm = (180 - np.degrees(orsO[path])) % 180

    # --- expected SPATIAL marginals: slide windows across the frame ---
    sf_sp, or_sp = [], []
    for f in fs[::130]:                                  # ~10 frames
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        for yc in range(hS, 400 - hS, hS):
            for xc in range(hS, 400 - hS, hS):
                sfv, _, _ = est(im[yc - hS:yc + hS, xc - hS:xc + hS], GES, GOS, orsS, sfsS)
                sf_sp.append(sfv)
        for yc in range(hO, 400 - hO, hO):
            for xc in range(hO, 400 - hO, hO):
                _, orv, _ = est(im[yc - hO:yc + hO, xc - hO:xc + hO], GEO, GOO, orsO, sfsO)
                or_sp.append(orv)
    sf_sp = np.array(sf_sp); or_sp = np.array(or_sp)

    print(f"SF  (win {SF_WIN}°): per-RF mean={sf_t.mean():.4f} std={sf_t.std():.4f} | "
          f"expected mean={sf_sp.mean():.4f} std={sf_sp.std():.4f} | token {SF0:.4f}, B_sf {BSF:.4f}")
    print(f"OR  (win {OR_WIN}°): per-RF raw mean={circ_mean(or_t):.1f} std={circ_std(or_t):.1f} | "
          f"smoothed mean={circ_mean(or_t_sm):.1f} std={circ_std(or_t_sm):.1f} | "
          f"expected mean={circ_mean(or_sp):.1f} std={circ_std(or_sp):.1f}")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].hist(sf_sp, bins=30, density=True, alpha=0.4, color="0.6",
               label=f"expected (spatial) μ={sf_sp.mean():.3f} σ={sf_sp.std():.3f}")
    ax[0].hist(sf_t, bins=20, density=True, histtype="step", lw=2, color="C0",
               label=f"per-RF temporal μ={sf_t.mean():.3f} σ={sf_t.std():.3f}")
    ax[0].axvline(SF0, ls="--", color="k", label=f"token {SF0:.3f}")
    ax[0].axvspan(SF0 - BSF, SF0 + BSF, color="green", alpha=0.1, label="±B_sf")
    ax[0].set(xlabel="SF (cpd)", ylabel="density", title=f"SF marginal — {SF_WIN}° window")
    ax[0].legend(fontsize=7)
    ax[1].hist(or_sp, bins=30, range=(0, 180), density=True, alpha=0.4, color="0.6",
               label=f"expected (spatial) σ={circ_std(or_sp):.0f}°")
    ax[1].hist(or_t, bins=20, range=(0, 180), density=True, histtype="step", lw=2, color="C0",
               label=f"per-RF raw σ={circ_std(or_t):.0f}°")
    ax[1].hist(or_t_sm, bins=20, range=(0, 180), density=True, histtype="step", lw=2, color="C3",
               label=f"per-RF smoothed σ={circ_std(or_t_sm):.0f}°")
    ax[1].axvline(0, ls="--", color="k"); ax[1].axvline(180, ls="--", color="k", label="token 0/180")
    ax[1].set(xlabel="OR (deg)", ylabel="density", title=f"OR marginal — {OR_WIN}° window")
    ax[1].legend(fontsize=7)
    fig.suptitle(f"Adaptive windows: SF {SF_WIN}° / OR {OR_WIN}° — do per-RF marginals match the "
                 f"stimulus spatial marginals?  ({CLOUD[:24]}…, cl90)", y=1.02, fontsize=10)
    fig.tight_layout()
    out = os.path.join(GOG, "_figs", "adaptive_window_diag.png")
    fig.savefig(out, dpi=140); print(f"[saved] {out}")


if __name__ == "__main__":
    main()
