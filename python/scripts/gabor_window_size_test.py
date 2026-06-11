"""Does a constant window (vs RF-sized) hurt SF? And what does it do to OR?

At one RF position (cl90), extract SF and OR across window sizes for the three
SF clouds (smallest-TF). Uses a WIDE Gabor SF band (NOT the cloud's own band) so
the recovery isn't forced toward the token. deg/px = 0.279.

Reads: SF accuracy = recovered vs token; OR locality = temporal std of OR.
"""
from __future__ import annotations
import glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import windowed_cloud_stat_recovery as W
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpi

GOG = os.path.expanduser("~/local_data/goggle_clouds")
DPP = 111.6 / 400.0
CX, CY = 267, 237
WINDOWS = [10, 15, 20, 30, 40]                       # deg (diameter)
CLOUDS = [("sf00p008", "VX0p382"), ("sf00p016", "VX0p191"), ("sf00p032", "VX0p095")]


def extract(cloud, half):
    fs = sorted(glob.glob(os.path.join(GOG, cloud, "*.png")))[::26]   # ~50 frames
    bank, ors, sfs = W.build_gabor_bank(2 * half, 2 * half, DPP, sf_cpd=(0.01, 0.25), n_sf=60)
    GE = np.stack([b[0] for b in bank]); GO = np.stack([b[1] for b in bank])
    sf_l, or_l = [], []
    for f in fs:
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        p = im[CY - half:CY + half, CX - half:CX + half]; pdc = p - p.mean()
        re = np.einsum("kij,ij->k", GE, pdc); iv = np.einsum("kij,ij->k", GO, pdc)
        en = (re * re + iv * iv).reshape(len(ors), len(sfs)); i, j = np.unravel_index(en.argmax(), en.shape)
        sf_l.append(sfs[j]); or_l.append((180 - np.degrees(ors[i])) % 180)
    sf_l = np.array(sf_l); orr = np.radians(2 * np.array(or_l))
    or_std = np.degrees(np.sqrt(-2 * np.log(np.abs(np.mean(np.exp(1j * orr)))))) / 2  # circ std (deg)
    return np.median(sf_l), or_std


def main():
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for sftok, vx in CLOUDS:
        tok_cpd = float(sftok.replace("sf", "").replace("p", ".")) / DPP   # token cpd
        cloud = f"theta0p000_Btheta0p785_{sftok}_Bsf0p005_{vx}_BV0p100"
        sfs_w, ors_w = [], []
        for wd in WINDOWS:
            half = max(8, int(round(wd / DPP / 2)))
            sf_m, or_s = extract(cloud, half)
            sfs_w.append(sf_m); ors_w.append(or_s)
            print(f"{sftok} (token {tok_cpd:.3f}) win {wd:>2}°: SF={sf_m:.4f}  OR temporal std={or_s:4.1f}° "
                  f"(cyc={tok_cpd*wd:.1f})")
        col = ax[0].plot(WINDOWS, sfs_w, "-o", label=f"{sftok} (tok {tok_cpd:.3f})")[0].get_color()
        ax[0].axhline(tok_cpd, ls="--", color=col, lw=0.8)
        ax[1].plot(WINDOWS, ors_w, "-o", color=col, label=f"{sftok}")
        print()
    ax[0].axvline(20, color="0.7", lw=1, ls=":")
    ax[0].set(xlabel="window diameter (deg)", ylabel="recovered SF (cpd)",
              title="SF accuracy vs window\n(dashed = token; bigger window → closer, esp. low SF)")
    ax[0].legend(fontsize=8)
    ax[1].axvline(20, color="0.7", lw=1, ls=":")
    ax[1].set(xlabel="window diameter (deg)", ylabel="OR temporal std (deg)",
              title="OR local fluctuation vs window\n(bigger window → washes it out)")
    ax[1].legend(fontsize=8)
    fig.suptitle("Constant vs RF-sized window — SF accuracy improves with size; OR locality shrinks "
                 "(cl90, 0.279°/px, wide SF band)", y=1.02)
    fig.tight_layout()
    out = os.path.join(GOG, "_figs", "window_size_test.png")
    fig.savefig(out, dpi=140); print(f"[saved] {out}")


if __name__ == "__main__":
    main()
