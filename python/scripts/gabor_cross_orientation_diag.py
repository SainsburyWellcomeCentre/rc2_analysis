"""Cross-orientation diagnosis for the smallest-TF goggle clouds.

Compares the four orientation tokens (−45/0/45/90°) of the smallest-TF group
(sf00p032, VX0p095) at one RF, to:
  (A) get a better SF estimate by pooling across orientations (SF should be
      orientation-independent → pooling tightens it);
  (B) diagnose the OR convention — token θ vs measured orientation (full-frame,
      where averaging is maximal), against the candidate mappings.
  (C) show the local (windowed) orientation fluctuation per cloud.

deg/px = 0.279 (wisecoco: 35.3mm/400px @12mm → FOV 111.6°; square=5° confirmed
in rc2_visual_stimuli/generate_rc2_sparse_noise.m).
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
DPP = 111.6 / 400.0                       # 0.279 deg/px (confirmed geometry)
CX, CY, HALF = 267, 237, 26               # cl90 RF: ~10° window centred (px)
SUFFIX = "sf00p032_Bsf0p005_VX0p095_BV0p100"   # smallest-TF, highest-SF
CLOUDS = [(-45, f"theta-0p785_Btheta0p785_{SUFFIX}"),
          (0,   f"theta0p000_Btheta0p785_{SUFFIX}"),
          (45,  f"theta0p785_Btheta0p785_{SUFFIX}"),
          (90,  f"theta1p571_Btheta0p785_{SUFFIX}")]

# SF band from the token (sf0 ± 3 B_sf), in cpd
sf0 = 0.032 / DPP; bsf = 0.005 / DPP
sf_lo, sf_hi = max(0.01, sf0 - 3 * bsf), sf0 + 3 * bsf
bank, ors, sfs = W.build_gabor_bank(2 * HALF, 2 * HALF, DPP, sf_cpd=(sf_lo, sf_hi), n_sf=30)
GE = np.stack([b[0] for b in bank]); GO = np.stack([b[1] for b in bank])
n_or, n_sf = len(ors), len(sfs)


def fullframe_orientation(fs):
    """Energy-weighted modulation orientation (axial) over the full frame."""
    acc = None
    for f in fs:
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        F = np.fft.fftshift(np.fft.fft2(im - im.mean())); P = np.abs(F) ** 2
        acc = P if acc is None else acc + P
    ny, nx = acc.shape
    yy, xx = np.indices(acc.shape); fx = xx - nx // 2; fy = yy - ny // 2
    r = np.sqrt(fx ** 2 + fy ** 2); m = (r > 3) & (r < 80)
    z = np.sum(acc[m] * np.exp(1j * 2 * np.arctan2(fy[m], fx[m]))) / np.sum(acc[m])
    return np.degrees(np.angle(z) / 2) % 180


def windowed(fs):
    """Per-frame Gabor SF (cpd) and raw modulation OR at the cl90 RF."""
    sf_l, or_l = [], []
    for f in fs:
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        p = im[CY - HALF:CY + HALF, CX - HALF:CX + HALF]; pdc = p - p.mean()
        re = np.einsum("kij,ij->k", GE, pdc); iv = np.einsum("kij,ij->k", GO, pdc)
        en = (re * re + iv * iv).reshape(n_or, n_sf); i, j = np.unravel_index(en.argmax(), en.shape)
        sf_l.append(sfs[j]); or_l.append(np.degrees(ors[i]))
    return np.array(sf_l), np.array(or_l)


def main():
    data = {}
    for tok, name in CLOUDS:
        fs = sorted(glob.glob(os.path.join(GOG, name, "*.png")))[::16]   # ~80 frames
        ff = fullframe_orientation(fs)
        sf_w, or_w = windowed(fs)
        data[tok] = dict(ff=ff, sf=sf_w, orr=or_w, n=len(fs))
        print(f"θ={tok:>3}°: full-frame OR={ff:6.1f}°  windowed SF med={np.median(sf_w):.4f}  "
              f"local-OR circmean={np.degrees(np.angle(np.mean(np.exp(1j*2*np.radians(or_w))))/2)%180:6.1f}°")

    all_sf = np.concatenate([data[t]["sf"] for t, _ in CLOUDS])
    pooled = np.median(all_sf)
    print(f"\nPOOLED SF across 4 orientations: median={pooled:.4f} cpd  (token {sf0:.3f})  n={len(all_sf)}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    # A: SF per orientation + pooled
    toks = [t for t, _ in CLOUDS]
    ax[0].boxplot([data[t]["sf"] for t in toks], positions=range(len(toks)), widths=0.6,
                  showfliers=False)
    ax[0].set_xticks(range(len(toks))); ax[0].set_xticklabels([f"{t}°" for t in toks])
    ax[0].axhline(sf0, ls="--", color="k", label=f"token {sf0:.3f}")
    ax[0].axhline(pooled, ls="-", color="C3", label=f"pooled median {pooled:.3f}")
    ax[0].set(xlabel="cloud orientation token", ylabel="windowed SF (cpd)",
              title=f"(A) SF across orientations\npooled is the better estimate (n={len(all_sf)})")
    ax[0].legend(fontsize=8)

    # B: token vs full-frame measured orientation + candidate conventions
    x = np.array(toks) % 180; y = np.array([data[t]["ff"] for t in toks])
    resid = np.degrees(np.abs(np.angle(np.exp(1j * 2 * np.radians(y - (180 - x)))) / 2))
    g = np.linspace(0, 180, 200)
    ax[1].plot(g, g, "C7--", lw=1, label="y=θ")
    ax[1].plot(g, (g + 90) % 180, "C0:", lw=1, label="y=θ+90")
    ax[1].plot(g, (180 - g) % 180, "C3-", lw=2.2, label="y=180−θ (= −θ)  ← FIT")
    ax[1].plot(x, y, "ks", ms=9)
    for xi, yi, t in zip(x, y, toks):
        ax[1].annotate(f"{t}°", (xi, yi), textcoords="offset points", xytext=(6, 4), fontsize=8)
    ax[1].set(xlabel="token θ (deg, mod 180)", ylabel="full-frame measured OR (deg)",
              title=f"(B) OR convention = REFLECTION (y-flip)\nstim OR=(180−measured); fit resid {resid.mean():.1f}°",
              xlim=(-5, 185), ylim=(-5, 185)); ax[1].legend(fontsize=8)
    print(f"OR convention fit (measured vs 180−token): mean residual {resid.mean():.2f}°")

    # C: local windowed-OR distributions per orientation
    for t in toks:
        ax[2].hist(data[t]["orr"], bins=24, range=(0, 180), histtype="step", lw=2, label=f"{t}°")
    ax[2].set(xlabel="windowed local OR (raw modulation, deg)", ylabel="count",
              title="(C) local OR fluctuation per cloud (cl90 RF)"); ax[2].legend(fontsize=8)

    fig.suptitle(f"Cross-orientation diagnosis — smallest-TF sf032 group, cl90 RF @({CX},{CY}), "
                 f"{DPP:.3f}°/px, SF band [{sf_lo:.3f},{sf_hi:.3f}]", y=1.03)
    fig.tight_layout()
    out = os.path.join(GOG, "_figs", "cross_orientation_diag.png")
    fig.savefig(out, dpi=140); print(f"[saved] {out}")


if __name__ == "__main__":
    main()
