"""SF×OR validation grid across the smallest-TF goggle clouds.

For the gain-1/30 batch (smallest TF), 3 SF tokens × 4 orientation tokens = 12
clouds, measure full-frame recovered SF and OR vs their nominal tokens — does the
extractor recover the right SF and OR across the whole parameter combination?

Full-frame (not RF-windowed) so every SF resolves (the small RF is sub-cycle for
the low SF). deg/px = 0.279 (wisecoco). OR uses the locked flip stim=(180−mod).
"""
from __future__ import annotations
import glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpi

GOG = os.path.expanduser("~/local_data/motion_clouds/saved_goggles")
DPP = 111.6 / 400.0
THETAS = [(-45, "theta-0p785"), (0, "theta0p000"), (45, "theta0p785"), (90, "theta1p571")]
# smallest-TF group (SF×VX const): per SF token its VX, + physical cpd (token/dpp)
SF_TOK = [("sf00p008", "VX0p382", 0.008 / DPP),
          ("sf00p016", "VX0p191", 0.016 / DPP),
          ("sf00p032", "VX0p095", 0.032 / DPP)]


def fullframe_sf_or(fs):
    """Full-frame radial-spectrum peak SF (cpd) and energy-weighted OR (flipped)."""
    acc = None
    for f in fs:
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        P = np.abs(np.fft.fftshift(np.fft.fft2(im - im.mean()))) ** 2
        acc = P if acc is None else acc + P
    ny, nx = acc.shape; cy, cx = ny // 2, nx // 2
    yy, xx = np.indices(acc.shape); fy = yy - cy; fx = xx - cx
    r = np.sqrt(fx ** 2 + fy ** 2)
    freqs = (r / nx) / DPP                       # cpd (square frame)
    nyq = 0.5 / DPP
    pw, edges = np.histogram(freqs.ravel(), bins=120, range=(0, nyq), weights=acc.ravel())
    cnt, _ = np.histogram(freqs.ravel(), bins=120, range=(0, nyq))
    prof = pw / np.maximum(cnt, 1); prof[0] = 0
    centers = 0.5 * (edges[:-1] + edges[1:])
    sf_pk = centers[np.argmax(prof)]
    m = (r > 3) & (r < 90)
    z = np.sum(acc[m] * np.exp(1j * 2 * np.arctan2(fy[m], fx[m]))) / np.sum(acc[m])
    or_meas = np.degrees(np.angle(z) / 2) % 180
    or_stim = (180 - or_meas) % 180              # locked flip
    return sf_pk, or_stim


def main():
    rows = []   # (sf_token_cpd, or_token, sf_meas, or_meas_stim)
    for sftok, vx, sf_cpd in SF_TOK:
        for orval, thetatok in THETAS:
            name = f"{thetatok}_Btheta0p785_{sftok}_Bsf0p005_{vx}_BV0p100"
            fs = sorted(glob.glob(os.path.join(GOG, name, "*.png")))[::24]   # ~55 frames
            if not fs:
                print(f"MISSING {name}"); continue
            sfm, orm = fullframe_sf_or(fs)
            rows.append((sf_cpd, orval % 180, sfm, orm))
            print(f"SF {sf_cpd:.3f} OR {orval:>3}°: measured SF={sfm:.4f}  OR={orm:6.1f}°")
    rows = np.array(rows)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    # SF recovery (x=token, y=measured), coloured by OR
    sc = ax[0].scatter(rows[:, 0], rows[:, 2], c=rows[:, 1], cmap="hsv", vmin=0, vmax=180,
                       s=80, edgecolor="k", zorder=3)
    lim = [0, max(rows[:, 0].max(), rows[:, 2].max()) * 1.1]
    ax[0].plot(lim, lim, "k--", lw=1, label="identity")
    for sf_cpd in np.unique(rows[:, 0]):
        ax[0].axvline(sf_cpd, color="0.85", lw=0.8, zorder=0)
    ax[0].set(xlabel="token SF (cpd)", ylabel="measured full-frame SF (cpd)", xlim=lim, ylim=lim,
              title="(A) SF recovery across SF×OR\n(colour = OR token; should sit on identity)")
    ax[0].legend(fontsize=8); plt.colorbar(sc, ax=ax[0], label="OR token (deg)")

    # OR recovery (x=token, y=measured stim OR), coloured by SF
    sc2 = ax[1].scatter(rows[:, 1], rows[:, 3], c=rows[:, 0], cmap="viridis",
                        s=80, edgecolor="k", zorder=3)
    ax[1].plot([0, 180], [0, 180], "k--", lw=1, label="identity")
    resid = np.degrees(np.abs(np.angle(np.exp(1j * 2 * np.radians(rows[:, 3] - rows[:, 1]))) / 2))
    ax[1].set(xlabel="token OR (deg, mod 180)", ylabel="measured stim OR = 180−mod (deg)",
              xlim=(-5, 185), ylim=(-5, 185),
              title=f"(B) OR recovery across SF×OR\nmean residual {resid.mean():.1f}° (locked flip)")
    ax[1].legend(fontsize=8); plt.colorbar(sc2, ax=ax[1], label="token SF (cpd)")

    fig.suptitle("SF×OR validation grid — smallest-TF group, full-frame, 0.279°/px", y=1.02)
    fig.tight_layout()
    out = os.path.join(GOG, "_figs", "sf_or_grid_diag.png")
    fig.savefig(out, dpi=140); print(f"\n[saved] {out}  | OR residual {resid.mean():.2f}°")


if __name__ == "__main__":
    main()
