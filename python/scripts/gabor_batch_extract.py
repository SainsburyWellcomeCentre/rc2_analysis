"""Batch per-RF SF(t)/OR(t) extraction with adaptive two-window recipe + validation.

Recipe (settled over the June-2026 sessions):
  SF  → large window (~25°), Gabor best-match, band = cloud's own sf0±k·B_sf.
        SF is spatially near-stationary, so a big window gives an accurate,
        marginal-matching estimate (and dodges the unreliable CSV RF sizes).
  OR  → small window (~10°), Gabor best-match → locked flip stim=(180−mod) →
        Viterbi temporal smoothing with λ AUTO-TUNED so the smoothed σ keeps
        ≥90% of the raw spread (rejects implausible jumps without flattening
        the real local-orientation drift).
  flag → Gabor energy-peak concentration (trust SF only above ~2 ≈ 1.5 cycles).

Validates by comparing per-RF TEMPORAL marginals to the stimulus SPATIAL
marginals (sliding the same windows across the frame = expected, by ergodicity),
then writes the per-(cloud,RF) time series (the GLM regressors) + a summary.

deg/px = 0.279 (wisecoco). Local data only. This is a VALIDATION batch (frames
subsampled); set STRIDE=1 for full-resolution production regressors.
"""
from __future__ import annotations
import csv, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import windowed_cloud_stat_recovery as W
from gabor_goggle_animation import viterbi_or
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpi

GOG = os.path.expanduser("~/local_data/goggle_clouds")
RFDIR = os.path.join(GOG, "_rfs"); OUT = os.path.join(GOG, "_extract"); os.makedirs(OUT, exist_ok=True)
DPP = 111.6 / 400.0
SF_WIN, OR_WIN = 25, 10
STRIDE = 10                                   # frame subsample for the validation batch
# smallest-TF group, θ=0, the 3 SF levels
CLOUDS = [("sf00p008", "VX0p382"), ("sf00p016", "VX0p191"), ("sf00p032", "VX0p095")]


def bank(win_deg, sf0, bsf):
    h = max(8, int(round(win_deg / DPP / 2)))
    lo, hi = max(0.01, sf0 - 4 * bsf), sf0 + 4 * bsf
    bk, ors, sfs = W.build_gabor_bank(2 * h, 2 * h, DPP, sf_cpd=(lo, hi), n_sf=40)
    return h, np.stack([b[0] for b in bk]), np.stack([b[1] for b in bk]), ors, sfs


def readgray(f):
    im = mpi.imread(f)
    return im[..., :3].mean(-1) if im.ndim == 3 else im


def est(patch, GE, GO):
    pdc = patch - patch.mean()
    re = np.einsum("kij,ij->k", GE, pdc); iv = np.einsum("kij,ij->k", GO, pdc)
    return (re * re + iv * iv)


def cstd(deg):
    z = np.mean(np.exp(1j * 2 * np.radians(deg)))
    return np.degrees(np.sqrt(max(-2 * np.log(abs(z)), 0))) / 2


def cmean(deg):
    return np.degrees(np.angle(np.mean(np.exp(1j * 2 * np.radians(deg)))) / 2) % 180


def auto_lambda(omarg, or_deg, raw_sf):
    """Largest λ keeping smoothed circular-σ ≥ 0.9·raw σ (preserve the real spread)."""
    raw_path = omarg.argmax(axis=1); raw_or = (180 - or_deg[raw_path]) % 180
    target = 0.9 * cstd(raw_or)
    best = 0.0
    for lam in [2e-4, 5e-4, 1e-3, 2e-3]:
        sm = (180 - or_deg[viterbi_or(omarg, or_deg, lam)]) % 180
        if cstd(sm) >= target:
            best = lam
    return best


def select_rfs(n=6):
    out = []
    for fp in glob.glob(os.path.join(RFDIR, "*_rf_metrics.csv")):
        probe = os.path.basename(fp).split("_")[0][-3:]
        for r in csv.DictReader(open(fp)):
            try:
                cx = float(r["centroid_azimuth_pixels"]); cy = float(r["centroid_elevation_pixels"])
                ms = min(float(r["size_azimuth_deg"]), float(r["size_elevation_deg"]))
            except (KeyError, ValueError):
                continue
            if 55 < cx < 345 and 55 < cy < 345 and 5 <= ms <= 20:
                out.append(dict(probe=probe, cl=r["cluster_id"], pol=r["rf_type"],
                                cx=int(round(cx)), cy=int(round(cy)), ms=ms))
    out = sorted(out, key=lambda d: d["ms"])
    return [out[i] for i in np.linspace(0, len(out) - 1, n).astype(int)]


def spatial_expected(frames, hS, GES, GOS, sfsS, hO, GEO, GOO, orsO):
    sf_sp, or_sp = [], []
    for im in frames:
        for yc in range(hS, 400 - hS, hS):
            for xc in range(hS, 400 - hS, hS):
                en = est(im[yc - hS:yc + hS, xc - hS:xc + hS], GES, GOS).reshape(-1, len(sfsS))
                sf_sp.append(sfsS[np.unravel_index(en.argmax(), en.shape)[1]])
        for yc in range(hO, 400 - hO, hO):
            for xc in range(hO, 400 - hO, hO):
                en = est(im[yc - hO:yc + hO, xc - hO:xc + hO], GEO, GOO).reshape(len(orsO), -1)
                or_sp.append((180 - np.degrees(orsO[np.unravel_index(en.argmax(), en.shape)[0]])) % 180)
    return np.array(sf_sp), np.array(or_sp)


def main():
    rfs = select_rfs(6)
    rows = []                                                       # summary
    ts_writer = csv.writer(open(os.path.join(OUT, "regressors.csv"), "w", newline=""))
    ts_writer.writerow(["cloud_sf", "cluster", "rf_type", "frame", "sf_cpd", "or_deg"])
    print(f"{len(CLOUDS)} clouds × {len(rfs)} RFs; SF win {SF_WIN}° OR win {OR_WIN}°, stride {STRIDE}")
    for sftok, vx in CLOUDS:
        sf0 = float(sftok.replace("sf", "").replace("p", ".")) / DPP; bsf = 0.005 / DPP
        cloud = f"theta0p000_Btheta0p785_{sftok}_Bsf0p005_{vx}_BV0p100"
        fs = sorted(glob.glob(os.path.join(GOG, cloud, "*.png")))[::STRIDE]
        hS, GES, GOS, orsS, sfsS = bank(SF_WIN, sf0, bsf)
        hO, GEO, GOO, orsO, sfsO = bank(OR_WIN, sf0, bsf)
        or_deg = np.degrees(orsO)
        frames = [readgray(f) for f in fs[::max(1, len(fs) // 8)]]
        exp_sf, exp_or = spatial_expected(frames, hS, GES, GOS, sfsS, hO, GEO, GOO, orsO)
        for rf in rfs:
            cx, cy = rf["cx"], rf["cy"]
            sf_t, omarg, conc = [], [], []
            for f in fs:
                im = readgray(f)
                enS = est(im[cy - hS:cy + hS, cx - hS:cx + hS], GES, GOS).reshape(len(orsS), len(sfsS))
                sf_t.append(sfsS[np.unravel_index(enS.argmax(), enS.shape)[1]])
                enO = est(im[cy - hO:cy + hO, cx - hO:cx + hO], GEO, GOO).reshape(len(orsO), len(sfsO))
                omarg.append(enO.sum(axis=1)); conc.append(enO.max() / enO.mean())
            sf_t = np.array(sf_t); omarg = np.array(omarg)
            lam = auto_lambda(omarg, or_deg, sf0)
            or_t = (180 - or_deg[viterbi_or(omarg, or_deg, lam)]) % 180
            for k, f in enumerate(fs):
                ts_writer.writerow([f"{sf0:.4f}", rf["cl"], rf["pol"], k, f"{sf_t[k]:.4f}", f"{or_t[k]:.1f}"])
            rows.append(dict(sf0=sf0, cl=rf["cl"], pol=rf["pol"], ms=rf["ms"],
                             sf_mean=sf_t.mean(), sf_std=sf_t.std(), or_mean=cmean(or_t), or_std=cstd(or_t),
                             exp_sf_mean=exp_sf.mean(), exp_or_std=cstd(exp_or),
                             conc=np.median(conc), lam=lam))
            print(f"  sf0={sf0:.3f} cl{rf['cl']:>3}{rf['pol'][0]} ms={rf['ms']:>2.0f}°: "
                  f"SF={sf_t.mean():.4f}(exp {exp_sf.mean():.4f}) OR σ={cstd(or_t):4.1f}(exp {cstd(exp_or):4.1f}) "
                  f"conc={np.median(conc):.1f} λ={lam:.0e}")

    # summary CSV + validation figure
    with open(os.path.join(OUT, "summary.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
    R = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    ax[0].scatter(R["exp_sf_mean"], R["sf_mean"], c=R["sf0"], cmap="viridis", s=70, edgecolor="k")
    lim = [0, R["sf_mean"].max() * 1.15]; ax[0].plot(lim, lim, "k--", lw=1)
    ax[0].set(xlabel="expected SF mean (cpd)", ylabel="per-RF SF mean", xlim=lim, ylim=lim,
              title="(A) SF mean: per-RF vs expected")
    ax[1].scatter(R["exp_or_std"], R["or_std"], c=R["sf0"], cmap="viridis", s=70, edgecolor="k")
    lim2 = [0, max(R["or_std"].max(), R["exp_or_std"].max()) * 1.15]
    ax[1].plot(lim2, lim2, "k--", lw=1); ax[1].set(xlabel="expected OR σ (deg)", ylabel="per-RF OR σ",
              xlim=lim2, ylim=lim2, title="(B) OR spread: per-RF vs expected")
    ax[2].scatter(R["ms"], R["conc"], c=R["sf0"], cmap="viridis", s=70, edgecolor="k")
    ax[2].axhline(2, color="r", ls=":", lw=1); ax[2].set(xlabel="RF min-size (deg)",
              ylabel="energy-peak concentration", title="(C) reliability flag vs RF size\n(red=trust threshold)")
    sm = plt.cm.ScalarMappable(cmap="viridis"); sm.set_array(R["sf0"]); plt.colorbar(sm, ax=ax[2], label="SF token (cpd)")
    fig.suptitle(f"Batch extraction validation — {len(CLOUDS)} clouds × {len(rfs)} RFs, "
                 f"SF {SF_WIN}°/OR {OR_WIN}°, auto-λ, 0.279°/px", y=1.02)
    fig.tight_layout(); out = os.path.join(OUT, "batch_validation.png"); fig.savefig(out, dpi=130)
    print(f"\n[saved] {out}\n[saved] {os.path.join(OUT,'regressors.csv')} + summary.csv ({len(rows)} cells)")


if __name__ == "__main__":
    main()
