"""Per-RF SF/OR extraction across ALL goggle clouds, with a per-cloud diagnostic.

For one receptive field (default cl90 / CAA-1124371, the ON RF at buffer px 268/236),
run the validated adaptive two-window Gabor recipe on every cloud folder:
    SF  -> 25deg window  (SF is spatially near-stationary; big window = accurate)
    OR  -> 10deg window  (small, to keep the local B_theta orientation spread)
           + Viterbi temporal smoothing, lambda auto-tuned to preserve >=90% raw sigma.
OR is reported in stimulus convention (180 - Gabor modulation = visual field); SF in cpd.

Writes the per-(cloud,frame) regressors (SF(t)/OR(t)) + a per-cloud diagnostic PNG:
  [SF marginal: observed(temporal) vs expected(spatial)] |
  [OR marginal: observed vs expected, offset from the cloud's theta token]    |
  [first frame, VISUAL-FIELD orientation: red=ON / blue=OFF RF edges +
   the SF (25deg) and OR (10deg) analysis windows overlaid]

Expected = slide the same windows across the frame (spatial marginal = what the RF
would see by ergodicity). deg/px = 0.279 (wisecoco). Local data only.
"""
from __future__ import annotations
import argparse, csv, glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import windowed_cloud_stat_recovery as W
from gabor_goggle_animation import viterbi_or
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpi
from matplotlib.patches import Circle

GOG = os.path.expanduser("~/local_data/motion_clouds/saved_goggles")
DPP = 111.6 / 400.0
SF_WIN, OR_WIN = 25, 10                         # deg diameters (validated recipe)
SF_MARGIN = 3                                   # Gabor SF band = sf0 +/- SF_MARGIN*B_sf


def bank(win_deg, sf0, bsf, n_sf=48, n_or=36):
    h = max(8, int(round(win_deg / DPP / 2)))
    lo, hi = max(0.01, sf0 - SF_MARGIN * bsf), sf0 + SF_MARGIN * bsf
    bk, ors, sfs = W.build_gabor_bank(2 * h, 2 * h, DPP, n_or=n_or, sf_cpd=(lo, hi), n_sf=n_sf)
    return h, np.stack([b[0] for b in bk]), np.stack([b[1] for b in bk]), ors, sfs, (lo, hi)


def readgray(f):
    im = mpi.imread(f)
    return im[..., :3].mean(-1) if im.ndim == 3 else im


def est(patch, GE, GO):
    pdc = patch - patch.mean()
    re = np.einsum("kij,ij->k", GE, pdc); iv = np.einsum("kij,ij->k", GO, pdc)
    return re * re + iv * iv


def cstd(deg):
    z = np.mean(np.exp(1j * 2 * np.radians(deg)))
    return np.degrees(np.sqrt(max(-2 * np.log(abs(z)), 0))) / 2


def cmean(deg):
    return np.degrees(np.angle(np.mean(np.exp(1j * 2 * np.radians(deg)))) / 2) % 180


def auto_lambda(omarg, or_deg):
    raw = (180 - or_deg[omarg.argmax(axis=1)]) % 180
    target = 0.9 * cstd(raw); best = 0.0
    for lam in [2e-4, 5e-4, 1e-3, 2e-3]:
        sm = (180 - or_deg[viterbi_or(omarg, or_deg, lam)]) % 180
        if cstd(sm) >= target:
            best = lam
    return best


def spatial_expected(frames, hS, GES, GOS, sfsS, hO, GEO, GOO, orsO):
    sf_sp, or_sp = [], []
    for im in frames:
        ny, nx = im.shape
        for yc in range(hS, ny - hS, hS):
            for xc in range(hS, nx - hS, hS):
                en = est(im[yc - hS:yc + hS, xc - hS:xc + hS], GES, GOS).reshape(-1, len(sfsS))
                sf_sp.append(sfsS[np.unravel_index(en.argmax(), en.shape)[1]])
        for yc in range(hO, ny - hO, hO):
            for xc in range(hO, nx - hO, hO):
                en = est(im[yc - hO:yc + hO, xc - hO:xc + hO], GEO, GOO).reshape(len(orsO), -1)
                or_sp.append((180 - np.degrees(orsO[np.unravel_index(en.argmax(), en.shape)[0]])) % 180)
    return np.array(sf_sp), np.array(or_sp)


def to_off(a, exp_or):
    return ((np.asarray(a, float) - exp_or + 90) % 180) - 90


def diagnostic_fig(out, cloud, frame0, on, off, cx, cy, sf_t, sf_exp, or_t, or_exp,
                   sf0, exp_or, band, conc):
    sf_lo, sf_hi = band
    fig = plt.figure(figsize=(15, 4.6))
    aSF = fig.add_subplot(1, 3, 1); aOR = fig.add_subplot(1, 3, 2); a0 = fig.add_subplot(1, 3, 3)
    # SF marginal (SF on y)
    aSF.hist(sf_exp, bins=26, range=(sf_lo, sf_hi), orientation="horizontal", density=True,
             color="0.7", alpha=0.55, label="expected (spatial)")
    aSF.hist(sf_t, bins=26, range=(sf_lo, sf_hi), orientation="horizontal", density=True,
             histtype="step", color="C2", lw=1.8, label="observed (temporal)")
    aSF.axhline(sf0, ls="--", color="k", lw=1, label=f"token {sf0:.3f}")
    aSF.set(ylim=(sf_lo, sf_hi), ylabel="SF (cpd)", xlabel="density")
    aSF.set_title(f"SF marginal  ({SF_WIN}° win)\nobs μ={sf_t.mean():.3f} σ={sf_t.std():.3f} | "
                  f"exp μ={sf_exp.mean():.3f}", fontsize=8)
    aSF.legend(fontsize=6, loc="upper right")
    # OR marginal (offset-from-token on y)
    aOR.hist(to_off(or_exp, exp_or), bins=26, range=(-90, 90), orientation="horizontal",
             density=True, color="0.7", alpha=0.55, label="expected")
    aOR.hist(to_off(or_t, exp_or), bins=26, range=(-90, 90), orientation="horizontal",
             density=True, histtype="step", color="C0", lw=1.8, label="observed")
    aOR.axhline(0, ls="--", color="k", lw=1)
    aOR.set(ylim=(-90, 90), ylabel=f"OR − token ({exp_or:.0f}°)", xlabel="density")
    aOR.set_title(f"OR marginal  ({OR_WIN}° win)\nobs σ={cstd(or_t):.0f}° | exp σ={cstd(or_exp):.0f}°",
                  fontsize=8)
    aOR.legend(fontsize=6, loc="upper right")
    # first frame + overlays, visual-field orientation
    a0.imshow(frame0, cmap="gray")
    if on is not None:
        a0.contour(on, levels=[0.5], colors="red", linewidths=1.6)
        a0.contour(off, levels=[0.5], colors="blue", linewidths=1.6)
    a0.add_patch(Circle((cx, cy), SF_WIN / DPP / 2, fill=False, ec="lime", lw=1.6, label="SF win"))
    a0.add_patch(Circle((cx, cy), OR_WIN / DPP / 2, fill=False, ec="orange", lw=1.6, label="OR win"))
    a0.set(xticks=[], yticks=[]); a0.invert_xaxis()
    a0.set_title(f"first frame (visual-field)\nred=ON blue=OFF | green=SF {SF_WIN}° orange=OR {OR_WIN}°\n"
                 f"conc={conc:.1f}", fontsize=8)
    fig.suptitle(f"{cloud}  ·  cl90 RF @({cx},{cy})  ·  SF token {sf0:.3f} cpd / OR token {exp_or:.0f}°  ·  "
                 f"0.279°/px", fontsize=9, y=1.02)
    fig.tight_layout(); fig.savefig(out, dpi=120, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cx", type=int, default=268); ap.add_argument("--cy", type=int, default=236)
    ap.add_argument("--rf-mask-npz", default=os.path.join(GOG, "_rfs/cl90_4371_rf_masks.npz"))
    ap.add_argument("--label", default="cl90")
    ap.add_argument("--stride", type=int, default=1, help="frame stride for regressors (1=full)")
    ap.add_argument("--one-cloud", default="", help="run a single cloud folder name (test)")
    ap.add_argument("--outdir", default=os.path.join(GOG, "_figs", "cl90_allclouds"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    extdir = os.path.join(GOG, "_extract"); os.makedirs(extdir, exist_ok=True)
    on = off = None
    if a.rf_mask_npz and os.path.exists(os.path.expanduser(a.rf_mask_npz)):
        m = np.load(os.path.expanduser(a.rf_mask_npz)); on, off = m["on"], m["off"]

    clouds = ([a.one_cloud] if a.one_cloud else
              sorted(os.path.basename(p) for p in glob.glob(os.path.join(GOG, "theta*"))
                     if os.path.isdir(p)))
    print(f"{a.label} @({a.cx},{a.cy}) over {len(clouds)} clouds; SF {SF_WIN}°/OR {OR_WIN}°, stride {a.stride}")
    reg = csv.writer(open(os.path.join(extdir, f"{a.label}_allclouds_regressors.csv"), "w", newline=""))
    reg.writerow(["cloud", "frame", "sf_cpd", "or_deg"])
    summ = []
    for cloud in clouds:
        msf = re.search(r"_sf(\d+p\d+)_", cloud); mbsf = re.search(r"Bsf(\d+p\d+)", cloud)
        mth = re.search(r"theta(-?\d+p\d+)", cloud)
        sf0 = float(msf.group(1).replace("p", ".")) / DPP
        bsf = float(mbsf.group(1).replace("p", ".")) / DPP
        exp_or = np.degrees(float(mth.group(1).replace("p", "."))) % 180
        hS, GES, GOS, orsS, sfsS, band = bank(SF_WIN, sf0, bsf)
        hO, GEO, GOO, orsO, sfsO, _ = bank(OR_WIN, sf0, bsf)
        or_deg = np.degrees(orsO)
        fs = sorted(glob.glob(os.path.join(GOG, cloud, "*.png")))
        fsub = fs[::a.stride]
        frame0 = readgray(fs[0])
        # measured temporal
        sf_t, omarg, conc = [], [], []
        for f in fsub:
            im = readgray(f)
            enS = est(im[a.cy - hS:a.cy + hS, a.cx - hS:a.cx + hS], GES, GOS).reshape(len(orsS), len(sfsS))
            sf_t.append(sfsS[np.unravel_index(enS.argmax(), enS.shape)[1]])
            enO = est(im[a.cy - hO:a.cy + hO, a.cx - hO:a.cx + hO], GEO, GOO).reshape(len(orsO), len(sfsO))
            omarg.append(enO.sum(axis=1)); conc.append(enO.max() / enO.mean())
        sf_t = np.array(sf_t); omarg = np.array(omarg)
        lam = auto_lambda(omarg, or_deg)
        or_t = (180 - or_deg[viterbi_or(omarg, or_deg, lam)]) % 180
        # expected spatial (8 frames spread)
        ef = [readgray(fs[i]) for i in np.linspace(0, len(fs) - 1, 8).astype(int)]
        sf_exp, or_exp = spatial_expected(ef, hS, GES, GOS, sfsS, hO, GEO, GOO, orsO)
        for k, f in enumerate(fsub):
            reg.writerow([cloud, k * a.stride, f"{sf_t[k]:.4f}", f"{or_t[k]:.1f}"])
        out = os.path.join(a.outdir, f"{cloud}.png")
        diagnostic_fig(out, cloud, frame0, on, off, a.cx, a.cy, sf_t, sf_exp, or_t, or_exp,
                       sf0, exp_or, band, float(np.median(conc)))
        summ.append(dict(cloud=cloud, sf_token=sf0, or_token=exp_or, sf_obs=sf_t.mean(),
                         sf_obs_std=sf_t.std(), sf_exp=sf_exp.mean(), or_obs_mean=cmean(or_t),
                         or_obs_std=cstd(or_t), or_exp_std=cstd(or_exp), conc=float(np.median(conc)), lam=lam))
        print(f"  {cloud[:46]:46} SF obs {sf_t.mean():.3f}(tok {sf0:.3f} exp {sf_exp.mean():.3f}) "
              f"OR obs σ{cstd(or_t):4.1f}(exp {cstd(or_exp):4.1f}) conc {np.median(conc):.1f}")
    with open(os.path.join(extdir, f"{a.label}_allclouds_summary.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(summ[0].keys())); wr.writeheader(); wr.writerows(summ)
    print(f"\n[saved] {a.outdir}/  ({len(clouds)} diagnostics)")
    print(f"[saved] {extdir}/{a.label}_allclouds_regressors.csv + _summary.csv")


if __name__ == "__main__":
    main()
