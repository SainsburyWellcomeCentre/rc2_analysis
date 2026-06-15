"""Gabor extractor animated over a FULL real goggle motion-cloud sequence at one RF.

Real-data analogue of `gabor_animation_goggles_rf15` (which ran on a synthetic
cloud). Reads every frame of one goggle cloud folder, windows at a real RF, runs
the optimised Gabor (fine SF bank) per frame, and renders a 3-panel GIF:
  clean patch | Gabor orientation×SF energy map (peak) | OR(t)/SF(t) + moving bar.

NB OR convention vs the nominal stimulus theta is NOT yet locked (pending the
4-orientation check), so OR is shown as the raw Gabor modulation angle (a
consistent local-orientation measure). SF is in cpd and validated on real frames.
Frames are read from a LOCAL copy (`~/local_data/motion_clouds/saved_goggles/`, not iCloud).
"""
from __future__ import annotations
import argparse, glob, os, re, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import windowed_cloud_stat_recovery as W

GOG = os.path.expanduser("~/local_data/motion_clouds/saved_goggles")
DEFAULT_CLOUD = (sorted(d for d in os.listdir(GOG) if d.startswith("theta"))[0]
                 if os.path.isdir(GOG) else "")


def viterbi_or(orient_marg, or_deg, lam):
    """Smoothest high-energy orientation path: max Σ logE − λ Σ (circular Δθ)².

    orient_marg: (T, K) per-frame orientation energy (marginalised over SF).
    Penalises implausible frame-to-frame OR jumps (circular, mod 180) so the
    track follows the energy but can't teleport. Returns bin indices per frame.
    """
    if lam <= 0:
        return orient_marg.argmax(axis=1)
    T, K = orient_marg.shape
    logE = np.log(orient_marg + 1e-12)
    logE = logE - logE.max(axis=1, keepdims=True)
    D = np.abs(or_deg[:, None] - or_deg[None, :]); D = np.minimum(D, 180 - D)
    trans = -lam * D ** 2                                    # (K_prev, K_cur)
    score = logE[0].copy(); back = np.zeros((T, K), int)
    for t in range(1, T):
        cand = score[:, None] + trans
        back[t] = cand.argmax(axis=0)
        score = cand.max(axis=0) + logE[t]
    path = np.zeros(T, int); path[-1] = int(score.argmax())
    for t in range(T - 1, 0, -1):
        path[t - 1] = back[t, path[t]]
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cloud", default=DEFAULT_CLOUD, help="cloud folder name under ~/local_data/motion_clouds/saved_goggles")
    p.add_argument("--cx", type=int, default=303); p.add_argument("--cy", type=int, default=232)
    p.add_argument("--half", type=int, default=27, help="RF half-size in px, azimuth (~15° at 0.275°/px)")
    p.add_argument("--half-y", type=int, default=0, help="RF half-size in px, elevation (0=use --half → square bbox)")
    p.add_argument("--mask", choices=["none", "gauss", "hard"], default="none",
                   help="window the patch by an elliptical RF mask (blob shape) before the Gabor")
    p.add_argument("--dpp", type=float, default=111.6 / 400.0,  # wisecoco: FOV 111.6°/400px
                   help="deg/px (0.279 = confirmed wisecoco geometry)")
    p.add_argument("--fps", type=float, default=60.0)
    p.add_argument("--gabor-nsf", type=int, default=40)
    p.add_argument("--gabor-nor", type=int, default=18, help="number of orientation bins in the Gabor energy")
    p.add_argument("--sf-margin", type=float, default=0.0,
                   help="if >0, restrict Gabor SF to sf0 ± margin*B_sf parsed from the cloud name")
    p.add_argument("--gauss-view", action="store_true",
                   help="circular RF + show patch×Gaussian (the Gabor's actual view) in panel 2")
    p.add_argument("--rf-deg", type=float, default=0.0, help="RF circle diameter (deg) overlay; 0=use 2*half")
    p.add_argument("--rf-mask-npz", default="", help="npz with on/off RF blob masks (buffer coords, from "
                   "the sparse-noise PDF); draws red=ON / blue=OFF subfield edges instead of the red circle")
    p.add_argument("--or-smooth", type=float, default=0.0002,
                   help="OR Viterbi λ (penalty per deg² jump). 2e-4 keeps the real σ≈25° while cutting "
                        "jumps ~1/3; 3e-3 over-smooths (σ→14°). Ideally auto-tuned so smoothed σ = expected σ.")
    p.add_argument("--rf-label", default="cl37 ON 15×15°")
    p.add_argument("--max-frames", type=int, default=0, help="limit frames (0=all) for quick layout checks")
    p.add_argument("--outdir", default=os.path.expanduser("~/local_data/motion_clouds/saved_goggles/_figs"))
    args = p.parse_args()

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpi
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.patches import Rectangle

    os.makedirs(args.outdir, exist_ok=True)
    cloud_dir = os.path.join(GOG, args.cloud)
    fs = sorted(glob.glob(cloud_dir + "/*.png"))
    if args.max_frames:
        fs = fs[:args.max_frames]
    print(f"cloud={args.cloud}  frames={len(fs)}  RF=({args.cx},{args.cy}) half={args.half}px  fps={args.fps}")

    from matplotlib.patches import Circle, Ellipse
    cx, cy = args.cx, args.cy
    hx = args.half; hy = args.half_y or args.half
    im0 = mpi.imread(fs[0]); im0 = im0[..., :3].mean(-1) if im0.ndim == 3 else im0

    # SF band: optionally restrict the Gabor search to the cloud's own SF band
    # (sf0 ± margin*B_sf), parsed from the cloud name — no point searching wide.
    sf_lo, sf_hi = 0.02, 0.25
    if args.sf_margin > 0:
        msf = re.search(r"_sf(\d+p\d+)_", args.cloud); mbsf = re.search(r"Bsf(\d+p\d+)", args.cloud)
        sf0 = float(msf.group(1).replace("p", ".")) / args.dpp
        bsf = float(mbsf.group(1).replace("p", ".")) / args.dpp
        sf_lo = max(0.01, sf0 - args.sf_margin * bsf); sf_hi = sf0 + args.sf_margin * bsf
        print(f"  cloud sf0={sf0:.3f} B_sf={bsf:.3f} cpd -> Gabor SF band [{sf_lo:.3f}, {sf_hi:.3f}]")
    bank, ors, sfs = W.build_gabor_bank(2 * hy, 2 * hx, args.dpp, n_or=args.gabor_nor,
                                        sf_cpd=(sf_lo, sf_hi), n_sf=args.gabor_nsf)
    GE = np.stack([b[0] for b in bank]); GO = np.stack([b[1] for b in bank])
    n_or, n_sf = len(ors), len(sfs)
    rfmask = (W.elliptical_mask(2 * hy, 2 * hx, 2 * hx, 2 * hy, args.mask)
              if args.mask != "none" else None)
    # gauss-view: visualise (panel 2) the Gabor's OWN circular Gaussian envelope
    # (analysis still uses the raw patch — the Gabor filters carry the envelope).
    env = None; sig = 0.35 * min(2 * hy, 2 * hx); rf_rad = (args.rf_deg / args.dpp / 2) if args.rf_deg else hx
    if args.gauss_view:
        yy, xx = np.mgrid[0:2 * hy, 0:2 * hx].astype(float); xx -= hx; yy -= hy
        env = np.exp(-(xx ** 2 + yy ** 2) / (2 * sig ** 2))

    patches, energies = [], []
    for f in fs:
        im = mpi.imread(f); im = im[..., :3].mean(-1) if im.ndim == 3 else im
        pa = im[cy - hy:cy + hy, cx - hx:cx + hx]
        if rfmask is not None:
            pa = pa * rfmask
        patches.append(pa)
        pdc = pa - pa.mean()
        re_ = np.einsum("kij,ij->k", GE, pdc); iv = np.einsum("kij,ij->k", GO, pdc)
        energies.append((re_ * re_ + iv * iv).reshape(n_or, n_sf))
    patches = np.array(patches); energies = np.array(energies)   # (T, n_or, n_sf)
    times = np.arange(len(fs)) / args.fps
    or_deg = np.degrees(ors)
    omarg = energies.sum(axis=2)                                 # orientation energy (T, n_or)
    # locked convention: stimulus OR = (180 − Gabor modulation) mod 180
    or_raw = (180.0 - or_deg[omarg.argmax(axis=1)]) % 180.0      # per-frame argmax (jumpy)
    path = viterbi_or(omarg, or_deg, args.or_smooth)             # temporal-smoothness prior
    or_t = (180.0 - or_deg[path]) % 180.0                        # smoothed stimulus OR
    sf_t = np.array([sfs[energies[t, path[t]].argmax()] for t in range(len(fs))])  # SF along smoothed OR
    sf_sm = W._lin_smooth_nan(sf_t, 9)
    emax = float(energies.max())
    j_raw = np.std((np.diff(or_raw) + 90) % 180 - 90); j_sm = np.std((np.diff(or_t) + 90) % 180 - 90)
    print(f"  SF over sequence: mean={np.mean(sf_t):.4f} median={np.median(sf_t):.4f} cpd | "
          f"OR frame-to-frame jump std: raw {j_raw:.1f}° -> Viterbi {j_sm:.1f}° (λ={args.or_smooth})")

    def disp(k):
        return patches[k] * env if env is not None else patches[k]
    ks = range(0, len(fs), max(1, len(fs) // 30))
    vlo = float(min(disp(k).min() for k in ks)); vhi = float(max(disp(k).max() for k in ks))

    samp = [mpi.imread(fs[i]) for i in np.linspace(0, len(fs) - 1, 40).astype(int)]
    samp = [s[..., :3].mean(-1) if s.ndim == 3 else s for s in samp]
    fvlo, fvhi = float(np.min(samp)), float(np.max(samp))
    fov_deg = im0.shape[1] * args.dpp

    # optional ON/OFF RF subfield masks (buffer coords) extracted from the sparse-noise PDF
    rf_on = rf_off = None
    if args.rf_mask_npz:
        _m = np.load(os.path.expanduser(args.rf_mask_npz))
        rf_on, rf_off = _m["on"], _m["off"]

    def add_rf_markers(ax):
        if rf_on is not None:
            # draw the actual RF subfield edges (red=ON, blue=OFF) instead of the circle;
            # keep the cyan Gabor-window outline so the analysis footprint stays visible
            ax.contour(rf_on, levels=[0.5], colors="red", linewidths=1.8)
            ax.contour(rf_off, levels=[0.5], colors="blue", linewidths=1.8)
            if args.gauss_view:
                ax.add_patch(Circle((cx, cy), 2 * sig, fill=False, ec="cyan", lw=1.2, ls="--"))
        elif args.gauss_view:
            ax.add_patch(Circle((cx, cy), rf_rad, fill=False, ec="red", lw=2))
            ax.add_patch(Circle((cx, cy), 2 * sig, fill=False, ec="cyan", lw=1.5, ls="--"))
        elif args.mask != "none":
            ax.add_patch(Ellipse((cx, cy), 2 * hx, 2 * hy, fill=False, ec="red", lw=2))
        else:
            ax.add_patch(Rectangle((cx - hx, cy - hy), 2 * hx, 2 * hy, fill=False, ec="red", lw=2))

    if args.gauss_view:
        # both labelled by DIAMETER so equal label ⟺ equal circle. Red circle radius
        # = rf_rad (Ø=2·rf_rad); cyan circle radius = 2σ (Ø=4σ). The Gabor 2σ window
        # (Ø~20°) is ~2× the RF (Ø~10°) by design — more cycles for a cleaner SF.
        rf_lbl = ("red=ON / blue=OFF RF edges (left=lateral, motion R→L)"
                  if rf_on is not None else
                  f"red RF Ø{2*rf_rad*args.dpp:.0f}° (left=lateral, motion R→L)")
        ftitle = (f"full frame {im0.shape[1]}px ≈ {fov_deg:.0f}° · visual-field view\n"
                  f"{rf_lbl}\n"
                  f"cyan Gabor 2σ Ø{4*sig*args.dpp:.0f}°")
        ptitle = "pixels inside the Gabor\n(patch × Gaussian, visual-field)"
    else:
        ftitle = (f"full frame {im0.shape[1]}px ≈ {fov_deg:.0f}° · visual-field view\n"
                  f"red RF {2*hx*args.dpp:.0f}×{2*hy*args.dpp:.0f}° (left=lateral, motion R→L)")
        ptitle = "RF patch (zoom, visual-field)"

    # --- expected SF/OR marginals: slide the SAME window across a few frames (spatial) ---
    mth = re.search(r"theta(-?\d+p\d+)", args.cloud)
    exp_or = (np.degrees(float(mth.group(1).replace("p", "."))) % 180) if mth else 0.0   # expected stim OR
    msf = re.search(r"_sf(\d+p\d+)_", args.cloud)
    sf0v = (float(msf.group(1).replace("p", ".")) / args.dpp) if msf else float(np.median(sf_t))
    ny0, nx0 = im0.shape
    sf_exp, or_exp = [], []
    for fi in np.linspace(0, len(fs) - 1, 8).astype(int):
        ime = mpi.imread(fs[fi]); ime = ime[..., :3].mean(-1) if ime.ndim == 3 else ime
        for yc in range(hy, ny0 - hy, hy):
            for xc in range(hx, nx0 - hx, hx):
                pp = ime[yc - hy:yc + hy, xc - hx:xc + hx]; pdc = pp - pp.mean()
                rr = np.einsum("kij,ij->k", GE, pdc); ii = np.einsum("kij,ij->k", GO, pdc)
                en = (rr * rr + ii * ii).reshape(n_or, n_sf); a_, b_ = np.unravel_index(en.argmax(), en.shape)
                sf_exp.append(sfs[b_]); or_exp.append((180 - np.degrees(ors[a_])) % 180)
    sf_exp = np.array(sf_exp); or_exp = np.array(or_exp)

    def to_off(a):  # OR offset from expected, wrapped to [-90, 90]
        return ((np.asarray(a, float) - exp_or + 90) % 180) - 90
    # energy-map centring: order orientation bins by their offset from expected so
    # the x-axis is a clean monotonic −90→+90 with expected (0) in the middle
    stim_off = to_off((180 - np.degrees(ors)) % 180)
    order = np.argsort(stim_off); off_sorted = stim_off[order]
    xt = [int(np.argmin(np.abs(off_sorted - t))) for t in (-90, -45, 0, 45, 90)]
    xtl = [f"{int(round(off_sorted[c])):+d}" for c in xt]
    ctr_col = int(np.argmin(np.abs(off_sorted)))

    fig = plt.figure(figsize=(23, 4.7))
    aSF = fig.add_subplot(1, 5, 1); aOR = fig.add_subplot(1, 5, 2)
    a0 = fig.add_subplot(1, 5, 3); a1 = fig.add_subplot(1, 5, 4); a2 = fig.add_subplot(1, 5, 5)
    # SF marginal (vertical: SF on y)
    aSF.hist(sf_exp, bins=24, range=(sf_lo, sf_hi), orientation="horizontal", density=True,
             color="0.7", alpha=0.55, label="expected")
    aSF.hist(sf_t, bins=24, range=(sf_lo, sf_hi), orientation="horizontal", density=True,
             histtype="step", color="C2", lw=1.6, label="measured")
    aSF.axhline(sf0v, ls="--", color="k", lw=0.8)
    sf_mark = aSF.axhline(sf_t[0], color="red", lw=2)
    aSF.set(ylim=(sf_lo, sf_hi), ylabel="SF (cpd)", xlabel="density"); aSF.set_title("SF marginal", fontsize=9)
    aSF.legend(fontsize=6, loc="upper right")
    # OR marginal (vertical: OR offset from expected on y)
    aOR.hist(to_off(or_exp), bins=24, range=(-90, 90), orientation="horizontal", density=True,
             color="0.7", alpha=0.55, label="expected")
    aOR.hist(to_off(or_t), bins=24, range=(-90, 90), orientation="horizontal", density=True,
             histtype="step", color="C0", lw=1.6, label="measured")
    aOR.axhline(0, ls="--", color="k", lw=0.8)
    or_mark = aOR.axhline(float(to_off(or_t[0])), color="red", lw=2)
    aOR.set(ylim=(-90, 90), ylabel=f"OR − expected ({exp_or:.0f}°)", xlabel="density")
    aOR.set_title("OR marginal", fontsize=9); aOR.legend(fontsize=6, loc="upper right")
    # full frame + patch
    imF = a0.imshow(im0, cmap="gray", vmin=fvlo, vmax=fvhi); add_rf_markers(a0)
    a0.set(xticks=[], yticks=[]); a0.set_title(ftitle, fontsize=8)
    im1 = a1.imshow(disp(0), cmap="gray", vmin=vlo, vmax=vhi)
    a1.set(xticks=[], yticks=[]); a1.set_title(ptitle, fontsize=8)
    # Display the two raw-image panels in VISUAL-FIELD orientation (mirror x), matching
    # the RF PDF's `xdir reverse`. Analysis/windowing stays in raw-buffer coords (markers
    # drawn at cx in data coords land on the left automatically); RF then reads on the left
    # = lateral and motion runs R→L, consistent with the energy map / OR marginal which are
    # already in stim convention (180−mod = visual field). See CHANGELOG 2026-06-11.
    a0.invert_xaxis(); a1.invert_xaxis()
    # 2D energy map centred on expected OR (columns sorted by offset → monotonic axis)
    st = np.arange(0, n_sf, max(1, n_sf // 6))

    def draw_emap(k):
        a2.clear()
        a2.imshow(energies[k][order].T, origin="lower", aspect="auto", cmap="magma", vmin=0, vmax=emax)
        pk = np.unravel_index(energies[k].argmax(), energies[k].shape)
        a2.plot(int(np.where(order == pk[0])[0][0]), pk[1], "co", ms=8, mfc="none", mew=2)
        a2.axvline(ctr_col, color="cyan", lw=0.8, ls=":")
        a2.set_xticks(xt); a2.set_xticklabels(xtl)
        a2.set_yticks(st); a2.set_yticklabels([f"{sfs[j]:.3f}" for j in st])
        a2.set_xlabel(f"OR − expected ({exp_or:.0f}°)"); a2.set_ylabel("SF (cpd)")
        a2.set_title(f"Gabor energy (centred)\nOR={or_t[k]:.0f}° SF={sf_t[k]:.3f}\nt={times[k]:.2f}s fr{k}",
                     fontsize=8)
    draw_emap(0)

    def update(k):
        full = mpi.imread(fs[k]); full = full[..., :3].mean(-1) if full.ndim == 3 else full
        imF.set_data(full); im1.set_data(disp(k))
        draw_emap(k)
        sf_mark.set_ydata([sf_t[k], sf_t[k]]); o = float(to_off(or_t[k])); or_mark.set_ydata([o, o])
        return [imF, im1]

    fig.suptitle(f"{args.cloud}  ·  {args.rf_label} @({cx},{cy})  ·  {len(fs)} fr @ {args.fps:g}Hz  ·  "
                 f"{args.dpp:.3f}°/px  ·  SF band [{sf_lo:.3f},{sf_hi:.3f}]  ·  OR=stimulus (180−mod, locked)",
                 y=1.06, fontsize=8)
    anim = FuncAnimation(fig, update, frames=len(fs), interval=1000 / args.fps, blit=False)
    fig.tight_layout(rect=[0, 0, 1, 0.84])
    tag = "".join(c if c.isalnum() else "_" for c in args.rf_label).strip("_")
    out = os.path.join(args.outdir, f"gabor_seq_fullframe_{args.cloud}__{tag}.gif")
    anim.save(out, writer=PillowWriter(fps=args.fps), dpi=70)
    print(f"[saved] {out}  ({os.path.getsize(out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
