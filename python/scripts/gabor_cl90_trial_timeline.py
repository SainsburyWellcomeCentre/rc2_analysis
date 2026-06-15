"""cl90 / 0.06 cpd cloud: peak-energy frame + Gabor map + OR(t)/SF(t) over a real trial.

Ties the per-frame Gabor extraction to a REAL trial's clock. The cloud frame index is
POSITION-locked: the MotionClouds generator renders 1300 frames over a 130 cm corridor =
**10 frames/cm** (motion_clouds_generator.py: "Ten frames per cm"). So the cloud advances
with distance travelled, held at frame 0 during the stationary prelude (baseline = static
first frame). Verified on CAA-1124371: the 0.06 cpd trials travel ~117–146 cm → ~1300
frames. The gain (1/30/2/30/4/30) selects WHICH cloud (VX/TF rung), not the advance rate;
tf = gain·speed falls out because advance is position-locked (ref speed 30 cm/s → 1/2/4 Hz).

    frame_index(t) = clip( round( 10 · ∫_onset^t |v(τ)| dτ ), 0, N-1 )   [v in cm/s]

Figure (PDF): [peak-energy frame, visual-field, ON/OFF + SF/OR windows] | [its Gabor energy
map] ; [SF(t)] ; [OR(t)] over the trial (baseline shaded, motion clear).

Run with the rc2 .venv (needs h5py + rc2 packages + numpy/scipy/matplotlib).
"""
from __future__ import annotations
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import windowed_cloud_stat_recovery as W
from gabor_goggle_animation import viterbi_or
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpi
from matplotlib.patches import Circle
from rc2_formatted_data_reader import StimulusLookup
from rc2_formatted_data_reader.trial_conditions import GOGGLES_STIMULUS
from rc2_glm.io import load_probe_data
from rc2_glm.config import GLMConfig

GOG = os.path.expanduser("~/local_data/motion_clouds/saved_goggles")
MC = os.path.expanduser("~/local_data/motion_clouds")
DPP = 111.6 / 400.0
SF_WIN, OR_WIN, SF_MARGIN = 25, 10, 3
CX, CY = 268, 236
CLOUD = "theta0p000_Btheta0p785_sf00p016_Bsf0p005_VX0p191_BV0p100"   # 0.06 cpd, gain 1/30, θ0
FRAMES_PER_CM = 10.0
TRIAL_ID = 96          # a VT trial showing this cloud (visual+translation coupled)


def bank(win_deg, sf0, bsf, n_sf=48, n_or=36):
    h = max(8, int(round(win_deg / DPP / 2)))
    lo, hi = max(0.01, sf0 - SF_MARGIN * bsf), sf0 + SF_MARGIN * bsf
    bk, ors, sfs = W.build_gabor_bank(2 * h, 2 * h, DPP, n_or=n_or, sf_cpd=(lo, hi), n_sf=n_sf)
    return h, np.stack([b[0] for b in bk]), np.stack([b[1] for b in bk]), ors, sfs, (lo, hi)


def readgray(f):
    im = mpi.imread(f); return im[..., :3].mean(-1) if im.ndim == 3 else im


def energy(patch, GE, GO):
    p = patch - patch.mean()
    re = np.einsum("kij,ij->k", GE, p); iv = np.einsum("kij,ij->k", GO, p)
    return re * re + iv * iv


def main():
    sf0 = 0.016 / DPP; bsf = 0.005 / DPP
    exp_or = 0.0
    fs = sorted(glob.glob(os.path.join(GOG, CLOUD, "*.png")))
    N = len(fs)
    hS, GES, GOS, orsS, sfsS, band = bank(SF_WIN, sf0, bsf)
    hO, GEO, GOO, orsO, sfsO, _ = bank(OR_WIN, sf0, bsf)
    or_deg = np.degrees(orsO)
    # per-frame: SF (25° win), OR marginal (10° win), and an SF-window peak energy
    sf_t, omarg, peak_en = [], [], []
    for f in fs:
        im = readgray(f)
        enS = energy(im[CY - hS:CY + hS, CX - hS:CX + hS], GES, GOS).reshape(len(orsS), len(sfsS))
        sf_t.append(sfsS[np.unravel_index(enS.argmax(), enS.shape)[1]]); peak_en.append(enS.max())
        enO = energy(im[CY - hO:CY + hO, CX - hO:CX + hO], GEO, GOO).reshape(len(orsO), len(sfsO))
        omarg.append(enO.sum(axis=1))
    sf_t = np.array(sf_t); omarg = np.array(omarg); peak_en = np.array(peak_en)
    # OR with auto-λ Viterbi
    raw = (180 - or_deg[omarg.argmax(1)]) % 180
    def cstd(d): z = np.mean(np.exp(1j * 2 * np.radians(d))); return np.degrees(np.sqrt(max(-2 * np.log(abs(z)), 0))) / 2
    target = 0.9 * cstd(raw); lam = 0.0
    for L in [2e-4, 5e-4, 1e-3, 2e-3]:
        if cstd((180 - or_deg[viterbi_or(omarg, or_deg, L)]) % 180) >= target: lam = L
    or_t = (180 - or_deg[viterbi_or(omarg, or_deg, lam)]) % 180

    # peak-energy frame: show the TWO analysis windows' OWN energy maps (SF 25° + OR 10°),
    # since SF and OR are read from different windows. Both on an absolute 0–180° OR axis.
    pk = int(peak_en.argmax())
    imP = readgray(fs[pk])
    enS_pk = energy(imP[CY - hS:CY + hS, CX - hS:CX + hS], GES, GOS).reshape(len(orsS), len(sfsS))
    enO_pk = energy(imP[CY - hO:CY + hO, CX - hO:CX + hO], GEO, GOO).reshape(len(orsO), len(sfsO))

    # ---- real trial clock ----
    lookup = StimulusLookup(os.path.join(MC, "motion_clouds_goggles_sequence_260420.mat"),
                            "/Volumes/margrie/mvelez/mateoData_mc/image_folders.mat",
                            stimulus_set=GOGGLES_STIMULUS)
    pdata = load_probe_data(os.path.join(MC, "formatted_data_goggles", "CAA-1124371_rec1_rec2_rec3.mat"),
                            config=GLMConfig(), stimulus_lookup=lookup, cluster_indices=[0])
    trial = next(t for t in pdata.trials if t.trial_id == TRIAL_ID)
    pt = np.asarray(trial.probe_t, float); v = np.abs(np.asarray(trial.velocity, float))
    mm = np.asarray(trial.motion_mask, bool)
    mi = np.flatnonzero(mm); t0 = pt[mi[0]]
    dt = np.gradient(pt)
    vgated = v * mm                                   # advance only during the motion period
    dist = np.cumsum(vgated * dt)                     # cm travelled (0 in prelude)
    frame_idx = np.clip(np.round(FRAMES_PER_CM * dist).astype(int), 0, N - 1)
    trel = pt - t0                                    # 0 at motion onset; baseline negative
    SFt = sf_t[frame_idx]; ORt = or_t[frame_idx]
    t_pk = float(trel[int(np.argmin(np.abs(frame_idx - pk)))])   # time the left panels (frame pk) refer to
    print(f"trial {TRIAL_ID} {trial.condition}: motion {trel[mi[-1]]:.1f}s, travel {dist[-1]:.0f}cm "
          f"-> max frame {frame_idx.max()}/{N-1}; peak-energy frame {pk}")

    # ---- figure ----
    vabs = np.abs(np.asarray(trial.velocity, float))
    tmin = max(trel[0], -2.5)                         # trim the long stationary prelude
    xlim = (tmin, trel[mi[-1]])
    sf_lo, sf_hi = sfsS[0], sfsS[-1]                  # SF-window map range (shared with SF(t))

    def draw_map(ax, en, ors, sfs, title):
        or_stim = (180 - np.degrees(ors)) % 180       # bank → stimulus OR, on [0,180)
        o = np.argsort(or_stim); oss = or_stim[o]
        ax.imshow(en[o].T, origin="lower", aspect="auto", cmap="magma")
        pi, pj = np.unravel_index(en.argmax(), en.shape)
        ax.plot(int(np.where(o == pi)[0][0]), pj, "o", mfc="none", mec="w", ms=9, mew=2)
        xt = [int(np.argmin(np.abs(oss - t))) for t in (0, 45, 90, 135, 180)]
        ax.set_xticks(xt); ax.set_xticklabels([0, 45, 90, 135, 180], fontsize=7)
        st = np.arange(0, len(sfs), max(1, len(sfs) // 5))
        ax.set_yticks(st); ax.set_yticklabels([f"{sfs[j]:.3f}" for j in st], fontsize=7)
        ax.set(xlabel="OR (deg)", ylabel="SF (cpd)"); ax.set_title(title, fontsize=8)

    def break_wrap(t, y, thr=90.0):                   # NaN-break the OR line across the 0/180 seam
        y = np.asarray(y, float); j = np.where(np.abs(np.diff(y)) > thr)[0]
        return np.insert(t, j + 1, np.nan), np.insert(y, j + 1, np.nan)

    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.7], hspace=0.4, wspace=0.2)
    on = off = None
    mpath = os.path.join(GOG, "_rfs/cl90_4371_rf_masks.npz")
    if os.path.exists(mpath):
        m = np.load(mpath); on, off = m["on"], m["off"]
    # (A) full frame  +  (A2) RF-area zoom (same overlays at alpha 0.4)
    def draw_frame(ax, alpha):
        ax.imshow(imP, cmap="gray")
        if on is not None:
            ax.contour(on, levels=[0.5], colors="red", linewidths=1.6, alpha=alpha)
            ax.contour(off, levels=[0.5], colors="blue", linewidths=1.6, alpha=alpha)
        ax.add_patch(Circle((CX, CY), SF_WIN / DPP / 2, fill=False, ec="lime", lw=1.6, alpha=alpha))
        ax.add_patch(Circle((CX, CY), OR_WIN / DPP / 2, fill=False, ec="orange", lw=1.6, alpha=alpha))
        ax.set(xticks=[], yticks=[])
    sub0 = gs[0, 0].subgridspec(1, 2, wspace=0.08)
    aF = fig.add_subplot(sub0[0, 0]); aZ = fig.add_subplot(sub0[0, 1])
    draw_frame(aF, 1.0); aF.invert_xaxis()
    aF.set_title(f"peak-energy frame {pk} (visual-field)\nred=ON blue=OFF · green=SF{SF_WIN}° orange=OR{OR_WIN}°",
                 fontsize=7)
    z = int(SF_WIN / DPP / 2) + 12                    # zoom half-window ~ SF radius + margin
    draw_frame(aZ, 0.4); aZ.set_xlim(CX - z, CX + z); aZ.set_ylim(CY + z, CY - z); aZ.invert_xaxis()
    aZ.set_title("RF area (zoom)\noverlays α0.4", fontsize=7)
    # (B1) SF-window energy map  (B2) OR-window energy map
    draw_map(fig.add_subplot(gs[1, 0]), enS_pk, orsS, sfsS,
             f"SF window ({SF_WIN}°) energy — frame {pk}\nSF={sf_t[pk]:.3f} cpd (sharp in SF)")
    draw_map(fig.add_subplot(gs[2, 0]), enO_pk, orsO, sfsO,
             f"OR window ({OR_WIN}°) energy — frame {pk}\nOR={or_t[pk]:.0f}° (sharp in OR)")
    # (C) SF(t) — same SF range as the SF-window map
    aS = fig.add_subplot(gs[0, 1])
    aS.axvspan(xlim[0], 0, color="0.9", label="baseline (static frame 0)")
    aS.axvline(t_pk, color="m", lw=1.3, label=f"left panels (frame {pk}, t={t_pk:.2f}s)")
    aS.plot(trel, SFt, color="C2", lw=0.9)
    aS.axhline(0.016 / DPP, ls="--", color="k", lw=0.8, label=f"token {0.016/DPP:.3f}")
    aS.set(ylabel="SF (cpd)", ylim=(sf_lo, sf_hi), xlim=xlim,
           title=f"cl90 · 0.06 cpd cloud · trial {TRIAL_ID} ({trial.condition})  —  "
           f"OR/SF over the trial (frame advance velocity-locked)")
    aS.legend(fontsize=7, loc="upper right"); aS.tick_params(labelbottom=True)
    # (D) OR(t) — centred on the token (offset, 180° span), line broken at the seam
    aO = fig.add_subplot(gs[1, 1], sharex=aS)
    ORoff = ((ORt - exp_or + 90) % 180) - 90
    aO.axvspan(xlim[0], 0, color="0.9")
    aO.axvline(t_pk, color="m", lw=1.3)
    aO.plot(*break_wrap(trel, ORoff), color="C0", lw=0.9)
    aO.axhline(0, ls="--", color="r", lw=0.8, label=f"token ({exp_or:.0f}°)")
    aO.set(ylabel=f"OR − token (°)", ylim=(-90, 90), title="", xlim=xlim)
    aO.set_yticks([-90, -45, 0, 45, 90])
    aO.legend(fontsize=7, loc="upper right"); aO.tick_params(labelbottom=True)
    # (E) velocity profile + frame index — makes the velocity-warp explicit
    aV = fig.add_subplot(gs[2, 1], sharex=aS)
    aV.axvspan(xlim[0], 0, color="0.9")
    aV.axvline(t_pk, color="m", lw=1.3)
    aV.plot(trel, vabs, color="0.35", lw=0.9, label="|velocity| (cm/s)")
    aV.set(xlabel="time from motion onset (s)", ylabel="|v| (cm/s)", xlim=xlim)
    aVt = aV.twinx()
    aVt.plot(trel, frame_idx, color="C3", lw=1.2, label="cloud frame = 10·dist")
    aVt.set_ylabel("cloud frame", color="C3"); aVt.tick_params(axis="y", labelcolor="C3")
    aVt.set_ylim(0, N)
    h1, l1 = aV.get_legend_handles_labels(); h2, l2 = aVt.get_legend_handles_labels()
    aV.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
    fig.suptitle(f"cl90 (CAA-1124371) · {CLOUD}  ·  frame=10·∫|v|dt (10 frames/cm, gain 1/30; "
                 f"motion {trel[mi[-1]]:.1f}s, {dist[-1]:.0f}cm→frame {frame_idx.max()})", fontsize=10, y=0.995)
    out = os.path.join(GOG, "_figs", "cl90_0p06cpd_trial_timeline.pdf")
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=130, bbox_inches="tight")
    print(f"[saved] {out}")


if __name__ == "__main__":
    main()
