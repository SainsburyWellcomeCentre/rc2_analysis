"""Windowed motion-cloud statistic recovery — how small a window can we trust?

WHY THIS EXISTS
---------------
The motion-clouds GLM wants SF / TF / orientation (OR) as *what the receptive
field actually saw*, extracted per frame from the portion of the cloud inside
each cluster's RF (goggles RF maps on ceph). The concern: a narrow RF window is
a poor spectral estimator. This script measures that degradation against KNOWN
ground truth, on a canvas and with windows that match the REAL goggle stimulus.

GOVERNING VARIABLE (the key insight)
------------------------------------
Windowed-extraction feasibility is set by DEGREES, not pixels. The number of
spatial cycles inside a window is

        cycles_in_window = SF[cpd] * window_extent[deg]

independent of render resolution. A window can't resolve a spatial frequency
below ~1 cycle across its extent. So we report SF recovery against
cycles_in_window, and the pixel grid only has to sample finely enough not to
alias the top SF (trivially satisfied at the real deg/px). TF lives on the
temporal axis (each pixel oscillates at TF = speed*SF) so it should be ~window
independent; OR is a local gradient so it degrades only as pixel-count noise.

REAL DISPLAY + RF GEOMETRY (confirmed from the code/data)
---------------------------------------------------------
- Screens motion clouds: 960 px @ 0.1024 deg/px (34.7 cm @ 15 cm -> 98.3 deg).
- Goggles: sparse-noise RF display ~400 px @ 0.275 deg/px (~110 deg FOV;
  degs_per_square=5, calculate_rf_mouse_goggles.m). Cloud token relation agrees:
  0.008 cpp / 0.03 cpd ~ 0.27 deg/px. The RFs (CAA-1124370/371) are GOGGLE RFs.
- RF sizes are real, ELLIPTICAL, quantised to 5 deg: azimuth 5-20 deg,
  elevation 10-30 deg, mean ~12 x 15 deg, with ON(white)/OFF(black) subfields.
  Loaded directly from *_rf_metrics.csv (size_azimuth_deg / size_elevation_deg).

UNITS NOTE (load-bearing — see project_motion_clouds_goggles_stimulus_units)
---------------------------------------------------------------------------
The GLM's SF tokens {0.003, 0.006, 0.012} are cycles-per-PIXEL, not cpd.
Physical SF = {0.03, 0.06, 0.12} cpd. We work in cpd throughout.

This is a methods simulation, not an rc2 analysis. Reads only the RF metrics CSV
(for real window sizes); writes only inside --outdir. Run in a fresh venv with
numpy / scipy / matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
from dataclasses import dataclass, field

import numpy as np

# ---------------------------------------------------------------------------
# Physical stimulus grid (cpd / Hz / rad).
# ---------------------------------------------------------------------------
SF_LEVELS_CPD = (0.03, 0.06, 0.12)
TF_LEVELS_HZ = (1.0, 2.0, 4.0)
THETA_LEVELS_RAD = (-np.pi / 4, 0.0, np.pi / 4, np.pi / 2)
B_SF_CPD = 0.02
B_THETA_RAD = np.pi / 4
B_V = 0.2

DEFAULT_RF_CSV_GLOB = (
    "/Volumes/margrie/laura/data transfer for laura/RFs googles/*_rf_metrics.csv"
)
# Empirical fallback (counts of az x el deg) measured from the two probes, used
# if the ceph CSV is unreachable. Keeps the sim runnable offline.
RF_FALLBACK_DEG = [
    (10, 15), (15, 15), (10, 10), (15, 20), (10, 20), (15, 10),
    (20, 20), (15, 30), (5, 10),
]


@dataclass
class DisplayProfile:
    name: str
    deg_per_px: float
    fov_w_deg: float      # canvas width (full-frame reference)
    fov_h_deg: float


SCREENS = DisplayProfile("screens", 0.1024, 98.3, 55.3)
GOGGLES = DisplayProfile("goggles", 0.275, 110.0, 80.0)
PROFILES = {"screens": SCREENS, "goggles": GOGGLES}


# ============================================================================
# 1. GENERATOR — envelope_gabor + random_cloud, in degrees at a real deg/px.
# ============================================================================
def make_cloud(
    nx: int, ny: int, n_frame: int, fps: float, deg_per_px: float,
    sf_0_cpd: float, b_sf_cpd: float, theta: float, b_theta: float,
    tf_hz: float, b_v: float, rng: np.random.Generator,
) -> np.ndarray:
    """Return a (n_frame, ny, nx) real motion-cloud movie in [0, 1].

    Gaussian envelope in 3-D Fourier space (radial SF band x orientation wedge x
    rigid-motion speed plane) filled with random phase. Power spectrum ~ env**2.
    Frequencies are in cycles/pixel; cpd is converted via deg_per_px.
    """
    fx = np.fft.fftfreq(nx)[None, None, :]
    fy = np.fft.fftfreq(ny)[None, :, None]
    ft = np.fft.fftfreq(n_frame)[:, None, None]

    f_radius = np.sqrt(fx**2 + fy**2)
    f_radius_safe = np.where(f_radius == 0, 1e-6, f_radius)

    sf_0 = sf_0_cpd * deg_per_px        # cpd -> cycles/pixel
    b_sf = b_sf_cpd * deg_per_px
    # Drift (px/frame) so dominant temporal freq at sf_0 is tf_hz:
    #   ft[cyc/frame] = V * sf_0[cyc/px] = tf_hz/fps  ->  V = tf/(fps*sf_0)
    v_x = (tf_hz / fps) / sf_0

    env_sf = np.exp(-0.5 * ((f_radius - sf_0) / b_sf) ** 2)
    angle = np.arctan2(fy, fx)
    dtheta = np.angle(np.exp(1j * 2 * (angle - theta))) / 2
    env_or = np.exp(-0.5 * (dtheta / b_theta) ** 2)
    env_speed = np.exp(-0.5 * ((ft + v_x * fx) / (b_v * f_radius_safe)) ** 2)

    envelope = env_sf * env_or * env_speed
    envelope = np.where(f_radius == 0, 0.0, envelope)  # kill DC (broadcasts over t)

    spectrum = envelope * np.exp(1j * rng.uniform(0, 2 * np.pi, envelope.shape))
    movie = np.real(np.fft.ifftn(spectrum))
    movie -= movie.min()
    if movie.max() > 0:
        movie /= movie.max()
    return movie


# ============================================================================
# 2. WINDOWS — elliptical (azimuth=x, elevation=y), real RF extents in deg.
# ============================================================================
def elliptical_mask(ny, nx, ax_px, ay_px, kind="gauss", cx=None, cy=None):
    """Elliptical window centred at (cx,cy) [default canvas centre].

    ax/ay are full extents (deg-derived) in px.
    """
    if cx is None:
        cx = nx / 2
    if cy is None:
        cy = ny / 2
    yy, xx = np.ogrid[:ny, :nx]
    if kind == "hard":
        m = ((xx - cx) / (ax_px / 2)) ** 2 + ((yy - cy) / (ay_px / 2)) ** 2 <= 1.0
        return m.astype(float)
    sx, sy = ax_px / 2.355, ay_px / 2.355   # treat extent as FWHM
    return np.exp(-(((xx - cx) ** 2) / (2 * sx**2) + ((yy - cy) ** 2) / (2 * sy**2)))


def apply_window2d(frame, mask):
    return frame * mask


# ============================================================================
# 3. ESTIMATORS
# ============================================================================
def recover_sf_cpd(frame2d, deg_per_px, n_bins=80):
    """Radially-averaged FFT magnitude -> PEAK spatial frequency (cpd).

    The peak of the radial profile tracks the stimulus SF; a power-weighted
    centroid is biased by the window's frequency support (compressed toward
    mid-band) and was misleading, so we report the peak.
    """
    f = np.fft.fftshift(np.fft.fft2(frame2d.astype(float)))
    mag = np.abs(f)
    ny, nx = frame2d.shape
    cy, cx = ny // 2, nx // 2
    yy, xx = np.indices(frame2d.shape)
    rr_cpp = np.sqrt(((xx - cx) / nx) ** 2 + ((yy - cy) / ny) ** 2)  # cycles/px
    freqs_cpd = (rr_cpp / deg_per_px).ravel()
    nyq_cpd = 0.5 / deg_per_px
    wsum, edges = np.histogram(freqs_cpd, bins=n_bins, range=(0, nyq_cpd),
                               weights=mag.ravel())
    cnt, _ = np.histogram(freqs_cpd, bins=n_bins, range=(0, nyq_cpd))
    prof = wsum / np.maximum(cnt, 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    prof[0] = 0.0  # drop the DC bin
    return float(centers[np.argmax(prof)]) if prof.max() > 0 else np.nan


def recover_or_rad(frame2d):
    from scipy.ndimage import convolve
    kx = np.array([[5, 8, 10, 8, 5], [4, 10, 20, 10, 4], [0, 0, 0, 0, 0],
                   [-4, -10, -20, -10, -4], [-5, -8, -10, -8, -5]], float)
    ky = kx.T[::-1].copy()
    img = frame2d.astype(float)
    gx = convolve(img, kx, mode="reflect")
    gy = convolve(img, ky, mode="reflect")
    edge = np.arctan2(-gy, gx) + np.pi / 2
    w = gx**2 + gy**2
    z = np.sum(w * np.exp(1j * 2 * edge)) / max(np.sum(w), 1e-12)
    return float(np.angle(z) / 2)


def build_gabor_bank(ny, nx, dpp, n_or=18, sf_cpd=(0.02, 0.25), n_sf=14,
                     sigma_frac=0.35):
    """Quadrature Gabor bank over (modulation orientation θ, spatial freq).

    Their (zebra-noise) approach: match the windowed patch to a library of Gabor
    wavelets. θ is the modulation direction (= SF-vector angle, generator's
    `theta` convention); bar orientation = θ + 90°. Even filters are zero-meaned
    so DC doesn't dominate. Returns the bank + axes for the energy map.
    """
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    xx -= cx; yy -= cy
    sigma = sigma_frac * min(ny, nx)
    env = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
    ors = np.linspace(0, np.pi, n_or, endpoint=False)
    sfs_cpd = np.geomspace(sf_cpd[0], sf_cpd[1], n_sf)
    bank = []
    for th in ors:
        xr = xx * np.cos(th) + yy * np.sin(th)              # along modulation
        for f in sfs_cpd * dpp:                             # cpd -> cyc/px
            ge = env * np.cos(2 * np.pi * f * xr)
            go = env * np.sin(2 * np.pi * f * xr)
            ge = ge - ge.mean()                             # kill DC of even filter
            bank.append((ge, go))
    return bank, ors, sfs_cpd


def gabor_response(patch, bank, ors, sfs_cpd):
    """Energy of each Gabor vs the (DC-removed) patch; return best OR/SF + map."""
    p = patch - patch.mean()
    n_or, n_sf = len(ors), len(sfs_cpd)
    energy = np.empty(len(bank))
    for k, (ge, go) in enumerate(bank):
        re = np.sum(p * ge); im = np.sum(p * go)
        energy[k] = re * re + im * im
    energy = energy.reshape(n_or, n_sf)
    i, j = np.unravel_index(np.argmax(energy), energy.shape)
    # bar orientation in the blob's y-up convention: bars ⊥ modulation, and the
    # image y-axis is flipped -> bar = 90 - θ (NOT θ+90; the difference is 2θ,
    # the mirror-about-45° artifact that also bit the gradient estimator).
    return dict(or_mod_deg=np.degrees(ors[i]),
                or_bar_deg=(90 - np.degrees(ors[i])) % 180,
                sf_cpd=float(sfs_cpd[j]), energy=energy, peak=(i, j))


def _fit_blobs(mask, min_area, sig, polarity):
    """Connected-component second-moment ellipses for a binary mask (y-up)."""
    from scipy import ndimage

    ny, nx = mask.shape
    cyc, cxc = ny / 2.0, nx / 2.0
    lbl, n = ndimage.label(mask)
    out = []
    for k in range(1, n + 1):
        ys, xs = np.where(lbl == k)
        if ys.size < min_area:
            continue
        cy_, cx_ = ys.mean(), xs.mean()
        xx = xs - cx_; yy = -(ys - cy_)                  # y-up
        cxx = (xx * xx).mean(); cyy = (yy * yy).mean(); cxy = (xx * yy).mean()
        th = 0.5 * np.arctan2(2 * cxy, cxx - cyy)
        tr, det = cxx + cyy, cxx * cyy - cxy * cxy
        disc = np.sqrt(max(tr * tr / 4 - det, 0.0))
        major = 2 * np.sqrt(max(tr / 2 + disc, 1e-9))
        minor = 2 * np.sqrt(max(tr / 2 - disc, 1e-9))
        d = np.hypot(cy_ - cyc, cx_ - cxc)
        w = np.exp(-0.5 * (d / sig) ** 2) * ys.size
        out.append(dict(cy=cy_, cx=cx_, theta=th, major=major, minor=minor,
                        area=ys.size, w=w, pol=polarity))
    return out


def blob_estimate(patch, dpp, thresh_q=0.65, min_area_frac=0.004,
                  center_sigma_frac=0.5):
    """Morphological OR/SF/bandwidth from BRIGHT and DARK blobs (Laura's idea).

    Binarise both polarities (bright = p > q-quantile, dark = p < (1-q)-quantile
    = the opposite-phase stripes), connected-component label, fit a second-moment
    ellipse per blob, weight by size x centrality. Then:
      OR  = weighted circular-mean of ALL (bright+dark) blob major-axis angles.
      bandwidth = circular spread of blob orientations (deg).
      SF (minor)   = 1 / (2 * weighted minor axis, deg)        [RAW, calib.]
      SF (spacing) = 1 / (2 * median bright<->nearest-dark centroid dist, deg)
                     — half-wavelength from opposite-phase spacing; works at
                     half a cycle (the real factor-2 win).
    """
    p = np.asarray(patch, float)
    ny, nx = p.shape
    bthr = np.quantile(p, thresh_q)
    dthr = np.quantile(p, 1 - thresh_q)
    bright_mask = p > bthr
    dark_mask = p < dthr
    min_area = max(4, int(min_area_frac * nx * ny))
    sig = center_sigma_frac * min(ny, nx) / 2.0
    bright = _fit_blobs(bright_mask, min_area, sig, +1)
    dark = _fit_blobs(dark_mask, min_area, sig, -1)
    blobs = bright + dark
    if not blobs:
        return dict(or_deg=np.nan, sf_minor=np.nan, sf_spacing=np.nan,
                    bw_deg=np.nan, n_blobs=0, blobs=[],
                    bright_mask=bright_mask, dark_mask=dark_mask,
                    bthr=bthr, dthr=dthr)
    w = np.array([b["w"] for b in blobs])
    th = np.array([b["theta"] for b in blobs])
    z = np.sum(w * np.exp(1j * 2 * th)) / np.sum(w)
    or_mean = np.angle(z) / 2
    bw = np.sqrt(max(-2 * np.log(max(abs(z), 1e-9)), 0.0)) / 2
    minor_deg = np.sum(w * np.array([b["minor"] for b in blobs])) / np.sum(w) * dpp
    sf_minor = 1.0 / (2 * minor_deg) if minor_deg > 0 else np.nan
    # opposite-phase bright<->dark nearest-neighbour vectors:
    #   |vector| -> half wavelength (SF);  direction -> modulation axis (OR).
    sf_spacing, or_vec_deg, vecs = np.nan, np.nan, []
    if bright and dark:
        bc = np.array([[b["cy"], b["cx"]] for b in bright])
        dc = np.array([[b["cy"], b["cx"]] for b in dark])
        d_bd = np.sqrt(((bc[:, None, :] - dc[None, :, :]) ** 2).sum(-1))
        bw_ = np.array([b["w"] for b in bright])
        dw_ = np.array([b["w"] for b in dark])
        ang_list, wt_list, dist_list = [], [], []
        for bi, dj in enumerate(d_bd.argmin(axis=1)):          # bright -> nearest dark
            dy = -(dc[dj, 0] - bc[bi, 0]); dx = dc[dj, 1] - bc[bi, 1]
            ang_list.append(np.arctan2(dy, dx)); wt_list.append(bw_[bi])
            dist_list.append(d_bd[bi, dj]); vecs.append((bc[bi, 1], bc[bi, 0], dx, -dy))
        for di, bj in enumerate(d_bd.argmin(axis=0)):          # dark -> nearest bright
            dy = -(bc[bj, 0] - dc[di, 0]); dx = bc[bj, 1] - dc[di, 1]
            ang_list.append(np.arctan2(dy, dx)); wt_list.append(dw_[di])
            dist_list.append(d_bd[bj, di])
        spacing_deg = np.median(dist_list) * dpp
        if spacing_deg > 0:
            sf_spacing = 1.0 / (2 * spacing_deg)
        wt = np.array(wt_list)
        zc = np.sum(wt * np.exp(1j * 2 * np.array(ang_list))) / np.sum(wt)
        or_vec_deg = np.degrees(np.angle(zc) / 2)              # modulation axis (≈ theta)
    return dict(or_deg=np.degrees(or_mean), or_vec_deg=or_vec_deg, sf_minor=sf_minor,
                sf_spacing=sf_spacing, bw_deg=np.degrees(bw), n_blobs=len(blobs),
                blobs=blobs, vecs=vecs, bright_mask=bright_mask, dark_mask=dark_mask,
                bthr=bthr, dthr=dthr)


def recover_tf_hz(patch_movie, fps):
    """Temporal FFT per pixel, power-averaged over the window -> peak Hz."""
    nf = patch_movie.shape[0]
    flat = patch_movie.reshape(nf, -1)
    flat = flat - flat.mean(axis=0, keepdims=True)
    power = (np.abs(np.fft.rfft(flat, axis=0)) ** 2).mean(axis=1)
    freqs = np.fft.rfftfreq(nf, d=1.0 / fps)
    power[0] = 0.0
    return float(freqs[np.argmax(power)])


# ============================================================================
# 4. RF TABLE
# ============================================================================
def load_rf_ellipses(csv_glob):
    """Return list of (size_az_deg, size_el_deg) from real RF metrics CSVs."""
    rows = []
    files = glob.glob(csv_glob)
    if not files:
        print(f"[warn] no RF CSV at {csv_glob}; using empirical fallback sizes")
        return [(float(a), float(e)) for a, e in RF_FALLBACK_DEG]
    for fp in files:
        with open(fp) as fh:
            for r in csv.DictReader(fh):
                try:
                    rows.append((float(r["size_azimuth_deg"]),
                                 float(r["size_elevation_deg"])))
                except (KeyError, ValueError):
                    continue
    print(f"[rf] loaded {len(rows)} real RF ellipses from {len(files)} file(s)")
    return rows


# ============================================================================
# 5. SWEEP
# ============================================================================
@dataclass
class Config:
    profile: DisplayProfile = field(default_factory=lambda: GOGGLES)
    n_frame: int = 64
    fps: float = 32.0
    n_real: int = 8
    # continuous square-window sweep, in degrees (full-frame ref = fov)
    sweep_deg: tuple = (90, 60, 40, 30, 20, 15, 10, 7, 5)
    seed: int = 0
    rf_csv_glob: str = DEFAULT_RF_CSV_GLOB
    records: list = field(default_factory=list)


def _canvas_px(cfg):
    nx = int(round(cfg.profile.fov_w_deg / cfg.profile.deg_per_px))
    ny = int(round(cfg.profile.fov_h_deg / cfg.profile.deg_per_px))
    return nx, ny


def run_sweep(cfg: Config):
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    rec = cfg.records

    def deg2px(d):
        return max(6, int(round(d / dpp)))

    # ---- (A) continuous square-window sweep: SF & TF vs window deg ----
    for sf in SF_LEVELS_CPD:
        for tf in TF_LEVELS_HZ:
            for _ in range(cfg.n_real):
                mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, sf, B_SF_CPD,
                                0.0, B_THETA_RAD, tf, B_V, rng)
                mid = mv[cfg.n_frame // 2]
                for wd in cfg.sweep_deg:
                    wpx = deg2px(wd)
                    m = elliptical_mask(ny, nx, wpx, wpx, "gauss")
                    cyc = sf * wd
                    rec.append(("SF_sq", sf, wd, cyc, recover_sf_cpd(mid * m, dpp)))
                    rec.append(("TF_sq", tf, wd, cyc,
                                recover_tf_hz(mv * m[None], cfg.fps)))

    # ---- (B) OR vs window deg ----
    for th in THETA_LEVELS_RAD:
        for _ in range(cfg.n_real):
            mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, 0.06, B_SF_CPD,
                            th, B_THETA_RAD, 2.0, B_V, rng)
            mid = mv[cfg.n_frame // 2]
            for wd in cfg.sweep_deg:
                wpx = deg2px(wd)
                m = elliptical_mask(ny, nx, wpx, wpx, "gauss")
                rec.append(("OR_sq", th, wd, np.nan, recover_or_rad(mid * m)))

    # ---- (C) REAL elliptical RF windows: SF & TF at each real (az,el) ----
    rfs = load_rf_ellipses(cfg.rf_csv_glob)
    uniq = sorted(set(rfs))
    for (az, el) in uniq:
        axpx, aypx = deg2px(az), deg2px(el)
        for sf in SF_LEVELS_CPD:
            tf = 2.0
            for _ in range(cfg.n_real):
                mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, sf, B_SF_CPD,
                                0.0, B_THETA_RAD, tf, B_V, rng)
                mid = mv[cfg.n_frame // 2]
                m = elliptical_mask(ny, nx, axpx, aypx, "gauss")
                cyc = sf * min(az, el)  # limiting axis
                rec.append(("SF_rf", sf, (az, el), cyc,
                            recover_sf_cpd(mid * m, dpp)))
                rec.append(("TF_rf", tf, (az, el), cyc,
                            recover_tf_hz(mv * m[None], cfg.fps)))
    return rec


# ============================================================================
# 6. PLOT + TABLE
# ============================================================================
def circ_err_deg(rec_rad, true_rad):
    d = np.angle(np.exp(1j * 2 * (rec_rad - true_rad))) / 2
    return abs(np.degrees(d))


def summarise_and_plot(rec, cfg, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    nx, ny = _canvas_px(cfg)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # SF vs window deg (square sweep)
    a = ax[0, 0]
    for sf in SF_LEVELS_CPD:
        rs = [r for r in rec if r[0] == "SF_sq" and r[1] == sf]
        wds = sorted({r[2] for r in rs})
        mu = [np.nanmean([r[4] for r in rs if r[2] == w]) for w in wds]
        sd = [np.nanstd([r[4] for r in rs if r[2] == w]) for w in wds]
        a.errorbar(wds, mu, yerr=sd, marker="o", capsize=3, label=f"{sf:.2f} cpd")
        a.axhline(sf, ls=":", lw=0.8, color="gray")
    a.set(xscale="log", xlabel="square window (deg)",
          ylabel="recovered SF (cpd)", title="SF vs window size")
    a.legend(fontsize=7)

    # SF recovery error vs cycles-in-window (the governing variable) + real RFs
    a = ax[0, 1]
    rs = [r for r in rec if r[0] == "SF_sq"]
    cyc = np.array([r[3] for r in rs])
    err = np.array([abs(r[4] - r[1]) / r[1] for r in rs])  # rel error
    order = np.argsort(cyc)
    a.plot(cyc[order], err[order], ".", alpha=0.3, color="gray", label="square sweep")
    # bin
    bins = np.logspace(np.log10(max(cyc.min(), 1e-2)), np.log10(cyc.max()), 12)
    idx = np.digitize(cyc, bins)
    bx = [cyc[idx == i].mean() for i in range(1, len(bins)) if np.any(idx == i)]
    by = [err[idx == i].mean() for i in range(1, len(bins)) if np.any(idx == i)]
    a.plot(bx, by, "-o", color="C3", label="binned mean")
    rf_cyc = sorted({r[3] for r in rec if r[0] == "SF_rf"})
    for c in rf_cyc:
        a.axvline(c, ls="--", lw=0.5, color="C0", alpha=0.5)
    a.axvline(1.0, color="k", lw=1.2)
    a.text(1.0, a.get_ylim()[1] * 0.9, " 1 cycle", fontsize=8)
    a.set(xscale="log", xlabel="cycles in window (SF_cpd x window_deg)",
          ylabel="relative SF error", title="SF error vs cycles-in-window\n(dashed = real RF sizes)")
    a.legend(fontsize=7)

    # TF vs window deg
    a = ax[1, 0]
    for tf in TF_LEVELS_HZ:
        rs = [r for r in rec if r[0] == "TF_sq" and r[1] == tf]
        wds = sorted({r[2] for r in rs})
        mu = [np.nanmean([r[4] for r in rs if r[2] == w]) for w in wds]
        sd = [np.nanstd([r[4] for r in rs if r[2] == w]) for w in wds]
        a.errorbar(wds, mu, yerr=sd, marker="s", capsize=3, label=f"{tf:.0f} Hz")
        a.axhline(tf, ls=":", lw=0.8, color="gray")
    a.set(xscale="log", xlabel="square window (deg)",
          ylabel="recovered TF (Hz)", title="TF vs window size (expect flat)")
    a.legend(fontsize=7)

    # OR error vs window deg
    a = ax[1, 1]
    rs = [r for r in rec if r[0] == "OR_sq"]
    wds = sorted({r[2] for r in rs})
    er = [np.nanmean([circ_err_deg(r[4], r[1]) for r in rs if r[2] == w]) for w in wds]
    a.plot(wds, er, "-^", color="purple")
    a.set(xscale="log", xlabel="square window (deg)",
          ylabel="orientation error (deg)", title="OR vs window size")

    fig.suptitle(f"Windowed cloud-stat recovery — {cfg.profile.name} "
                 f"({cfg.profile.deg_per_px:.3f} deg/px, canvas {nx}x{ny}px)")
    fig.tight_layout()
    out = os.path.join(outdir, f"windowed_recovery_{cfg.profile.name}.png")
    fig.savefig(out, dpi=140)
    print(f"[saved] {out}")

    # real-RF table
    print("\n=== SF recovery at REAL RF sizes (rel error, mean over reals) ===")
    rf_rows = sorted({r[2] for r in rec if r[0] == "SF_rf"})
    print(f"{'RF az x el':>12} | " + " | ".join(f"{s:.2f}cpd" for s in SF_LEVELS_CPD))
    for rf in rf_rows:
        cells = []
        for sf in SF_LEVELS_CPD:
            vals = [abs(r[4] - sf) / sf for r in rec
                    if r[0] == "SF_rf" and r[2] == rf and r[1] == sf]
            cells.append(f"{np.nanmean(vals):6.2f}")
        cyc = [r[3] for r in rec if r[0] == "SF_rf" and r[2] == rf][0]
        print(f"{str(rf):>12} | " + " | ".join(cells) + f"   (min-axis cyc@0.12={0.12*min(rf):.2f})")


def load_rf_full(csv_glob):
    """Return list of dicts with az/el size (deg) and centroid az/el (deg)."""
    out = []
    for fp in glob.glob(csv_glob):
        with open(fp) as fh:
            for r in csv.DictReader(fh):
                try:
                    out.append(dict(
                        az=float(r["size_azimuth_deg"]), el=float(r["size_elevation_deg"]),
                        caz=float(r["centroid_azimuth_deg"]), cel=float(r["centroid_elevation_deg"]),
                        rf_type=r.get("rf_type", "")))
                except (KeyError, ValueError):
                    continue
    return out


def _crop(arr, cy, cx, hy, hx):
    """Crop arr[..., y, x] of half-extents (hy,hx) about (cy,cx), clamped."""
    y0, y1 = max(0, cy - hy), min(arr.shape[-2], cy + hy)
    x0, x1 = max(0, cx - hx), min(arr.shape[-1], cx + hx)
    return arr[..., y0:y1, x0:x1]


def plot_examples(cfg, outdir, rf_scale=1.0):
    """Diagnostic: example cloud + RF mask + local-property fluctuation.

    Reframes around Laura's model: a local patch is a realisation with its OWN
    dominant OR/SF/TF fluctuating around the global envelope means. We show the
    cloud, the RF window, a spatial map of local orientation, and the spread of
    locally-measured OR/SF/TF at real RF sizes (spread = fluctuation, not error).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_theta, g_sf, g_tf = np.pi / 4, 0.06, 2.0  # global params of the example

    # one example cloud
    mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                    g_theta, B_THETA_RAD, g_tf, B_V, rng)
    mid = mv[cfg.n_frame // 2]

    # a representative real RF (median-ish 15x15) at its real centroid
    rfs = load_rf_full(cfg.rf_csv_glob)
    rep = next((r for r in rfs if r["az"] == 15 and r["el"] == 15), None)
    if rep is None:
        rep = dict(az=15, el=15, caz=cfg.profile.fov_w_deg / 2,
                   cel=cfg.profile.fov_h_deg / 2)
    cx = int(round(rep["caz"] / dpp)); cy = int(round(rep["cel"] / dpp))
    cx = min(max(cx, 0), nx - 1); cy = min(max(cy, 0), ny - 1)
    az_deg = rep["az"] * rf_scale; el_deg = rep["el"] * rf_scale  # scaled extents
    ax_px = az_deg / dpp; ay_px = el_deg / dpp

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # A: example cloud + RF ellipse
    a = ax[0, 0]
    a.imshow(mid, cmap="gray", extent=[0, nx * dpp, ny * dpp, 0])
    a.add_patch(Ellipse((rep["caz"], rep["cel"]), az_deg, el_deg,
                        fill=False, ec="red", lw=2))
    a.set(title=f"example cloud frame (theta={np.degrees(g_theta):.0f}°, "
          f"SF={g_sf} cpd)\nred = RF {az_deg:.0f}x{el_deg:.0f}° (x{rf_scale:g}) @"
          f"({rep['caz']:.0f},{rep['cel']:.0f})°", xlabel="azimuth (deg)",
          ylabel="elevation (deg)")

    # B: windowed patch (cloud x elliptical gaussian at the RF centroid), zoomed
    a = ax[0, 1]
    mask = elliptical_mask(ny, nx, ax_px, ay_px, "gauss", cx=cx, cy=cy)
    patch = _crop(mid * mask, cy, cx, int(ay_px), int(ax_px))  # ~2x RF for context
    a.imshow(patch, cmap="gray")
    a.set(title="what the RF sees (cloud x Gaussian RF mask)", xticks=[], yticks=[])

    # C: local-orientation map (sliding RF-min-axis window)
    a = ax[0, 2]
    win = max(8, int(round(min(az_deg, el_deg) / dpp)))
    step = max(4, win // 2)
    ys = range(win, ny - win, step); xs = range(win, nx - win, step)
    ormap = np.full((len(list(ys)), len(list(xs))), np.nan)
    for iy, yc in enumerate(range(win, ny - win, step)):
        for ix, xc in enumerate(range(win, nx - win, step)):
            p = _crop(mid, yc, xc, win // 2, win // 2)
            ormap[iy, ix] = np.degrees(recover_or_rad(p))
    im = a.imshow(ormap, cmap="twilight", vmin=-90, vmax=90,
                  extent=[0, nx * dpp, ny * dpp, 0], aspect="auto")
    plt.colorbar(im, ax=a, label="local OR (deg)")
    a.set(title=f"local orientation map ({min(az_deg,el_deg):.0f}° window)\n"
          f"global theta={np.degrees(g_theta):.0f}° — note fluctuation",
          xlabel="azimuth (deg)", ylabel="elevation (deg)")

    # sample local estimates at two RF sizes across several realisations
    def sample_local(stat, az_deg, el_deg, n_clouds=3):
        vals = []
        for _ in range(n_clouds):
            m2 = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                            g_theta, B_THETA_RAD, g_tf, B_V, rng)
            f2 = m2[cfg.n_frame // 2]
            hx, hy = int(az_deg / dpp / 2), int(el_deg / dpp / 2)
            for yc in range(hy + 2, ny - hy - 2, max(6, hy)):
                for xc in range(hx + 2, nx - hx - 2, max(6, hx)):
                    if stat == "OR":
                        vals.append(np.degrees(recover_or_rad(_crop(f2, yc, xc, hy, hx))))
                    elif stat == "SF":
                        vals.append(recover_sf_cpd(_crop(f2, yc, xc, hy, hx), dpp))
                    else:  # TF
                        vals.append(recover_tf_hz(_crop(m2, yc, xc, hy, hx), cfg.fps))
        return np.array(vals)

    small = (5.0 * rf_scale, 10.0 * rf_scale); med = (15.0 * rf_scale, 15.0 * rf_scale)
    # D: local OR distribution
    a = ax[1, 0]
    a.hist(sample_local("OR", *small), bins=30, alpha=0.5, density=True,
           label=f"{small[0]:.0f}x{small[1]:.0f}°")
    a.hist(sample_local("OR", *med), bins=30, alpha=0.5, density=True,
           label=f"{med[0]:.0f}x{med[1]:.0f}°")
    a.axvline(np.degrees(g_theta), color="k", lw=2, label="global theta")
    a.axvspan(np.degrees(g_theta - B_THETA_RAD), np.degrees(g_theta + B_THETA_RAD),
              color="gray", alpha=0.15, label="±B_theta")
    a.set(title="local OR distribution (spread = fluctuation)",
          xlabel="local OR (deg)", ylabel="density"); a.legend(fontsize=7)

    # E: local SF distribution
    a = ax[1, 1]
    a.hist(sample_local("SF", *small), bins=30, alpha=0.5, density=True,
           label=f"{small[0]:.0f}x{small[1]:.0f}°")
    a.hist(sample_local("SF", *med), bins=30, alpha=0.5, density=True,
           label=f"{med[0]:.0f}x{med[1]:.0f}°")
    a.axvline(g_sf, color="k", lw=2, label="global SF")
    a.set(title="local SF distribution (peak)\ncollapses to window floor at small RF",
          xlabel="local SF (cpd)", ylabel="density"); a.legend(fontsize=7)

    # F: local TF distribution
    a = ax[1, 2]
    a.hist(sample_local("TF", *small), bins=20, alpha=0.5, density=True,
           label=f"{small[0]:.0f}x{small[1]:.0f}°")
    a.hist(sample_local("TF", *med), bins=20, alpha=0.5, density=True,
           label=f"{med[0]:.0f}x{med[1]:.0f}°")
    a.axvline(g_tf, color="k", lw=2, label="global TF")
    a.set(title="local TF distribution (tight around global)",
          xlabel="local TF (Hz)", ylabel="density"); a.legend(fontsize=7)

    fig.suptitle(f"Cloud / RF / local-fluctuation diagnostic — {cfg.profile.name} "
                 f"({dpp:.3f} deg/px, RF x{rf_scale:g})")
    fig.tight_layout()
    suffix = "" if rf_scale == 1.0 else f"_rfx{rf_scale:g}"
    out = os.path.join(outdir, f"cloud_rf_local_diagnostic_{cfg.profile.name}{suffix}.png")
    fig.savefig(out, dpi=140)
    print(f"[saved] {out}")


def plot_timecourse(cfg, outdir, duration_s=1.0, rf_scales=(1.0, 2.0)):
    """Per-frame local SF/OR/TF over ~duration_s, at a fixed RF, for two radii.

    Dashed line = the global (design) mean of each parameter. Solid lines = the
    local windowed estimate per frame for each RF radius — i.e. the fluctuation
    that would become the offset regressor. TF(t) uses a sliding temporal window.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_theta, g_sf, g_tf = np.pi / 4, 0.06, 2.0

    mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                    g_theta, B_THETA_RAD, g_tf, B_V, rng)

    rfs = load_rf_full(cfg.rf_csv_glob)
    rep = next((r for r in rfs if r["az"] == 15 and r["el"] == 15), None)
    if rep is None:
        rep = dict(az=15, el=15, caz=cfg.profile.fov_w_deg / 2,
                   cel=cfg.profile.fov_h_deg / 2)
    cx = min(max(int(round(rep["caz"] / dpp)), 0), nx - 1)
    cy = min(max(int(round(rep["cel"] / dpp)), 0), ny - 1)

    n_show = min(cfg.n_frame, int(round(duration_s * cfg.fps)))
    times = np.arange(n_show) / cfg.fps
    tf_half = min(cfg.n_frame, 32) // 2   # ~1 s sliding window for TF(t)
    start = min(tf_half, max(0, cfg.n_frame - n_show))  # center: full TF windows

    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    colors = {1.0: "C0", 2.0: "C1"}
    for scale in rf_scales:
        az_px = rep["az"] * scale / dpp; el_px = rep["el"] * scale / dpp
        mask = elliptical_mask(ny, nx, az_px, el_px, "gauss", cx=cx, cy=cy)
        mvm = mv * mask                      # masked movie (nf, ny, nx)
        hy, hx = int(el_px), int(az_px)      # ~2x RF bbox for context
        sf_t, or_t, tf_t = [], [], []
        for k in range(n_show):
            i = start + k
            patch = _crop(mvm[i], cy, cx, hy, hx)
            sf_t.append(recover_sf_cpd(patch, dpp))
            or_t.append(np.degrees(recover_or_rad(patch)))
            t0, t1 = max(0, i - tf_half), min(cfg.n_frame, i + tf_half)
            tf_t.append(recover_tf_hz(_crop(mvm[t0:t1], cy, cx, hy, hx), cfg.fps))
        lbl = f"RF {rep['az']*scale:.0f}x{rep['el']*scale:.0f}° (x{scale:g})"
        c = colors.get(scale, None)
        ax[0].plot(times, sf_t, "-o", ms=3, color=c, label=lbl)
        ax[1].plot(times, or_t, "-o", ms=3, color=c, label=lbl)
        ax[2].plot(times, tf_t, "-o", ms=3, color=c, label=lbl)

    ax[0].axhline(g_sf, ls="--", color="k", label="global mean (design)")
    ax[0].set(ylabel="local SF (cpd)", title="SF(t)")
    ax[1].axhline(np.degrees(g_theta), ls="--", color="k")
    ax[1].set(ylabel="local OR (deg)", title="OR(t)")
    ax[2].axhline(g_tf, ls="--", color="k")
    ax[2].set(ylabel="local TF (Hz)", xlabel="time (s)",
              title=f"TF(t) [sliding {2*tf_half} frames ≈ {2*tf_half/cfg.fps:.1f}s; "
              f"flat = constant sim drift]")
    ax[0].legend(fontsize=8, ncol=3)
    fig.suptitle(f"Per-frame local SF/OR/TF — {cfg.profile.name} "
                 f"({dpp:.3f} deg/px), one RF @({rep['caz']:.0f},{rep['cel']:.0f})°, "
                 f"{duration_s:g}s @{cfg.fps:g}Hz")
    fig.tight_layout()
    out = os.path.join(outdir, f"local_timecourse_{cfg.profile.name}.png")
    fig.savefig(out, dpi=140)
    print(f"[saved] {out}")


def test_one_frame(cfg, outdir, rf_scales=(1.0, 2.0), thresh_q=0.65):
    """Single-frame head-to-head: blob estimator vs FFT/gradient vs truth.

    Hard-crops the RF (no soft mask) so the blob method isn't fighting the
    aperture's low-freq blob; centrality weighting plays the 'within-RF' role.
    Saves a figure with detected blob ellipses overlaid so we can eyeball it.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_theta, g_sf, g_tf = np.pi / 4, 0.06, 2.0

    mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                    g_theta, B_THETA_RAD, g_tf, B_V, rng)
    frame = mv[cfg.n_frame // 2]
    cx, cy = nx // 2, ny // 2

    print(f"\n=== one-frame test | truth: SF={g_sf} cpd, theta={np.degrees(g_theta):.0f}° ===")
    print(f"{'RF':>10} | {'FFT SF':>8} {'blob SFsp':>9} | {'grad OR':>8} "
          f"{'blob OR':>8} | {'blob BW':>8} {'nblob':>6}")

    fig, axes = plt.subplots(len(rf_scales), 2, figsize=(11, 5 * len(rf_scales)))
    if len(rf_scales) == 1:
        axes = axes[None, :]
    for ri, scale in enumerate(rf_scales):
        az = 15.0 * scale; el = 15.0 * scale
        hx = max(6, int(az / dpp / 2)); hy = max(6, int(el / dpp / 2))
        patch = _crop(frame, cy, cx, hy, hx)
        fft_sf = recover_sf_cpd(patch, dpp)
        grad_or = np.degrees(recover_or_rad(patch))
        b = blob_estimate(patch, dpp, thresh_q=thresh_q)
        print(f"{az:>5.0f}x{el:<4.0f} | {fft_sf:>8.3f} {b['sf_spacing']:>9.3f} | "
              f"{grad_or:>8.1f} {b['or_deg']:>8.1f} | "
              f"{b['bw_deg']:>8.1f} {b['n_blobs']:>6d}")

        # left: patch + bright(cyan)/dark(orange) ellipses with major-axis lines
        a = axes[ri, 0]
        a.imshow(patch, cmap="gray", origin="upper")
        for bl in b["blobs"]:
            ec = "cyan" if bl["pol"] > 0 else "orange"
            a.add_patch(Ellipse((bl["cx"], bl["cy"]), width=2 * bl["major"],
                                height=2 * bl["minor"], angle=-np.degrees(bl["theta"]),
                                fill=False, ec=ec, lw=1.3))
            a.plot([bl["cx"], bl["cx"] + bl["major"] * np.cos(bl["theta"])],
                   [bl["cy"], bl["cy"] - bl["major"] * np.sin(bl["theta"])],
                   color=ec, lw=1)
        a.set(title=f"RF {az:.0f}x{el:.0f}° + ellipses (cyan=bright, orange=dark)\n"
              f"blob OR(major-axis avg)={b['or_deg']:.0f}°  grad OR={grad_or:.0f}°  "
              f"truth {np.degrees(g_theta):.0f}°", xticks=[], yticks=[])
        # right: bright/dark binarised composite
        a = axes[ri, 1]
        comp = np.zeros((*patch.shape, 3))
        comp[b["bright_mask"]] = [0.2, 0.8, 0.9]   # cyan = bright stripes
        comp[b["dark_mask"]] = [0.95, 0.6, 0.1]    # orange = dark stripes
        a.imshow(comp, origin="upper")
        a.set(title=f"bright/dark masks (q={thresh_q}) {b['n_blobs']} blobs\n"
              f"SF spacing={b['sf_spacing']:.3f}  FFT={fft_sf:.3f}  truth {g_sf}",
              xticks=[], yticks=[])

    fig.suptitle(f"One-frame blob(bright+dark) vs FFT/gradient — {cfg.profile.name} "
                 f"({dpp:.3f} deg/px). OR = major-axis avg; SF = opposite-phase spacing.")
    fig.tight_layout()
    out = os.path.join(outdir, f"oneframe_blob_test_{cfg.profile.name}.png")
    fig.savefig(out, dpi=140)
    print(f"[saved] {out}")


def _circ_smooth(deg, win):
    """Circular moving-average of an orientation series (axial, mod 180°)."""
    if win < 2:
        return np.asarray(deg, float)
    z = np.exp(1j * 2 * np.radians(deg))
    k = np.ones(win) / win
    zs = np.convolve(z, k, mode="same")
    return np.degrees(np.angle(zs) / 2)


def _lin_smooth_nan(x, win):
    """Centred moving-average ignoring NaNs (short series, loop is fine)."""
    x = np.asarray(x, float); n = len(x); h = win // 2
    out = np.full(n, np.nan)
    for i in range(n):
        seg = x[max(0, i - h):i + h + 1]
        seg = seg[~np.isnan(seg)]
        if seg.size:
            out[i] = seg.mean()
    return out


def _circ_smooth_nan(deg, win):
    """Centred circular (axial, mod 180°) moving-average ignoring NaNs."""
    z = np.exp(1j * 2 * np.radians(np.asarray(deg, float)))
    n = len(z); h = win // 2
    out = np.full(n, np.nan)
    for i in range(n):
        seg = z[max(0, i - h):i + h + 1]
        seg = seg[~np.isnan(seg)]
        if seg.size:
            out[i] = np.degrees(np.angle(seg.mean()) / 2)
    return out


def compare_or_methods(cfg, outdir, duration_s=1.5, rf_scales=(1.0, 2.0),
                       thresh_q=0.65, smooth_win=5):
    """OR(t): blob major-axis-average (raw + time-smoothed) vs gradient, two radii.

    Tests Laura's claim that the blob major-axis OR is the better measure and
    just needs temporal smoothing — head-to-head with the gradient/structure OR,
    against the design mean (dashed). New file; nothing overwritten.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_theta, g_sf, g_tf = np.pi / 4, 0.06, 2.0
    mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                    g_theta, B_THETA_RAD, g_tf, B_V, rng)
    cx, cy = nx // 2, ny // 2
    n_show = min(cfg.n_frame, int(round(duration_s * cfg.fps)))
    times = np.arange(n_show) / cfg.fps

    fig, axes = plt.subplots(len(rf_scales), 1, figsize=(11, 4.2 * len(rf_scales)),
                             sharex=True, squeeze=False)
    print(f"\n=== OR method comparison | truth theta={np.degrees(g_theta):.0f}° "
          f"(circ-err deg, mean over {n_show} frames) ===")
    for ri, scale in enumerate(rf_scales):
        hx = max(6, int(15.0 * scale / dpp / 2)); hy = hx
        blob_or, grad_or = [], []
        for i in range(n_show):
            patch = _crop(mv[i], cy, cx, hy, hx)
            blob_or.append(blob_estimate(patch, dpp, thresh_q=thresh_q)["or_deg"])
            grad_or.append(np.degrees(recover_or_rad(patch)))
        blob_or = np.array(blob_or); grad_or = np.array(grad_or)
        blob_sm = _circ_smooth(blob_or, smooth_win)
        truth = np.degrees(g_theta)

        def cerr(x):
            d = np.angle(np.exp(1j * 2 * (np.radians(x) - g_theta))) / 2
            return np.nanmean(np.abs(np.degrees(d)))
        print(f"RF {15*scale:.0f}°: blob raw={cerr(blob_or):4.1f}  "
              f"blob smooth={cerr(blob_sm):4.1f}  gradient={cerr(grad_or):4.1f}")

        a = axes[ri, 0]
        a.plot(times, blob_or, "-", color="C0", alpha=0.35, label="blob major-axis (raw)")
        a.plot(times, blob_sm, "-o", color="C0", ms=3, label=f"blob smoothed ({smooth_win}f)")
        a.plot(times, grad_or, "-s", color="C1", ms=3, alpha=0.8, label="gradient/structure")
        a.axhline(truth, ls="--", color="k", label="design theta")
        a.set(ylabel="local OR (deg)",
              title=f"RF {15*scale:.0f}×{15*scale:.0f}°  (truth {truth:.0f}°)")
        a.legend(fontsize=8, ncol=2)
    axes[-1, 0].set_xlabel("time (s)")
    fig.suptitle(f"OR estimator comparison: blob major-axis vs gradient — "
                 f"{cfg.profile.name} ({dpp:.3f} deg/px), {duration_s:g}s")
    fig.tight_layout()
    out = os.path.join(outdir, f"cloud_rf_local_diagnostic_{cfg.profile.name}_ORcompare.png")
    fig.savefig(out, dpi=140)
    print(f"[saved] {out}")


def animate_patch(cfg, outdir, rf_scale=1.0, duration_s=1.0, thresh_q=0.65,
                  fps_play=6, smooth_win=5, sf_smooth_win=9):
    """Animated GIF: clean patch | patch+fitted ellipses | SF/OR time courses.

    For judging whether the per-frame fluctuations are REAL: watch whether the
    fitted ellipses track the drifting cloud structure (real) or jump around
    (noise). A moving red bar marks the current frame on the time courses.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    from matplotlib.animation import FuncAnimation, PillowWriter

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_theta, g_sf, g_tf = np.pi / 4, 0.06, 2.0
    mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                    g_theta, B_THETA_RAD, g_tf, B_V, rng)
    cx, cy = nx // 2, ny // 2
    n = min(cfg.n_frame, int(round(duration_s * cfg.fps)))
    times = np.arange(n) / cfg.fps
    hx = hy = max(6, int(15.0 * rf_scale / dpp / 2))

    patches, blobs, or_b, or_g, sf_sp = [], [], [], [], []
    for i in range(n):
        p = _crop(mv[i], cy, cx, hy, hx)
        b = blob_estimate(p, dpp, thresh_q=thresh_q)
        patches.append(p); blobs.append(b)
        or_b.append(b["or_deg"]); sf_sp.append(b["sf_spacing"])
        or_g.append(np.degrees(recover_or_rad(p)))

    pstack = np.stack(patches)
    vlo, vhi = float(pstack.min()), float(pstack.max())   # fixed contrast
    ph, pw = patches[0].shape
    or_b_sm = _circ_smooth_nan(or_b, smooth_win)          # smoothed blob OR
    sf_sp_sm = _lin_smooth_nan(sf_sp, sf_smooth_win)      # smoothed blob SF (heavier)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.6))
    im1 = ax1.imshow(patches[0], cmap="gray", vmin=vlo, vmax=vhi)
    ax1.set(xticks=[], yticks=[])

    # static time courses on ax3 (+ twin for SF), moving bar updated each frame
    ax3.plot(times, or_b, "-", color="C0", alpha=0.3)
    ax3.plot(times, or_b_sm, "-o", color="C0", ms=3, label="blob OR (smoothed)")
    ax3.axhline(np.degrees(g_theta), ls="--", color="k", lw=1, label="θ (design mean)")
    ax3.set(xlabel="time (s)", ylabel="OR (deg)")
    ax3.legend(loc="upper left", fontsize=7)
    ax3b = ax3.twinx()
    ax3b.plot(times, sf_sp, "-", color="C2", alpha=0.3)
    ax3b.plot(times, sf_sp_sm, "-^", color="C2", ms=3, label="blob SF (smoothed)")
    ax3b.axhline(g_sf, ls=":", color="C2", lw=1)
    ax3b.set_ylabel("SF (cpd)", color="C2"); ax3b.tick_params(axis="y", colors="C2")
    ax3b.legend(loc="upper right", fontsize=7)
    vline = ax3.axvline(times[0], color="red", lw=2)

    def update(i):
        im1.set_data(patches[i])
        ax1.set_title(f"cloud patch (RF {15*rf_scale:.0f}°)  t={times[i]:.2f}s")
        ax2.clear()
        comp = np.full(patches[i].shape, 0.5)          # neither = mid gray
        comp[blobs[i]["bright_mask"]] = 1.0            # bright blob = white
        comp[blobs[i]["dark_mask"]] = 0.0             # dark blob   = black
        ax2.imshow(comp, cmap="gray", vmin=0, vmax=1)
        for bl in blobs[i]["blobs"]:
            ec = "cyan" if bl["pol"] > 0 else "orange"
            ax2.add_patch(Ellipse((bl["cx"], bl["cy"]), 2 * bl["major"],
                                  2 * bl["minor"], angle=-np.degrees(bl["theta"]),
                                  fill=False, ec=ec, lw=1.4))
            ax2.plot([bl["cx"], bl["cx"] + bl["major"] * np.cos(bl["theta"])],
                     [bl["cy"], bl["cy"] - bl["major"] * np.sin(bl["theta"])],
                     color=ec, lw=1)
        # pin view to the patch so it does not jump as ellipses spill over
        ax2.set_xlim(-0.5, pw - 0.5); ax2.set_ylim(ph - 0.5, -0.5)
        ax2.set_aspect("equal"); ax2.set_autoscale_on(False)
        ax2.set(xticks=[], yticks=[],
                title=f"binary (white=bright, black=dark) + ellipses\n"
                      f"OR={or_b_sm[i]:.0f}°  SF={sf_sp_sm[i]:.3f} cpd")
        vline.set_xdata([times[i], times[i]])
        return [im1]

    anim = FuncAnimation(fig, update, frames=n, interval=1000 / fps_play, blit=False)
    fig.suptitle(f"Cloud patch + blob ellipses over {duration_s:g}s — "
                 f"{cfg.profile.name} ({dpp:.3f} deg/px). truth: SF={g_sf}, θ=45°",
                 y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = os.path.join(outdir, f"patch_animation_{cfg.profile.name}_rf{15*rf_scale:.0f}.gif")
    anim.save(out, writer=PillowWriter(fps=fps_play), dpi=90)
    plt.close(fig)
    print(f"[saved] {out}  ({n} frames @ {fps_play} fps)")


def animate_gabor(cfg, outdir, rf_scale=1.0, duration_s=1.0, thresh_q=0.65,
                  fps_play=6, smooth_win=5, sf_smooth_win=9):
    """Animated GIF using THEIR (Gabor-wavelet) estimator on the same patch.

    Panels: clean patch | orientation×SF Gabor-energy map (peak = best wavelet) |
    OR(t)/SF(t) from Gabor (bold) with the blob estimates overlaid faint for
    direct comparison, design means dashed, moving time bar. Same seed/RF as the
    blob GIF so the patches are identical.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_theta, g_sf, g_tf = np.pi / 4, 0.06, 2.0
    mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                    g_theta, B_THETA_RAD, g_tf, B_V, rng)
    cx, cy = nx // 2, ny // 2
    n = min(cfg.n_frame, int(round(duration_s * cfg.fps)))
    times = np.arange(n) / cfg.fps
    hx = hy = max(6, int(15.0 * rf_scale / dpp / 2))
    pshape = _crop(mv[0], cy, cx, hy, hx).shape
    bank, ors, sfs = build_gabor_bank(pshape[0], pshape[1], dpp)

    patches, energies, or_g, sf_g, or_b, sf_b = [], [], [], [], [], []
    for i in range(n):
        p = _crop(mv[i], cy, cx, hy, hx)
        patches.append(p)
        gb = gabor_response(p, bank, ors, sfs)
        energies.append(gb["energy"])
        or_g.append(gb["or_bar_deg"]); sf_g.append(gb["sf_cpd"])
        bb = blob_estimate(p, dpp, thresh_q=thresh_q)
        or_b.append(bb["or_deg"]); sf_b.append(bb["sf_spacing"])
    vlo, vhi = float(np.min(patches)), float(np.max(patches))
    emax = float(np.max(energies))
    ph, pw = pshape
    or_g_sm = _circ_smooth_nan(or_g, smooth_win); sf_g_sm = _lin_smooth_nan(sf_g, sf_smooth_win)
    or_b_sm = _circ_smooth_nan(or_b, smooth_win); sf_b_sm = _lin_smooth_nan(sf_b, sf_smooth_win)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.6))
    im1 = ax1.imshow(patches[0], cmap="gray", vmin=vlo, vmax=vhi)
    ax1.set(xticks=[], yticks=[])

    ax3.plot(times, or_g_sm, "-o", color="C3", ms=3, label="Gabor OR")
    ax3.plot(times, or_b_sm, "--", color="C0", alpha=0.7, label="blob OR")
    ax3.axhline(np.degrees(g_theta), ls="--", color="k", lw=1, label="θ (design mean)")
    ax3.set(xlabel="time (s)", ylabel="OR (deg)")
    ax3.legend(loc="upper left", fontsize=7)
    ax3b = ax3.twinx()
    ax3b.plot(times, sf_g_sm, "-^", color="C2", ms=3, label="Gabor SF")
    ax3b.plot(times, sf_b_sm, ":", color="C2", alpha=0.6, label="blob SF")
    ax3b.axhline(g_sf, ls=":", color="C2", lw=1)
    ax3b.set_ylabel("SF (cpd)", color="C2"); ax3b.tick_params(axis="y", colors="C2")
    ax3b.legend(loc="upper right", fontsize=7)
    vline = ax3.axvline(times[0], color="red", lw=2)

    or_ticks = np.arange(0, len(ors), max(1, len(ors) // 6))
    sf_ticks = np.arange(0, len(sfs), max(1, len(sfs) // 6))

    def update(i):
        im1.set_data(patches[i])
        ax1.set_title(f"cloud patch (RF {15*rf_scale:.0f}°)  t={times[i]:.2f}s")
        ax2.clear()
        ax2.imshow(energies[i].T, origin="lower", aspect="auto", cmap="magma",
                   vmin=0, vmax=emax)
        pk = np.unravel_index(np.argmax(energies[i]), energies[i].shape)
        ax2.plot(pk[0], pk[1], "co", ms=9, mfc="none", mew=2)
        ax2.set_xticks(or_ticks); ax2.set_xticklabels([f"{np.degrees(ors[k]):.0f}" for k in or_ticks])
        ax2.set_yticks(sf_ticks); ax2.set_yticklabels([f"{sfs[k]:.3f}" for k in sf_ticks])
        ax2.set(xlabel="modulation θ (deg)", ylabel="SF (cpd)",
                title=f"Gabor energy map — best wavelet\nOR(bar)={or_g[i]:.0f}°  SF={sf_g[i]:.3f} cpd")
        vline.set_xdata([times[i], times[i]])
        return [im1]

    anim = FuncAnimation(fig, update, frames=n, interval=1000 / fps_play, blit=False)
    fig.suptitle(f"Gabor-wavelet estimator (Skriabine/Shinn) on the same patch — "
                 f"{cfg.profile.name} ({dpp:.3f} deg/px). truth: SF={g_sf}, θ=45°",
                 y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = os.path.join(outdir, f"gabor_animation_{cfg.profile.name}_rf{15*rf_scale:.0f}.gif")
    anim.save(out, writer=PillowWriter(fps=fps_play), dpi=90)
    plt.close(fig)
    print(f"[saved] {out}  ({n} frames @ {fps_play} fps)")


def _fullframe_sf_density(frame, dpp, edges):
    """Radial power spectrum of a frame, binned into SF edges (cpd)."""
    f = np.fft.fftshift(np.fft.fft2(frame - frame.mean()))
    mag2 = np.abs(f) ** 2
    ny, nx = frame.shape
    yy, xx = np.indices(frame.shape)
    rr = np.sqrt(((xx - nx // 2) / nx) ** 2 + ((yy - ny // 2) / ny) ** 2)
    freqs_cpd = (rr / dpp).ravel()
    pw, _ = np.histogram(freqs_cpd, bins=edges, weights=mag2.ravel())
    return pw


def compare_sf_distributions(cfg, outdir, rf_scale=1.0, n_real=4, n_frame=96,
                             grid=3, thresh_q=0.65, sf_max=0.25, gabor_nsf=14):
    """Pool ellipse-SF and Gabor-SF over many patches; compare to the stimulus
    SF distribution (full-frame radial power spectrum). Wasserstein = closeness.
    """
    from scipy.stats import wasserstein_distance

    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_sf = 0.06
    hx = hy = max(6, int(15.0 * rf_scale / dpp / 2))
    pshape = (2 * hy, 2 * hx)
    bank, ors, sfs = build_gabor_bank(pshape[0], pshape[1], dpp, n_sf=gabor_nsf)
    GE = np.stack([b[0] for b in bank]); GO = np.stack([b[1] for b in bank])
    n_or, n_sf = len(ors), len(sfs)

    # patch-centre grid (avoid edges)
    ycs = np.linspace(hy + 2, ny - hy - 2, grid).astype(int)
    xcs = np.linspace(hx + 2, nx - hx - 2, grid).astype(int)
    edges = np.linspace(0, sf_max, 80)
    centers = 0.5 * (edges[:-1] + edges[1:])

    ell_sf, gab_sf, ref_pw = [], [], np.zeros(len(centers))
    for _ in range(n_real):
        mv = make_cloud(nx, ny, n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                        np.pi / 4, B_THETA_RAD, 2.0, B_V, rng)
        for fi in range(n_frame):
            frame = mv[fi]
            ref_pw += _fullframe_sf_density(frame, dpp, edges)
            for yc in ycs:
                for xc in xcs:
                    p = _crop(frame, yc, xc, hy, hx)
                    s = blob_estimate(p, dpp, thresh_q=thresh_q)["sf_spacing"]
                    if np.isfinite(s):
                        ell_sf.append(s)
                    pdc = p - p.mean()
                    re = np.einsum("kij,ij->k", GE, pdc)
                    im = np.einsum("kij,ij->k", GO, pdc)
                    en = (re * re + im * im).reshape(n_or, n_sf)
                    gab_sf.append(float(sfs[np.argmax(en) % n_sf]))

    ell_sf = np.array(ell_sf); gab_sf = np.array(gab_sf)
    ref_density = ref_pw / ref_pw.sum()
    # expected samples ∝ stimulus power, restricted to the analysed SF range
    exp_samp = rng.choice(centers, size=20000, p=ref_density)
    ref_mean = float(np.sum(centers * ref_density))
    w_ell = wasserstein_distance(ell_sf, exp_samp)
    w_gab = wasserstein_distance(gab_sf, exp_samp)
    winner = "ELLIPSE" if w_ell < w_gab else "GABOR"

    print(f"\n=== SF distribution vs stimulus (full-frame) | RF {15*rf_scale:.0f}°, "
          f"{len(ell_sf)} patches ===")
    print(f"expected (power-wtd) mean SF = {ref_mean:.4f} cpd  (design sf_0={g_sf})")
    print(f"ellipse: mean={ell_sf.mean():.4f} median={np.median(ell_sf):.4f}  "
          f"Wasserstein={w_ell:.4f}")
    print(f"gabor  : mean={gab_sf.mean():.4f} median={np.median(gab_sf):.4f}  "
          f"Wasserstein={w_gab:.4f}")
    print(f"closest to expected: {winner}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.fill_between(centers, ref_density / np.diff(edges)[0], color="0.8",
                    label="stimulus SF power (full-frame, expected)")
    ax.hist(ell_sf, bins=edges, density=True, histtype="step", lw=2, color="C0",
            label=f"ellipse (spacing)  W={w_ell:.4f}")
    ax.hist(gab_sf, bins=edges, density=True, histtype="step", lw=2, color="C3",
            label=f"Gabor (best wavelet)  W={w_gab:.4f}")
    ax.axvline(g_sf, ls="--", color="k", label=f"design sf₀={g_sf}")
    ax.axvline(ref_mean, ls=":", color="0.4", label=f"expected mean={ref_mean:.3f}")
    ax.set(xlabel="spatial frequency (cpd)", ylabel="density", xlim=(0, sf_max),
           title=f"Pooled local-SF vs stimulus SF distribution — RF {15*rf_scale:.0f}° "
           f"({len(ell_sf)} patches). Closest: {winner}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(outdir, f"sf_distribution_compare_{cfg.profile.name}_rf{15*rf_scale:.0f}.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"[saved] {out}")


def plot_examples_gabor(cfg, outdir, rf_scale=1.0, gabor_nsf=40):
    """Local-fluctuation diagnostic using the OPTIMISED Gabor (fine SF bank,
    corrected OR convention) — the Gabor analogue of cloud_rf_local_diagnostic.

    Panels: example cloud + RF | windowed patch | local-OR map (Gabor) |
    local-SF map (Gabor) | OR distribution vs design | SF distribution vs design.
    New file; the gradient/blob diagnostic is left untouched.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    nx, ny = _canvas_px(cfg)
    dpp = cfg.profile.deg_per_px
    g_theta, g_sf = np.pi / 4, 0.06
    g_or_bar = (90 - np.degrees(g_theta)) % 180          # design bar orientation
    B_theta_deg = np.degrees(B_THETA_RAD)

    bank_cache = {}

    def get_bank(h, w):
        if (h, w) not in bank_cache:
            bk, ors, sfs = build_gabor_bank(h, w, dpp, n_sf=gabor_nsf)
            GE = np.stack([b[0] for b in bk]); GO = np.stack([b[1] for b in bk])
            bank_cache[(h, w)] = (GE, GO, ors, sfs)
        return bank_cache[(h, w)]

    def gabor_best(patch):
        GE, GO, ors, sfs = get_bank(*patch.shape)
        pdc = patch - patch.mean()
        re = np.einsum("kij,ij->k", GE, pdc); im = np.einsum("kij,ij->k", GO, pdc)
        en = (re * re + im * im).reshape(len(ors), len(sfs))
        i, j = np.unravel_index(en.argmax(), en.shape)
        return (90 - np.degrees(ors[i])) % 180, float(sfs[j])

    mv = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                    g_theta, B_THETA_RAD, 2.0, B_V, rng)
    mid = mv[cfg.n_frame // 2]
    rfs = load_rf_full(cfg.rf_csv_glob)
    rep = next((r for r in rfs if r["az"] == 15 and r["el"] == 15), None) or \
        dict(az=15, el=15, caz=cfg.profile.fov_w_deg / 2, cel=cfg.profile.fov_h_deg / 2)
    az_deg, el_deg = rep["az"] * rf_scale, rep["el"] * rf_scale
    cx = min(max(int(round(rep["caz"] / dpp)), 0), nx - 1)
    cy = min(max(int(round(rep["cel"] / dpp)), 0), ny - 1)
    ax_px, ay_px = az_deg / dpp, el_deg / dpp

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # A: cloud + RF
    a = ax[0, 0]
    a.imshow(mid, cmap="gray", extent=[0, nx * dpp, ny * dpp, 0])
    a.add_patch(Ellipse((rep["caz"], rep["cel"]), az_deg, el_deg, fill=False,
                        ec="red", lw=2))
    a.set(title=f"example cloud (θ={np.degrees(g_theta):.0f}°→bar {g_or_bar:.0f}°, "
          f"SF={g_sf})\nRF {az_deg:.0f}×{el_deg:.0f}°", xlabel="azimuth (deg)",
          ylabel="elevation (deg)")

    # B: windowed patch
    a = ax[0, 1]
    mask = elliptical_mask(ny, nx, ax_px, ay_px, "gauss", cx=cx, cy=cy)
    a.imshow(_crop(mid * mask, cy, cx, int(ay_px), int(ax_px)), cmap="gray")
    a.set(title="what the RF sees (cloud × Gaussian RF)", xticks=[], yticks=[])

    # C/D: local OR & SF maps (sliding Gabor window = min RF axis)
    win = max(10, int(round(min(az_deg, el_deg) / dpp)))
    step = max(5, win // 2)
    ys = list(range(win, ny - win, step)); xs = list(range(win, nx - win, step))
    ormap = np.full((len(ys), len(xs)), np.nan)
    sfmap = np.full((len(ys), len(xs)), np.nan)
    for iy, yc in enumerate(ys):
        for ix, xc in enumerate(xs):
            o, s = gabor_best(_crop(mid, yc, xc, win // 2, win // 2))
            ormap[iy, ix] = o; sfmap[iy, ix] = s
    a = ax[0, 2]
    im = a.imshow(ormap, cmap="twilight", vmin=0, vmax=180,
                  extent=[0, nx * dpp, ny * dpp, 0], aspect="auto")
    plt.colorbar(im, ax=a, label="local OR (deg)")
    a.set(title=f"Gabor local-OR map ({min(az_deg,el_deg):.0f}° win)\n"
          f"design bar OR={g_or_bar:.0f}°", xlabel="azimuth (deg)", ylabel="elev (deg)")
    a = ax[1, 0]
    im = a.imshow(sfmap, cmap="viridis", vmin=0.03, vmax=0.12,
                  extent=[0, nx * dpp, ny * dpp, 0], aspect="auto")
    plt.colorbar(im, ax=a, label="local SF (cpd)")
    a.set(title=f"Gabor local-SF map ({min(az_deg,el_deg):.0f}° win)\n"
          f"design SF={g_sf} cpd", xlabel="azimuth (deg)", ylabel="elev (deg)")

    # sample local Gabor estimates at two RF sizes
    def sample(az, el, n_clouds=3):
        ors_, sfs_ = [], []
        for _ in range(n_clouds):
            m2 = make_cloud(nx, ny, cfg.n_frame, cfg.fps, dpp, g_sf, B_SF_CPD,
                            g_theta, B_THETA_RAD, 2.0, B_V, rng)
            f2 = m2[cfg.n_frame // 2]
            hx, hy = int(az / dpp / 2), int(el / dpp / 2)
            for yc in range(hy + 2, ny - hy - 2, max(8, hy)):
                for xc in range(hx + 2, nx - hx - 2, max(8, hx)):
                    o, s = gabor_best(_crop(f2, yc, xc, hy, hx))
                    ors_.append(o); sfs_.append(s)
        return np.array(ors_), np.array(sfs_)

    small = (5.0 * rf_scale, 10.0 * rf_scale); med = (15.0 * rf_scale, 15.0 * rf_scale)
    os_s, sf_s = sample(*small); os_m, sf_m = sample(*med)

    # E: OR distribution
    a = ax[1, 1]
    a.hist(os_s, bins=30, range=(0, 180), alpha=0.5, density=True,
           label=f"{small[0]:.0f}×{small[1]:.0f}°")
    a.hist(os_m, bins=30, range=(0, 180), alpha=0.5, density=True,
           label=f"{med[0]:.0f}×{med[1]:.0f}°")
    a.axvline(g_or_bar, color="k", lw=2, label="design bar OR")
    a.axvspan(g_or_bar - B_theta_deg, g_or_bar + B_theta_deg, color="gray",
              alpha=0.15, label="±B_theta")
    a.set(title="Gabor local OR distribution", xlabel="local OR (deg)",
          ylabel="density"); a.legend(fontsize=7)

    # F: SF distribution
    a = ax[1, 2]
    a.hist(sf_s, bins=30, range=(0.02, 0.16), alpha=0.5, density=True,
           label=f"{small[0]:.0f}×{small[1]:.0f}°")
    a.hist(sf_m, bins=30, range=(0.02, 0.16), alpha=0.5, density=True,
           label=f"{med[0]:.0f}×{med[1]:.0f}°")
    a.axvline(g_sf, color="k", lw=2, label="design SF")
    a.axvspan(g_sf - B_SF_CPD, g_sf + B_SF_CPD, color="gray", alpha=0.15,
              label="±B_sf")
    a.set(title="Gabor local SF distribution", xlabel="local SF (cpd)",
          ylabel="density"); a.legend(fontsize=7)

    fig.suptitle(f"OPTIMISED GABOR local diagnostic — {cfg.profile.name} "
                 f"({dpp:.3f} deg/px, {gabor_nsf} SF levels). truth SF={g_sf}, bar OR={g_or_bar:.0f}°")
    fig.tight_layout()
    out = os.path.join(outdir, f"cloud_rf_local_diagnostic_{cfg.profile.name}_gabor.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"[saved] {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile", choices=list(PROFILES), default="goggles")
    p.add_argument("--outdir", default="./windowed_cloud_sim")
    p.add_argument("--rf-csv", default=DEFAULT_RF_CSV_GLOB)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--examples", action="store_true",
                   help="only the cloud/RF/local-fluctuation diagnostic (fast)")
    p.add_argument("--rf-scale", type=float, default=1.0,
                   help="multiply RF extents in the diagnostic (e.g. 2 = double radius)")
    p.add_argument("--timecourse", action="store_true",
                   help="per-frame local SF/OR/TF over --duration s, two RF radii")
    p.add_argument("--duration", type=float, default=1.0, help="timecourse length (s)")
    p.add_argument("--oneframe", action="store_true",
                   help="single-frame blob-estimator head-to-head vs FFT/gradient")
    p.add_argument("--thresh-q", type=float, default=0.65,
                   help="binarisation quantile for the blob estimator")
    p.add_argument("--compare-or", action="store_true",
                   help="OR(t): blob major-axis (raw+smoothed) vs gradient, two radii")
    p.add_argument("--animate", action="store_true",
                   help="animated GIF: patch | patch+ellipses | SF/OR time courses")
    p.add_argument("--animate-gabor", action="store_true",
                   help="animated GIF using the Gabor-wavelet estimator (their method)")
    p.add_argument("--sf-dist", action="store_true",
                   help="pool ellipse-SF vs Gabor-SF, compare to stimulus SF distribution")
    p.add_argument("--gabor-nsf", type=int, default=14,
                   help="number of SF levels in the Gabor bank (finer = less quantized)")
    p.add_argument("--examples-gabor", action="store_true",
                   help="local-fluctuation diagnostic using the optimised Gabor estimator")
    args = p.parse_args()

    cfg = Config(profile=PROFILES[args.profile], rf_csv_glob=args.rf_csv)
    if args.quick:
        cfg.n_frame, cfg.n_real = 48, 3
        cfg.profile = DisplayProfile(cfg.profile.name + "-quick",
                                     cfg.profile.deg_per_px * 1.6,  # coarser->fewer px
                                     cfg.profile.fov_w_deg, cfg.profile.fov_h_deg)
    nx, ny = _canvas_px(cfg)
    print(f"profile={cfg.profile.name} deg/px={cfg.profile.deg_per_px:.4f} "
          f"canvas={nx}x{ny}px ({cfg.profile.fov_w_deg:.0f}x{cfg.profile.fov_h_deg:.0f} deg) "
          f"frames={cfg.n_frame}@{cfg.fps}Hz n_real={cfg.n_real}")
    if args.animate:
        animate_patch(cfg, args.outdir, rf_scale=args.rf_scale,
                      duration_s=args.duration, thresh_q=args.thresh_q)
        return
    if args.animate_gabor:
        animate_gabor(cfg, args.outdir, rf_scale=args.rf_scale,
                      duration_s=args.duration, thresh_q=args.thresh_q)
        return
    if args.examples_gabor:
        nsf = args.gabor_nsf if args.gabor_nsf != 14 else 40
        plot_examples_gabor(cfg, args.outdir, rf_scale=args.rf_scale, gabor_nsf=nsf)
        return
    if args.sf_dist:
        compare_sf_distributions(cfg, args.outdir, rf_scale=args.rf_scale,
                                 thresh_q=args.thresh_q, gabor_nsf=args.gabor_nsf)
        return
    if args.compare_or:
        compare_or_methods(cfg, args.outdir, thresh_q=args.thresh_q)
        return
    if args.oneframe:
        test_one_frame(cfg, args.outdir, thresh_q=args.thresh_q)
        return
    if args.timecourse:
        plot_timecourse(cfg, args.outdir, duration_s=args.duration)
        return
    if args.examples:
        plot_examples(cfg, args.outdir, rf_scale=args.rf_scale)
        return
    rec = run_sweep(cfg)
    summarise_and_plot(rec, cfg, args.outdir)
    plot_examples(cfg, args.outdir, rf_scale=args.rf_scale)


if __name__ == "__main__":
    main()
