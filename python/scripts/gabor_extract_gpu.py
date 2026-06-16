"""GPU Gabor extraction of per-frame SF(t)/OR(t) for ALL RFs of ONE cloud.

Stimulus-only (probe/trial-independent): windows each RF on every cloud frame and
reads local SF (25° window) + OR (10° window, Viterbi-smoothed) with the validated
recipe. Designed as a SLURM job-array task — one cloud per `$SLURM_ARRAY_TASK_ID` —
streaming frames straight from ceph. torch on GPU (falls back to CPU); the einsum is
trivial, the bottleneck is PNG I/O, so the array parallelises that across nodes.

Output: one parquet per cloud, columns
    probe, cluster, rf_type, cx, cy, frame, sf_cpd, or_deg, concentration, edge_frac
The trial-clock velocity-warp (frame = 10·∫|v|dt) stays a separate cheap local step.
"""
from __future__ import annotations
import argparse, csv, glob, os, re
import numpy as np
from PIL import Image

DPP = 111.6 / 400.0
SF_WIN, OR_WIN, SF_MARGIN = 25, 10, 3
GABOR_NOR, GABOR_NSF = 36, 48
CLOUD_ROOT = "/ceph/margrie/mvelez/mateoData_mc/saved_goggles"
RF_DIR = "/ceph/margrie/laura/data transfer for laura/RFs googles"


def build_bank(win_deg, sf0, bsf):
    h = max(8, int(round(win_deg / DPP / 2))); n = 2 * h
    lo, hi = max(0.01, sf0 - SF_MARGIN * bsf), sf0 + SF_MARGIN * bsf
    cy = cx = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n].astype(np.float64); xx -= cx; yy -= cy
    env = np.exp(-(xx ** 2 + yy ** 2) / (2 * (0.35 * n) ** 2))
    ors = np.linspace(0, np.pi, GABOR_NOR, endpoint=False)
    sfs = np.geomspace(lo, hi, GABOR_NSF)
    GE, GO = [], []
    for th in ors:
        xr = xx * np.cos(th) + yy * np.sin(th)
        for f in sfs * DPP:
            ge = env * np.cos(2 * np.pi * f * xr); go = env * np.sin(2 * np.pi * f * xr)
            GE.append(ge - ge.mean()); GO.append(go)
    return h, np.asarray(GE, np.float32), np.asarray(GO, np.float32), ors, sfs


def viterbi_or(omarg, or_deg, lam):
    if lam <= 0:
        return omarg.argmax(1)
    T, K = omarg.shape
    logE = np.log(omarg + 1e-12); logE -= logE.max(1, keepdims=True)
    D = np.abs(or_deg[:, None] - or_deg[None, :]); D = np.minimum(D, 180 - D)
    trans = -lam * D ** 2
    score = logE[0].copy(); back = np.zeros((T, K), int)
    for t in range(1, T):
        cand = score[:, None] + trans; back[t] = cand.argmax(0); score = cand.max(0) + logE[t]
    path = np.zeros(T, int); path[-1] = int(score.argmax())
    for t in range(T - 1, 0, -1):
        path[t - 1] = back[t, path[t]]
    return path


def cstd(deg):
    z = np.mean(np.exp(1j * 2 * np.radians(deg)))
    return np.degrees(np.sqrt(max(-2 * np.log(abs(z)), 0))) / 2


def auto_lambda(omarg, or_deg):
    raw = (180 - or_deg[omarg.argmax(1)]) % 180
    target = 0.9 * cstd(raw); best = 0.0
    for lam in (2e-4, 5e-4, 1e-3, 2e-3):
        if cstd((180 - or_deg[viterbi_or(omarg, or_deg, lam)]) % 180) >= target:
            best = lam
    return best


def load_rfs(rf_dir=RF_DIR):
    rfs = []
    for fp in sorted(glob.glob(os.path.join(os.path.expanduser(rf_dir), "*_rf_metrics.csv"))):
        probe = os.path.basename(fp).split("_")[0]
        for r in csv.DictReader(open(fp)):
            try:
                cx = float(r["centroid_azimuth_pixels"]); cy = float(r["centroid_elevation_pixels"])
            except (KeyError, ValueError):
                continue
            if not (np.isfinite(cx) and np.isfinite(cy)):
                continue  # nofit / un-imputable rows carry NaN centroids
            rfs.append(dict(probe=probe, cluster=r["cluster_id"], rf_type=r["rf_type"],
                            cx=int(round(cx)), cy=int(round(cy))))
    return rfs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cloud-index", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--outdir", default="/ceph/margrie/laura/goggle_gabor/out")
    ap.add_argument("--rf-dir", default=RF_DIR,
                    help="dir of <probe>_rf_metrics.csv RF centres (ceph default; "
                         "point at the merged dir for the 1.5-sigma + imputed cohort)")
    ap.add_argument("--cloud-root", default=CLOUD_ROOT,
                    help="dir of theta* cloud-frame folders (ceph default; "
                         "local saved_goggles works too)")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--cpu", action="store_true")
    a = ap.parse_args()
    cloud_root = os.path.expanduser(a.cloud_root)
    os.makedirs(a.outdir, exist_ok=True)
    import torch
    dev = torch.device("cpu" if a.cpu or not torch.cuda.is_available() else "cuda")

    clouds = sorted(os.path.basename(p) for p in glob.glob(os.path.join(cloud_root, "theta*"))
                    if os.path.isdir(p))
    cloud = clouds[a.cloud_index]
    sftok = re.search(r"_sf(\d+p\d+)_", cloud).group(1).replace("p", ".")
    sf0 = float(sftok) / DPP; bsf = 0.005 / DPP
    print(f"[{a.cloud_index}] {cloud}  dev={dev}", flush=True)

    fs = sorted(glob.glob(os.path.join(cloud_root, cloud, "*.png")))[::a.stride]
    frames = np.stack([np.asarray(Image.open(f).convert("L"), np.float32) / 255.0 for f in fs])
    T, H, W = frames.shape
    hS, GES, GOS, orsS, sfsS = build_bank(SF_WIN, sf0, bsf)
    hO, GEO, GOO, orsO, sfsO = build_bank(OR_WIN, sf0, bsf)
    or_degO = np.degrees(orsO)
    PAD = max(hS, hO) + 1
    fr = torch.from_numpy(np.pad(frames, ((0, 0), (PAD, PAD), (PAD, PAD)), mode="reflect")).to(dev)
    tens = {k: torch.from_numpy(v).to(dev) for k, v in
            dict(GES=GES, GOS=GOS, GEO=GEO, GOO=GOO).items()}

    def energy(cy, cx, h, GE, GO):
        c = cy + PAD; d = cx + PAD
        patch = fr[:, c - h:c + h, d - h:d + h]
        patch = patch - patch.mean(dim=(1, 2), keepdim=True)
        re = torch.einsum("kij,tij->tk", GE, patch); iv = torch.einsum("kij,tij->tk", GO, patch)
        return re * re + iv * iv

    rfs = load_rfs(a.rf_dir)
    rows = []
    for rf in rfs:
        cx, cy = rf["cx"], rf["cy"]
        edge = (cx - hS < 0) or (cx + hS > W) or (cy - hS < 0) or (cy + hS > H)
        enS = energy(cy, cx, hS, tens["GES"], tens["GOS"]).reshape(T, len(orsS), len(sfsS))
        sf_t = sfsS[(enS.amax(dim=1).argmax(dim=1)).cpu().numpy()]
        conc = float(np.median((enS.amax(dim=(1, 2)) / enS.mean(dim=(1, 2))).cpu().numpy()))
        enO = energy(cy, cx, hO, tens["GEO"], tens["GOO"]).reshape(T, len(orsO), len(sfsO))
        omarg = enO.sum(dim=2).cpu().numpy()
        lam = auto_lambda(omarg, or_degO)
        or_t = (180 - or_degO[viterbi_or(omarg, or_degO, lam)]) % 180
        for k in range(T):
            rows.append((rf["probe"], rf["cluster"], rf["rf_type"], cx, cy,
                         k * a.stride, float(sf_t[k]), float(or_t[k]), conc, float(edge)))
    import pandas as pd
    df = pd.DataFrame(rows, columns=["probe", "cluster", "rf_type", "cx", "cy",
                                     "frame", "sf_cpd", "or_deg", "concentration", "edge"])
    df["cloud"] = cloud
    out = os.path.join(a.outdir, f"{cloud}.parquet")
    df.to_parquet(out, index=False)
    print(f"[{a.cloud_index}] wrote {out}  ({len(rfs)} RFs × {T} frames = {len(rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
