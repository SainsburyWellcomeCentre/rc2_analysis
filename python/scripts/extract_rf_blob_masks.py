"""Extract ON/OFF receptive-field subfield masks from the sparse-noise combined PDF.

The goggle RF maps are STA blobs (`RF_raw = STA > mean+2.5·std`, bwareaopen) computed
on-the-fly by `calculate_compare_rf.m` and are NOT saved as arrays — the only durable
artifact is the per-cluster combined PDF (`*_rf_by_cluster_goggles_combined.pdf`). This
script reads a cluster's page, thresholds the red (ON/white) and blue (OFF/black) blobs
in the two "proc" panels, and maps them into the 400px stimulus-buffer coordinate frame
(the same frame the motion-cloud PNGs live in), so the blob edges can be overlaid on the
cloud at the RF location.

Coordinate convention (validated against the rf_metrics centroid, 2026-06-11):
  The proc panels are drawn `imagesc(RF_g); set(gca,'xdir','reverse')` — i.e. VISUAL-FIELD
  orientation (high azimuth on the left). `centroid_azimuth_pixels` is the RAW BUFFER column.
  Panel fraction (fx from panel-left, fy from panel-top) maps to buffer pixels by:
      az_px = (1 - fx) * W        (x is mirrored: buffer right = panel left)
      el_px =      fy  * W        (y is direct: top = elevation 0)
  Calibrated on cl90/CAA-1124371: ON panel (fx0.331, fy0.595) -> (az0.669, el0.595) vs
  CSV (0.670, 0.591); OFF (fx0.338,fy0.597) -> (0.662,0.597) vs (0.664,0.593). <1% both axes.

Output: npz with bool arrays `on`, `off` (shape W×W, buffer coords). Feed to
`gabor_goggle_animation.py --rf-mask-npz`.

Usage:
  python extract_rf_blob_masks.py --pdf <combined.pdf> --cluster 90 --out <masks.npz>
  (page is auto-located by the "ClusterNN" text; panel frames auto-detected from the render.)
"""
from __future__ import annotations
import argparse, os, subprocess, tempfile
import numpy as np
from scipy import ndimage as ndi
import matplotlib.image as mpi

W = 400  # stimulus-buffer / cloud-frame size in px (wisecoco, 0.279°/px)


def find_page(pdf, cluster):
    """Return the 1-based page index whose text contains 'Cluster<NN>'."""
    import re
    info = subprocess.run(["pdfinfo", pdf], capture_output=True, text=True).stdout
    npages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
    for p in range(1, npages + 1):
        txt = subprocess.run(["pdftotext", "-f", str(p), "-l", str(p), pdf, "-"],
                             capture_output=True, text=True).stdout
        if re.search(rf"Cluster\s*{cluster}\b", txt.replace(" ", "")) or \
           re.search(rf"Cluster\s*{cluster}\b", txt):
            return p
    raise SystemExit(f"cluster {cluster} not found in {pdf}")


def render_page(pdf, page, dpi=300):
    d = tempfile.mkdtemp()
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(dpi),
                    "-png", pdf, os.path.join(d, "pg")], check=True)
    f = [x for x in os.listdir(d) if x.endswith(".png")][0]
    im = mpi.imread(os.path.join(d, f))
    return im[..., :3]


def detect_panel_frames(rgb):
    """Find the proc-panel black rectangle borders (right column). Returns (xs, ys) groups."""
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    H, Wp, _ = rgb.shape
    dark = (r < 0.25) & (g < 0.25) & (b < 0.25)
    rows = [y for y in range(H) if dark[y].sum() > 0.13 * Wp
            and (np.ptp(np.where(dark[y])[0]) if dark[y].any() else 0) > 0.15 * Wp]
    cols = [x for x in range(Wp) if dark[:, x].sum() > 0.10 * H
            and (np.ptp(np.where(dark[:, x])[0]) if dark[:, x].any() else 0) > 0.10 * H]

    def group(a, gap=10):
        a = sorted(a); out, cur = [], [a[0]]
        for v in a[1:]:
            (cur.append(v) if v - cur[-1] <= gap else (out.append(int(np.mean(cur))), cur.clear(), cur.append(v)))
        out.append(int(np.mean(cur))); return out
    return group(cols), group(rows)


def blob_to_buffer(mask, xL, xR, yT, yB):
    az = np.arange(W); el = np.arange(W)
    AZ, EL = np.meshgrid(az, el)
    fx = 1.0 - AZ / W                      # az_px = (1 - fx) * W  (x mirrored)
    fy = EL / W                            # el_px = fy * W        (y direct)
    px = np.clip((xL + fx * (xR - xL)).astype(int), 0, mask.shape[1] - 1)
    py = np.clip((yT + fy * (yB - yT)).astype(int), 0, mask.shape[0] - 1)
    Mb = ndi.binary_closing(mask[py, px], iterations=2)
    lab, n = ndi.label(Mb)
    if n > 1:                              # keep the largest connected component
        sizes = ndi.sum(Mb, lab, range(1, n + 1)); Mb = lab == (1 + int(np.argmax(sizes)))
    return Mb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--cluster", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()
    page = find_page(a.pdf, a.cluster)
    rgb = render_page(a.pdf, page, a.dpi)
    xs, ys = detect_panel_frames(rgb)
    # right column = last vertical-border pair; ON = top row, OFF = bottom row
    xL, xR = xs[-2], xs[-1]
    (yT_on, yB_on), (yT_off, yB_off) = (ys[0], ys[1]), (ys[2], ys[3])
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    red = (r > 0.45) & (g < 0.55) & (b < 0.55)
    blue = (b > 0.45) & (r < 0.55) & (g < 0.55)
    on = blob_to_buffer(red, xL, xR, yT_on, yB_on)
    off = blob_to_buffer(blue, xL, xR, yT_off, yB_off)
    for nm, M in [("ON", on), ("OFF", off)]:
        yy, xx = np.where(M)
        print(f"{nm}: n={M.sum()} centroid az={xx.mean():.0f} el={yy.mean():.0f} "
              f"(page {page}, panel x[{xL},{xR}])")
    np.savez_compressed(a.out, on=on, off=off)
    print(f"[saved] {a.out}")


if __name__ == "__main__":
    main()
