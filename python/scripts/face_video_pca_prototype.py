"""Face-video PCA prototype (Stringer-2019 style) — short slice OR full session.

Per `/Users/laura/.claude/plans/motion-clouds-project-you-misty-owl.md`.

Two modes:
  - Default: in-memory slice driven by `--max-seconds` (default 60 s, full res 480x640).
  - `--streaming`: IncrementalPCA over the whole H.264 stream at `--downsample` 4
    (default), two passes (fit + project). Memory bounded by `--batch-size`.

Figure titles include the actual duration processed and the per-PC threshold
quantile used for the heatmap overlay.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import IncrementalPCA, TruncatedSVD

RC2_PY_SRC = Path("/Users/laura/source/github/SainsburyWellcomeCentre/rc2_analysis/python/src")
if str(RC2_PY_SRC) not in sys.path:
    sys.path.insert(0, str(RC2_PY_SRC))

from rc2_formatted_data_reader.reader import FormattedDataReader  # noqa: E402


DEFAULTS = dict(
    video_path=Path("/Volumes/margrie/mvelez/mateoData_cameras/CAA-1123243_rec1/camera0.mp4"),
    mat_path=Path("~/local_data/motion_clouds/formatted_data/CAA-1123243_rec1.mat").expanduser(),
    fig_dir=Path("~/local_data/motion_clouds/figures/glm/exploration/face_video_pca").expanduser(),
    array_dir=Path("~/local_data/motion_clouds/face_video_pca").expanduser(),
    probe="CAA-1123243_rec1",
    start_seconds=0.0,
    max_seconds=60.0,
    frame_stride=5,
    downsample=1,
    streaming_downsample=4,
    batch_size=2000,
    n_components=20,
    roi=None,
    pc_top_quantile=0.10,
)


def _to_gray_downsampled(frame_rgb: np.ndarray, downsample: int) -> np.ndarray:
    g = frame_rgb.mean(axis=2).astype(np.float32)
    if downsample <= 1:
        return g
    h, w = g.shape
    h2 = (h // downsample) * downsample
    w2 = (w // downsample) * downsample
    g = g[:h2, :w2]
    return g.reshape(h2 // downsample, downsample, w2 // downsample, downsample).mean(axis=(1, 3))


def load_slice(video_path, start_seconds, max_seconds, frame_stride, downsample, roi):
    import imageio.v2 as imageio
    reader = imageio.get_reader(str(video_path), format="ffmpeg")
    meta = reader.get_meta_data()
    fps = float(meta.get("fps", 60.0))
    start_idx = int(round(start_seconds * fps))
    end_idx = start_idx + int(round(max_seconds * fps))
    frames = []
    first_native_idx = None
    try:
        for i, frame in enumerate(reader):
            if i < start_idx:
                continue
            if i >= end_idx:
                break
            if (i - start_idx) % frame_stride != 0:
                continue
            g = _to_gray_downsampled(frame, downsample)
            if roi is not None:
                y0, y1, x0, x1 = roi
                g = g[y0:y1, x0:x1]
            if first_native_idx is None:
                first_native_idx = i
            frames.append(g)
    finally:
        reader.close()
    arr = np.stack(frames, axis=0)
    return arr, dict(
        n_frames_kept=arr.shape[0],
        height=arr.shape[1],
        width=arr.shape[2],
        native_fps=fps,
        effective_fps=fps / frame_stride,
        first_native_idx=first_native_idx,
    )


def iter_video_chunks(video_path, frame_stride, downsample, batch_size, roi):
    """Yield (chunk, batch_index). Each chunk has shape (B+1, H, W) when a
    previous frame is available, otherwise (B, H, W); the consumer takes its
    L1 frame-difference to get B (or B-1 for the very first chunk) ME frames."""
    import imageio.v2 as imageio
    reader = imageio.get_reader(str(video_path), format="ffmpeg")
    meta = reader.get_meta_data()
    fps = float(meta.get("fps", 60.0))
    buf: list[np.ndarray] = []
    prev_last: np.ndarray | None = None
    batch_idx = 0
    try:
        for i, frame in enumerate(reader):
            if i % frame_stride != 0:
                continue
            g = _to_gray_downsampled(frame, downsample)
            if roi is not None:
                y0, y1, x0, x1 = roi
                g = g[y0:y1, x0:x1]
            buf.append(g)
            if len(buf) >= batch_size:
                if prev_last is not None:
                    chunk = np.stack([prev_last] + buf, axis=0)
                else:
                    chunk = np.stack(buf, axis=0)
                prev_last = buf[-1]
                yield chunk, batch_idx, fps
                batch_idx += 1
                buf = []
        if buf:
            if prev_last is not None:
                chunk = np.stack([prev_last] + buf, axis=0)
            else:
                chunk = np.stack(buf, axis=0)
            yield chunk, batch_idx, fps
    finally:
        reader.close()


def fit_streaming(video_path, frame_stride, downsample, batch_size, n_components, roi):
    ipca = IncrementalPCA(n_components=n_components)
    H = W = -1
    representative_frame = None
    n_diffs_total = 0
    fps = 60.0
    for chunk, batch_idx, fps in iter_video_chunks(
        video_path, frame_stride, downsample, batch_size, roi
    ):
        if H < 0:
            _, H, W = chunk.shape
        me_chunk = np.abs(np.diff(chunk, axis=0)).astype(np.float32)
        me_flat = me_chunk.reshape(me_chunk.shape[0], H * W)
        if me_flat.shape[0] < n_components:
            continue
        ipca.partial_fit(me_flat)
        n_diffs_total += me_flat.shape[0]
        if representative_frame is None:
            representative_frame = chunk[len(chunk) // 2].copy()
        print(f"  fit chunk {batch_idx + 1}: +{me_flat.shape[0]} frames "
              f"(running total {n_diffs_total})", flush=True)
    return ipca, (H, W), representative_frame, n_diffs_total, fps


def project_streaming(video_path, frame_stride, downsample, batch_size, ipca, roi, n_total):
    pc1 = np.empty(n_total, dtype=np.float64)
    cursor = 0
    for chunk, batch_idx, _ in iter_video_chunks(
        video_path, frame_stride, downsample, batch_size, roi
    ):
        me_chunk = np.abs(np.diff(chunk, axis=0)).astype(np.float32)
        if me_chunk.shape[0] == 0:
            continue
        H, W = me_chunk.shape[1], me_chunk.shape[2]
        me_flat = me_chunk.reshape(me_chunk.shape[0], H * W)
        proj = ipca.transform(me_flat)
        end = min(cursor + proj.shape[0], n_total)
        proj = proj[: end - cursor]
        pc1[cursor:end] = proj[:, 0]
        cursor = end
        print(f"  project chunk {batch_idx + 1}: cursor={cursor}/{n_total}", flush=True)
        if cursor >= n_total:
            break
    return pc1[:cursor]


def run_svd_in_memory(me: np.ndarray, n_components: int) -> dict:
    T, H, W = me.shape
    M = me.reshape(T, H * W)
    pixel_mean = M.mean(axis=0, keepdims=True)
    M -= pixel_mean  # in-place; caller must not reuse `me`
    svd = TruncatedSVD(n_components=n_components, algorithm="randomized", random_state=0)
    U_scaled = svd.fit_transform(M)
    V = svd.components_.reshape(n_components, H, W)
    return dict(
        U=U_scaled,
        V=V,
        singular_values=svd.singular_values_,
        explained_variance_ratio=svd.explained_variance_ratio_,
        pixel_mean=pixel_mean.reshape(H, W),
    )


def _suptitle(probe: str, duration_min: float, extra: str = "") -> str:
    base = f"Face motion-energy PCA — {probe}, {duration_min:.1f} min"
    return f"{base}  ({extra})" if extra else base


def plot_scree(evr: np.ndarray, out_path: Path, probe: str, duration_min: float) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 3.7), dpi=140)
    ax.plot(np.arange(1, len(evr) + 1), evr, "o-", lw=1, ms=4)
    ax.set_xlabel("PC index")
    ax.set_ylabel("Explained variance ratio")
    ax.set_yscale("log")
    ax.set_title(_suptitle(probe, duration_min, "scree"))
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"))
    fig.savefig(out_path.with_suffix(".png"))
    plt.close(fig)


def plot_spatial_pcs(V, background, evr, out_path, probe, duration_min,
                     n_cols=5, top_quantile=0.10):
    K = V.shape[0]
    n_rows = int(np.ceil(K / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.6 * n_cols, 2.6 * n_rows), dpi=140)
    axes = np.atleast_1d(axes).flatten()
    bg_norm = (background - background.min()) / max(background.max() - background.min(), 1e-9)
    for k in range(K):
        ax = axes[k]
        pc = V[k]
        abs_pc = np.abs(pc)
        thresh = float(np.quantile(abs_pc, 1.0 - top_quantile))
        vmax = float(np.percentile(abs_pc, 99.5))
        alpha_mask = (abs_pc >= thresh).astype(np.float32)
        ax.imshow(bg_norm, cmap="gray", interpolation="nearest")
        ax.imshow(pc, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                  alpha=alpha_mask, interpolation="nearest")
        ax.set_title(f"PC{k + 1}  ({evr[k] * 100:.1f}%)", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    for k in range(K, len(axes)):
        axes[k].set_visible(False)
    fig.suptitle(_suptitle(probe, duration_min,
                           f"showing |loading| top {top_quantile * 100:.0f}% per PC"),
                 y=1.01)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def sanity_check_pc1_vs_scalar_me(pc1_trace, first_native_idx, frame_stride,
                                  mat_path, out_path, probe, duration_min):
    reader = FormattedDataReader(str(mat_path))
    me_scalar = reader.camera0()
    if me_scalar is None:
        raise RuntimeError(f"camera0 not present in {mat_path}")
    indices = first_native_idx + (np.arange(len(pc1_trace)) + 1) * frame_stride
    indices = indices[indices < len(me_scalar)]
    pc1_sub = pc1_trace[: len(indices)]
    me_sub = np.asarray(me_scalar[indices], dtype=np.float64)
    valid = np.isfinite(pc1_sub) & np.isfinite(me_sub)
    pc1_sub = pc1_sub[valid]
    me_sub = me_sub[valid]
    r = float(np.corrcoef(pc1_sub, me_sub)[0, 1])
    fig, axes = plt.subplots(2, 1, figsize=(9, 4.5), dpi=140, sharex=True)
    axes[0].plot(me_sub, lw=0.4, color="k")
    axes[0].set_ylabel("session.camera0\n(scalar ME)")
    axes[1].plot(pc1_sub, lw=0.4, color="C0")
    axes[1].set_ylabel("PC1 temporal trace")
    axes[1].set_xlabel("Subsampled frame index")
    axes[0].set_title(_suptitle(probe, duration_min,
                                f"PC1 vs scalar ME — Pearson r = {r:.3f}, n={len(pc1_sub)}"))
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"))
    fig.savefig(out_path.with_suffix(".png"))
    plt.close(fig)
    return r


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--video", type=Path, default=DEFAULTS["video_path"])
    parser.add_argument("--mat", type=Path, default=DEFAULTS["mat_path"])
    parser.add_argument("--fig-dir", type=Path, default=DEFAULTS["fig_dir"])
    parser.add_argument("--array-dir", type=Path, default=DEFAULTS["array_dir"])
    parser.add_argument("--probe", default=DEFAULTS["probe"])
    parser.add_argument("--start-seconds", type=float, default=DEFAULTS["start_seconds"])
    parser.add_argument("--max-seconds", type=float, default=DEFAULTS["max_seconds"])
    parser.add_argument("--frame-stride", type=int, default=DEFAULTS["frame_stride"])
    parser.add_argument("--downsample", type=int, default=DEFAULTS["downsample"])
    parser.add_argument("--batch-size", type=int, default=DEFAULTS["batch_size"])
    parser.add_argument("--n-components", type=int, default=DEFAULTS["n_components"])
    parser.add_argument("--roi", type=int, nargs=4, default=None, metavar=("Y0", "Y1", "X0", "X1"))
    parser.add_argument("--pc-top-quantile", type=float, default=DEFAULTS["pc_top_quantile"])
    parser.add_argument("--streaming", action="store_true",
                        help="Stream the full video via IncrementalPCA (ignores --max-seconds).")
    parser.add_argument("--streaming-downsample", type=int,
                        default=DEFAULTS["streaming_downsample"],
                        help="Spatial downsample factor used in --streaming mode (default 4).")
    parser.add_argument("--replot-from-npz", action="store_true")
    args = parser.parse_args()

    args.fig_dir.mkdir(parents=True, exist_ok=True)
    args.array_dir.mkdir(parents=True, exist_ok=True)
    roi = tuple(args.roi) if args.roi is not None else None
    npz_path = args.array_dir / f"{args.probe}_face_pcs.npz"

    if args.replot_from_npz:
        print(f"[replot] loading {npz_path}")
        data = np.load(npz_path)
        duration_min = float(data.get("duration_seconds", 0.0)) / 60.0
        plot_scree(data["explained_variance_ratio"], args.fig_dir / "scree",
                   args.probe, duration_min)
        plot_spatial_pcs(
            data["spatial_components"], data["representative_frame"],
            data["explained_variance_ratio"], args.fig_dir / "spatial_pcs_grid",
            args.probe, duration_min, top_quantile=args.pc_top_quantile,
        )
        print(f"[replot] re-rendered ({duration_min:.1f} min) "
              f"with top_quantile={args.pc_top_quantile}")
        return 0

    if args.streaming:
        ds = args.streaming_downsample
        print(f"[streaming 1/3] IncrementalPCA fit (stride={args.frame_stride}, "
              f"downsample={ds}, batch={args.batch_size}, k={args.n_components})",
              flush=True)
        ipca, (H, W), representative_frame, n_diffs, fps = fit_streaming(
            args.video, args.frame_stride, ds, args.batch_size, args.n_components, roi,
        )
        duration_seconds = (n_diffs + 1) * args.frame_stride / fps
        duration_min = duration_seconds / 60.0
        print(f"      fit done: H={H} W={W}  n_diffs={n_diffs}  "
              f"duration={duration_min:.2f} min", flush=True)
        V = ipca.components_.reshape(args.n_components, H, W)
        evr = ipca.explained_variance_ratio_
        print(f"      explained variance top 5: {evr[:5]}", flush=True)
        print(f"      cumulative top-{args.n_components}: {evr.sum():.3f}", flush=True)

        print("[streaming 2/3] projecting onto PCs (second pass)", flush=True)
        pc1_trace = project_streaming(
            args.video, args.frame_stride, ds, args.batch_size, ipca, roi, n_diffs,
        )

        print("[streaming 3/3] writing figures + arrays", flush=True)
        first_native_idx = 0
        result = dict(
            U=None, V=V, singular_values=ipca.singular_values_,
            explained_variance_ratio=evr, pixel_mean=ipca.mean_.reshape(H, W),
        )
        plot_scree(evr, args.fig_dir / "scree", args.probe, duration_min)
        plot_spatial_pcs(V, representative_frame, evr,
                         args.fig_dir / "spatial_pcs_grid",
                         args.probe, duration_min,
                         top_quantile=args.pc_top_quantile)
        r = sanity_check_pc1_vs_scalar_me(
            pc1_trace, first_native_idx, args.frame_stride, args.mat,
            args.fig_dir / "pc1_vs_scalar_me", args.probe, duration_min,
        )
        np.savez_compressed(
            npz_path,
            spatial_components=V,
            pc1_trace=pc1_trace,
            singular_values=ipca.singular_values_,
            explained_variance_ratio=evr,
            pixel_mean=result["pixel_mean"],
            representative_frame=representative_frame,
            frame_stride=args.frame_stride,
            downsample=ds,
            start_seconds=0.0,
            max_seconds=duration_seconds,
            duration_seconds=duration_seconds,
            first_native_idx=first_native_idx,
            roi=np.array(roi) if roi is not None else np.array([]),
            video_path=str(args.video),
            mode="streaming",
        )
        print(f"      {npz_path}")
        print(f"\nDone (streaming). Duration {duration_min:.2f} min, "
              f"PC1 vs scalar-ME Pearson r = {r:.3f}")
        return 0

    print(f"[1/4] loading {args.max_seconds:g}s slice from t={args.start_seconds:g}s "
          f"(stride={args.frame_stride}, downsample={args.downsample})")
    frames, info = load_slice(
        args.video, args.start_seconds, args.max_seconds,
        args.frame_stride, args.downsample, roi,
    )
    duration_seconds = info["n_frames_kept"] * args.frame_stride / info["native_fps"]
    duration_min = duration_seconds / 60.0
    print(f"      kept {info['n_frames_kept']} frames at {info['effective_fps']:.2f} Hz, "
          f"shape {info['height']}x{info['width']}, duration {duration_min:.2f} min")
    print(f"      memory footprint: {frames.nbytes / 1e6:.0f} MB")

    print("[2/4] computing motion-energy movie")
    me = np.abs(np.diff(frames, axis=0)).astype(np.float32)
    representative_frame = frames[len(frames) // 2].copy()
    del frames

    print(f"[3/4] TruncatedSVD (k={args.n_components}, n_samples={me.shape[0]}, "
          f"n_features={me.shape[1] * me.shape[2]})")
    result = run_svd_in_memory(me, args.n_components)
    print(f"      explained variance ratio top 5: {result['explained_variance_ratio'][:5]}")
    print(f"      cumulative top-{args.n_components}: "
          f"{result['explained_variance_ratio'].sum():.3f}")

    print("[4/4] writing figures + arrays")
    plot_scree(result["explained_variance_ratio"], args.fig_dir / "scree",
               args.probe, duration_min)
    plot_spatial_pcs(
        result["V"], representative_frame, result["explained_variance_ratio"],
        args.fig_dir / "spatial_pcs_grid", args.probe, duration_min,
        top_quantile=args.pc_top_quantile,
    )
    r = sanity_check_pc1_vs_scalar_me(
        result["U"][:, 0], info["first_native_idx"], args.frame_stride,
        args.mat, args.fig_dir / "pc1_vs_scalar_me",
        args.probe, duration_min,
    )

    np.savez_compressed(
        npz_path,
        spatial_components=result["V"],
        temporal_projections=result["U"],
        singular_values=result["singular_values"],
        explained_variance_ratio=result["explained_variance_ratio"],
        pixel_mean=result["pixel_mean"],
        representative_frame=representative_frame,
        frame_stride=args.frame_stride,
        downsample=args.downsample,
        start_seconds=args.start_seconds,
        max_seconds=args.max_seconds,
        duration_seconds=duration_seconds,
        first_native_idx=info["first_native_idx"],
        roi=np.array(roi) if roi is not None else np.array([]),
        video_path=str(args.video),
        mode="in_memory",
    )
    print(f"      {npz_path}")
    print(f"\nDone. Duration {duration_min:.2f} min, "
          f"PC1 vs scalar-ME Pearson r = {r:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
