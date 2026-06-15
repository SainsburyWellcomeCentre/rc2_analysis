"""Speed × face-ME firing-rate heatmap on probe 243.

Two modes:
  - `--mode representatives` (default): one cluster per category — Speed-only,
    Speed+ME (additive), Speed+ME+(MExSpeed) (interaction), ME-only. 4 panels.
  - `--mode all-me-speed-combos`: every cluster whose GLM-selected vars are a
    subset of {Speed, ME_face, ME_face_x_Speed}, grouped by category, sorted by
    cv_bps. Renders a tall grid.

Per cluster, restricts to motion bins (condition != 'stationary'), bins speed
and me_face_raw at 5% quantiles each (20 bins / axis from the pooled motion-bin
distribution so panels are comparable), and renders mean firing rate per cell.
Cells with < 10 samples are masked.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RC2_PY_SRC = Path("/Users/laura/source/github/SainsburyWellcomeCentre/rc2_analysis/python/src")
if str(RC2_PY_SRC) not in sys.path:
    sys.path.insert(0, str(RC2_PY_SRC))

from rc2_glm.io import load_probe_data  # noqa: E402
from rc2_glm.config import GLMConfig  # noqa: E402
from rc2_glm.time_binning import bin_probe  # noqa: E402

PROBE_ID = "CAA-1123243_rec1"
MAT_PATH = Path("~/local_data/motion_clouds/formatted_data/CAA-1123243_rec1.mat").expanduser()
GLM_CMP_CSV = Path(
    "~/local_data/motion_clouds/figures/glm/current_with_ME_3_probes/glm_model_comparison.csv"
).expanduser()
OUT_DIR = Path(
    "~/local_data/motion_clouds/figures/glm/exploration/speed_x_me_heatmap"
).expanduser()
N_BINS = 20
MIN_SAMPLES_PER_CELL = 10


CATEGORY_MASKS: dict[str, list[str]] = {
    "Speed": ["Speed"],
    "Speed + ME_face": ["Speed+ME_face", "ME_face+Speed"],
    "Speed + ME_face + ME×Speed": [
        "Speed+ME_face+ME_face_x_Speed", "ME_face+Speed+ME_face_x_Speed"
    ],
    "ME_face": ["ME_face"],
}
ALL_ALLOWED_SELECTED_VARS = sorted({v for vs in CATEGORY_MASKS.values() for v in vs})


def _category_for(selected_vars: str) -> str | None:
    for label, options in CATEGORY_MASKS.items():
        if selected_vars in options:
            return label
    return None


def pick_representatives(cmp_df: pd.DataFrame) -> list[tuple[str, int, float]]:
    """One cluster per category, sorted to best cv_bps."""
    df = cmp_df[cmp_df["probe_id"] == PROBE_ID].copy()
    picks = []
    for label, options in CATEGORY_MASKS.items():
        sub = df[df["time_selected_vars"].isin(options)]
        if sub.empty:
            continue
        best = sub.sort_values("time_Selected_cv_bps", ascending=False).iloc[0]
        picks.append((label, int(best.cluster_id), float(best.time_Selected_cv_bps)))
    return picks


def pick_all_me_speed_combos(cmp_df: pd.DataFrame) -> list[tuple[str, int, float]]:
    """All probe-243 clusters whose selected_vars is in CATEGORY_MASKS.

    Ordered by category (Speed → Speed+ME → triple → ME-only), and within
    each category by cv_bps descending.
    """
    df = cmp_df[cmp_df["probe_id"] == PROBE_ID].copy()
    picks = []
    for label, options in CATEGORY_MASKS.items():
        sub = df[df["time_selected_vars"].isin(options)].sort_values(
            "time_Selected_cv_bps", ascending=False,
        )
        for _, row in sub.iterrows():
            picks.append((label, int(row.cluster_id), float(row.time_Selected_cv_bps)))
    return picks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=["representatives", "all-me-speed-combos"],
        default="representatives",
    )
    parser.add_argument("--cols", type=int, default=0,
                        help="Override grid columns (0 = auto).")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] reading {GLM_CMP_CSV.name} + picking clusters (mode={args.mode})")
    cmp_df = pd.read_csv(GLM_CMP_CSV)
    if args.mode == "representatives":
        picks = pick_representatives(cmp_df)
    else:
        picks = pick_all_me_speed_combos(cmp_df)
    for label, cid, cv in picks:
        print(f"      {label:30s}  cluster {cid:4d}  cv_bps={cv:.3f}")

    print(f"[2/4] loading {PROBE_ID} + bin_probe (100 ms)")
    config = GLMConfig()
    probe = load_probe_data(str(MAT_PATH), config=config)
    df = bin_probe(probe)
    print(f"      bin_probe shape: {df.shape}")

    motion = df[df["condition"] != "stationary"].copy()
    motion = motion.dropna(subset=["me_face_raw", "speed"])
    print(f"      motion-only rows: {len(motion):,}")

    # Pooled quantile edges so all panels are on the same axes
    speed_edges = np.quantile(motion["speed"].to_numpy(), np.linspace(0, 1, N_BINS + 1))
    me_edges = np.quantile(motion["me_face_raw"].to_numpy(), np.linspace(0, 1, N_BINS + 1))
    # nudge to avoid duplicate edges (low-speed often has many zero values)
    speed_edges = np.unique(speed_edges)
    me_edges = np.unique(me_edges)
    if len(speed_edges) < 2 or len(me_edges) < 2:
        raise RuntimeError("Degenerate quantile edges")
    n_speed_bins = len(speed_edges) - 1
    n_me_bins = len(me_edges) - 1
    print(f"      speed bins: {n_speed_bins}  ({speed_edges[0]:.2f} ... {speed_edges[-1]:.2f} cm/s)")
    print(f"      me bins:    {n_me_bins}  ({me_edges[0]:.3f} ... {me_edges[-1]:.3f})")

    motion["speed_bin"] = np.clip(
        np.digitize(motion["speed"], speed_edges) - 1, 0, n_speed_bins - 1
    )
    motion["me_bin"] = np.clip(
        np.digitize(motion["me_face_raw"], me_edges) - 1, 0, n_me_bins - 1
    )

    print(f"[3/4] computing per-cluster heatmaps (FR = spike_count / {config.time_bin_width}s)")
    panels = []
    for label, cid, cv in picks:
        sub = motion[motion["cluster_id"] == cid]
        if sub.empty:
            print(f"      WARNING: cluster {cid} not found in binned data, skipping")
            continue
        cell_mean = sub.groupby(["me_bin", "speed_bin"])["spike_count"].mean().unstack()
        cell_count = sub.groupby(["me_bin", "speed_bin"])["spike_count"].size().unstack()
        # Reindex to full grid, NaN where missing
        cell_mean = cell_mean.reindex(index=range(n_me_bins), columns=range(n_speed_bins))
        cell_count = cell_count.reindex(index=range(n_me_bins), columns=range(n_speed_bins))
        fr = cell_mean.to_numpy(dtype=np.float64) / config.time_bin_width  # Hz
        count = cell_count.to_numpy(dtype=np.float64)
        fr = np.where(count >= MIN_SAMPLES_PER_CELL, fr, np.nan)
        panels.append(dict(label=label, cid=cid, cv=cv, fr=fr, count=count, n_bins=len(sub)))
        n_valid = np.isfinite(fr).sum()
        print(f"      cluster {cid:4d}  ({label[:36]:36s})  "
              f"valid cells {n_valid}/{n_speed_bins * n_me_bins}  "
              f"FR range [{np.nanmin(fr):.1f}, {np.nanmax(fr):.1f}] Hz")

    print(f"[4/4] rendering figure")
    n_panels = len(panels)
    if args.cols > 0:
        n_cols = args.cols
    elif n_panels <= 4:
        n_cols = 2
    elif n_panels <= 12:
        n_cols = 4
    else:
        n_cols = 5
    n_rows = int(np.ceil(n_panels / n_cols))
    panel_w, panel_h = (5.5, 4.5) if n_panels <= 4 else (3.4, 3.0)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(panel_w * n_cols, panel_h * n_rows), dpi=140,
        squeeze=False,
    )
    axes = axes.flatten()

    title_fs = 10 if n_panels <= 4 else 8
    label_fs = 9 if n_panels <= 4 else 7

    for ax, panel in zip(axes, panels):
        fr = panel["fr"]
        vmin = np.nanpercentile(fr, 2)
        vmax = np.nanpercentile(fr, 98)
        im = ax.imshow(
            fr,
            origin="lower",
            aspect="auto",
            cmap="viridis",
            vmin=vmin, vmax=vmax,
            extent=(speed_edges[0], speed_edges[-1], me_edges[0], me_edges[-1]),
        )
        ax.set_xlabel("Speed (cm/s)", fontsize=label_fs)
        ax.set_ylabel("ME_face_raw", fontsize=label_fs)
        ax.tick_params(labelsize=label_fs - 1)
        ax.set_title(
            f"cluster {panel['cid']}  [{panel['label']}]\n"
            f"cv_bps={panel['cv']:.2f}",
            fontsize=title_fs,
        )
        cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cbar.set_label("FR (Hz)", fontsize=label_fs)
        cbar.ax.tick_params(labelsize=label_fs - 1)

    for ax in axes[n_panels:]:
        ax.set_visible(False)

    fig.suptitle(
        f"Speed × face-ME joint tuning — {PROBE_ID}, 100 ms motion bins, "
        f"{N_BINS}×{N_BINS} 5%-quantile cells (≥{MIN_SAMPLES_PER_CELL} samples)  "
        f"[mode={args.mode}, n_clusters={n_panels}]",
        y=1.002, fontsize=11,
    )
    fig.tight_layout()
    out_base = OUT_DIR / f"speed_x_me_heatmap_{PROBE_ID}_{args.mode}"
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    print(f"      wrote {out_base.with_suffix('.png')}")
    print(f"      wrote {out_base.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
