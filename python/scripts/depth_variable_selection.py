"""Cortical depth × GLM-variable-selection — motion-clouds GLM runs.

For each retained cluster in a GLM run, the Hardcastle forward selection keeps
a subset of {Speed, TF, SF, OR, ME_face}. This script attaches each cluster's
continuous cortical depth (read from the formatted .mat) and renders a
2×(N+1) faceted figure: a leftmost cluster-spread column (every cluster as a
depth-positioned dot), then one column per selected variable, with
top row = raw count of clusters per depth bin that selected the variable and
bottom row = that count as a fraction of clusters in the bin.

ME_face is included only when the run actually fit it (the no-ME `current/`
run carries the column but it is all-False, so the figure drops to 4 columns).

Layer boundaries (VISp2/3 .. VISp6b) are *data-derived* — midpoints of each
layer's mean depth among the labelled clusters — and drawn as dashed lines,
MATLAB MIDepthPlot-style. They are not anatomical ground truth.

Colours reuse the forward-selection-summary palette from
rc2_analysis/python/src/rc2_glm/plots.py:599-614.

Output: ~/local_data/motion_clouds/figures/glm/exploration/depth_variable_selection/
  depth_variable_selection_<label>.{pdf,png}
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import blended_transform_factory


def _rc2_python_src() -> Path:
    """rc2_analysis python/src — metonymy uses ~/Source, synecdoche ~/source."""
    home = Path.home()
    for src_seg in ("Source", "source"):
        cand = (
            home / src_seg / "github" / "SainsburyWellcomeCentre"
            / "rc2_analysis" / "python" / "src"
        )
        if cand.exists():
            return cand
    raise RuntimeError(f"rc2_analysis python/src not found under {home}")


RC2_PY_SRC = _rc2_python_src()
if str(RC2_PY_SRC) not in sys.path:
    sys.path.insert(0, str(RC2_PY_SRC))

from rc2_formatted_data_reader import FormattedDataReader  # noqa: E402

DEFAULT_RUN_DIR = Path(
    "~/local_data/motion_clouds/figures/glm/current_with_ME_3_probes"
).expanduser()
FORMATTED_DIR = Path("~/local_data/motion_clouds/formatted_data").expanduser()
OUT_DIR = Path(
    "~/local_data/motion_clouds/figures/glm/exploration/depth_variable_selection"
).expanduser()

# Forward-selection-summary palette — copied verbatim from
# rc2_analysis/python/src/rc2_glm/plots.py:599-614. Kept inline so the script
# stays self-contained (the source names are underscore-private).
VAR_COLORS: dict[str, tuple[float, float, float]] = {
    "Speed": (0.17, 0.63, 0.17),    # green
    "TF": (1.00, 0.50, 0.05),       # orange
    "SF": (0.95, 0.85, 0.10),       # yellow
    "OR": (0.84, 0.15, 0.16),       # red
    "ME_face": (0.40, 0.20, 0.55),  # purple
}
# variable -> selection-flag column in glm_model_comparison.csv
VAR_FLAG: dict[str, str] = {
    "Speed": "time_is_speed_tuned",
    "TF": "time_is_tf_tuned",
    "SF": "time_is_sf_tuned",
    "OR": "time_is_or_tuned",
    "ME_face": "time_is_me_face_tuned",
}
# Speed/TF/SF/OR are always columns; ME_face is included only if the run
# actually selected it for at least one cluster.
STIMULUS_VARS = ["Speed", "TF", "SF", "OR"]
ALL_VARIABLES = [*STIMULUS_VARS, "ME_face"]
CANONICAL_LAYERS = ["VISp1", "VISp2/3", "VISp4", "VISp5", "VISp6a", "VISp6b"]


def _as_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1"))


def _probes_str(probe_ids: list[str]) -> str:
    """'CAA-1123243_rec1' ... -> 'CAA-1123243/244/466/467'."""
    nums = sorted({p.split("CAA-1123")[-1].split("_")[0] for p in probe_ids})
    return "CAA-1123" + "/".join(nums)


def load_clusters(run_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    """Join glm_model_comparison + prefilter region; detect selected variables."""
    cmp_df = pd.read_csv(run_dir / "glm_model_comparison.csv")
    pre_df = pd.read_csv(run_dir / "prefilter_decision_tree.csv")
    df = cmp_df.merge(
        pre_df[["probe_id", "cluster_id", "region"]],
        on=["probe_id", "cluster_id"], how="inner", validate="one_to_one",
    )
    variables: list[str] = []
    for var in ALL_VARIABLES:
        col = VAR_FLAG[var]
        if col not in df.columns:
            continue
        df[col] = _as_bool(df[col])
        # stimulus vars always shown; ME_face only if the run fit it
        if var in STIMULUS_VARS or df[col].any():
            variables.append(var)
    keep = ["probe_id", "cluster_id", "region", *(VAR_FLAG[v] for v in variables)]
    return df[keep], variables


def attach_depth(df: pd.DataFrame) -> pd.DataFrame:
    """Read continuous cortical depth from each probe's formatted .mat.

    Cross-checks the reader's region against the prefilter CSV region — a
    mismatch means the cluster_id -> array-index mapping is wrong.
    """
    records: list[tuple[str, int, float]] = []
    for probe_id, grp in df.groupby("probe_id"):
        mat_path = FORMATTED_DIR / f"{probe_id}.mat"
        if not mat_path.exists():
            raise FileNotFoundError(f"formatted data missing: {mat_path}")
        with FormattedDataReader(mat_path) as reader:
            id_to_idx = {int(c): i for i, c in enumerate(reader.cluster_ids())}
            for _, row in grp.iterrows():
                cid = int(row["cluster_id"])
                if cid not in id_to_idx:
                    raise RuntimeError(
                        f"{probe_id}: cluster_id {cid} absent from formatted data"
                    )
                idx = id_to_idx[cid]
                region_reader = reader.cluster_region(idx)
                if region_reader != row["region"]:
                    raise RuntimeError(
                        f"{probe_id} cluster {cid}: region mismatch — "
                        f"CSV={row['region']!r} reader={region_reader!r}"
                    )
                records.append((probe_id, cid, reader.cluster_depth(idx)))
    depth_df = pd.DataFrame(records, columns=["probe_id", "cluster_id", "depth_um"])
    return df.merge(depth_df, on=["probe_id", "cluster_id"], how="left")


def layer_geometry(
    df: pd.DataFrame,
) -> tuple[list[str], dict[str, float], list[float], bool]:
    """Per-layer mean depth, inter-layer boundaries, and depth orientation."""
    present = [L for L in CANONICAL_LAYERS if L in set(df["region"])]
    mean_depth = {
        L: float(df.loc[df["region"] == L, "depth_um"].mean()) for L in present
    }
    boundaries = [
        (mean_depth[u] + mean_depth[d]) / 2.0
        for u, d in zip(present[:-1], present[1:])
    ]
    # ascending == depth grows from pia into cortex (superficial layer shallowest)
    ascending = all(
        mean_depth[u] <= mean_depth[d]
        for u, d in zip(present[:-1], present[1:])
    )
    return present, mean_depth, boundaries, ascending


def bin_and_count(
    df: pd.DataFrame, bin_um: float, variables: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Fixed-width depth bins; per (bin, variable) selection count + bin totals."""
    d = df["depth_um"].to_numpy()
    lo = np.floor(d.min() / bin_um) * bin_um
    hi = np.ceil(d.max() / bin_um) * bin_um
    edges = np.arange(lo, hi + bin_um, bin_um)
    centers = (edges[:-1] + edges[1:]) / 2.0
    bin_idx = np.clip(np.digitize(d, edges) - 1, 0, len(centers) - 1)
    df = df.assign(_bin=bin_idx)
    n_per_bin = np.array(
        [(bin_idx == b).sum() for b in range(len(centers))], dtype=float
    )
    counts = {
        var: np.array(
            [df.loc[df["_bin"] == b, VAR_FLAG[var]].sum() for b in range(len(centers))],
            dtype=float,
        )
        for var in variables
    }
    return edges, centers, n_per_bin, counts


def render(
    edges: np.ndarray,
    centers: np.ndarray,
    n_per_bin: np.ndarray,
    counts: dict[str, np.ndarray],
    depths: np.ndarray,
    variables: list[str],
    present: list[str],
    mean_depth: dict[str, float],
    boundaries: list[float],
    ascending: bool,
    bin_um: float,
    n_clusters: int,
    label: str,
    run_name: str,
    probes_str: str,
) -> Path:
    bin_h = float(edges[1] - edges[0])
    depth_lo, depth_hi = float(edges[0]), float(edges[-1])
    fractions = {
        var: np.divide(
            counts[var], n_per_bin,
            out=np.zeros_like(counts[var]), where=n_per_bin > 0,
        )
        for var in variables
    }
    count_xmax = max((c.max() for c in counts.values()), default=1.0)

    # alternating faint bands between consecutive layer boundaries
    band_edges = sorted([depth_lo, *boundaries, depth_hi])
    spread_color = (0.25, 0.25, 0.27)

    def _draw_layers(ax) -> None:
        for i in range(len(band_edges) - 1):
            if i % 2 == 0:
                ax.axhspan(
                    band_edges[i], band_edges[i + 1],
                    color="black", alpha=0.05, zorder=0,
                )
        for b in boundaries:
            ax.axhline(b, color="black", ls="--", lw=0.6, zorder=1)

    n_vars = len(variables)
    fig = plt.figure(figsize=(3.06 * (0.72 + n_vars), 8.5), dpi=140)
    mosaic = [
        ["spread", *[f"{v}|c" for v in variables]],
        ["spread", *[f"{v}|f" for v in variables]],
    ]
    axd = fig.subplot_mosaic(
        mosaic, sharey=True, width_ratios=[0.72, *([1] * n_vars)],
    )

    # --- leftmost column: every cluster as a depth-positioned dot ---
    ax_s = axd["spread"]
    _draw_layers(ax_s)
    jitter = np.random.default_rng(0).uniform(0.24, 0.92, size=len(depths))
    ax_s.scatter(
        jitter, depths, s=18, color=spread_color,
        edgecolor="white", linewidth=0.4, alpha=0.85, zorder=3,
    )
    ax_s.set_xlim(0, 1)
    ax_s.set_xticks([])
    ax_s.set_title("clusters", fontsize=12, fontweight="bold",
                   color=spread_color, pad=8)
    ax_s.text(0.94, 0.96, f"n={n_clusters}", transform=ax_s.transAxes,
              ha="right", va="top", fontsize=8.5, fontweight="bold",
              color=spread_color)
    ax_s.set_ylabel("Cortical depth (µm)", fontsize=10)
    ax_s.tick_params(labelsize=8)

    # --- per-variable count (top) + fraction (bottom) panels ---
    for var in variables:
        for suffix, data, is_count in (
            ("c", counts, True), ("f", fractions, False),
        ):
            ax = axd[f"{var}|{suffix}"]
            _draw_layers(ax)
            ax.barh(
                centers, data[var], height=bin_h * 0.92,
                color=VAR_COLORS[var], edgecolor="white", linewidth=0.4,
                zorder=3,
            )
            ax.tick_params(labelsize=8, labelleft=False)
            if is_count:
                total = int(counts[var].sum())
                ax.set_title(
                    var, color=VAR_COLORS[var], fontsize=12,
                    fontweight="bold", pad=8,
                )
                ax.text(
                    0.96, 0.96, f"n={total}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8.5,
                    color=VAR_COLORS[var], fontweight="bold",
                )
                ax.set_xlim(0, count_xmax * 1.12)
                ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
                ax.set_xlabel("clusters selecting", fontsize=9)
            else:
                ax.set_xlim(0, 1.0)
                ax.set_xlabel("fraction of depth bin", fontsize=9)

    # depth orientation: superficial layer at the top of every panel (sharey)
    if ascending:
        ax_s.set_ylim(depth_hi, depth_lo)
    else:
        ax_s.set_ylim(depth_lo, depth_hi)

    # layer labels in each band, on the cluster-spread panel
    trans = blended_transform_factory(ax_s.transAxes, ax_s.transData)
    for layer in present:
        ax_s.text(
            0.04, mean_depth[layer], layer, transform=trans,
            fontsize=8, fontweight="semibold", va="center", ha="left",
            zorder=4,
            bbox=dict(boxstyle="round,pad=0.18", fc="white",
                      ec="none", alpha=0.8),
        )

    fig.suptitle(
        "Cortical depth × GLM variable selection — "
        f"{run_name} ({probes_str}), {n_clusters} clusters, "
        f"{bin_um:g} µm depth bins",
        fontsize=13, y=0.985,
    )
    fig.text(
        0.5, 0.012,
        "Layer boundaries are midpoints of per-layer mean depth "
        "(data-derived from the labelled clusters, not anatomical ground truth).",
        ha="center", fontsize=8, style="italic", color="0.4",
    )
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.955), w_pad=2.2)

    # row labels, centred in the gap between the spread panel and the bars
    first = variables[0]
    gap_x = (axd["spread"].get_position().x1
             + axd[f"{first}|c"].get_position().x0) / 2.0
    for key, text in ((f"{first}|c", "RAW COUNTS"),
                      (f"{first}|f", "FRACTION OF DEPTH BIN")):
        pos = axd[key].get_position()
        fig.text(
            gap_x, pos.y0 + pos.height / 2.0, text, rotation=90,
            va="center", ha="center", fontsize=10.5, fontweight="bold",
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_base = OUT_DIR / f"depth_variable_selection_{label}"
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    return out_base


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir", type=Path, default=DEFAULT_RUN_DIR,
        help="GLM run directory (holds glm_model_comparison.csv + "
             "prefilter_decision_tree.csv). Default: the 3-probe ME-aware run.",
    )
    parser.add_argument(
        "--label", default="3probe",
        help="Output filename slug: depth_variable_selection_<label>.{pdf,png}.",
    )
    parser.add_argument(
        "--run-name", default="3-probe ME-aware run",
        help="Human-readable run name for the figure title.",
    )
    parser.add_argument(
        "--depth-bin-um", type=float, default=50.0,
        help="Depth-bin width in microns (default 50).",
    )
    args = parser.parse_args()
    run_dir = args.run_dir.expanduser()

    print(f"[1/5] loading GLM CSVs from {run_dir}")
    df, variables = load_clusters(run_dir)
    probes = sorted(df["probe_id"].unique())
    print(
        f"      {len(df)} clusters across {len(probes)} probes: {probes}"
    )
    print(f"      variables in this run: {variables}")

    print("[2/5] attaching continuous cortical depth from formatted .mat")
    df = attach_depth(df)
    n_nan = int(df["depth_um"].isna().sum())
    if n_nan:
        print(f"      WARNING: {n_nan} clusters have NaN depth — dropping them")
        df = df.dropna(subset=["depth_um"]).reset_index(drop=True)
    print(
        f"      depth range: {df['depth_um'].min():.1f} – "
        f"{df['depth_um'].max():.1f} µm  (region cross-check passed)"
    )
    print("      per-variable selection totals:")
    for var in variables:
        print(f"        {var:8s} {int(df[VAR_FLAG[var]].sum()):3d}/{len(df)}")

    print("[3/5] deriving layer geometry")
    present, mean_depth, boundaries, ascending = layer_geometry(df)
    orient = "from pia (deeper = larger)" if ascending else "from probe tip"
    print(f"      depth orientation: {orient}")
    for layer in present:
        n_layer = int((df["region"] == layer).sum())
        print(f"        {layer:9s} mean depth {mean_depth[layer]:7.1f} µm  (n={n_layer})")
    print(f"      inter-layer boundaries (µm): {[round(b, 1) for b in boundaries]}")

    print(f"[4/5] binning depth at {args.depth_bin_um:g} µm")
    edges, centers, n_per_bin, counts = bin_and_count(
        df, args.depth_bin_um, variables
    )
    print(f"      {len(centers)} depth bins, {int(n_per_bin.sum())} clusters placed")

    print(f"[5/5] rendering 2×{len(variables) + 1} faceted figure "
          f"(cluster-spread + {len(variables)} variables)")
    out_base = render(
        edges, centers, n_per_bin, counts, df["depth_um"].to_numpy(),
        variables, present, mean_depth, boundaries, ascending,
        args.depth_bin_um, len(df), args.label, args.run_name,
        _probes_str(probes),
    )
    print(f"      wrote {out_base.with_suffix('.png')}")
    print(f"      wrote {out_base.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
