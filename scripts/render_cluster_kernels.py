"""Render per-cluster kernel-shape plots for an existing GLM run.

For each (probe, cluster) in the run's ``glm_coefficients.csv``, builds
a 4×6 panel (Null/Selected/Additive/FullInteraction × Speed/TF/SF/OR/
History/Spd×TF) showing the basis × β kernel that contributes to log λ.
Saves to ``<run_dir>/_runs/<probe>/figs/cluster_<id>_kernels.{pdf,png}``.

Usage::

    python scripts/render_cluster_kernels.py
        [--run-dir ~/local_data/motion_clouds/figures/glm/current]
        [--cluster CLUSTER_ID]   # render just one cluster across probes
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python" / "src"))
from rc2_glm.config import GLMConfig
from rc2_glm.plots import plot_cluster_kernels

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("kernels")

DEFAULT_RUN = Path(
    "/Users/lauraporta/local_data/motion_clouds/figures/glm/current"
)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, default=DEFAULT_RUN,
                   help="Top-level run dir (must contain _runs/<probe>/glm_coefficients.csv)")
    p.add_argument("--cluster", type=int, default=None,
                   help="Restrict to one cluster ID (across all probes)")
    p.add_argument("--bin-width", type=float, default=None,
                   help="Override config.time_bin_width (in seconds). "
                        "Use 0.02 for the 20ms run output.")
    args = p.parse_args()

    runs_root = args.run_dir / "_runs"
    if not runs_root.is_dir():
        log.error("no _runs/ dir under %s", args.run_dir)
        return 1

    config = GLMConfig()
    if args.bin_width is not None:
        # GLMConfig is frozen — work with a replacement
        from dataclasses import replace
        config = replace(config, time_bin_width=args.bin_width)
        log.info("config.time_bin_width overridden → %.3f s", args.bin_width)

    n_total = 0
    for probe_dir in sorted(runs_root.iterdir()):
        if not probe_dir.is_dir():
            continue
        coef_path = probe_dir / "glm_coefficients.csv"
        if not coef_path.is_file():
            continue
        coefs = pd.read_csv(coef_path)
        coefs = coefs[coefs["glm_type"] == "time"]
        figs_dir = probe_dir / "figs"
        figs_dir.mkdir(exist_ok=True)
        clusters = (
            [args.cluster] if args.cluster is not None
            else sorted(coefs["cluster_id"].unique())
        )
        for cid in clusters:
            sub = coefs[coefs["cluster_id"] == cid]
            if sub.empty:
                continue
            try:
                fig = plot_cluster_kernels(
                    probe_dir.name, int(cid), sub, config,
                )
                out = figs_dir / f"cluster_{int(cid)}_kernels.pdf"
                fig.savefig(out)
                fig.savefig(out.with_suffix(".png"), dpi=150)
                plt.close(fig)
                n_total += 1
            except Exception as exc:
                log.warning("cluster %d on %s failed: %s", cid, probe_dir.name, exc)
        log.info("probe %s: rendered %d kernel plots",
                 probe_dir.name, len(clusters))

    log.info("done — wrote %d cluster_<id>_kernels.pdf files under %s",
             n_total, runs_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
