"""Prompt-06 (2026-04-30): face motion energy as a Phase-1 GLM candidate.

Runs ``rc2-glm`` on the 4 passive probes with the new defaults:

  - ``--no-history``       (history → stand-by)
  - ``--with-onset-kernel`` (onset → restored)
  - ME_face on by default  (no flag needed; new default since 2026-04-30)

Each probe lands in
  ``~/local_data/motion_clouds/figures/glm/exploration/face_me_onset_no_history/_runs/<probe>_rec1/``

Per-probe outputs (model_comparison, selection_history, coefficients,
diagnostics, run.log) come from the pipeline; this driver also emits a
single aggregate ``camera_availability.csv`` summarising which probes
have a real camera trace vs the placeholder.

After all four probes finish, it pulls per-cluster Selected CV-bps and
selected_vars into one ``glm_model_comparison_aggregated.csv``. The
3-way comparison vs ``glm/current/`` and ``glm/legacy_with_onset/``
lives in a separate script (``compare_face_me_vs_current.py``).

Usage::

    /Users/laura/mambaforge3/envs/rc2_analysis/bin/python -m \\
        scripts.run_face_me_onset_no_history

Both ``synecdoche`` (HOME=/Users/laura) and ``metonymy`` (HOME=/Users/lauraporta)
work — HOME-relative paths throughout, no machine-specific hardcodes.
"""
from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("face_me_onset_no_history")

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
HOME = Path.home()
ROOT = HOME / "local_data" / "motion_clouds"
FORMATTED_DIR = ROOT / "formatted_data"
CLUSTER_FILTER = (
    ROOT / "figures" / "matlab_reference" / "prefilter_decision_tree.csv"
)
OUT_ROOT = (
    ROOT / "figures" / "glm" / "exploration" / "face_me_onset_no_history"
)
PROBES = (
    "CAA-1123243_rec1",
    "CAA-1123244_rec1",
    "CAA-1123466_rec1",
    "CAA-1123467_rec1",
)
# Conda env named per ``prompts/rc2-analysis/_git-workflow.md``. Resolve
# the interpreter via ``$CONDA_PREFIX/envs/...`` when running inside a
# conda activation; else fall back to the per-machine canonical path.
PYTHON = (
    Path(os.environ["CONDA_PREFIX"]).parent.parent  # base envs dir
    / "envs" / "rc2_analysis" / "bin" / "python"
    if os.environ.get("CONDA_PREFIX")
    else HOME / "mambaforge3" / "envs" / "rc2_analysis" / "bin" / "python"
)


# ----------------------------------------------------------------------
# Camera availability matrix
# ----------------------------------------------------------------------
def emit_camera_availability_matrix() -> pd.DataFrame:
    """Inspect each probe's session.camera0 / camera1 / camera_t.

    Real camera traces are (1, ~3.7e5) at ~60 Hz. Sessions without a
    recorded camera write a (2,) all-zero placeholder (observed
    2026-04-30 on CAA-1123244_rec1). Surface this BEFORE fitting so the
    user knows which probes will produce no ME_face columns.
    """
    rows: list[dict] = []
    for probe in PROBES:
        mat = FORMATTED_DIR / f"{probe}.mat"
        with h5py.File(mat, "r") as f:
            sess = f["sessions"]
            for key in ("camera0", "camera1"):
                if key not in sess:
                    rows.append({
                        "probe_id": probe, "camera": key,
                        "n_samples": 0, "is_placeholder": True,
                        "first_5": "", "is_real_camera": False,
                    })
                    continue
                arr = np.asarray(sess[key][()]).ravel()
                # (2,) all-zero is the observed placeholder pattern.
                is_placeholder = (arr.size < 1000) or bool(np.all(arr == 0.0))
                rows.append({
                    "probe_id": probe, "camera": key,
                    "n_samples": int(arr.size),
                    "is_placeholder": is_placeholder,
                    "first_5": ",".join(f"{x:.4f}" for x in arr[:5]),
                    "is_real_camera": not is_placeholder,
                })
            cam_t = (
                np.asarray(sess["camera_t"][()]).ravel()
                if "camera_t" in sess else np.array([])
            )
            rows.append({
                "probe_id": probe, "camera": "camera_t",
                "n_samples": int(cam_t.size),
                "is_placeholder": cam_t.size < 1000,
                "first_5": (
                    ",".join(f"{x:.4f}" for x in cam_t[:5])
                    if cam_t.size else ""
                ),
                "is_real_camera": cam_t.size >= 1000,
            })
    df = pd.DataFrame(rows)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out_csv = OUT_ROOT / "camera_availability.csv"
    df.to_csv(out_csv, index=False)
    log.info("wrote %s", out_csv)

    # Pretty summary in the log.
    log.info("camera availability matrix:")
    pivot = df.pivot(index="probe_id", columns="camera", values="is_real_camera")
    log.info("\n%s", pivot.to_string())
    return df


# ----------------------------------------------------------------------
# Per-probe pipeline run
# ----------------------------------------------------------------------
def run_probe(probe: str) -> Path:
    """Invoke rc2-glm on one probe with the prompt-06 flag set."""
    out_dir = OUT_ROOT / "_runs" / probe
    cmp_path = out_dir / "glm_model_comparison.csv"
    if cmp_path.exists():
        log.info("skipping %s — output already at %s", probe, cmp_path)
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    mat = FORMATTED_DIR / f"{probe}.mat"
    # Stimulus-lookup files are needed so V/VT trials get a real
    # batch_gain (otherwise tf = NaN * mean_speed = NaN, B_tf has NaN
    # values, and the variance-pruning bypass in prediction mode lets
    # the NaN propagate into the design matrix → Additive and
    # FullInteraction cv_bps land as NaN. Verified 2026-04-30 on
    # CAA-1123243).
    mc_sequence = ROOT / "motion_cloud_sequence_250414.mat"
    mc_folders = ROOT / "image_folders.mat"
    cmd = [
        str(PYTHON), "-m", "rc2_glm.pipeline",
        str(mat), str(out_dir),
        "--cluster-filter-csv", str(CLUSTER_FILTER),
        "--mc-sequence", str(mc_sequence),
        "--mc-folders", str(mc_folders),
        "--no-history",          # explicit even though default; survives a re-flip
        "--with-onset-kernel",   # explicit, same reason
        "--n-jobs", "4",
        # Per-cluster plots ARE rendered (no --plot-clusters 0). The
        # cluster_<id>_overview.pdf, cluster_<id>_tuning.pdf, and
        # cluster_<id>_kernels.pdf panels carry the per-cell tuning
        # detail and ME_face main-effect / × Speed interaction shape.
    ]
    log.info("running %s", probe)
    log.info("  cmd: %s", " ".join(cmd))
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        log.error("%s pipeline FAILED (exit %d)", probe, proc.returncode)
        raise RuntimeError(f"{probe} pipeline run failed")
    log.info("  done %s", probe)
    return out_dir


# ----------------------------------------------------------------------
# Aggregate
# ----------------------------------------------------------------------
def aggregate(out_dirs: dict[str, Path]) -> pd.DataFrame:
    """Concatenate per-probe glm_model_comparison.csv into one frame."""
    frames = []
    for probe, out_dir in out_dirs.items():
        cmp_path = out_dir / "glm_model_comparison.csv"
        if not cmp_path.exists():
            log.warning("missing %s — skipping in aggregation", cmp_path)
            continue
        df = pd.read_csv(cmp_path)
        df["probe_id"] = probe
        frames.append(df)
    if not frames:
        raise RuntimeError("no per-probe output to aggregate")
    out = pd.concat(frames, ignore_index=True)

    # Convenience flags surfacing ME_face selection. Robust to absence of
    # the column (e.g., on a probe where ME never made it into Selected,
    # 'time_selected_vars' values may be 'Speed' etc. with no 'ME_face').
    sv = out["time_selected_vars"].fillna("Null")
    out["me_face_selected"] = sv.str.contains("ME_face", regex=False, na=False)
    out["me_face_x_speed_selected"] = sv.str.contains(
        "ME_face_x_Speed", regex=False, na=False
    )
    return out


def render_aggregate_forward_selection(df: pd.DataFrame, out_dir: Path) -> None:
    """4-probe-combined forward-selection summary, mirroring per-probe.

    Calls ``rc2_glm.plots.plot_forward_selection_summary`` on the
    concatenated comparison_df. ME_face / ME_face × Speed bars are
    derived from the selected_vars string by the plot function (works on
    legacy CSVs that don't have the dedicated `time_has_me_face_x_speed`
    column).
    """
    from rc2_glm.plots import plot_forward_selection_summary, save_figure
    fig = plot_forward_selection_summary(df)
    fig.suptitle(
        f"Forward Selection Results — face_me_onset_no_history "
        f"(4 probes, n={len(df)} clusters)",
        fontsize=12, fontweight="bold",
    )
    paths = save_figure(fig, out_dir / "figs" / "forward_selection_summary",
                        fmt="pdf")
    for p in paths:
        log.info("wrote %s", p)


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUT_ROOT / "figs").mkdir(parents=True, exist_ok=True)
    log.info("OUT_ROOT: %s", OUT_ROOT)
    log.info("PYTHON: %s", PYTHON)

    cam_df = emit_camera_availability_matrix()
    cam_face = cam_df[cam_df["camera"] == "camera0"]
    n_real = int(cam_face["is_real_camera"].sum())
    log.info("face cameras with real data: %d/4 probes", n_real)
    if n_real == 0:
        log.error("no probe has real camera0 data — aborting before fits")
        return 1

    out_dirs: dict[str, Path] = {}
    for probe in PROBES:
        out_dirs[probe] = run_probe(probe)

    df = aggregate(out_dirs)
    csv_path = OUT_ROOT / "glm_model_comparison_aggregated.csv"
    df.to_csv(csv_path, index=False)
    log.info("wrote %s (%d rows across %d probes)",
             csv_path, len(df), df["probe_id"].nunique())

    # Aggregate forward-selection summary across all probes (the
    # grouped-across-probes panel the per-probe plots can't show).
    render_aggregate_forward_selection(df, OUT_ROOT)

    # Selection-pattern preview
    print()
    print("=== Selection pattern preview ===")
    sv = df["time_selected_vars"].fillna("Null")
    print(sv.value_counts().head(20).to_string())
    print()
    n_total = len(df)
    n_me = int(df["me_face_selected"].sum())
    n_mexspd = int(df["me_face_x_speed_selected"].sum())
    print(f"clusters total: {n_total}")
    print(f"clusters with ME_face selected: {n_me} ({100*n_me/n_total:.1f}%)")
    print(f"clusters with ME_face_x_Speed selected: {n_mexspd} ({100*n_mexspd/n_total:.1f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
