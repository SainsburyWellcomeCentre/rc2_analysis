"""Merge first-pass + 2nd-pass goggle RF metrics and impute the selected-cohort gaps.

The goggle analysis cohort is the formatted file's curated ``selected_clusters``
(anatomy-independent — goggle anatomy is largely ``unknownLocation``). For the
RF-local GLM every selected cluster needs an SF/OR, i.e. an RF centre. This
script builds, per probe, one merged ``<probe>_rf_metrics.csv`` (the name the
cohort Gabor extraction globs, ``gabor_extract_gpu.py``) that:

- keeps every MEASURED RF (first pass 2.5σ + 2nd pass 1.5σ, de-duped, 2nd-pass
  wins) — these are real STA blobs;
- for every cluster in ``selected_clusters`` that has NO measured RF, IMPUTES a
  centre = the mean centroid of measured RFs on the **same probe** in the **same
  VISp layer** (``region``). Imputed rows are flagged ``rf_type=imputed`` /
  ``source=imputed`` so the layer-mean stand-in stays distinguishable from a
  measured RF all the way into the parquet (which carries ``rf_type``).

A selected cluster whose layer has no measured RF on its probe cannot be imputed
— it is reported and left without an RF (it would drop out of the rf_local
cohort) rather than given a fabricated centre.

Usage:
  python merge_rf_metrics.py \
      --rf-dir        ~/local_data/motion_clouds/saved_goggles/_rfs \
      --formatted-dir ~/local_data/motion_clouds/formatted_data_goggles \
      --out-dir       ~/local_data/motion_clouds/saved_goggles/_rfs/merged
  (2nd-pass CSVs read from <rf-dir>/2ndpass/ by default.)
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

from rc2_formatted_data_reader.reader import FormattedDataReader

# Columns averaged when imputing a gap cluster from its layer cohort.
_AVG_COLS = [
    "centroid_azimuth_pixels", "centroid_elevation_pixels",
    "centroid_azimuth_deg", "centroid_elevation_deg",
]
_SIZE_COLS = [
    "size_azimuth_pixels", "size_elevation_pixels",
    "size_azimuth_deg", "size_elevation_deg",
]
_OUT_COLS = ["cluster_id", "region", "rf_type", "source", *_AVG_COLS, *_SIZE_COLS]


def _read_csv(path: str, default_source: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "source" not in df.columns:
        df["source"] = default_source
    return df


def _selected_regions(formatted_dir: str, probe: str) -> dict[int, str]:
    """``{selected_cluster_id: region_str}`` from the probe's formatted .mat."""
    with FormattedDataReader(os.path.join(formatted_dir, f"{probe}.mat")) as r:
        ids = r.cluster_ids()
        out: dict[int, str] = {}
        for idx in r.selected_cluster_indices():
            out[int(ids[idx])] = r.cluster_region(int(idx))
    return out


def merge_probe(rf_dir: str, secondpass_dir: str, formatted_dir: str,
                probe_stem: str) -> pd.DataFrame:
    """Merged + selected-gap-imputed RF table for one probe."""
    first_fp = os.path.join(rf_dir, f"{probe_stem}_rf_metrics.csv")
    second_fp = os.path.join(secondpass_dir, f"{probe_stem}_rf_metrics_2ndpass.csv")

    parts = []
    if os.path.exists(first_fp):
        parts.append(_read_csv(first_fp, "firstpass"))
    if os.path.exists(second_fp):
        parts.append(_read_csv(second_fp, "2ndpass"))
    if not parts:
        raise FileNotFoundError(f"no RF CSVs for {probe_stem} in {rf_dir}")
    allrows = pd.concat(parts, ignore_index=True)

    # Measured rows = a real blob (white/black) from either pass; drop the
    # MATLAB 'nofit'/'none' placeholders (NaN centroid).
    is_placeholder = (allrows["source"].astype(str) == "nofit") | \
                     (allrows["rf_type"].astype(str) == "none")
    measured = allrows[~is_placeholder].dropna(subset=["centroid_azimuth_pixels"]).copy()
    # De-dup a (cluster, rf_type) seen in both passes -> keep the relaxed 2nd-pass.
    measured["_pref"] = (measured["source"] == "2ndpass").astype(int)
    measured = (measured.sort_values("_pref")
                .drop_duplicates(["cluster_id", "rf_type"], keep="last")
                .drop(columns="_pref"))

    # Per-layer (region) mean centroid over measured RFs on THIS probe.
    layer_mean = measured.groupby("region")[_AVG_COLS].mean()
    measured_clusters = set(measured["cluster_id"])

    # Impute every SELECTED cluster that has no measured RF.
    sel_regions = _selected_regions(formatted_dir, probe_stem)
    imputed_rows, uncovered = [], []
    for cid, region in sorted(sel_regions.items()):
        if cid in measured_clusters:
            continue
        if region not in layer_mean.index:
            uncovered.append((cid, region))
            continue
        row = {"cluster_id": cid, "region": region, "rf_type": "imputed",
               "source": "imputed"}
        for col in _AVG_COLS:
            row[col] = float(layer_mean.loc[region, col])
        for col in _SIZE_COLS:
            row[col] = np.nan
        imputed_rows.append(row)

    out = pd.concat(
        [measured[_OUT_COLS], pd.DataFrame(imputed_rows, columns=_OUT_COLS)],
        ignore_index=True,
    )

    sel = set(sel_regions)
    covered = len(sel & set(out["cluster_id"]))
    print(f"  {probe_stem}: selected={len(sel)} | measured RFs total={len(measured_clusters)} "
          f"| imputed gaps={len(imputed_rows)} | selected covered={covered}/{len(sel)}")
    for cid, region in uncovered:
        print(f"     [WARN] selected cluster {cid} ({region}): no measured RF in this "
              f"layer on this probe -> NOT imputed (drops from rf_local cohort)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rf-dir",
                    default="~/local_data/motion_clouds/saved_goggles/_rfs",
                    help="dir with canonical first-pass <probe>_rf_metrics.csv")
    ap.add_argument("--secondpass-dir", default="",
                    help="dir with <probe>_rf_metrics_2ndpass.csv (default <rf-dir>/2ndpass)")
    ap.add_argument("--formatted-dir",
                    default="~/local_data/motion_clouds/formatted_data_goggles",
                    help="dir with <probe>.mat (read selected_clusters + regions)")
    ap.add_argument("--out-dir",
                    default="~/local_data/motion_clouds/saved_goggles/_rfs/merged",
                    help="output dir for merged <probe>_rf_metrics.csv")
    ap.add_argument("--probes", nargs="*",
                    default=["CAA-1124370_rec1_rec2_rec3", "CAA-1124371_rec1_rec2_rec3"])
    a = ap.parse_args()

    rf_dir = os.path.expanduser(a.rf_dir)
    second_dir = os.path.expanduser(a.secondpass_dir) if a.secondpass_dir \
        else os.path.join(rf_dir, "2ndpass")
    formatted_dir = os.path.expanduser(a.formatted_dir)
    out_dir = os.path.expanduser(a.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    for probe in a.probes:
        print(f"[{probe}]")
        merged = merge_probe(rf_dir, second_dir, formatted_dir, probe)
        out_fp = os.path.join(out_dir, f"{probe}_rf_metrics.csv")
        merged.to_csv(out_fp, index=False)
        print(f"  [saved] {out_fp}\n")


if __name__ == "__main__":
    main()
