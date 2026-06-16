"""Merge first-pass + 2nd-pass goggle RF metrics and impute the no-blob targets.

The goggle RF pipeline now has two STA passes:

- **first pass** (2.5 sigma): the canonical ``<probe>_rf_metrics.csv`` (no
  ``source`` column) -- the clusters that got an RF at the strict threshold.
- **2nd pass** (1.5 sigma, ``calculate_rf_mouse_goggles_2ndpass.m``): a relaxed
  remap of V-responsive no-RF target clusters -> ``<probe>_rf_metrics_2ndpass.csv``,
  carrying a ``source`` column (``2ndpass`` for a measured blob, ``nofit`` for a
  target that produced no blob even at 1.5 sigma).

This script bridges the two into one merged ``<probe>_rf_metrics.csv`` (the exact
name the cohort Gabor extraction globs -- ``gabor_extract_gpu.py``), and for every
``nofit`` target IMPUTES an RF centre as the mean centroid of the MEASURED RFs on
the **same probe** in the **same VISp layer** (the ``region`` column already
resolves layer: VISp4 / VISp5 / VISp6a / ...). Imputed rows are flagged
``rf_type='imputed'`` / ``source='imputed'`` so they stay distinguishable from
measured RFs all the way into the parquet (which carries ``rf_type``) and the GLM.

A target whose layer has NO measured RF on its probe cannot be imputed -- it is
dropped with a loud warning rather than given a fabricated centre.

Usage:
  python merge_rf_metrics.py \
      --rf-dir   ~/local_data/motion_clouds/saved_goggles/_rfs \
      --out-dir  ~/local_data/motion_clouds/saved_goggles/_rfs/merged
  (2nd-pass CSVs are read from <rf-dir>/2ndpass/ by default.)
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

# Columns averaged when imputing a no-blob target from its layer cohort.
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


def merge_probe(rf_dir: str, secondpass_dir: str, probe_stem: str) -> pd.DataFrame:
    """Merged + imputed RF table for one probe (``probe_stem`` = file stem)."""
    first_fp = os.path.join(rf_dir, f"{probe_stem}_rf_metrics.csv")
    second_fp = os.path.join(secondpass_dir, f"{probe_stem}_rf_metrics_2ndpass.csv")

    parts = []
    if os.path.exists(first_fp):
        parts.append(_read_csv(first_fp, "firstpass"))
    second = _read_csv(second_fp, "2ndpass") if os.path.exists(second_fp) else None
    if second is not None:
        parts.append(second)
    if not parts:
        raise FileNotFoundError(f"no RF CSVs for {probe_stem} in {rf_dir}")

    allrows = pd.concat(parts, ignore_index=True)

    # Measured rows = a real blob (white/black) from either pass. nofit rows
    # (rf_type 'none'/source 'nofit', NaN centroid) are the imputation targets.
    is_nofit = (allrows["source"].astype(str) == "nofit") | \
               (allrows["rf_type"].astype(str) == "none")
    measured = allrows[~is_nofit].copy()
    measured = measured.dropna(subset=["centroid_azimuth_pixels"])
    # De-dup a (cluster, rf_type) seen in both passes -> keep the relaxed 2nd-pass.
    measured["_pref"] = (measured["source"] == "2ndpass").astype(int)
    measured = (measured.sort_values("_pref")
                .drop_duplicates(["cluster_id", "rf_type"], keep="last")
                .drop(columns="_pref"))

    # Per-layer (region) mean centroid over measured RFs on THIS probe.
    layer_mean = measured.groupby("region")[_AVG_COLS].mean()

    imputed_rows = []
    nofit = allrows[is_nofit]
    for _, r in nofit.iterrows():
        region = r["region"]
        cid = r["cluster_id"]
        if region not in layer_mean.index:
            print(f"  [WARN] {probe_stem} cluster {cid}: layer '{region}' has no "
                  f"measured RF on this probe -> cannot impute, DROPPED")
            continue
        row = {"cluster_id": cid, "region": region, "rf_type": "imputed",
               "source": "imputed"}
        for col in _AVG_COLS:
            row[col] = float(layer_mean.loc[region, col])
        for col in _SIZE_COLS:
            row[col] = np.nan
        imputed_rows.append(row)
        print(f"  [impute] {probe_stem} cluster {cid} ({region}): "
              f"az={row['centroid_azimuth_pixels']:.1f}px "
              f"el={row['centroid_elevation_pixels']:.1f}px "
              f"(layer mean of {int((measured['region'] == region).sum())} measured RFs)")

    out = pd.concat(
        [measured[_OUT_COLS], pd.DataFrame(imputed_rows, columns=_OUT_COLS)],
        ignore_index=True,
    )
    n_meas = (out["source"] != "imputed").sum()
    n_imp = (out["source"] == "imputed").sum()
    print(f"  {probe_stem}: {n_meas} measured rows + {n_imp} imputed rows "
          f"({out['cluster_id'].nunique()} clusters)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rf-dir",
                    default="~/local_data/motion_clouds/saved_goggles/_rfs",
                    help="dir with canonical first-pass <probe>_rf_metrics.csv")
    ap.add_argument("--secondpass-dir", default="",
                    help="dir with <probe>_rf_metrics_2ndpass.csv (default <rf-dir>/2ndpass)")
    ap.add_argument("--out-dir",
                    default="~/local_data/motion_clouds/saved_goggles/_rfs/merged",
                    help="output dir for merged <probe>_rf_metrics.csv")
    ap.add_argument("--probes", nargs="*",
                    default=["CAA-1124370_rec1_rec2_rec3", "CAA-1124371_rec1_rec2_rec3"])
    a = ap.parse_args()

    rf_dir = os.path.expanduser(a.rf_dir)
    second_dir = os.path.expanduser(a.secondpass_dir) if a.secondpass_dir \
        else os.path.join(rf_dir, "2ndpass")
    out_dir = os.path.expanduser(a.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    for probe in a.probes:
        print(f"[{probe}]")
        merged = merge_probe(rf_dir, second_dir, probe)
        out_fp = os.path.join(out_dir, f"{probe}_rf_metrics.csv")
        merged.to_csv(out_fp, index=False)
        print(f"  [saved] {out_fp}\n")


if __name__ == "__main__":
    main()
