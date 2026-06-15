"""RF-local SF/OR regressor lookup (goggles, opt-in).

Bridges the Gabor extraction (``gabor_extract_gpu.py`` → one parquet per cloud,
cols ``probe, cluster, rf_type, cx, cy, frame, sf_cpd, or_deg, concentration,
edge, cloud``) to the GLM. For one probe it builds, per *cluster* and per
*cloud*, the per-frame local spatial-frequency ``SF(frame)`` (cpd) and
orientation ``OR(frame)`` (deg, 0–180) that the cluster's receptive field sees
on that motion cloud.

Two design facts baked in here:

- **ON/OFF pooling.** A cluster can carry both a white (ON) and a black (OFF)
  RF. Their local *SF* is essentially identical (cohort ΔSF ~1e-4 cpd), but
  *OR* can diverge (median ~16°, up to ~84°). Per Laura's convention we pool
  the polarities by the **circular mean** for OR (plain mean for SF, which is
  the same either way). Single-polarity clusters use their one RF.

- **Cloud join on (theta, sf, VX) tokens, not the folder string.** The rendered
  frame folders the extraction ran on carry a ``BV0p100`` velocity-bandwidth
  token; the presentation metadata (``image_folders_goggles.mat``,
  ``StimulusLookup.cloud_name``) carries ``BV0p200``. BV is the cloud's
  velocity bandwidth — it does not change the spatial SF/OR at an RF — and the
  two name sets are a clean bijection on (theta, sf, VX), so we key on those
  three tokens.

The per-frame → per-bin mapping (the velocity-locked frame clock,
``frame = 10·∫|v|dt``) lives in ``time_binning``, which holds the sample-level
velocity; this module only serves the per-cloud frame arrays.
"""

from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# (theta, sf, VX) token regexes — the cloud identity used to join a trial's
# cloud to its extracted SF/OR, robust to the BV-token discrepancy.
_RE_THETA = re.compile(r"theta(-?\d+p\d+)")
_RE_SF = re.compile(r"_sf(\d+p\d+)_")
_RE_VX = re.compile(r"VX(\d+p\d+)")

CloudKey = tuple[str, str, str]


def cloud_key(name: str) -> CloudKey | None:
    """(theta, sf, VX) token key for a motion-cloud folder name, or None."""
    mt, ms, mv = _RE_THETA.search(name), _RE_SF.search(name), _RE_VX.search(name)
    if mt is None or ms is None or mv is None:
        return None
    return (mt.group(1), ms.group(1), mv.group(1))


@dataclass
class RFLookup:
    """Per-(cluster, cloud) RF-local SF(frame)/OR(frame) arrays for one probe."""

    probe_id: str
    n_frames: int
    # (cluster_id, cloud_key) -> (sf_cpd[n_frames], or_deg[n_frames])
    _store: dict[tuple[int, CloudKey], tuple[np.ndarray, np.ndarray]]
    # cloud_key -> (sf_cpd_const, or_deg_const): the per-cloud cohort NOMINAL
    # (mean SF, circular-mean OR over all RF clusters/frames) — the constant
    # stand-in SF/OR for clusters WITHOUT an RF under the "_all" fallback
    # (config.rf_sf_or_nominal_fallback). Empty when the loader didn't build it.
    _nominal: dict[CloudKey, tuple[float, float]] = field(default_factory=dict)

    @property
    def clusters(self) -> set[int]:
        """Cluster ids that have an RF (and so RF-local SF/OR) on this probe."""
        return {cid for (cid, _) in self._store}

    def nominal(self, cloud_name: str | None):
        """Per-cloud cohort ``(sf_cpd, or_deg)`` for a no-RF cluster, or None.

        The "_all"-mode stand-in: a cluster with no identifiable RF is fit with
        this cloud's cohort-typical SF/OR (constant over the trial) instead of
        being dropped — the continuous analogue of the categorical token level.
        """
        if cloud_name is None:
            return None
        key = cloud_key(cloud_name)
        if key is None:
            return None
        return self._nominal.get(key)

    def get(self, cluster_id: int, cloud_name: str | None):
        """``(sf_arr, or_arr)`` for this cluster+cloud, or ``None`` if absent.

        ``None`` when the cluster has no RF, the trial has no cloud, or the
        cloud's tokens aren't in the extraction — the caller then leaves that
        trial's SF/OR undefined (NaN), exactly like a grey-screen bin.
        """
        if cloud_name is None:
            return None
        key = cloud_key(cloud_name)
        if key is None:
            return None
        return self._store.get((int(cluster_id), key))


def load_rf_sf_or(
    parquet_dir: str, probe_id: str, min_concentration: float = 0.0
) -> RFLookup:
    """Build an :class:`RFLookup` for ``probe_id`` from the cohort parquets.

    ``probe_id`` is the GLM/formatted form (``CAA-1124370_rec1_rec2_rec3``);
    the parquet ``probe`` column is the short form (``CAA-1124370``), matched on
    the leading token. Fails loudly (Pattern 10) when no parquet is found or no
    RF survives for the probe — a silent empty lookup would drop the whole
    cohort downstream.
    """
    pdir = os.path.expanduser(parquet_dir)
    files = sorted(glob.glob(os.path.join(pdir, "*.parquet")))
    if not files:
        raise FileNotFoundError(
            f"RF SF/OR parquet dir has no *.parquet: {pdir}. Run "
            "gabor_extract_gpu.py (cohort) and mirror its output here."
        )
    probe_short = probe_id.split("_")[0]

    store: dict[tuple[int, CloudKey], tuple[np.ndarray, np.ndarray]] = {}
    n_frames_seen: set[int] = set()
    for fp in files:
        cloud = os.path.basename(fp)[: -len(".parquet")]
        key = cloud_key(cloud)
        if key is None:
            continue
        d = pd.read_parquet(
            fp, columns=["probe", "cluster", "rf_type", "frame", "sf_cpd",
                         "or_deg", "concentration"]
        )
        d = d[d["probe"].astype(str) == probe_short]
        if min_concentration > 0.0:
            d = d[d["concentration"] >= min_concentration]
        if d.empty:
            continue
        # Numeric cluster ids only (RF metrics are numeric); drop the rest.
        cid = pd.to_numeric(d["cluster"], errors="coerce")
        d = d.assign(cluster_i=cid).dropna(subset=["cluster_i"])
        d["cluster_i"] = d["cluster_i"].astype(int)
        # Circular components for the ON/OFF OR pool (orientation is
        # π-periodic, hence the doubled angle).
        two_or = np.radians(2.0 * d["or_deg"].to_numpy(dtype=np.float64))
        d = d.assign(_c2=np.cos(two_or), _s2=np.sin(two_or))
        # Pool polarities per (cluster, frame): mean SF, vector-mean OR.
        agg = (
            d.groupby(["cluster_i", "frame"], sort=True)
            .agg(sf=("sf_cpd", "mean"), c2=("_c2", "mean"), s2=("_s2", "mean"))
            .reset_index()
        )
        agg["or_deg"] = np.degrees(np.arctan2(agg["s2"], agg["c2"]) / 2.0) % 180.0
        for cluster_i, g in agg.groupby("cluster_i", sort=False):
            frames = g["frame"].to_numpy(dtype=int)
            n = int(frames.max()) + 1
            n_frames_seen.add(n)
            sf_arr = np.full(n, np.nan, dtype=np.float64)
            or_arr = np.full(n, np.nan, dtype=np.float64)
            sf_arr[frames] = g["sf"].to_numpy(dtype=np.float64)
            or_arr[frames] = g["or_deg"].to_numpy(dtype=np.float64)
            store[(int(cluster_i), key)] = (sf_arr, or_arr)

    if not store:
        raise ValueError(
            f"no RF SF/OR rows for probe {probe_short} in {pdir} "
            f"(min_concentration={min_concentration}). Check the probe id / "
            "concentration gate — an empty lookup would drop the cohort."
        )
    n_frames = max(n_frames_seen) if n_frames_seen else 0

    # Per-cloud cohort NOMINAL (SF mean, OR circular-mean over all RF
    # clusters+frames) — the constant stand-in for no-RF clusters in "_all"
    # mode. SF barely varies across RFs so this recovers the ~3 token SF levels;
    # OR is vector-averaged (π-periodic) to the ~4 token orientations.
    by_cloud: dict[CloudKey, list[tuple[np.ndarray, np.ndarray]]] = {}
    for (_, key), arrs in store.items():
        by_cloud.setdefault(key, []).append(arrs)
    nominal: dict[CloudKey, tuple[float, float]] = {}
    for key, arrs in by_cloud.items():
        sf_all = np.concatenate([a[0] for a in arrs])
        or_all = np.concatenate([a[1] for a in arrs])
        sf_fin = sf_all[np.isfinite(sf_all)]
        or_fin = or_all[np.isfinite(or_all)]
        if sf_fin.size == 0 or or_fin.size == 0:
            continue
        two = np.radians(2.0 * or_fin)
        or_deg = float(np.degrees(
            np.arctan2(np.sin(two).mean(), np.cos(two).mean()) / 2.0) % 180.0)
        nominal[key] = (float(sf_fin.mean()), or_deg)

    return RFLookup(probe_id=probe_id, n_frames=n_frames, _store=store,
                    _nominal=nominal)
