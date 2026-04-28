"""Loaders for the MATLAB-precomputed Speed / TF tuning-curve cache.

The MATLAB pipeline writes a per-recording cache beside each formatted .mat:

    <formatted_data_dir>/
        <probe_id>.mat
        csvs/
            tuning_curves/<probe_id>.mat     # Speed cache (3 trial groups)
            tf_tuning_curves/<probe_id>.mat  # TF cache    (3 trial groups)

Each cache file is MATLAB v7.3 (HDF5). Top-level layout::

    trial_groups   (3,1) char-array refs   ['VT', 'V', 'T_Vstatic']
    tuning_curves  (3,1) group     refs   one group per trial group, each
                                            holding per-cluster bin edges,
                                            bin centres, per-trial-per-bin
                                            firing rates, stationary baseline.

Per-trial-group group fields:
    bin_edges       (n_clusters, 1) refs → each (n_bins+1,) float
    bin_centres     (n_clusters, 1) refs → each (n_bins,) float
    cluster_id      (n_clusters, 1) refs → each (1,1) integer
    prc_per_bin     (n_clusters, 1) refs → each (1,1); asserted == 5.0
    tuning          (n_clusters, 1) refs → each (n_trials × n_bins,) float
                                                (reshape (n_trials, n_bins))
                                                — per-trial-per-bin firing rate (Hz)
    trial_ids       (n_clusters, 1) refs → each (n_trials,) int
    stationary_fr   (n_clusters, 1) refs → each (n_trials,) float
    stationary_time (n_clusters, 1) refs → each (n_trials,) float
    timing          (n_clusters, 1) refs → each (n_trials × n_bins,) float
                                                — per-trial-per-bin sample
                                                duration (used by MATLAB as
                                                the FR denominator)

This module exposes both the bin edges (used by the GLM time-binning) and
the per-trial-per-bin firing rates (used by the Observed-row plot to match
MATLAB's tuning-curve PDFs exactly — Laura's call 2026-04-27, see
results/motion-clouds/2026-04-27-per-trial-tuning-uncertainty.md).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np

logger = logging.getLogger("rc2_glm.precomputed_bins")

_TRIAL_GROUPS: tuple[str, ...] = ("VT", "V", "T_Vstatic")


@dataclass
class PerConditionTuning:
    """One trial-group's per-cluster tuning cache (Speed *or* TF, one condition).

    All per-cluster arrays are indexed by ``cluster_id`` (the MATLAB cluster
    ID, not a 0-based array index). ``bin_edges``/``bin_centres`` are
    per-recording quantities (identical across clusters), so they live as
    plain arrays.
    """
    bin_edges: np.ndarray            # (n_bins+1,)
    bin_centres: np.ndarray          # (n_bins,)
    cluster_ids: np.ndarray          # (n_clusters,) sorted
    tuning: dict[int, np.ndarray]            # cluster_id → (n_trials, n_bins)
    trial_ids: dict[int, np.ndarray]         # cluster_id → (n_trials,)
    stationary_fr: dict[int, np.ndarray]     # cluster_id → (n_trials,)
    stationary_time: dict[int, np.ndarray]   # cluster_id → (n_trials,)


@dataclass
class PrecomputedBinEdges:
    """Per-recording tuning cache for Speed and TF, one entry per condition.

    Backwards-compatible name (used to be edges-only); now also exposes the
    per-trial-per-bin firing-rate cache that MATLAB plots from. Methods
    ``speed_edges`` / ``tf_edges`` keep working unchanged.
    """
    speed_by_group: dict[str, PerConditionTuning] = field(default_factory=dict)
    tf_by_group: dict[str, PerConditionTuning] = field(default_factory=dict)

    # --- Edge accessors (unchanged signature) ---

    def speed_edges(self, condition: str) -> np.ndarray | None:
        c = self.speed_by_group.get(condition)
        return c.bin_edges if c is not None else None

    def tf_edges(self, condition: str) -> np.ndarray | None:
        c = self.tf_by_group.get(condition)
        return c.bin_edges if c is not None else None

    # --- Centres (new) ---

    def speed_centres(self, condition: str) -> np.ndarray | None:
        c = self.speed_by_group.get(condition)
        return c.bin_centres if c is not None else None

    def tf_centres(self, condition: str) -> np.ndarray | None:
        c = self.tf_by_group.get(condition)
        return c.bin_centres if c is not None else None

    # --- Per-cluster tuning data (new) ---

    def speed_tuning(self, condition: str, cluster_id: int) -> np.ndarray | None:
        """Per-trial-per-bin firing rate (Hz), shape (n_trials, n_bins). None if missing."""
        c = self.speed_by_group.get(condition)
        return None if c is None else c.tuning.get(int(cluster_id))

    def tf_tuning(self, condition: str, cluster_id: int) -> np.ndarray | None:
        c = self.tf_by_group.get(condition)
        return None if c is None else c.tuning.get(int(cluster_id))

    def speed_trial_ids(self, condition: str, cluster_id: int) -> np.ndarray | None:
        c = self.speed_by_group.get(condition)
        return None if c is None else c.trial_ids.get(int(cluster_id))

    def tf_trial_ids(self, condition: str, cluster_id: int) -> np.ndarray | None:
        c = self.tf_by_group.get(condition)
        return None if c is None else c.trial_ids.get(int(cluster_id))

    def stationary_fr(
        self, variable: str, condition: str, cluster_id: int,
    ) -> np.ndarray | None:
        """Per-trial stationary firing rate for ``variable`` ('speed'|'tf')."""
        store = self.speed_by_group if variable == "speed" else self.tf_by_group
        c = store.get(condition)
        return None if c is None else c.stationary_fr.get(int(cluster_id))


def load_precomputed_bin_edges(
    formatted_mat_path: Path | str,
) -> PrecomputedBinEdges | None:
    """Locate ``csvs/{tuning_curves,tf_tuning_curves}/<probe>.mat`` and load.

    Returns None (and logs a warning) if either file is missing — callers
    should fall back to a linspace default so the pipeline still runs on
    recordings where the MATLAB preprocessing step hasn't been executed.
    """
    formatted_mat_path = Path(formatted_mat_path)
    probe_stem = formatted_mat_path.stem
    csvs_dir = formatted_mat_path.parent / "csvs"
    speed_path = csvs_dir / "tuning_curves" / f"{probe_stem}.mat"
    tf_path = csvs_dir / "tf_tuning_curves" / f"{probe_stem}.mat"

    if not speed_path.exists() or not tf_path.exists():
        logger.warning(
            "precomputed bin caches not found for probe %r "
            "(expected %s and %s) — falling back to linspace bins, "
            "Observed-row plots will recompute from raw spikes",
            probe_stem, speed_path, tf_path,
        )
        return None

    try:
        speed_by_group = _read_cache_file(speed_path)
        tf_by_group = _read_cache_file(tf_path)
    except (OSError, KeyError, ValueError) as exc:
        logger.warning(
            "failed to read precomputed cache for %r (%s) — "
            "falling back to linspace bins / recomputed Observed",
            probe_stem, exc,
        )
        return None

    return PrecomputedBinEdges(
        speed_by_group=speed_by_group,
        tf_by_group=tf_by_group,
    )


def _read_cache_file(path: Path) -> dict[str, PerConditionTuning]:
    """Return ``{'VT': PerConditionTuning, ...}`` from one cache file."""
    out: dict[str, PerConditionTuning] = {}
    with h5py.File(path, "r") as f:
        trial_groups = f["trial_groups"]
        tuning_curves = f["tuning_curves"]
        n_groups = trial_groups.shape[0]
        for i in range(n_groups):
            name = _deref_char_array(f, trial_groups[i, 0])
            grp = f[tuning_curves[i, 0]]

            prc_refs = grp["prc_per_bin"]
            prc_val = float(np.array(f[prc_refs[0, 0]]).flatten()[0])
            if not np.isclose(prc_val, 5.0):
                raise ValueError(
                    f"{path} trial_group={name!r}: prc_per_bin={prc_val} "
                    f"(expected 5.0). The MATLAB bin convention has "
                    f"changed — align Python tuning-curve bins before "
                    f"reusing these edges."
                )

            # Bin edges + centres are per-recording (same across clusters).
            edges_refs = grp["bin_edges"]
            centres_refs = grp["bin_centers"]  # MATLAB spelling
            first_edges = np.array(f[edges_refs[0, 0]]).flatten().astype(np.float64)
            first_centres = np.array(f[centres_refs[0, 0]]).flatten().astype(np.float64)

            # Per-cluster fields.
            cid_refs = grp["cluster_id"]
            tuning_refs = grp["tuning"]
            trial_ids_refs = grp["trial_ids"]
            stat_fr_refs = grp["stationary_fr"]
            stat_time_refs = grp["stationary_time"]

            n_clusters = cid_refs.shape[0]
            cluster_ids: list[int] = []
            tuning: dict[int, np.ndarray] = {}
            trial_ids: dict[int, np.ndarray] = {}
            stat_fr: dict[int, np.ndarray] = {}
            stat_time: dict[int, np.ndarray] = {}

            n_bins = first_centres.size
            for j in range(n_clusters):
                cid = int(np.array(f[cid_refs[j, 0]]).flatten()[0])
                cluster_ids.append(cid)

                trial_id_arr = np.array(
                    f[trial_ids_refs[j, 0]]
                ).flatten().astype(int)
                tuning_flat = np.array(
                    f[tuning_refs[j, 0]]
                ).flatten().astype(np.float64)

                # MATLAB stores tuning as (n_trials × n_bins,) flat; the
                # cluster's n_trials = len(trial_id_arr).
                n_trials = trial_id_arr.size
                expected = n_trials * n_bins
                if tuning_flat.size != expected:
                    logger.warning(
                        "%s %s cluster %d: tuning length %d != "
                        "n_trials*n_bins (%d*%d=%d) — dropping cluster",
                        path.name, name, cid, tuning_flat.size,
                        n_trials, n_bins, expected,
                    )
                    cluster_ids.pop()
                    continue
                tuning[cid] = tuning_flat.reshape(n_trials, n_bins)
                trial_ids[cid] = trial_id_arr
                stat_fr[cid] = np.array(
                    f[stat_fr_refs[j, 0]]
                ).flatten().astype(np.float64)
                stat_time[cid] = np.array(
                    f[stat_time_refs[j, 0]]
                ).flatten().astype(np.float64)

            out[name] = PerConditionTuning(
                bin_edges=first_edges,
                bin_centres=first_centres,
                cluster_ids=np.asarray(sorted(cluster_ids), dtype=int),
                tuning=tuning,
                trial_ids=trial_ids,
                stationary_fr=stat_fr,
                stationary_time=stat_time,
            )
    missing = [g for g in _TRIAL_GROUPS if g not in out]
    if missing:
        raise KeyError(
            f"{path}: expected trial groups {_TRIAL_GROUPS}, missing {missing}"
        )
    return out


def _deref_char_array(f: h5py.File, ref) -> str:
    """MATLAB char arrays are uint16 refs; convert to a Python string."""
    arr = np.array(f[ref]).flatten()
    return "".join(chr(int(x)) for x in arr)
