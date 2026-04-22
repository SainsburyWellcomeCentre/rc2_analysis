"""Loaders for MATLAB-precomputed 5%-quantile tuning-curve bin edges.

The MATLAB pipeline caches per-trial-group velocity/TF bin edges beside each
formatted recording so Speed and TF tuning curves use bins that each hold
~5% of motion-masked samples (MATLAB ``lib/classes/tuning/VelocityBins.m``,
``prc_per_bin = 5``).

Layout beside the formatted .mat::

    <formatted_data_dir>/
        <probe_id>.mat
        csvs/
            tuning_curves/<probe_id>.mat     # Speed edges, 3 trial groups
            tf_tuning_curves/<probe_id>.mat  # TF edges, 3 trial groups

Each file is MATLAB v7.3 (HDF5): ``trial_groups`` (3×1 char-array refs) and
``tuning_curves`` (3×1 group refs, one per trial group). Each group holds
``bin_edges`` (96×1 object refs — one per cluster, but identical within a
trial-group/recording since edges are per-recording quantiles).

``prc_per_bin`` is asserted at 5.0 when loading — if the MATLAB pipeline
changes the convention this module will raise so Python stays in sync.
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
class PrecomputedBinEdges:
    """5% quantile bin edges per trial group for Speed and TF panels."""

    speed_by_group: dict[str, np.ndarray] = field(default_factory=dict)
    tf_by_group: dict[str, np.ndarray] = field(default_factory=dict)

    def speed_edges(self, condition: str) -> np.ndarray | None:
        """Return edges for the Speed panel for ``condition`` (``'VT'``/``'T_Vstatic'``)."""
        return self.speed_by_group.get(condition)

    def tf_edges(self, condition: str) -> np.ndarray | None:
        """Return edges for the TF panel for ``condition`` (``'VT'``/``'V'``)."""
        return self.tf_by_group.get(condition)


def load_precomputed_bin_edges(
    formatted_mat_path: Path | str,
) -> PrecomputedBinEdges | None:
    """Locate ``csvs/{tuning_curves,tf_tuning_curves}/<probe>.mat`` and load edges.

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
            "precomputed bin edges not found for probe %r "
            "(expected %s and %s) — falling back to linspace bins",
            probe_stem, speed_path, tf_path,
        )
        return None

    try:
        speed_by_group = _read_edges_file(speed_path)
        tf_by_group = _read_edges_file(tf_path)
    except (OSError, KeyError, ValueError) as exc:
        logger.warning(
            "failed to read precomputed bin edges for %r (%s) — "
            "falling back to linspace bins", probe_stem, exc,
        )
        return None

    return PrecomputedBinEdges(
        speed_by_group=speed_by_group,
        tf_by_group=tf_by_group,
    )


def _read_edges_file(path: Path) -> dict[str, np.ndarray]:
    """Return ``{'VT': edges, 'V': edges, 'T_Vstatic': edges}`` from one file."""
    out: dict[str, np.ndarray] = {}
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

            edges_refs = grp["bin_edges"]
            first_edges = np.array(f[edges_refs[0, 0]]).flatten().astype(np.float64)
            out[name] = first_edges
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
