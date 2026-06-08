"""Trial condition mapping and motion-cloud stimulus parsing.

The trial-condition mapping (`trial_sequence` → label) comes from
`PassiveProtocolAlwaysVisSession.m`.

Cloud-name parsing matches `parse_cloud_name` in
`scripts/glm_single_cluster_analysis.m` (lines 5477–5503), with the
SF / orientation / batch_gain dictionaries defined at lines 150–164.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

CONDITION_MAP: dict[int, str] = {1: "VT", 2: "V", 3: "T_Vstatic"}

# Orientation (theta) is shared between the screens and goggles stimulus sets —
# the goggles port (Mateo, master e2a65b5) left the OR branch unchanged.
OR_VALUES: dict[str, float] = {
    "theta0p000": 0.0,
    "theta0p785": np.pi / 4,
    "theta1p571": np.pi / 2,
    "theta-0p785": -np.pi / 4,
}


@dataclass(frozen=True)
class StimulusSet:
    """The motion-cloud parameter grid for one experimental modality.

    Two modalities exist: the original ``screens`` set and the ``goggles``
    set (head-mounted VR, from 2026-03). They share the orientation grid and
    the three-rung gain ladder, but differ in the SF levels and the cloud-name
    patterns that map to each gain rung — see ``SCREENS_STIMULUS`` /
    ``GOGGLES_STIMULUS`` and ``scripts/preprocess/create_tables_by_batch.m``.
    """

    sf_values: dict[str, float]
    batch_patterns: tuple[tuple[str, ...], ...]
    batch_gains: tuple[float, ...]
    exclude_patterns: tuple[str, ...] = ()
    or_values: dict[str, float] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.or_values is None:
            object.__setattr__(self, "or_values", OR_VALUES)


# Original screen-displayed motion clouds (CAA-1123xxx probes).
SCREENS_STIMULUS = StimulusSet(
    sf_values={"sf00p003": 0.003, "sf00p006": 0.006, "sf00p012": 0.012},
    batch_patterns=(
        ("sf00p003_Bsf0p002_VX1p002", "sf00p006_Bsf0p002_VX0p501", "sf00p012_Bsf0p002_VX0p250"),
        ("sf00p003_Bsf0p002_VX2p003", "sf00p006_Bsf0p002_VX1p002", "sf00p012_Bsf0p002_VX0p501"),
        ("sf00p003_Bsf0p002_VX4p006", "sf00p006_Bsf0p002_VX2p003", "sf00p012_Bsf0p002_VX1p002"),
    ),
    batch_gains=(1 / 30, 2 / 30, 4 / 30),
    exclude_patterns=("theta0p000_Btheta3p142_sf00p006_Bsf0p004_VX0p000_BV2p000",),
)

# Goggles motion clouds (CAA-1124370/371, mouse VR goggles, 2026-03).
# Ported from master e2a65b5 (Mateo). New SF grid (008/016/032) and VX values,
# but each batch reproduces the same SF×VX product as the matching screens
# batch (≈0.003 / 0.006 / 0.012), so the gain ladder is preserved.
# NOTE: batch_gains are carried over from screens on that design-equivalence
# argument; confirm with Mateo that the goggles speed→TF coupling gains are the
# same three rungs and not rescaled for the goggles display geometry.
GOGGLES_STIMULUS = StimulusSet(
    sf_values={"sf00p008": 0.008, "sf00p016": 0.016, "sf00p032": 0.032},
    batch_patterns=(
        ("sf00p008_Bsf0p005_VX0p382", "sf00p016_Bsf0p005_VX0p191", "sf00p032_Bsf0p005_VX0p095"),
        ("sf00p008_Bsf0p005_VX0p764", "sf00p016_Bsf0p005_VX0p382", "sf00p032_Bsf0p005_VX0p191"),
        ("sf00p008_Bsf0p005_VX1p528", "sf00p016_Bsf0p005_VX0p764", "sf00p032_Bsf0p005_VX0p382"),
    ),
    batch_gains=(1 / 30, 2 / 30, 4 / 30),
    exclude_patterns=(),  # the screens control stim is screens-only
)

# Backward-compatible module-level aliases (screens). Existing callers and
# tests import these / call parse_cloud_name with no stimulus set → screens.
SF_VALUES: dict[str, float] = SCREENS_STIMULUS.sf_values
BATCH_PATTERNS = [list(p) for p in SCREENS_STIMULUS.batch_patterns]
BATCH_GAINS: list[float] = list(SCREENS_STIMULUS.batch_gains)
EXCLUDE_PATTERNS: list[str] = list(SCREENS_STIMULUS.exclude_patterns)

_VX_RE = re.compile(r"VX(\d+p\d+)")


def parse_cloud_name(
    name: str, stimulus_set: StimulusSet = SCREENS_STIMULUS
) -> dict[str, float | bool]:
    """Parse SF, orientation, VX, batch gain from a motion-cloud name.

    Returns a dict with keys: sf, orientation, vx, batch_gain, excluded.
    Values are NaN for SF / orientation / VX / batch_gain when unmatched.
    ``stimulus_set`` selects the screens (default) or goggles parameter grid.
    """
    if any(p in name for p in stimulus_set.exclude_patterns):
        return {"sf": np.nan, "orientation": np.nan, "vx": np.nan,
                "batch_gain": np.nan, "excluded": True}

    sf = next((v for k, v in stimulus_set.sf_values.items() if k in name), np.nan)
    orientation = next(
        (v for k, v in stimulus_set.or_values.items() if name.startswith(k)), np.nan
    )

    vx = np.nan
    m = _VX_RE.search(name)
    if m:
        vx = float(m.group(1).replace("p", "."))

    batch_gain = np.nan
    for patterns, gain in zip(stimulus_set.batch_patterns, stimulus_set.batch_gains):
        if any(p in name for p in patterns):
            batch_gain = gain
            break

    return {"sf": sf, "orientation": orientation, "vx": vx,
            "batch_gain": batch_gain, "excluded": False}


class StimulusLookup:
    """Map per-trial stimulus parameters from motion-cloud metadata.

    Loads `motion_cloud_sequence_*.mat` and `image_folders.mat` (both v5
    MATLAB files, read with `scipy.io.loadmat`) and exposes per-trial
    stimulus parameters via `stimulus_params(trial_id)`.

    `trial_id` is 1-based (matches the trial_id stored in the formatted
    .mat file).
    """

    def __init__(
        self,
        mc_sequence_path: str | Path,
        mc_folders_path: str | Path,
        stimulus_set: StimulusSet = SCREENS_STIMULUS,
    ):
        from scipy.io import loadmat

        self._stimulus_set = stimulus_set
        seq = loadmat(str(mc_sequence_path))
        self._presentation: np.ndarray = np.asarray(
            seq["presentation_sequence"]
        ).flatten().astype(int)

        folders = loadmat(str(mc_folders_path))
        image_folders = folders["image_folders"]
        names: list[str] = []
        for row in image_folders:
            entry = row[0]
            name_arr = entry["name"] if isinstance(entry, np.void) else entry[0]
            names.append(str(np.asarray(name_arr).flatten()[0]))
        self._cloud_names: list[str] = names

    @staticmethod
    def condition(trial_sequence: int) -> str:
        return CONDITION_MAP[int(trial_sequence)]

    @property
    def n_trials(self) -> int:
        return len(self._presentation)

    @property
    def n_clouds(self) -> int:
        return len(self._cloud_names)

    def cloud_name(self, trial_id: int) -> str | None:
        idx = self._presentation[trial_id - 1]
        if idx < 1 or idx > len(self._cloud_names):
            return None
        return self._cloud_names[idx - 1]

    def trial_profile_id(self, trial_id: int) -> int:
        """Speed-profile identifier (1 or 2) for the given trial.

        The rc2 motion-clouds stimulus presents the same set of clouds
        twice, once under each of two canonical speed profiles (the
        two reproduced velocity trajectories Laura refers to). The
        ``presentation_sequence`` array is assembled such that trial
        IDs in the first half correspond to profile 1, and the second
        half to profile 2. Mirrors the MATLAB convention in
        ``scripts/glm_single_cluster_analysis.m:2255-2293``
        (``sp_midpoint = floor(n_total_seq / 2);
        sp_fold(trial_id > sp_midpoint) = 2``).

        ``trial_id`` is 1-based.
        """
        mid = len(self._presentation) // 2
        return 1 if int(trial_id) <= mid else 2

    def stimulus_params(self, trial_id: int) -> dict[str, Any]:
        name = self.cloud_name(trial_id)
        if name is None:
            return {"sf": np.nan, "orientation": np.nan, "vx": np.nan,
                    "batch_gain": np.nan, "excluded": True, "cloud_name": None}
        params = parse_cloud_name(name, self._stimulus_set)
        params["cloud_name"] = name
        return params
