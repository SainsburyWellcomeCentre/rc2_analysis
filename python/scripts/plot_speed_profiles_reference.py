"""Reference figure: the two passive replay speed profiles.

Loads one representative trial per (profile 1, profile 2) from probe
``CAA-1123243_rec1`` — the T_Vstatic condition, so the velocity channel
is the motorised stage (``stage``). Velocity-vs-time is aligned at
motion onset (``trial.motion_mask`` first True). Saves a small 2-panel
PDF as a reference for the fr_me_corr session figures (which separate
trials by profile but don't show the profile shape itself).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rc2_formatted_data_reader import StimulusLookup
from rc2_glm.config import GLMConfig
from rc2_glm.io import load_probe_data

DATA_ROOT = Path("~/local_data/motion_clouds").expanduser()
MAT_PATH = DATA_ROOT / "formatted_data_3probe" / "CAA-1123243_rec1.mat"
SEQ_PATH = DATA_ROOT / "motion_cloud_sequence_250414.mat"
FOL_PATH = DATA_ROOT / "image_folders.mat"
OUT_PDF = (
    DATA_ROOT / "figures" / "glm" / "exploration" / "fr_me_corr"
    / "speed_profiles_reference.pdf"
)
XLIM_S = (-4.0, 4.0)


def _trial_t_relative(trial) -> tuple[np.ndarray, np.ndarray]:
    """Velocity vs time, aligned at motion onset (t=0)."""
    motion_idx = np.flatnonzero(trial.motion_mask)
    if motion_idx.size == 0:
        return np.zeros(0), np.zeros(0)
    t_motion_start = float(trial.probe_t[int(motion_idx[0])])
    t_rel = trial.probe_t - t_motion_start
    return t_rel, np.asarray(trial.velocity, dtype=np.float64)


def main() -> None:
    lookup = StimulusLookup(SEQ_PATH, FOL_PATH)
    probe = load_probe_data(
        MAT_PATH, config=GLMConfig(),
        stimulus_lookup=lookup, visp_only=True,
    )
    # Pick the first T_Vstatic trial per profile (stage-only, no replay,
    # so the velocity is the raw motorised-stage trajectory).
    picks: dict[int, int] = {}
    for i, t in enumerate(probe.trials):
        if t.excluded or t.condition != "T_Vstatic":
            continue
        if t.profile_id in (1, 2) and t.profile_id not in picks:
            picks[t.profile_id] = i
        if len(picks) == 2:
            break
    if len(picks) != 2:
        raise SystemExit(f"could not find both profiles, picked: {picks}")

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 2.6), sharey=True)
    for ax, profile in zip(axes, (1, 2)):
        trial = probe.trials[picks[profile]]
        t_rel, vel = _trial_t_relative(trial)
        ax.plot(t_rel, vel, color="black", lw=1.0)
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="black", lw=0.5, ls="--")
        ax.axvline(0.1, color="black", lw=0.7, alpha=0.7)
        ax.set_xlim(*XLIM_S)
        ax.set_xlabel("Time from motion onset (s)")
        ax.set_title(
            f"profile {profile}  (trial_id={trial.trial_id}, "
            f"channel={trial.velocity_channel})",
            fontsize=9,
        )
    axes[0].set_ylabel("Velocity (cm/s)")
    fig.suptitle(
        f"Passive replay speed profiles — {probe.probe_id}, T_Vstatic example trial each",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    plt.close(fig)
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
