"""Test fixtures.

Resolves a real formatted .mat file from one of:
  1. RC2_FORMATTED_DATA_PATH env var (full path to a .mat file), or
  2. RC2_FORMATTED_DATA_DIR env var (directory; first *.mat is used), or
  3. ~/local_data/motion_clouds/formatted_data/  (default for laptops)
Tests skip cleanly if no file is found.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

DEFAULT_DATA_DIRS = [
    Path.home() / "local_data" / "motion_clouds" / "formatted_data",
]

DEFAULT_MC_ROOT = Path.home() / "local_data" / "motion_clouds"


def _find_mat_file() -> Path | None:
    explicit = os.environ.get("RC2_FORMATTED_DATA_PATH")
    if explicit:
        p = Path(explicit)
        if p.is_file():
            return p

    dir_override = os.environ.get("RC2_FORMATTED_DATA_DIR")
    candidates: list[Path] = []
    if dir_override:
        candidates.append(Path(dir_override))
    candidates.extend(DEFAULT_DATA_DIRS)

    for d in candidates:
        if not d.is_dir():
            continue
        mats = sorted(d.glob("*.mat"))
        if mats:
            return mats[0]
    return None


@pytest.fixture(scope="session")
def formatted_mat_path() -> Path:
    p = _find_mat_file()
    if p is None:
        pytest.skip(
            "No formatted .mat file found. Set RC2_FORMATTED_DATA_PATH "
            "or RC2_FORMATTED_DATA_DIR, or place files under "
            f"{DEFAULT_DATA_DIRS[0]}"
        )
    return p


@pytest.fixture(scope="session")
def mc_sequence_path() -> Path:
    p = Path(os.environ.get(
        "RC2_MC_SEQUENCE_PATH",
        DEFAULT_MC_ROOT / "motion_cloud_sequence_250414.mat",
    ))
    if not p.is_file():
        pytest.skip(f"Motion-cloud sequence file not found: {p}")
    return p


@pytest.fixture(scope="session")
def mc_folders_path() -> Path:
    p = Path(os.environ.get(
        "RC2_MC_FOLDERS_PATH",
        DEFAULT_MC_ROOT / "image_folders.mat",
    ))
    if not p.is_file():
        pytest.skip(f"Motion-cloud folders file not found: {p}")
    return p
