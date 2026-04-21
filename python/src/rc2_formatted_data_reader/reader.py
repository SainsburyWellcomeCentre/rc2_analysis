"""Lazy reader for MATLAB v7.3 (HDF5) formatted data files."""

from __future__ import annotations

import re
import csv
from pathlib import Path

import h5py
import numpy as np

from rc2_formatted_data_reader import masks
from rc2_formatted_data_reader.trial_conditions import CONDITION_MAP


def _deref_scalar(f: h5py.File, ref) -> float:
    return float(np.asarray(f[ref][()]).flat[0])


def _deref_string(f: h5py.File, ref) -> str:
    raw = np.asarray(f[ref][()]).flatten()
    return "".join(chr(int(c)) for c in raw)


def _deref_array(f: h5py.File, ref) -> np.ndarray:
    return np.asarray(f[ref][()]).squeeze()


def _hdf5_string(ds: h5py.Dataset) -> str:
    raw = np.asarray(ds[()]).flatten()
    return "".join(chr(int(c)) for c in raw)


class FormattedDataReader:
    """Lazy reader for MATLAB formatted data files (.mat v7.3 / HDF5).

    Opens the file once, provides indexed access to clusters, trials,
    and session-level arrays without loading everything into memory.

    Methods take 0-based array indices into the clusters / trials
    datasets (not MATLAB cluster IDs). Use `cluster_ids()` for the
    mapping to MATLAB IDs.
    """

    def __init__(self, mat_path: str | Path):
        self._path = Path(mat_path)
        self._f: h5py.File | None = h5py.File(self._path, "r")
        self._replay_offsets = self._load_replay_offsets()

    def _load_replay_offsets(self) -> dict[int, int]:
        offsets: dict[int, int] = {}
        csv_path = self._path.parent / "csvs" / "trial_matched_offsets" / self._path.with_suffix(".csv").name
        if not csv_path.exists():
            return offsets
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    replay_id = int(row["replay_trial_id"])
                    offset = int(float(row["sample_offset"]))
                    offsets[replay_id] = offset
                except (ValueError, KeyError):
                    continue
        return offsets

    # --- Context manager / lifecycle ---

    def __enter__(self) -> "FormattedDataReader":
        return self

    def __exit__(self, *args) -> None:
        self.close()

    def close(self) -> None:
        if self._f is not None:
            self._f.close()
            self._f = None

    @property
    def _file(self) -> h5py.File:
        if self._f is None:
            raise ValueError("FormattedDataReader is closed")
        return self._f

    # --- Top-level metadata ---

    @property
    def probe_id(self) -> str:
        return _hdf5_string(self._file["probe_id"])

    @property
    def fs(self) -> float:
        return float(self._file["sessions"]["fs"][()].flat[0])

    @property
    def n_clusters(self) -> int:
        return int(self._file["clusters"]["id"].shape[0])

    @property
    def n_trials(self) -> int:
        return int(self._file["sessions"]["trials"]["trial_id"].shape[0])

    @property
    def n_samples(self) -> int:
        return int(self._file["sessions"]["probe_t"].shape[1])

    # --- Cluster-level access ---

    def cluster_ids(self) -> np.ndarray:
        f = self._file
        refs = f["clusters"]["id"][:, 0]
        return np.array([_deref_scalar(f, r) for r in refs], dtype=np.int64)

    def cluster_region(self, cluster_idx: int) -> str:
        f = self._file
        return _deref_string(f, f["clusters"]["region_str"][cluster_idx, 0])

    def cluster_regions(self) -> list[str]:
        f = self._file
        refs = f["clusters"]["region_str"][:, 0]
        return [_deref_string(f, r) for r in refs]

    def visp_cluster_indices(self) -> np.ndarray:
        regions = self.cluster_regions()
        return np.array(
            [i for i, r in enumerate(regions) if re.match(r"VISp", r)],
            dtype=np.int64,
        )

    def visp_cluster_ids(self) -> np.ndarray:
        idx = self.visp_cluster_indices()
        if idx.size == 0:
            return np.array([], dtype=np.int64)
        return self.cluster_ids()[idx]

    def spike_times(self, cluster_idx: int) -> np.ndarray:
        f = self._file
        arr = _deref_array(f, f["clusters"]["spike_times"][cluster_idx, 0])
        return np.atleast_1d(arr).astype(np.float64)

    def cluster_class(self, cluster_idx: int) -> str:
        f = self._file
        return _deref_string(f, f["clusters"]["class"][cluster_idx, 0])

    def cluster_depth(self, cluster_idx: int) -> float:
        f = self._file
        return _deref_scalar(f, f["clusters"]["depth"][cluster_idx, 0])

    # --- Trial-level access ---

    def trial_bounds(self, trial_idx: int) -> tuple[int, int]:
        f = self._file
        start = int(_deref_scalar(f, f["sessions"]["trials"]["start_idx"][trial_idx, 0]))
        end = int(_deref_scalar(f, f["sessions"]["trials"]["end_idx"][trial_idx, 0]))
        trial_id = self.trial_id(trial_idx)
        if trial_id in self._replay_offsets:
            offset = self._replay_offsets[trial_id]
            start += offset
            end += offset
        return start, end

    def trial_id(self, trial_idx: int) -> int:
        f = self._file
        return int(_deref_scalar(f, f["sessions"]["trials"]["trial_id"][trial_idx, 0]))

    def trial_protocol(self, trial_idx: int) -> str:
        f = self._file
        return _deref_string(f, f["sessions"]["trials"]["protocol"][trial_idx, 0])

    def _trial_sequence(self, trial_idx: int) -> int:
        f = self._file
        cfg_ref = f["sessions"]["trials"]["config"][trial_idx, 0]
        cfg = f[cfg_ref]
        return int(np.asarray(cfg["trial_sequence"][()]).flat[0])

    def trial_condition(self, trial_idx: int) -> str:
        return CONDITION_MAP[self._trial_sequence(trial_idx)]

    def _trial_config_scalar(self, trial_idx: int, key: str, default: float = 0.0) -> float:
        f = self._file
        cfg_ref = f["sessions"]["trials"]["config"][trial_idx, 0]
        cfg = f[cfg_ref]
        if key in cfg:
            val = cfg[key][0, 0]
            if type(val).__name__ == 'Reference':
                return _deref_scalar(f, val)
            return float(val)
        return default

    def is_replay(self, trial_idx: int) -> bool:
        """Returns True if the trial is a replay protocol."""
        return self.trial_protocol(trial_idx) in ("StageOnly", "ReplayOnly")

    def trial_analysis_mask(self, trial_idx: int) -> np.ndarray:
        """Returns the boolean analysis window mask replicating AlignedTrial logic."""
        fs = self.fs
        start, end = self.trial_bounds(trial_idx)
        solenoid = self.solenoid(start, end)
        protocol = self.trial_protocol(trial_idx)

        remove_after_solenoid_start = 0.2
        remove_at_end = 0.55
        add_before_solenoid_low = 2.0

        mask = np.zeros(len(solenoid), dtype=bool)

        if protocol in ("Coupled", "EncoderOnly"):
            idx = np.where(solenoid < 2.5)[0]
            if len(idx) > 0:
                idx = idx[:-1] # MATLAB removes last idx
                p1 = np.arange(int(idx[0] - add_before_solenoid_low * fs), idx[0])
                p2 = np.arange(int(idx[0] + remove_after_solenoid_start * fs), int(idx[-1] - remove_at_end * fs) + 1)
                idx_new = np.concatenate([p1, p2])
                idx_new = idx_new[(idx_new >= 0) & (idx_new < len(solenoid))]
                mask[idx_new] = True

        elif protocol in ("CoupledMismatch", "EncoderOnlyMismatch"):
            after_mm_to_remove = 3.0
            idx = np.where(solenoid < 2.5)[0]
            if len(idx) > 0:
                idx = idx[:-1]
                p1 = np.arange(int(idx[0] - add_before_solenoid_low * fs), idx[0])
                p2 = np.arange(int(idx[0] + remove_after_solenoid_start * fs), int(idx[-1] - remove_at_end * fs) + 1)
                idx_new = np.concatenate([p1, p2])

                sess = self._file["sessions"]
                teensy_gain = sess["minidaq_ao0"][0, start:end] if "minidaq_ao0" in sess else sess["teensy_gain"][0, start:end]
                is_high = (teensy_gain > 2.5).astype(int)
                mm_onsets = np.where(np.diff(is_high) == 1)[0]
                if len(mm_onsets) > 0:
                    mm_onset_idx = mm_onsets[0] + 1
                    p1_new = np.arange(idx_new[0] if len(idx_new) > 0 else 0, mm_onset_idx + 1)
                    p2_new = np.arange(int(mm_onset_idx + after_mm_to_remove * fs), (idx_new[-1] + 1) if len(idx_new) > 0 else len(solenoid))
                    idx_new = np.concatenate([p1_new, p2_new])

                idx_new = idx_new[(idx_new >= 0) & (idx_new < len(solenoid))]
                mask[idx_new] = True

        elif protocol in ("StageOnly", "ReplayOnly"):
            sol_high = np.where(solenoid > 2.5)[0]
            start_off = 0
            if len(sol_high) > 0:
                start_off = int(sol_high[0] + remove_after_solenoid_start * fs)

            vel = self.velocity(start, end)
            position_tr = np.cumsum(vel) / fs

            c_start_pos = self._trial_config_scalar(trial_idx, "start_pos")
            c_fwd_limit = self._trial_config_scalar(trial_idx, "forward_limit")
            thresh = 0.98 * (c_start_pos - c_fwd_limit) / 10.0

            pos_high = np.where(position_tr > thresh)[0]
            end_off = len(solenoid) - 1
            if len(pos_high) > 0:
                end_off = pos_high[0]

            idx_new = np.arange(start_off, end_off + 1)
            idx_new = idx_new[(idx_new >= 0) & (idx_new < len(solenoid))]
            mask[idx_new] = True

        rc2_t = np.arange(len(solenoid)) / fs
        baseline = np.zeros(len(solenoid), dtype=bool)
        if protocol in ("Coupled", "EncoderOnly", "CoupledMismatch", "EncoderOnlyMismatch"):
            sol_down = np.where(np.diff((solenoid > 2.5).astype(int)) == -1)[0]
            if len(sol_down) > 0:
                t_ref = rc2_t[sol_down[0]]
                baseline = (rc2_t > (t_ref - 3.0)) & (rc2_t < (t_ref - 1.0))
        elif protocol in ("StageOnly", "ReplayOnly"):
            sol_up = np.where(np.diff((solenoid > 2.5).astype(int)) == 1)[0]
            if len(sol_up) > 0:
                t_ref = rc2_t[sol_up[0]]
                baseline = (rc2_t > (t_ref + 2.0)) & (rc2_t < (t_ref + 4.0))

        return mask | baseline

    # --- Session-level arrays (lazy slicing) ---

    def probe_t(self, start: int, end: int) -> np.ndarray:
        return self._file["sessions"]["probe_t"][0, start:end]

    def velocity(self, start: int, end: int) -> np.ndarray:
        return self._file["sessions"]["filtered_teensy"][0, start:end]

    def solenoid(self, start: int, end: int) -> np.ndarray:
        return self._file["sessions"]["solenoid"][0, start:end]

    def stage(self, start: int, end: int) -> np.ndarray:
        return self._file["sessions"]["stage"][0, start:end]

    def photodiode(self, start: int, end: int) -> np.ndarray:
        return self._file["sessions"]["photodiode"][0, start:end]

    def _camera(
        self, key: str, start: int | None, end: int | None
    ) -> np.ndarray | None:
        sess = self._file["sessions"]
        if key not in sess:
            return None
        cam = sess[key]
        if start is None and end is None:
            return cam[()].squeeze()
        s = 0 if start is None else start
        e = cam.shape[-1] if end is None else end
        return cam[0, s:e]

    def camera0(self, start: int | None = None, end: int | None = None) -> np.ndarray | None:
        return self._camera("camera0", start, end)

    def camera1(self, start: int | None = None, end: int | None = None) -> np.ndarray | None:
        return self._camera("camera1", start, end)

    def camera_t(self) -> np.ndarray | None:
        sess = self._file["sessions"]
        if "camera_t" not in sess:
            return None
        return np.asarray(sess["camera_t"][()]).squeeze()

    # --- Motion / stationary masks (simple thresholding) ---

    def motion_mask(
        self, trial_idx: int, start: int, end: int, threshold: float = 1.0
    ) -> np.ndarray:
        base_mask = masks.motion_mask(self.velocity(start, end), threshold)
        return base_mask & self.trial_analysis_mask(trial_idx)

    def stationary_mask(
        self, trial_idx: int, start: int, end: int, threshold: float = 1.0
    ) -> np.ndarray:
        base_mask = masks.stationary_mask(self.velocity(start, end), threshold)
        return base_mask & self.trial_analysis_mask(trial_idx)
