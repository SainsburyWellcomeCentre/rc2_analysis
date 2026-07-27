#!/usr/bin/env python3
"""
Extract UnitMatch / DeepUnitMatch RawWaveforms for sorted NP2 sessions.

This is the shared input-preparation step for cross-session unit tracking
(Option A). It is run ONLY as part of the (optional, independent) unit-match
workflow -- never during per-session sorting -- so the sort stays lean.

For each session it writes:
    <ks_dir>/RawWaveforms/Unit{id}_RawSpikes.npy   (one per unit, shape (T, 384, 2))
    <ks_dir>/cluster_group.tsv                      (good/mua/noise, for unit selection)

These are exactly the files UnitMatchPy.utils.paths_from_KS / load_good_waveforms
(and hence DeepUnitMatch) look for.

WHY extract from the CatGT AP bin (not the destriped recording):
  UnitMatch/DeepUnitMatch match units on their *spatial* waveform footprint
  across channels. Destriping (highpass_spatial_filter / CAR) deliberately
  removes shared spatial structure, which would corrupt that footprint. The
  CatGT *.ap.bin in this pipeline is temporally band-passed only (the spatial
  step was moved into SpikeInterface for sorting), so it preserves the spatial
  footprint and is the correct source for raw waveforms. Spike sample indices
  from KS4 align 1:1 with the CatGT bin (destriping does not move samples).

Mirrors the official UnitMatch SpikeInterface demo
(UnitMatchPy/Demo Notebooks/UMPy_spike_interface_demo.ipynb):
  read_spikeglx -> phase_shift -> split into 2 CV halves
  -> SortingAnalyzer per half -> templates -> erd.save_avg_waveforms.

Requires (in the spikeinterface conda env):
    spikeinterface[full], and UnitMatchPy on sys.path (see add_unitmatch_to_path).
"""

import os
import re
import sys
import shutil
import argparse
import numpy as np


# Default clone of github.com/EnnyvanBeest/UnitMatch (override with --unitmatch-repo).
DEFAULT_UNITMATCH_REPO = r'C:\Users\Lab\SWC\UnitMatch'

# ---- ms windows for the extracted waveform (matches UnitMatch SI demo) ----
MS_BEFORE = 1.0
MS_AFTER = 2.0
MAX_SPIKES_PER_UNIT = 500   # spikes averaged per CV half
N_AP_CHANNELS = 384         # NP2.0 neural channels (excludes the SY sync channel)


def add_unitmatch_to_path(unitmatch_repo_dir):
    """Put UnitMatchPy (and the DeepUnitMatch package) on sys.path.

    unitmatch_repo_dir is the clone of github.com/EnnyvanBeest/UnitMatch.
    UnitMatchPy lives at <repo>/UnitMatchPy and is imported as `UnitMatchPy`;
    DeepUnitMatch lives at <repo>/UnitMatchPy/DeepUnitMatch.
    """
    umpy_parent = os.path.join(unitmatch_repo_dir, 'UnitMatchPy')
    for p in (umpy_parent, unitmatch_repo_dir):
        if p not in sys.path:
            sys.path.insert(0, p)


def discover_sessions(recordings_root):
    """Find every sorted session under a single folder.

    A session is any directory named imec<N>_ks4 that contains KS4 output
    (spike_times.npy), searched recursively under recordings_root. Returns the
    list sorted by path -- which is chronological when the session/probe folder
    names are date- or sequence-ordered. The discovered order is the order units
    are matched in, so ALWAYS check the printed list before trusting results.

    Point this at, e.g., an animal's output folder containing several
    catgt_<run>_g0/<run>_g0_imec0/imec0_ks4 sessions.
    """
    root = os.path.abspath(recordings_root)
    if not os.path.isdir(root):
        raise NotADirectoryError(f'Not a folder: {root}')
    sessions = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if re.match(r'imec\d+_ks4$', os.path.basename(dirpath)) and \
                'spike_times.npy' in filenames:
            sessions.append(dirpath)
    return sorted(sessions)


def _ensure_cluster_group_tsv(ks_dir):
    """Make sure <ks_dir>/cluster_group.tsv exists for unit selection.

    The sorting pipeline writes phy/cluster_group.tsv (Bombcell labels) and
    cluster_groups.csv. UnitMatch's paths_from_KS / load_good_waveforms expect
    cluster_group.tsv at the KS directory root, so place a copy there.
    """
    dst = os.path.join(ks_dir, 'cluster_group.tsv')
    if os.path.isfile(dst):
        return dst

    phy_tsv = os.path.join(ks_dir, 'phy', 'cluster_group.tsv')
    if os.path.isfile(phy_tsv):
        shutil.copy2(phy_tsv, dst)
        return dst

    # Fall back to deriving it from cluster_groups.csv (cluster_id,group)
    csv_path = os.path.join(ks_dir, 'cluster_groups.csv')
    if os.path.isfile(csv_path):
        import pandas as pd
        df = pd.read_csv(csv_path)
        id_col = 'cluster_id' if 'cluster_id' in df.columns else df.columns[0]
        grp_col = 'group' if 'group' in df.columns else df.columns[-1]
        df[[id_col, grp_col]].rename(
            columns={id_col: 'cluster_id', grp_col: 'group'}
        ).to_csv(dst, sep='\t', index=False)
        return dst

    raise FileNotFoundError(
        f'No cluster labels found for {ks_dir} '
        f'(looked for phy/cluster_group.tsv and cluster_groups.csv). '
        f'Run the sorting pipeline with Bombcell enabled first.')


def extract_session(ks_dir, overwrite=False):
    """Extract and save RawWaveforms for one sorted session.

    Parameters
    ----------
    ks_dir : str
        Path to an imec{prb}_ks4 directory (KS4 output: spike_times.npy,
        spike_clusters.npy, params.py, channel_positions.npy).
        The CatGT *.ap.bin is read from its parent directory.
    overwrite : bool
        If False and a non-empty RawWaveforms/ already exists, skip.
    """
    import spikeinterface.full as si
    import spikeinterface.extractors as se
    import spikeinterface.preprocessing as spre
    import UnitMatchPy.extract_raw_data as erd

    raw_dir = os.path.join(ks_dir, 'RawWaveforms')
    if (not overwrite) and os.path.isdir(raw_dir) and \
            any(f.endswith('_RawSpikes.npy') for f in os.listdir(raw_dir)):
        print(f'  RawWaveforms already present, skipping: {ks_dir}')
        _ensure_cluster_group_tsv(ks_dir)
        return raw_dir

    # --- recording: CatGT AP bin (temporal band-pass only, spatial footprint intact) ---
    catgt_dir = os.path.dirname(ks_dir)
    recording = se.read_spikeglx(catgt_dir, stream_id='imec0.ap', load_sync_channel=False)
    recording = spre.phase_shift(recording)   # correct per-channel ADC sample delay
    n_samples = recording.get_num_samples()

    # --- sorting: KS4 native output in ks_dir ---
    sorting = se.read_kilosort(ks_dir)
    all_unit_ids = np.array(sorting.get_unit_ids(), dtype=int)
    print(f'  {len(all_unit_ids)} units, {recording.get_num_channels()} channels, '
          f'{n_samples / recording.get_sampling_frequency():.1f} s')

    # --- split into 2 cross-validation halves (first/second half of the recording) ---
    mid = n_samples // 2
    rec_halves = [recording.frame_slice(0, mid),
                  recording.frame_slice(mid, n_samples)]
    sort_halves = [sorting.frame_slice(0, mid),
                   sorting.frame_slice(mid, n_samples)]

    # --- templates per half via SortingAnalyzer (sparse=False -> full 384-ch footprint) ---
    t_halves = []
    for h in (0, 1):
        ana = si.create_sorting_analyzer(sort_halves[h], rec_halves[h], sparse=False)
        ana.compute('random_spikes', method='uniform', max_spikes_per_unit=MAX_SPIKES_PER_UNIT)
        ana.compute('waveforms', ms_before=MS_BEFORE, ms_after=MS_AFTER, dtype='float32')
        ana.compute('templates')
        t_halves.append(ana.get_extension('templates').get_data())  # (n_units, n_samples, n_ch)

    # (n_units, spike_width, n_channels, 2)  -- the UnitMatch RawWaveforms layout
    avg_waves = np.stack((t_halves[0], t_halves[1]), axis=-1)

    # Save ALL units (named by their KS id); good-unit selection happens later
    # from cluster_group.tsv inside load_good_waveforms.
    erd.save_avg_waveforms(avg_waves, ks_dir, all_unit_ids,
                           good_units=None, extract_good_units_only=False)

    _ensure_cluster_group_tsv(ks_dir)
    print(f'  Saved RawWaveforms -> {raw_dir}')
    return raw_dir


def ensure_raw_waveforms(sessions, unitmatch_repo_dir=None, overwrite=False):
    """Extract RawWaveforms for every session that does not already have them.

    Parameters
    ----------
    sessions : list of str
        imec{prb}_ks4 directories, one per session.
    unitmatch_repo_dir : str, optional
        Clone of the UnitMatch repo; added to sys.path so UnitMatchPy imports.
    overwrite : bool
        Re-extract even if RawWaveforms/ already exists.
    """
    if unitmatch_repo_dir is not None:
        add_unitmatch_to_path(unitmatch_repo_dir)

    print(f'Extracting RawWaveforms for {len(sessions)} session(s)...')
    for i, ks_dir in enumerate(sessions):
        print(f'[{i + 1}/{len(sessions)}] {ks_dir}')
        extract_session(ks_dir, overwrite=overwrite)
    print('RawWaveforms extraction complete.\n')


# ============================================================
# Command-line use: point at ONE folder containing all sessions
#   python extract_raw_waveforms.py  D:\path\to\recordings_root
# ============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract UnitMatch RawWaveforms for every sorted session '
                    'found under a single folder.')
    parser.add_argument('recordings_root',
                        help='Folder containing the sorted sessions '
                             '(searched recursively for imec*_ks4).')
    parser.add_argument('--unitmatch-repo', default=DEFAULT_UNITMATCH_REPO,
                        help=f'Clone of github.com/EnnyvanBeest/UnitMatch '
                             f'(default: {DEFAULT_UNITMATCH_REPO}).')
    parser.add_argument('--overwrite', action='store_true',
                        help='Re-extract even if RawWaveforms/ already exists.')
    args = parser.parse_args()

    sessions = discover_sessions(args.recordings_root)
    print(f'Discovered {len(sessions)} session(s) under {args.recordings_root}:')
    for i, s in enumerate(sessions):
        print(f'  [{i + 1}] {s}')
    if not sessions:
        raise SystemExit('No imec*_ks4 sessions found.')

    ensure_raw_waveforms(sessions, unitmatch_repo_dir=args.unitmatch_repo,
                         overwrite=args.overwrite)
