#!/usr/bin/env python3
"""
Interactive visual review of a completed match_sessions.py run, using
UnitMatchPy's own GUI (GUI.py) -- per-pair waveform overlay, multi-channel
raw-waveform footprint, score histograms, spatial trajectory, and ACG, for
whichever candidate pair you select.

match_sessions.py itself cannot open this GUI: it runs headless (matplotlib
Agg backend, invoked as a MATLAB subprocess with no display attached), and
the GUI is a blocking Tkinter window meant for a human to sit at. Run this
script separately, after match_sessions.py has finished, from a machine with
a display.

Reloads what match_sessions.py already saved (MatchProb.npy, ClusInfo.pickle,
UMparam.pickle, 'UM Scores.npz') instead of re-running the Naive Bayes
classifier. The raw per-unit waveforms and channel geometry are NOT saved by
UnitMatchPy's save_to_output, so those are rebuilt from Session*/ (seconds,
not the ~10-15 min the full matching run takes).

Usage:
    conda activate spikeinterface
    python browse_matches.py <unit_match_save_dir>

Controls (in the GUI window): q/m = confirm match, e/n = confirm non-match,
arrow keys = navigate pairs. Closing the window returns the manually curated
match/non-match lists (printed, not auto-saved -- see UnitMatchPy's own
util.curate_matches + save_to_output(matches_curated=...) if you want to
persist manual curation back into MatchTable.csv).
"""

import os
import sys
import json
import pickle
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from match_sessions import DEFAULT_UNITMATCH_REPO, add_unitmatch_to_path  # noqa: E402


def browse_matches(save_dir, unitmatch_repo=DEFAULT_UNITMATCH_REPO):
    add_unitmatch_to_path(unitmatch_repo)
    import UnitMatchPy.default_params as default_params
    import UnitMatchPy.utils as util
    import UnitMatchPy.overlord as ov
    import UnitMatchPy.GUI as gui

    print('Loading saved matching results...')
    with open(os.path.join(save_dir, 'ClusInfo.pickle'), 'rb') as f:
        clus_info = pickle.load(f)
    with open(os.path.join(save_dir, 'UMparam.pickle'), 'rb') as f:
        param = pickle.load(f)
    output_prob_matrix = np.load(os.path.join(save_dir, 'MatchProb.npy'))
    scores_to_include = dict(np.load(os.path.join(save_dir, 'UM Scores.npz')))

    n_sessions = param['n_sessions']
    wave_paths = [os.path.join(save_dir, f'Session{i}') for i in range(n_sessions)]
    unit_label_paths = [os.path.join(wp, 'bombcell_labels.tsv') for wp in wave_paths]
    channel_pos = [np.load(os.path.join(wp, 'channel_locations.npy')) for wp in wave_paths]

    print('Rebuilding raw waveforms (not saved by UnitMatchPy -- quick, not a re-match)...')
    waveform, session_id, session_switch, within_session, good_units, param = \
        util.load_good_waveforms(wave_paths, unit_label_paths, param, good_units_only=True)
    extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)
    waveform = extracted_wave_properties.get('waveform', waveform)

    # total_score isn't saved by save_to_output (only used to build
    # MatchTable.csv's TotalScore column, ordered by matrix POSITION via
    # np.meshgrid over clus_info['original_ids'] -- see make_match_table).
    # Reconstruct it the same way: MatchTable.csv is one row per (row, col)
    # position in row-major order, so a plain reshape recovers the matrix.
    import pandas as pd
    match_table = pd.read_csv(os.path.join(save_dir, 'MatchTable.csv'))
    n_units = output_prob_matrix.shape[0]
    total_score = match_table['TotalScore'].values.reshape(n_units, n_units)

    print('Launching UnitMatchPy GUI (close the window when done)...')
    gui.process_info_for_GUI(
        output_prob_matrix,
        param['match_threshold'],
        scores_to_include,
        total_score,
        extracted_wave_properties['amplitude'],
        extracted_wave_properties['spatial_decay'],
        extracted_wave_properties['avg_centroid'],
        extracted_wave_properties['avg_waveform'],
        extracted_wave_properties['avg_waveform_per_tp'],
        extracted_wave_properties['good_wave_idxs'],
        extracted_wave_properties['max_site'],
        extracted_wave_properties['max_site_mean'],
        waveform,
        within_session,
        channel_pos,
        clus_info,
        param,
    )
    is_match, not_match, matches_gui = gui.run_GUI()

    print(f'\nManually confirmed matches   : {len(is_match)}')
    print(f'Manually confirmed non-matches: {len(not_match)}')
    print('(not auto-saved -- see UnitMatchPy.utils.curate_matches + '
          'save_to_output(matches_curated=...) to persist these into MatchTable.csv)')
    return is_match, not_match, matches_gui


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Browse a completed match_sessions.py run in the UnitMatchPy GUI.')
    parser.add_argument('save_dir', help='unit_match save directory (from match_sessions.py).')
    parser.add_argument('--unitmatch-repo', default=DEFAULT_UNITMATCH_REPO,
                        help=f'Clone of EnnyvanBeest/UnitMatch (default: {DEFAULT_UNITMATCH_REPO}).')
    args = parser.parse_args()
    browse_matches(args.save_dir, unitmatch_repo=args.unitmatch_repo)
