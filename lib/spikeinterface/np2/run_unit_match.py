#!/usr/bin/env python3
"""
UnitMatch script for tracking units across chronic recording sessions.

Matches units (neurons) across multiple sessions from the same animal using
waveform similarity. Uses UnitMatchPy with Bombcell output from
spikeGLX_pipeline_np2.py.

Prerequisites:
  - All sessions must already be sorted with spikeGLX_pipeline_np2.py
    (including the Bombcell step, which produces RawWaveforms)
  - UnitMatchPy installed: pip install UnitMatchPy
  - Same conda environment as the sorting pipeline (spikeinterface)

Usage:
  1. Edit the "User input" section below
  2. conda activate spikeinterface
  3. python run_unit_match.py

Output (saved to save_dir):
  MatchTable.csv         : unit pairs with match probability and similarity scores
  UniqueIDConversion.mat : cluster IDs with unique IDs shared across sessions
  + figures and additional match metrics
"""

import os
import numpy as np

# Default clone of github.com/EnnyvanBeest/UnitMatch (override with --unitmatch-repo).
DEFAULT_UNITMATCH_REPO = r'C:\Users\Lab\SWC\UnitMatch'


def check_sessions(sessions):
    """Verify that each session directory has the expected Bombcell output."""
    all_ok = True
    for s in sessions:
        if not os.path.isdir(s):
            print(f'  [ERROR] directory not found: {s}')
            all_ok = False
            continue

        # check for RawWaveforms (bombcell or qMetrics subfolder)
        raw_waveforms_found = any(
            os.path.isdir(os.path.join(s, sub, 'RawWaveforms'))
            for sub in ['bombcell', 'qMetrics', '']
        )
        if not raw_waveforms_found:
            print(f'  [WARNING] no RawWaveforms folder found in: {s}')
            print(f'            Run spikeGLX_pipeline_np2.py with Bombcell enabled first.')
            all_ok = False

        # check for cluster labels
        label_files = ['cluster_bc_unitType.tsv', 'cluster_group.tsv']
        label_found = any(os.path.isfile(os.path.join(s, f)) for f in label_files)
        # also check inside phy/ subfolder
        if not label_found:
            label_found = any(
                os.path.isfile(os.path.join(s, 'phy', f)) for f in label_files
            )
        if not label_found:
            print(f'  [WARNING] no cluster label file found in: {s}')
            all_ok = False

        # check channel_positions.npy
        if not os.path.isfile(os.path.join(s, 'channel_positions.npy')):
            print(f'  [WARNING] channel_positions.npy not found in: {s}')
            all_ok = False

    return all_ok


def run_unit_match(sessions, save_dir, unitmatch_repo=DEFAULT_UNITMATCH_REPO,
                   match_threshold=0.75, good_units_only=True):
    """
    Run UnitMatchPy across a list of sorted session directories.

    Parameters
    ----------
    sessions : list of str
        Paths to imec0_ks4 directories (one per session, chronological order).
    save_dir : str
        Directory where results are saved.
    unitmatch_repo : str
        Clone of EnnyvanBeest/UnitMatch (provides UnitMatchPy on sys.path).
    match_threshold : float
        Probability threshold above which a unit pair is called a match.
    good_units_only : bool
        If True, only match units labelled 'good' by Bombcell.
    """
    from extract_raw_waveforms import add_unitmatch_to_path
    add_unitmatch_to_path(unitmatch_repo)
    try:
        import UnitMatchPy.utils as util
        import UnitMatchPy.overlord as ov
        import UnitMatchPy.bayes_functions as bf
        import UnitMatchPy.save_utils as su
        import UnitMatchPy.assign_unique_id as aid
        import UnitMatchPy.default_params as default_params
    except ImportError:
        raise ImportError(
            'UnitMatchPy is not installed.\n'
            'Install it with: pip install UnitMatchPy'
        )

    print('\n' + '='*60)
    print('UnitMatch — cross-session unit matching')
    print('='*60)
    print(f'Sessions ({len(sessions)}):')
    for i, s in enumerate(sessions):
        print(f'  [{i+1}] {s}')
    print(f'Save dir : {save_dir}')
    print(f'Threshold: {match_threshold}')
    print()

    # --- prepare inputs: extract RawWaveforms for any session missing them ---
    # The sorting pipeline runs Bombcell in label-only mode (no RawWaveforms),
    # so they are generated here, on demand, from the CatGT AP bin.
    from extract_raw_waveforms import ensure_raw_waveforms
    ensure_raw_waveforms(sessions, unitmatch_repo_dir=unitmatch_repo)

    # --- verify inputs ---
    print('[0] Checking session directories...')
    ok = check_sessions(sessions)
    if not ok:
        raise RuntimeError(
            'Some session directories are missing required files. '
            'See warnings above.'
        )
    print('    All sessions OK.\n')

    os.makedirs(save_dir, exist_ok=True)

    # --- initialise parameters ---
    print('[1] Initialising parameters...')
    param = {}
    param['KS_dirs'] = sessions
    param = default_params.get_default_param(param=param)

    # --- discover waveform and label paths ---
    print('[2] Discovering waveform and label paths...')
    wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(
        sessions, param=param
    )
    param = util.get_probe_geometry(channel_pos[0], param)
    print(f'    Waveform paths found: {[str(p) for p in wave_paths]}')

    # --- load waveforms ---
    print('[3] Loading waveforms...')
    waveform, session_id, session_switch, within_session, good_units, param = \
        util.load_good_waveforms(
            wave_paths, unit_label_paths, param,
            good_units_only=good_units_only
        )

    n_units = sum(len(g) for g in good_units)
    print(f'    Loaded {n_units} units across {len(sessions)} sessions.')

    clus_info = {
        'good_units':   good_units,
        'session_switch': session_switch,
        'session_id':   session_id,
        'original_ids': np.concatenate(good_units),
    }

    # --- extract waveform properties ---
    print('[4] Extracting waveform properties...')
    extracted_wave_properties = ov.extract_parameters(
        waveform, channel_pos, clus_info, param
    )

    # --- compute matching scores ---
    print('[5] Computing matching scores...')
    total_score, candidate_pairs, scores_to_include, predictors = \
        ov.extract_metric_scores(
            extracted_wave_properties, session_switch,
            within_session, param, niter=2
        )

    # --- Naive Bayes classification ---
    print('[6] Running Naive Bayes classifier...')
    prior_match = 1 - (param['n_expected_matches'] / param['n_units'] ** 2)
    priors = np.array([prior_match, 1 - prior_match])

    labels = candidate_pairs.astype(int)
    cond   = np.unique(labels)

    parameter_kernels = bf.get_parameter_kernels(
        scores_to_include, labels, cond, param, add_one=1
    )
    probability = bf.apply_naive_bayes(
        parameter_kernels, priors, predictors, param, cond
    )

    output_prob_matrix = probability[:, 1].reshape(
        param['n_units'], param['n_units']
    )

    # --- apply threshold ---
    output_threshold = (output_prob_matrix > match_threshold).astype(float)
    matches = np.argwhere(output_threshold == 1)
    # exclude within-session self-matches
    cross_session_matches = [
        m for m in matches
        if session_id[m[0]] != session_id[m[1]]
    ]
    print(f'    Found {len(cross_session_matches)} cross-session unit matches '
          f'(threshold = {match_threshold}).')

    # --- assign unique IDs ---
    print('[7] Assigning unique IDs across sessions...')
    UIDs = aid.assign_unique_id(output_prob_matrix, param, clus_info)

    # --- save results ---
    print('[8] Saving results...')
    avg_centroid       = extracted_wave_properties['avg_centroid']
    avg_waveform       = extracted_wave_properties['avg_waveform']
    avg_waveform_per_tp = extracted_wave_properties['avg_waveform_per_tp']
    max_site           = extracted_wave_properties['max_site']

    su.save_to_output(
        save_dir,
        scores_to_include, matches, output_prob_matrix,
        avg_centroid, avg_waveform, avg_waveform_per_tp, max_site,
        total_score, output_threshold, clus_info, param,
        UIDs=UIDs, matches_curated=None, save_match_table=True
    )

    print('\n' + '='*60)
    print(f'UnitMatch complete.')
    print(f'  Cross-session matches : {len(cross_session_matches)}')
    print(f'  Results saved to      : {save_dir}')
    print('='*60 + '\n')

    return output_prob_matrix, matches, UIDs


if __name__ == '__main__':
    import argparse
    from extract_raw_waveforms import discover_sessions

    parser = argparse.ArgumentParser(
        description='Classic UnitMatchPy cross-session matching. Point at ONE '
                    'folder containing all the recordings you want to match.')
    parser.add_argument('recordings_root',
                        help='Folder with the sorted sessions '
                             '(searched recursively for imec*_ks4).')
    parser.add_argument('--save-dir', default=None,
                        help='Output folder (default: <recordings_root>/unit_match).')
    parser.add_argument('--unitmatch-repo', default=DEFAULT_UNITMATCH_REPO,
                        help=f'Clone of EnnyvanBeest/UnitMatch (default: {DEFAULT_UNITMATCH_REPO}).')
    parser.add_argument('--threshold', type=float, default=0.75,
                        help='Match-probability threshold (default 0.75).')
    parser.add_argument('--include-mua', action='store_true',
                        help='Also match mua units (default: good units only).')
    args = parser.parse_args()

    sessions = discover_sessions(args.recordings_root)
    print(f'Discovered {len(sessions)} session(s) under {args.recordings_root}:')
    for i, s in enumerate(sessions):
        print(f'  [{i + 1}] {s}')
    if len(sessions) < 2:
        raise SystemExit(f'Need at least 2 sessions to match; found {len(sessions)}.')

    save_dir = args.save_dir or os.path.join(
        os.path.abspath(args.recordings_root), 'unit_match')

    run_unit_match(sessions, save_dir,
                   unitmatch_repo=args.unitmatch_repo,
                   match_threshold=args.threshold,
                   good_units_only=not args.include_mua)
