#!/usr/bin/env python3
"""
Cross-session unit matching for chronic Neuropixels recordings (UnitMatchPy).

Tracks the same neurons across N sorted sessions (N >= 2) by matching their
waveform spatial footprint. Uses the official SpikeInterface/UnitMatchPy
integration (make_UnitMatch_folder_from_sorting_analyzers), which reads
directly from a SortingAnalyzer -- no manual RawWaveforms extraction needed.

Each session's SortingAnalyzer is built on the bandpass_only recording
(saved by the sorting pipeline, see spikeGLX_pipeline_np2.py) rather than the
destriped one used for sorting: destriping deliberately removes the
cross-channel spatial structure that UnitMatch matches on. It is built
sparse (radius_um=150), matching UnitMatch's own internal channel radius
(see UnitMatchPy default_params.py, channel_radius=150) -- there is no
accuracy gain from a denser footprint, only extra computation.

"good" units are taken from each session's selected_clusters.txt (the
pipeline's actual curation -- ClusterSelectionTable's automated decision plus
any manual keep/discard edit), NOT recomputed from Bombcell's own default
thresholds -- see apply_selected_clusters_labels().

Usage:
    conda activate spikeinterface
    python match_sessions.py <imec0_ks4_dir_1> <imec0_ks4_dir_2> [...] [--save-dir DIR]

Prerequisites: each session already sorted (bandpass_only.ap/,
sorting_analyzer/, cluster_groups.csv, selected_clusters.txt present under
its imec{prb}_ks4 dir).

Output (saved to save_dir, default <first session's parent>/unit_match):
    MatchTable.csv       : unit pairs with match probability and similarity scores
    MatchingOverview.png : total-score / probability / final-match matrices
    ClusInfo.pickle, UMparam.pickle, MatchProb.npy, Matches.npy,
    UM Scores.npz, WaveformInfo.npz : intermediate data (see UnitMatchPy's
    save_utils.save_to_output for details)
"""

import os
import io
import re
import json
import shutil
import argparse
import contextlib
import numpy as np
import pandas as pd

# Default clone of github.com/EnnyvanBeest/UnitMatch (override with --unitmatch-repo).
DEFAULT_UNITMATCH_REPO = r'C:\Users\Lab\SWC\UnitMatch'

MATCH_THRESHOLD = 0.75


def add_unitmatch_to_path(unitmatch_repo_dir):
    """Put UnitMatchPy on sys.path. unitmatch_repo_dir is the clone of
    github.com/EnnyvanBeest/UnitMatch; UnitMatchPy lives at <repo>/UnitMatchPy."""
    import sys
    umpy_parent = os.path.join(unitmatch_repo_dir, 'UnitMatchPy')
    if umpy_parent not in sys.path:
        sys.path.insert(0, umpy_parent)


def load_bandpass_analyzer(ks4_dir):
    """Build a sparse SortingAnalyzer (radius_um=150) on the session's
    bandpass_only recording + its existing Kilosort4 sort, with the
    extensions make_UnitMatch_folder_from_sorting_analyzers requires.

    Matching is occasional (only for chronic-tracking animals), unlike
    sorting which always runs -- so this analyzer is built here, on first
    use, rather than by the sorting pipeline itself. It is cached to disk
    (bandpass_sorting_analyzer/, next to sorting_analyzer/) so re-running
    match_sessions on the same sessions (e.g. a different --threshold)
    does not repeat this computation.
    """
    import spikeinterface.full as si
    import spikeinterface.extractors as se

    # progress_bar=True (SI default) emits \r-updated tqdm lines; MATLAB's
    # Command Window (this script's usual caller, via RC2Preprocess.match_
    # sessions) renders each update as a new line and truncates long output
    # well before the run finishes. Disabled here so MATLAB logs stay readable.
    si.set_global_job_kwargs(n_jobs=12, chunk_duration='1s', progress_bar=False)

    data_dir = os.path.dirname(ks4_dir)
    bandpass_dir = os.path.join(data_dir, 'bandpass_only.ap')
    if not os.path.isdir(bandpass_dir):
        raise FileNotFoundError(
            f'No bandpass_only recording found at {bandpass_dir} -- '
            f're-run the sorting pipeline for this session first.'
        )

    analyzer_folder = os.path.join(ks4_dir, 'bandpass_sorting_analyzer')
    if os.path.isdir(analyzer_folder):
        return si.load_sorting_analyzer(analyzer_folder)

    recording = si.load(bandpass_dir)
    sorting = se.read_kilosort(ks4_dir)

    analyzer = si.create_sorting_analyzer(
        sorting=sorting, recording=recording, sparse=True,
        method='radius', radius_um=150.0,
        format='binary_folder', folder=analyzer_folder,
    )
    analyzer.compute('random_spikes', method='uniform', max_spikes_per_unit=500)
    analyzer.compute('waveforms', ms_before=1.5, ms_after=2.5, dtype='float32')
    analyzer.compute('templates', operators=['average'])
    analyzer.compute('spike_amplitudes')  # required for 'amplitude_cutoff'/'amplitude_median'
    analyzer.compute('spike_locations')   # required for the 'drift' quality metric
    analyzer.compute('noise_levels')      # required for the 'snr' quality metric
    analyzer.compute('template_metrics')
    analyzer.compute('principal_components', n_components=5, mode='by_channel_local')
    # required for 'mahalanobis'/'nearest_neighbor'/'d_prime'/'silhouette'

    # Same metric_names as spikeGLX_pipeline_np2.py's SortingAnalyzer step:
    # bombcell_label_units (called inside make_UnitMatch_folder_from_
    # sorting_analyzers) aborts if any of its required metrics are missing.
    analyzer.compute('quality_metrics', metric_names=[
        'num_spikes', 'firing_rate', 'presence_ratio', 'snr',
        'isi_violation', 'rp_violation', 'amplitude_cutoff',
        'amplitude_median', 'drift', 'mahalanobis', 'd_prime',
        'nearest_neighbor', 'silhouette',
    ])
    return analyzer


def check_sessions(ks4_dirs):
    """Verify every session directory has a completed sort."""
    all_ok = True
    for d in ks4_dirs:
        if not os.path.isdir(d):
            print(f'  [ERROR] directory not found: {d}')
            all_ok = False
            continue
        if not os.path.isdir(os.path.join(d, 'sorting_analyzer')):
            print(f'  [ERROR] no sorting_analyzer/ found in: {d}')
            all_ok = False
        if not os.path.isfile(os.path.join(d, 'cluster_groups.csv')):
            print(f'  [ERROR] no cluster_groups.csv found in: {d}')
            all_ok = False
        if not os.path.isfile(os.path.join(d, 'selected_clusters.txt')):
            print(f'  [ERROR] no selected_clusters.txt found in: {d}')
            all_ok = False
        data_dir = os.path.dirname(d)
        if not os.path.isdir(os.path.join(data_dir, 'bandpass_only.ap')):
            print(f'  [ERROR] no bandpass_only.ap/ found next to: {d}')
            all_ok = False
    return all_ok


def apply_selected_clusters_labels(ks4_dir, session_dir):
    """Overwrite session_dir/bombcell_labels.tsv so 'good' means "kept in
    this session's selected_clusters.txt" instead of the default Bombcell
    thresholds make_UnitMatch_folder_from_sorting_analyzers used internally.

    selected_clusters.txt is the pipeline's actual curation result --
    ClusterSelectionTable's automated decision, plus any manual keep/discard
    edit a researcher made in clusters_to_check.csv (see
    RC2Preprocess.create_check_clusters_csv/create_selected_clusters_txt).
    Matching should track that, not silently redo Bombcell with its own
    default thresholds on units that were never actually kept.
    """
    # Same read as UnitMatchPy.utils.load_tsv (plain pd.read_csv, no
    # index_col) -- keep this consistent with how the file is read back.
    labels_path = os.path.join(session_dir, 'bombcell_labels.tsv')
    labels_df = pd.read_csv(labels_path, sep='\t')
    id_col, label_col = labels_df.columns[0], labels_df.columns[-1]

    selected_path = os.path.join(ks4_dir, 'selected_clusters.txt')
    selected_ids = set(np.loadtxt(selected_path, dtype=int).reshape(-1))

    labels_df[label_col] = [
        'good' if unit_id in selected_ids else 'discard'
        for unit_id in labels_df[id_col]
    ]
    labels_df.to_csv(labels_path, sep='\t', index=False)


def save_match_evaluation(save_dir, output_prob_matrix, param, within_session,
                           session_switch, match_threshold, util):
    """Run UnitMatchPy's own evaluate_output() (self-match / false-negative /
    false-positive rates -- see UnitMatchPy demo notebooks) and save its
    printed summary to MatchEvaluation.txt instead of letting it disappear
    into stdout, which is where the upstream function leaves it."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        util.evaluate_output(
            output_prob_matrix, param, within_session, session_switch,
            match_threshold=match_threshold,
        )
    text = buf.getvalue()
    print(text)
    with open(os.path.join(save_dir, 'MatchEvaluation.txt'), 'w') as f:
        f.write(text)


def save_matched_pair_figures(save_dir, cross_session_matches, output_prob_matrix,
                               extracted_wave_properties, clus_info, session_id,
                               ks4_dirs, max_figures=200):
    """Save one waveform-overlay figure per cross-session matched pair, so a
    researcher can visually sanity-check matches without the interactive
    UnitMatchPy GUI (GUI.py's plot_avg_waveforms, reimplemented headless here
    -- the GUI itself is Tkinter-only and can't run from a batch pipeline).

    Each figure shows both units' CV-averaged waveform (at their own peak
    channel) side by side, labelled with their session and TRUE Kilosort
    cluster_id (via clus_info['good_units'], not the internal 0..n_units-1
    index used everywhere else in UnitMatchPy) and the pair's match probability.

    max_figures caps output for large matched sets (e.g. a 10-session chronic
    run) -- pairs are sorted by probability descending, so the cap keeps the
    most confident matches, which are also the most useful to spot-check.
    """
    import matplotlib.pyplot as plt

    if not cross_session_matches:
        print('    No cross-session matches to plot.')
        return

    out_dir = os.path.join(save_dir, 'MatchedPairs')
    os.makedirs(out_dir, exist_ok=True)

    avg_waveform = extracted_wave_properties['avg_waveform']  # (spike_width, n_units, cv)
    good_units = clus_info['good_units']  # per-session arrays of true cluster_id
    session_names = [os.path.basename(os.path.dirname(d)) for d in ks4_dirs]

    session_switch = clus_info['session_switch']  # cumulative unit-count boundaries per session

    def true_cluster_id(unit_idx):
        sid = session_id[unit_idx]
        offset = unit_idx - session_switch[sid]
        return int(good_units[sid][offset].squeeze()), sid

    # output_threshold (and therefore cross_session_matches) is symmetric --
    # (a, b) and (b, a) both appear for the same matched pair. Keep only one
    # ordered copy per unordered pair before plotting/counting.
    unique_pairs = {tuple(sorted((int(a), int(b)))) for a, b in cross_session_matches}

    pairs_sorted = sorted(
        unique_pairs,
        key=lambda m: output_prob_matrix[m[0], m[1]],
        reverse=True,
    )[:max_figures]

    for a, b in pairs_sorted:
        cid_a, sid_a = true_cluster_id(a)
        cid_b, sid_b = true_cluster_id(b)
        prob = output_prob_matrix[a, b]

        fig, ax = plt.subplots(figsize=(5, 4))
        # linewidth/alpha differ so a near-perfect overlap (a good match) still
        # shows both curves instead of the second one fully hiding the first.
        ax.plot(avg_waveform[:, a, :].mean(axis=-1), color='g', linewidth=2.5, alpha=0.6,
                label=f'{session_names[sid_a]}  cluster {cid_a}')
        ax.plot(avg_waveform[:, b, :].mean(axis=-1), color='b', linewidth=1.2, alpha=0.9,
                label=f'{session_names[sid_b]}  cluster {cid_b}')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Amplitude (a.u.)')
        ax.set_title(f'p(match) = {prob:.3f}')
        ax.legend(fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(
            out_dir, f'{session_names[sid_a]}_c{cid_a}__{session_names[sid_b]}_c{cid_b}.png'
        ), dpi=120)
        plt.close(fig)

    print(f'    Saved {len(pairs_sorted)} matched-pair figures -> {out_dir}')


def match_sessions(ks4_dirs, save_dir, unitmatch_repo=DEFAULT_UNITMATCH_REPO,
                    match_threshold=MATCH_THRESHOLD):
    """Run UnitMatchPy across a list of sorted session directories.

    Parameters
    ----------
    ks4_dirs : list of str
        Paths to imec{prb}_ks4 directories (one per session, chronological order).
    save_dir : str
        Directory where results are saved.
    unitmatch_repo : str
        Clone of EnnyvanBeest/UnitMatch (provides UnitMatchPy on sys.path).
    match_threshold : float
        Probability threshold above which a unit pair is called a match.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    add_unitmatch_to_path(unitmatch_repo)
    try:
        from UnitMatchPy.save_utils import make_UnitMatch_folder_from_sorting_analyzers
        import UnitMatchPy.default_params as default_params
        import UnitMatchPy.utils as util
        import UnitMatchPy.overlord as ov
        import UnitMatchPy.bayes_functions as bf
        import UnitMatchPy.save_utils as su
        import UnitMatchPy.assign_unique_id as aid
    except ImportError as e:
        raise ImportError(
            f'UnitMatchPy import failed ({e}).\n'
            'Install/update it with: pip install -U UnitMatchPy'
        )

    print('\n' + '=' * 60)
    print('UnitMatch -- cross-session unit matching')
    print('=' * 60)
    print(f'Sessions ({len(ks4_dirs)}):')
    for i, d in enumerate(ks4_dirs):
        print(f'  [{i + 1}] {d}')
    print(f'Save dir : {save_dir}')
    print(f'Threshold: {match_threshold}')
    print()

    print('[0] Checking session directories...')
    if not check_sessions(ks4_dirs):
        raise RuntimeError('Some session directories are missing required files (see above).')
    print('    All sessions OK.\n')

    os.makedirs(save_dir, exist_ok=True)

    # make_UnitMatch_folder_from_sorting_analyzers (UnitMatchPy) does a plain
    # Session{i}.mkdir() with no exist_ok -- it aborts with FileExistsError if
    # save_dir was already used by a previous run (success or failure). Clear
    # any leftover SessionN/ folders from a prior run of THIS save_dir first
    # so re-running match_sessions on the same sessions doesn't require the
    # caller to manually delete the folder.
    for name in os.listdir(save_dir):
        if re.fullmatch(r'Session\d+', name):
            shutil.rmtree(os.path.join(save_dir, name))

    print('[1] Building bandpass-only SortingAnalyzers + UnitMatch folders...')
    analyzers = [load_bandpass_analyzer(d) for d in ks4_dirs]
    make_UnitMatch_folder_from_sorting_analyzers(analyzers=analyzers, save_dir=save_dir)

    # 'good' now means "kept in this session's selected_clusters.txt" (the
    # pipeline's actual curation, incl. any manual edit) instead of the
    # default Bombcell thresholds used above.
    for i, ks4_dir in enumerate(ks4_dirs):
        apply_selected_clusters_labels(ks4_dir, os.path.join(save_dir, f'Session{i}'))

    print('[2] Loading waveforms and labels...')
    n_sessions = len(ks4_dirs)
    wave_paths = [os.path.join(save_dir, f'Session{i}') for i in range(n_sessions)]
    unit_label_paths = [os.path.join(wp, 'bombcell_labels.tsv') for wp in wave_paths]
    channel_pos = [np.load(os.path.join(wp, 'channel_locations.npy')) for wp in wave_paths]

    with open(os.path.join(wave_paths[0], 'waveform_params.json')) as f:
        waveform_params = json.load(f)

    param = default_params.get_default_param()
    param.update(waveform_params)
    param = util.get_probe_geometry(channel_pos[0], param)

    waveform, session_id, session_switch, within_session, good_units, param = \
        util.load_good_waveforms(wave_paths, unit_label_paths, param, good_units_only=True)

    n_units = sum(len(g) for g in good_units)
    print(f'    Loaded {n_units} units across {n_sessions} sessions.\n')

    clus_info = {
        'good_units':     good_units,
        'session_switch': session_switch,
        'session_id':     session_id,
        'original_ids':   np.concatenate(good_units),
    }

    print('[3] Extracting waveform properties and matching scores...')
    extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)
    param['n_units'] = n_units

    total_score, candidate_pairs, scores_to_include, predictors = ov.extract_metric_scores(
        extracted_wave_properties, session_switch, within_session, param, niter=2
    )

    print('[4] Running Naive Bayes classifier...')
    prior_match = 1 - (param['n_expected_matches'] / param['n_units'] ** 2)
    priors = np.array([prior_match, 1 - prior_match])

    labels = candidate_pairs.astype(int)
    cond = np.unique(labels)
    parameter_kernels = bf.get_parameter_kernels(scores_to_include, labels, cond, param, add_one=1)
    probability = bf.apply_naive_bayes(parameter_kernels, priors, predictors, param, cond)
    output_prob_matrix = probability[:, 1].reshape(n_units, n_units)

    output_threshold = (output_prob_matrix > match_threshold).astype(float)
    matches = np.argwhere(output_threshold == 1)
    # output_threshold is symmetric: (a, b) and (b, a) both appear for the
    # same matched pair, so len(cross_session_matches) alone double-counts
    # every match -- n_unique_matches below is the real "N units tracked
    # across sessions" figure.
    cross_session_matches = [m for m in matches if session_id[m[0]] != session_id[m[1]]]
    n_unique_matches = len({tuple(sorted((int(a), int(b)))) for a, b in cross_session_matches})
    print(f'    Found {n_unique_matches} cross-session unit matches '
          f'(threshold = {match_threshold}).')

    print('[5] Assigning unique IDs across sessions...')
    UIDs = aid.assign_unique_id(output_prob_matrix, param, clus_info)

    # Recorded so browse_matches.py can reopen the UnitMatchPy GUI on this
    # exact run later without re-running the classifier -- save_to_output
    # doesn't persist match_threshold itself, only the resulting matches.
    param['match_threshold'] = match_threshold

    print('[6] Saving results...')
    su.save_to_output(
        save_dir,
        scores_to_include, matches, output_prob_matrix,
        extracted_wave_properties['avg_centroid'],
        extracted_wave_properties['avg_waveform'],
        extracted_wave_properties['avg_waveform_per_tp'],
        extracted_wave_properties['max_site'],
        total_score, output_threshold, clus_info, param,
        UIDs=UIDs, matches_curated=None, save_match_table=True,
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, mat, title in zip(
            axes, (total_score, output_prob_matrix, output_threshold),
            ('Total score', 'Match probability',
             f'Final matches (n={n_unique_matches} unique pairs)')):
        im = ax.imshow(mat, cmap='viridis', aspect='auto')
        ax.set_title(title)
        ax.set_xlabel('Unit')
        ax.set_ylabel('Unit')
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'MatchingOverview.png'), dpi=150)
    plt.close(fig)

    print('[7] Evaluating match quality...')
    save_match_evaluation(
        save_dir, output_prob_matrix, param, within_session, session_switch,
        match_threshold, util,
    )

    print('[8] Saving per-pair figures for visual review...')
    save_matched_pair_figures(
        save_dir, cross_session_matches, output_prob_matrix,
        extracted_wave_properties, clus_info, session_id, ks4_dirs,
    )

    print('\n' + '=' * 60)
    print('UnitMatch complete.')
    print(f'  Cross-session matches : {n_unique_matches} unique pairs')
    print(f'  Results saved to      : {save_dir}')
    print('=' * 60 + '\n')

    return output_prob_matrix, matches, UIDs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Cross-session unit matching across sorted Neuropixels sessions.')
    parser.add_argument('ks4_dirs', nargs='+',
                        help='imec{prb}_ks4 directories, one per session, in chronological order.')
    parser.add_argument('--save-dir', default=None,
                        help='Output folder (default: <parent of first session>/unit_match).')
    parser.add_argument('--unitmatch-repo', default=DEFAULT_UNITMATCH_REPO,
                        help=f'Clone of EnnyvanBeest/UnitMatch (default: {DEFAULT_UNITMATCH_REPO}).')
    parser.add_argument('--threshold', type=float, default=MATCH_THRESHOLD,
                        help=f'Match-probability threshold (default {MATCH_THRESHOLD}).')
    args = parser.parse_args()

    if len(args.ks4_dirs) < 2:
        raise SystemExit(f'Need at least 2 sessions to match; got {len(args.ks4_dirs)}.')

    save_dir = args.save_dir or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(args.ks4_dirs[0]))), 'unit_match')

    match_sessions(args.ks4_dirs, save_dir,
                   unitmatch_repo=args.unitmatch_repo,
                   match_threshold=args.threshold)
