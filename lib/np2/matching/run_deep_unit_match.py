#!/usr/bin/env python3
"""
DeepUnitMatch cross-session unit tracking for Neuropixels 2.0 chronic recordings.

Tracks the same neurons across sessions/days using the DEEP-NEURAL-NETWORK
variant of UnitMatch (van Beest et al., Nature Methods).

Probe scope: NP2.0 only, but ANY shank count (this lab uses 2- and 4-shank).
The pretrained model shipped with the repo was trained on "Npix 2.0 4-shank"
data; the actual probe geometry is read per session from channel_positions.npy
(util.get_probe_geometry), so 2-shank NP2.0 recordings are handled too. The
model generalises across NP2.0 layouts because it matches on the local spatial
waveform footprint rather than a fixed global shank map -- but if you match
2-shank data, sanity-check the results (the model was not trained on it).

This is an OPTIONAL, INDEPENDENT step. It is never run by the sorting pipeline;
you launch it by hand after all sessions of an animal are sorted, by pointing
it at ONE folder that contains every recording you want to match:

    conda activate spikeinterface
    python run_deep_unit_match.py  D:\path\to\recordings_root

All imec*_ks4 sessions under that folder are discovered automatically (sorted by
path = chronological when folders are date/sequence-named). No path editing.

    Options:
      --save-dir DIR       output folder (default <root>/unit_match_deep)
      --threshold 0.5      match-probability threshold
      --dist-thresh 50     max drift-corrected distance (um)
      --include-mua        match mua units too (default: good only)
      --re-extract         force-rebuild RawWaveforms
      --unitmatch-repo DIR clone of EnnyvanBeest/UnitMatch (has default)

Pipeline (per the UnitMatch repo's PaperAnalyses/run_deepunitmatch_batch.py):
    0. ensure RawWaveforms/ exist for every session   (extract_raw_waveforms.py)
    1. paths_from_KS + load_good_waveforms             (UnitMatchPy I/O)
    2. DeepUnitMatch model -> snippets -> NN inference -> similarity matrix
    3. drift correction + Naive-Bayes refinement       (shared with UMPy)
    4. directional threshold -> matches -> unique IDs
    5. save MatchTable.csv / UniqueIDConversion + overview figure

Requires (in the spikeinterface conda env):
    pip install UnitMatchPy        (I/O + Naive Bayes + assign IDs)
    torch                          (already installed for KS4; GPU recommended)
    a clone of github.com/EnnyvanBeest/UnitMatch  (DeepUnitMatch package + model)

Output (saved to the save dir, default <recordings_root>\\unit_match_deep):
    MatchTable.csv          : unit pairs with match probability
    UniqueIDConversion.*    : cluster IDs with cross-session unique IDs
    MatchingOverview.png    : similarity / probability / final-match matrices
"""

import os
import sys
import numpy as np

# Default clone of github.com/EnnyvanBeest/UnitMatch (override with --unitmatch-repo).
# Contains UnitMatchPy + DeepUnitMatch + the pretrained NP2.0-4shank model.
DEFAULT_UNITMATCH_REPO = r'C:\Users\Lab\SWC\UnitMatch'


def run_deep_unit_match(sessions, save_dir, unitmatch_repo=DEFAULT_UNITMATCH_REPO,
                        good_units_only=True, thresh=0.5, dist_thresh=50,
                        re_extract=False):
    """Run DeepUnitMatch across a list of sorted imec0_ks4 directories."""

    # local prep step (Option A): put UnitMatch on the path and make sure every
    # session has RawWaveforms/ (extracted from its CatGT bin + KS4 sorting)
    from extract_raw_waveforms import ensure_raw_waveforms, add_unitmatch_to_path
    add_unitmatch_to_path(unitmatch_repo)
    ensure_raw_waveforms(sessions, unitmatch_repo_dir=unitmatch_repo,
                         overwrite=re_extract)

    import copy
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    import UnitMatchPy.default_params as default_params
    import UnitMatchPy.utils as util
    import UnitMatchPy.overlord as ov
    import UnitMatchPy.bayes_functions as bf
    import UnitMatchPy.assign_unique_id as aid
    import UnitMatchPy.save_utils as su
    import UnitMatchPy.metric_functions as mf
    from DeepUnitMatch.utils import param_fun
    from DeepUnitMatch.testing import test
    from DeepUnitMatch.utils import helpers

    device = 'cuda' if test.torch.cuda.is_available() else 'cpu'
    print(f'\n{"=" * 60}\nDeepUnitMatch -- cross-session unit tracking')
    print(f'Device   : {device}')
    print(f'Sessions : {len(sessions)}')
    for i, s in enumerate(sessions):
        print(f'  [{i + 1}] {s}')
    print(f'Save dir : {save_dir}\n{"=" * 60}\n')

    os.makedirs(save_dir, exist_ok=True)
    tmp_path = os.path.join(save_dir, 'tmp_waveforms')
    os.makedirs(tmp_path, exist_ok=True)

    # --- [1] load waveforms + labels (UnitMatchPy I/O) ---
    print('[1] Loading waveforms and labels...')
    wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(sessions)
    param = {'KS_dirs': sessions}
    param = default_params.get_default_param(param=param)
    param = util.get_probe_geometry(channel_pos[0], param)

    waveform, session_id, session_switch, within_session, good_units, param = \
        util.load_good_waveforms(wave_paths, unit_label_paths, param,
                                 good_units_only=good_units_only)
    param['good_units'] = good_units
    n_total = waveform.shape[0]
    print(f'    {n_total} units across {param["n_sessions"]} session(s)\n')

    clus_info = {
        'good_units':     param['good_units'],
        'session_switch': session_switch,
        'session_id':     session_id,
        'original_ids':   np.concatenate(param['good_units']),
    }

    # --- [2] DeepUnitMatch: model -> snippets -> NN similarity ---
    print('[2] DeepUnitMatch inference...')
    model = test.load_trained_model(device=device)

    unit_ids = np.concatenate(param['good_units']).squeeze()
    _, _, kept_idx = param_fun.get_snippets(
        waveform, channel_pos, session_id,
        save_path=tmp_path, unit_ids=unit_ids, param=param)

    # re-sync arrays if get_snippets rejected any units
    if len(kept_idx) < len(waveform):
        waveform, session_id, session_switch, _, good_units, param = \
            util.filter_units_by_index(
                waveform, session_id, session_switch, good_units, kept_idx, param)
        param['good_units'] = good_units
        clus_info = {
            'good_units':     param['good_units'],
            'session_switch': session_switch,
            'session_id':     session_id,
            'original_ids':   np.concatenate(param['good_units']),
        }

    data_dir = os.path.join(tmp_path, 'processed_waveforms')
    sim_matrix = test.inference(model, data_dir)
    print(f'    similarity matrix: {sim_matrix.shape}\n')

    # --- [3] drift correction + Naive-Bayes refinement (same as UMPy) ---
    print('[3] Drift correction + Naive-Bayes refinement...')
    extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)
    sessions_idx = np.unique(session_id)

    # pre-pass: collect NN match labels across every session pair for shared drift correction
    labels_full = np.eye(n_total)
    pair_matches_cache = {}
    for r1 in sessions_idx:
        for r2 in sessions_idx:
            if r1 >= r2:
                continue
            mask = np.isin(session_id, [r1, r2])
            sim_mat = sim_matrix[mask][:, mask]
            indices = np.where(mask)[0]
            n = int(np.sum(mask))

            df = helpers.create_dataframe(
                [param['good_units'][r1], param['good_units'][r2]],
                sim_mat, session_list=[r1, r2])
            matches = test.get_matches(df, sim_mat, session_id[indices], data_dir,
                                       dist_thresh=dist_thresh)
            pair_matches_cache[(r1, r2)] = matches

            subsessionid = np.array(
                [r1] * len(param['good_units'][r1]) +
                [r2] * len(param['good_units'][r2]))
            labels_pair = np.eye(n)
            for (recses1, recses2), group in matches.groupby(by=['RecSes1', 'RecSes2']):
                asmatrix = group['match'].values.reshape(
                    len(param['good_units'][recses1]),
                    len(param['good_units'][recses2])).astype(int)
                labels_pair[np.ix_(subsessionid == recses1, subsessionid == recses2)] = asmatrix
            labels_full[np.ix_(indices, indices)] = labels_pair

    avg_centroid        = extracted_wave_properties['avg_centroid'].copy()
    avg_waveform_per_tp = extracted_wave_properties['avg_waveform_per_tp'].copy()
    _, avg_centroid, avg_waveform_per_tp = mf.drift_n_sessions(
        labels_full.astype(bool), session_switch, avg_centroid, avg_waveform_per_tp,
        sim_matrix, param)

    probs           = np.zeros(sim_matrix.shape)
    distance_matrix = np.zeros(sim_matrix.shape)
    for r1 in sessions_idx:
        for r2 in sessions_idx:
            if r1 >= r2:
                continue
            mask    = np.isin(session_id, [r1, r2])
            sim_mat = sim_matrix[mask][:, mask]
            n       = int(np.sum(mask))
            matches = pair_matches_cache[(r1, r2)]

            labels       = np.eye(sim_mat.shape[0])
            subsessionid = np.array(
                [r1] * len(param['good_units'][r1]) +
                [r2] * len(param['good_units'][r2]))
            for (recses1, recses2), group in matches.groupby(by=['RecSes1', 'RecSes2']):
                asmatrix = group['match'].values.reshape(
                    len(param['good_units'][recses1]),
                    len(param['good_units'][recses2])).astype(int)
                labels[np.ix_(subsessionid == recses1, subsessionid == recses2)] = asmatrix

            avg_waveform_per_tp_pair = avg_waveform_per_tp[:, mask, :, :]
            avg_waveform_per_tp_flip = mf.flip_dim(avg_waveform_per_tp_pair, param, n)
            euclid_dist              = mf.get_Euclidean_dist(avg_waveform_per_tp_flip, param, n)
            centroid_dist, _         = mf.centroid_metrics(euclid_dist, param)

            scores_to_incl    = {'similarity': sim_mat, 'distance': centroid_dist}
            n_units           = int(np.sqrt(len(matches)))
            priors            = np.array([1 - 2 / n_units, 2 / n_units])
            parameter_kernels = bf.get_parameter_kernels(
                scores_to_incl, labels, np.unique(labels), param)
            predictors  = np.stack(list(scores_to_incl.values()), axis=2)
            probability = bf.apply_naive_bayes(
                parameter_kernels, priors, predictors, param, np.unique(labels))
            prob_matrix = probability[:, 1].reshape(n_units, n_units)

            probs[np.ix_(mask, mask)]           = prob_matrix
            distance_matrix[np.ix_(mask, mask)] = centroid_dist

    # --- [4] final matches + unique IDs ---
    final_matches = test.directional_filter(probs, session_id, thresh)
    n_matches = int(np.sum(final_matches)) // 2
    print(f'    {n_matches} cross-session matches (threshold={thresh})\n')

    UIDs = aid.assign_unique_id(probs, param, clus_info)

    # --- [5] save ---
    print('[5] Saving results...')
    su.save_to_output(
        save_dir,
        {'distance': distance_matrix},
        np.argwhere(final_matches),
        probs,
        extracted_wave_properties['avg_centroid'],
        extracted_wave_properties['avg_waveform'],
        extracted_wave_properties['avg_waveform_per_tp'],
        extracted_wave_properties['max_site'],
        distance_matrix,
        final_matches,
        clus_info,
        param,
        UIDs=UIDs,
        matches_curated=None,
        save_match_table=True,
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, mat, title in zip(
            axes, (sim_matrix, probs, final_matches),
            ('Similarity', 'Match probability', f'Final matches (n={n_matches})')):
        im = ax.imshow(mat, cmap='viridis', aspect='auto')
        ax.set_title(title)
        ax.set_xlabel('Unit')
        ax.set_ylabel('Unit')
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'MatchingOverview.png'), dpi=150)
    plt.close(fig)

    print(f'\n{"=" * 60}\nDeepUnitMatch complete.')
    print(f'  Cross-session matches : {n_matches}')
    print(f'  Results saved to      : {save_dir}\n{"=" * 60}\n')

    return probs, final_matches, UIDs


if __name__ == '__main__':
    import argparse
    from extract_raw_waveforms import discover_sessions

    parser = argparse.ArgumentParser(
        description='DeepUnitMatch cross-session unit tracking. Point at ONE '
                    'folder containing all the recordings you want to match.')
    parser.add_argument('recordings_root',
                        help='Folder with the sorted sessions '
                             '(searched recursively for imec*_ks4).')
    parser.add_argument('--save-dir', default=None,
                        help='Output folder (default: <recordings_root>/unit_match_deep).')
    parser.add_argument('--unitmatch-repo', default=DEFAULT_UNITMATCH_REPO,
                        help=f'Clone of EnnyvanBeest/UnitMatch (default: {DEFAULT_UNITMATCH_REPO}).')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Match-probability threshold (default 0.5).')
    parser.add_argument('--dist-thresh', type=float, default=50,
                        help='Max drift-corrected centroid distance in um (default 50).')
    parser.add_argument('--include-mua', action='store_true',
                        help='Also match mua units (default: good units only).')
    parser.add_argument('--re-extract', action='store_true',
                        help='Force re-extraction of RawWaveforms.')
    args = parser.parse_args()

    sessions = discover_sessions(args.recordings_root)
    print(f'Discovered {len(sessions)} session(s) under {args.recordings_root}:')
    for i, s in enumerate(sessions):
        print(f'  [{i + 1}] {s}')
    if len(sessions) < 2:
        raise SystemExit('Need at least 2 sessions to match; found '
                         f'{len(sessions)}.')

    save_dir = args.save_dir or os.path.join(
        os.path.abspath(args.recordings_root), 'unit_match_deep')

    run_deep_unit_match(sessions, save_dir,
                        unitmatch_repo=args.unitmatch_repo,
                        good_units_only=not args.include_mua,
                        thresh=args.threshold, dist_thresh=args.dist_thresh,
                        re_extract=args.re_extract)
