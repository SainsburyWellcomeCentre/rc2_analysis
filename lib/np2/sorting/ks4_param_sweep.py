#!/usr/bin/env python3
"""
Kilosort4 parameter sweep -- TEST ONLY, does not change the pipeline.

Re-sorts ONE already-processed session with a few KS4 clustering/merging
parameters changed one at a time, to check whether KS4 overmerges units, and
compares each result against the session's existing sort (the "baseline").

  python ks4_param_sweep.py                  # run every condition in run_conditions
  python ks4_param_sweep.py ccg020 ccg015    # run only these conditions
  python ks4_param_sweep.py --check          # integrity check only (see below)

Parameters tested (Kilosort 4.1.7 behaviour, checked in its source):
  ccg_threshold (0.25)       final merge step only (template_matching.
                             merging_function): two clusters merge if template
                             similarity > 0.5 AND their CCG is refractory
                             (R12 < ccg_threshold). Lower -> fewer merges.
                             Does NOT affect splits (swarmsplitter hard-codes 0.25).
  cluster_neighbors (10)     k of the kNN graph used for clustering. Lower ->
                             sparser graph -> clusters split more easily.
  cluster_downsampling (20)  every Nth spike is a graph landmark, capped at
                             max_cluster_subset (25000) per clustering center --
                             so it is tested with the cap removed.

The pipeline is NEVER edited: the test values only exist in this script's own
run_sorter() call. spikeGLX_pipeline_np2.py is imported (read only) to reuse
its destriping, KS4 settings and Bombcell thresholds, so every condition
differs from the pipeline ONLY by the parameters listed in CONDITIONS.

Destriped recording: detect_bad_channels picks random chunks (seed=None), so
re-destriping per condition could flag different bad channels each time. The
destriped recording is therefore written ONCE to out_root/_destriped_recording
(~ the size of bandpass_only.ap) and every condition sorts that same file.
KS4 reads it directly (binary compatible), so no per-run recording.dat copy.

Output (out_root/):
  _destriped_recording/            shared KS4 input (delete when the sweep is over)
  <condition>/
    ks4/                           SpikeInterface KS4 folder (sorter_output/ inside)
    sorting_analyzer/              small analyzer (no PCA)
    bombcell/unit_labels.csv       Bombcell labels + failed thresholds
    csv/clusters_to_check.csv      same format as ClusterSelectionTable.curation_table
    csv/overmerge_splits.csv       baseline clusters this condition split (see below)
    csv/new_merges.csv             baseline clusters this condition merged together
    condition.json                 parameters + summary numbers
  sweep_summary.csv                one row per condition (+ the baseline)
  sweep_log.csv                    status/runtime per step
  integrity_before.json / integrity_after.json

Also written next to the baseline (NEW files only -- nothing existing is touched):
  <baseline>/csv/clusters_to_check_<condition>.csv
  <baseline>/csv/overmerge_splits_<condition>.csv

overmerge_splits: a new unit is a "piece" of a baseline cluster if >= 50% of
its spikes coincide (+-0.4 ms) with that baseline cluster and it holds >= 5%
of the baseline cluster's spikes. A baseline cluster with >= 2 pieces was
split by this condition -- a candidate overmerge in the baseline. Cluster ids
in baseline_cluster_id are the baseline's own (open <baseline>/phy in Phy).

Integrity check (--check, and automatically at the end of every run):
SHA-256 of the pipeline, kilosort/parameters.py and the baseline's curation
files is recorded before the first run and compared afterwards; the KS4
defaults and the pipeline's run_sorter() arguments are also checked.
"""

import os
import sys
import json
import re
import hashlib
import argparse
import traceback
from datetime import datetime

import numpy as np
import pandas as pd


# ============================================================
# User input
# ============================================================

session_dir = r'D:\mvelez\mateoData_probe\ilie_pipeline\CAA-1124720\output\preprocessed_CAA-1124720_rec1_rec2_rec3_g0\CAA-1124720_rec1_rec2_rec3_g0_imec0'
baseline_ks4_dir = os.path.join(session_dir, 'imec0_ks4')
out_root = r'D:\mvelez\mateoData_probe\ilie_pipeline\CAA-1124720\ks4_sweep'

# Run order: re-run of the defaults first (run-to-run noise), then
# ccg_threshold (priority), then the clustering-stage parameters.
run_conditions = ['base_rep', 'ccg020', 'ccg015', 'ccg010', 'nn7', 'nn5', 'ds10_nocap', 'ds5_nocap']

# Optional: baseline cluster ids already seen as overmerged in Phy. The
# summary reports which of them each condition split.
known_overmerged = []

# Only these KS4 parameters differ from the pipeline, per condition.
CONDITIONS = {
    'base_rep':   {},
    'ccg020':     {'ccg_threshold': 0.20},
    'ccg015':     {'ccg_threshold': 0.15},
    'ccg010':     {'ccg_threshold': 0.10},
    'nn7':        {'cluster_neighbors': 7},
    'nn5':        {'cluster_neighbors': 5},
    'ds10_nocap': {'cluster_downsampling': 10, 'max_cluster_subset': None},
    'ds5_nocap':  {'cluster_downsampling': 5,  'max_cluster_subset': None},
}

# Split/merge matching (see docstring).
PIECE_MIN_FRAC_OF_NEW = 0.5
PIECE_MIN_FRAC_OF_BASELINE = 0.05
MATCH_DELTA_MS = 0.4

# ============================================================
# End of user input
# ============================================================

HERE = os.path.dirname(os.path.abspath(__file__))
PIPELINE_PY = os.path.join(HERE, 'spikeGLX_pipeline_np2.py')
DESTRIPED_DIR = os.path.join(out_root, '_destriped_recording')
LOG_CSV = os.path.join(out_root, 'sweep_log.csv')
SUMMARY_CSV = os.path.join(out_root, 'sweep_summary.csv')

# KS4 values the pipeline must still be on after the sweep (4.1.7 defaults).
EXPECTED_KS4_DEFAULTS = {'ccg_threshold': 0.25, 'cluster_neighbors': 10,
                         'cluster_downsampling': 20, 'max_cluster_subset': 25000,
                         'acg_threshold': 0.2}
EXPECTED_PIPELINE_RUN_SORTER_KWARGS = {'sorter_name', 'recording', 'folder', 'do_CAR', 'nblocks',
                                       'Th_universal', 'Th_learned', 'remove_existing_folder'}

sys.path.insert(0, HERE)
import spikeGLX_pipeline_np2 as pipe  # noqa: E402  (definitions only; main() is guarded)


def log(condition, step, status, message=''):
    os.makedirs(out_root, exist_ok=True)
    new = not os.path.isfile(LOG_CSV)
    with open(LOG_CSV, 'a') as f:
        if new:
            f.write('timestamp,condition,step,status,message\n')
        msg = str(message).replace('\n', ' ').replace(',', ';')
        f.write(f'{datetime.now():%Y-%m-%d %H:%M:%S},{condition},{step},{status},{msg}\n')
    print(f'[{condition}] {step}: {status} {message}')


# ---- Integrity (the pipeline and baseline must be untouched) ----

def integrity_files():
    import kilosort
    files = {
        'spikeGLX_pipeline_np2.py': PIPELINE_PY,
        'kilosort/parameters.py': os.path.join(os.path.dirname(kilosort.__file__), 'parameters.py'),
    }
    for rel in ['csv/clusters_to_check.csv', 'csv/metrics.csv', 'csv/waveform_metrics.csv',
                'cluster_groups.csv', 'bombcell/unit_labels.csv', 'spike_clusters.npy',
                'sorter_output/spike_clusters.npy']:
        files[f'baseline/{rel}'] = os.path.join(baseline_ks4_dir, *rel.split('/'))
    return files


def sha256(path):
    if not os.path.isfile(path):
        return None
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def integrity_snapshot():
    from kilosort.parameters import DEFAULT_SETTINGS
    with open(PIPELINE_PY) as f:
        src = f.read()
    call = re.search(r"ss\.run_sorter\((.*?)\n\s*\)", src, re.S)
    kwargs = sorted(set(re.findall(r'^\s*(\w+)\s*=', call.group(1), re.M))) if call else None
    return {
        'timestamp': f'{datetime.now():%Y-%m-%d %H:%M:%S}',
        'sha256': {k: sha256(p) for k, p in integrity_files().items()},
        'ks4_defaults': {k: DEFAULT_SETTINGS.get(k) for k in EXPECTED_KS4_DEFAULTS},
        'pipeline_run_sorter_kwargs': kwargs,
        'pipeline_ks_settings': {'nblocks': pipe.ks_nblocks, 'Th_universal': pipe.ks_Th_universal,
                                 'Th_learned': pipe.ks_Th_learned},
    }


def record_integrity_before():
    path = os.path.join(out_root, 'integrity_before.json')
    if not os.path.isfile(path):
        os.makedirs(out_root, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(integrity_snapshot(), f, indent=2)
        print(f'Recorded integrity snapshot -> {path}')


def check_integrity():
    before_path = os.path.join(out_root, 'integrity_before.json')
    after = integrity_snapshot()
    with open(os.path.join(out_root, 'integrity_after.json'), 'w') as f:
        json.dump(after, f, indent=2)

    problems = []
    if os.path.isfile(before_path):
        with open(before_path) as f:
            before = json.load(f)
        for name, h in before['sha256'].items():
            if after['sha256'].get(name) != h:
                problems.append(f'{name} changed since {before["timestamp"]}')
    else:
        problems.append('integrity_before.json missing -- nothing to compare against')
    for k, v in EXPECTED_KS4_DEFAULTS.items():
        if after['ks4_defaults'].get(k) != v:
            problems.append(f'KS4 default {k} = {after["ks4_defaults"].get(k)} (expected {v})')
    if set(after['pipeline_run_sorter_kwargs'] or []) != EXPECTED_PIPELINE_RUN_SORTER_KWARGS:
        problems.append(f'pipeline run_sorter kwargs = {after["pipeline_run_sorter_kwargs"]}')

    print('\n' + '=' * 60)
    if problems:
        print('INTEGRITY CHECK: PROBLEMS FOUND')
        for p in problems:
            print('  - ' + p)
    else:
        print('INTEGRITY CHECK: OK -- pipeline, KS4 defaults and baseline files unchanged')
        print(f'  KS4 defaults: {after["ks4_defaults"]}')
        print(f'  pipeline run_sorter: {after["pipeline_run_sorter_kwargs"]}')
    print('=' * 60)
    log('-', 'integrity_check', 'ok' if not problems else 'PROBLEMS', '; '.join(problems))
    return not problems


# ---- Recording ----

def get_destriped_recording():
    import spikeinterface.full as si
    marker = os.path.join(DESTRIPED_DIR, '_complete.json')
    if os.path.isfile(marker):
        print(f'Loading shared destriped recording: {DESTRIPED_DIR}')
        return si.load(DESTRIPED_DIR)

    bandpass_dir = os.path.join(session_dir, 'bandpass_only.ap')
    print(f'Destriping {bandpass_dir} once (shared by every condition)...')
    t0 = datetime.now()
    rec, bad_ids = pipe.ibl_destripe_by_shank(si.load(bandpass_dir))
    rec = rec.save(folder=DESTRIPED_DIR, overwrite=True)
    elapsed = (datetime.now() - t0).total_seconds()
    with open(marker, 'w') as f:
        json.dump({'bad_channel_ids': [str(c) for c in bad_ids], 'elapsed_s': elapsed}, f, indent=2)
    log('-', 'destripe', 'done', f'{len(bad_ids)} bad channels; {elapsed:.0f}s')
    return rec


# ---- Per condition ----

def run_ks4(name, rec):
    import spikeinterface.sorters as ss
    cond_dir = os.path.join(out_root, name)
    ks_dir = os.path.join(cond_dir, 'ks4')
    marker = os.path.join(cond_dir, '_ks4_done.json')
    if os.path.isfile(marker):
        with open(marker) as f:
            return json.load(f)['elapsed_s']

    os.makedirs(cond_dir, exist_ok=True)
    params = CONDITIONS[name]
    log(name, 'kilosort4', 'start', json.dumps(params))
    t0 = datetime.now()
    ss.run_sorter(
        sorter_name='kilosort4',
        recording=rec,
        folder=ks_dir,
        do_CAR=False,
        nblocks=pipe.ks_nblocks,
        Th_universal=pipe.ks_Th_universal,
        Th_learned=pipe.ks_Th_learned,
        remove_existing_folder=True,
        **params,
    )
    elapsed = (datetime.now() - t0).total_seconds()
    # pc_features.npy (~4 GB) is not used by anything below.
    pcf = os.path.join(ks_dir, 'sorter_output', 'pc_features.npy')
    if os.path.isfile(pcf):
        os.remove(pcf)
    with open(marker, 'w') as f:
        json.dump({'elapsed_s': elapsed, 'params': params}, f, indent=2)
    log(name, 'kilosort4', 'done', f'{elapsed:.0f}s')
    return elapsed


def analyze(name, rec):
    """Small SortingAnalyzer (what Bombcell needs + template similarity), then Bombcell."""
    import spikeinterface.full as si
    import spikeinterface.curation as sc

    cond_dir = os.path.join(out_root, name)
    sorting = si.read_sorter_folder(os.path.join(cond_dir, 'ks4'))
    analyzer = si.create_sorting_analyzer(
        sorting=sorting, recording=rec, format='binary_folder',
        folder=os.path.join(cond_dir, 'sorting_analyzer'), overwrite=True,
    )
    analyzer.compute('random_spikes', method='uniform', max_spikes_per_unit=500)
    analyzer.compute('waveforms', ms_before=1.5, ms_after=2.5, dtype='float32')
    analyzer.compute('templates', operators=['average', 'std'])
    analyzer.compute('unit_locations', method='monopolar_triangulation')
    analyzer.compute('spike_amplitudes')
    analyzer.compute('spike_locations')
    analyzer.compute('noise_levels')
    analyzer.compute('template_similarity')
    analyzer.compute('template_metrics', include_multi_channel_metrics=True)
    analyzer.compute('quality_metrics', metric_names=[
        'num_spikes', 'firing_rate', 'presence_ratio', 'snr', 'isi_violation',
        'rp_violation', 'amplitude_cutoff', 'amplitude_median', 'drift',
    ])

    labels_df = sc.bombcell_label_units(
        analyzer, thresholds=pipe.bombcell_thresholds, split_non_somatic_good_mua=True
    )
    metrics = analyzer.get_metrics_extension_data()
    metrics.index.name = 'unit_id'
    labels_df = labels_df.join(pipe.bombcell_failed_thresholds(metrics, pipe.bombcell_thresholds))
    bc_dir = os.path.join(cond_dir, 'bombcell')
    os.makedirs(bc_dir, exist_ok=True)
    labels_df.to_csv(os.path.join(bc_dir, 'unit_labels.csv'))
    metrics.to_csv(os.path.join(bc_dir, 'all_metrics.csv'))
    return sorting, analyzer, labels_df['bombcell_label'], metrics


def curation_table(unit_ids, labels):
    """Python copy of ClusterSelectionTable.curation_table (MATLAB)."""
    group = [str(labels.get(u, '')) for u in unit_ids]
    is_non_somatic = [int(g.startswith('non_soma_')) for g in group]
    keep = [int(g == 'good') for g in group]  # somatic 'good' only
    return pd.DataFrame({'cluster_id': unit_ids, 'bombcell_group': group,
                         'is_non_somatic': is_non_somatic, 'keep_pipeline': keep, 'keep': keep})


def compare_to_baseline(name, sorting, analyzer, labels, metrics, base):
    from spikeinterface.comparison import compare_two_sorters

    cmp = compare_two_sorters(base['sorting'], sorting, sorting1_name='baseline',
                              sorting2_name=name, delta_time=MATCH_DELTA_MS, verbose=False)
    counts = cmp.match_event_count.astype(float)       # rows: baseline, cols: new
    n_base = cmp.event_counts1.astype(float)
    n_new = cmp.event_counts2.astype(float)
    frac_of_new = counts.div(n_new.replace(0, np.nan), axis=1).fillna(0)
    frac_of_base = counts.div(n_base.replace(0, np.nan), axis=0).fillna(0)

    sim = analyzer.get_extension('template_similarity').get_data()
    sim_idx = {u: i for i, u in enumerate(analyzer.unit_ids)}

    def fmt(values, f='{}'):
        return ';'.join(f.format(v) for v in values)

    # Splits: baseline cluster -> >= 2 new pieces
    pieces = {}
    for j in counts.columns:
        i = frac_of_new[j].idxmax()
        if frac_of_new.at[i, j] >= PIECE_MIN_FRAC_OF_NEW and frac_of_base.at[i, j] >= PIECE_MIN_FRAC_OF_BASELINE:
            pieces.setdefault(i, []).append(j)
    split_rows = []
    for i, js in pieces.items():
        if len(js) < 2:
            continue
        js = sorted(js, key=lambda j: -frac_of_base.at[i, j])
        parent_rp = base['metrics'].at[i, 'rp_contamination'] if i in base['metrics'].index else np.nan
        piece_rp = [metrics.at[j, 'rp_contamination'] for j in js]
        pair_sims = [sim[sim_idx[a], sim_idx[b]] for k, a in enumerate(js) for b in js[k + 1:]]
        split_rows.append({
            'baseline_cluster_id': i,
            'n_spikes': int(n_base[i]),
            'bombcell_label': base['labels'].get(i, ''),
            'rp_contamination': parent_rp,
            'n_pieces': len(js),
            'new_cluster_ids': fmt(js),
            'new_n_spikes': fmt([int(n_new[j]) for j in js]),
            'frac_of_parent': fmt([frac_of_base.at[i, j] for j in js], '{:.2f}'),
            'new_rp_contamination': fmt(piece_rp, '{:.3f}'),
            'new_bombcell_labels': fmt([labels.get(j, '') for j in js]),
            'pieces_cleaner': bool(np.all(np.array(piece_rp) < parent_rp)),
            'max_template_similarity_between_pieces': max(pair_sims),
            'known_overmerged': i in known_overmerged,
        })
    splits = pd.DataFrame(split_rows, columns=[
        'baseline_cluster_id', 'n_spikes', 'bombcell_label', 'rp_contamination', 'n_pieces',
        'new_cluster_ids', 'new_n_spikes', 'frac_of_parent', 'new_rp_contamination',
        'new_bombcell_labels', 'pieces_cleaner', 'max_template_similarity_between_pieces',
        'known_overmerged'])
    splits = splits.sort_values(['pieces_cleaner', 'n_spikes'], ascending=[False, False])

    # Merges: >= 2 baseline clusters -> the same new unit
    absorbed = {}
    for i in counts.index:
        j = frac_of_base.loc[i].idxmax()
        if frac_of_base.at[i, j] >= PIECE_MIN_FRAC_OF_NEW and frac_of_new.at[i, j] >= PIECE_MIN_FRAC_OF_BASELINE:
            absorbed.setdefault(j, []).append(i)
    merges = pd.DataFrame([
        {'new_cluster_id': j, 'baseline_cluster_ids': fmt(sorted(iis)),
         'new_n_spikes': int(n_new[j]), 'new_rp_contamination': metrics.at[j, 'rp_contamination'],
         'new_bombcell_label': labels.get(j, '')}
        for j, iis in absorbed.items() if len(iis) >= 2
    ], columns=['new_cluster_id', 'baseline_cluster_ids', 'new_n_spikes',
                'new_rp_contamination', 'new_bombcell_label'])
    return splits, merges


def label_stats(unit_ids, labels, metrics, ks_label_tsv):
    labels = pd.Series(labels).reindex(unit_ids)
    counts = labels.value_counts()
    rp = metrics['rp_contamination'].reindex(unit_ids)
    in_use = labels.isin(['good', 'mua', 'non_soma_good', 'non_soma_mua'])
    n_ks_good = np.nan
    if os.path.isfile(ks_label_tsv):
        n_ks_good = int((pd.read_csv(ks_label_tsv, sep='\t')['KSLabel'] == 'good').sum())
    return {
        'n_units': len(unit_ids),
        'n_ks_good': n_ks_good,
        'bc_good': int(counts.get('good', 0)),
        'bc_mua': int(counts.get('mua', 0)),
        'bc_noise': int(counts.get('noise', 0)),
        'bc_non_soma_good': int(counts.get('non_soma_good', 0)),
        'bc_non_soma_mua': int(counts.get('non_soma_mua', 0)),
        'pct_rp_at_cap': round(100 * float((rp >= 0.999).mean()), 1),
        'median_rp_good_mua': round(float(rp[in_use].median()), 3),
        'n_units_lt300_spikes': int((metrics['num_spikes'].reindex(unit_ids) < 300).sum()),
    }


def load_baseline():
    import spikeinterface.full as si
    labels = pd.read_csv(os.path.join(baseline_ks4_dir, 'bombcell', 'unit_labels.csv'),
                         index_col=0)['bombcell_label']
    metrics = pd.read_csv(os.path.join(baseline_ks4_dir, 'bombcell', 'all_metrics.csv'), index_col='unit_id')
    sorting = si.read_sorter_folder(baseline_ks4_dir)
    stats = label_stats(list(sorting.unit_ids), labels, metrics,
                        os.path.join(baseline_ks4_dir, 'sorter_output', 'cluster_KSLabel.tsv'))
    return {'sorting': sorting, 'labels': labels, 'metrics': metrics, 'stats': stats}


def process_condition(name, rec, base):
    cond_dir = os.path.join(out_root, name)
    done = os.path.join(cond_dir, 'condition.json')
    if os.path.isfile(done):
        print(f'[{name}] already finished -- skipping (delete {done} to redo the analysis)')
        return

    ks_elapsed = run_ks4(name, rec)

    t0 = datetime.now()
    log(name, 'analysis', 'start')
    sorting, analyzer, labels, metrics = analyze(name, rec)
    unit_ids = list(sorting.unit_ids)

    csv_dir = os.path.join(cond_dir, 'csv')
    os.makedirs(csv_dir, exist_ok=True)
    curation = curation_table(unit_ids, labels)
    curation.to_csv(os.path.join(csv_dir, 'clusters_to_check.csv'), index=False)

    splits, merges = compare_to_baseline(name, sorting, analyzer, labels, metrics, base)
    splits.to_csv(os.path.join(csv_dir, 'overmerge_splits.csv'), index=False)
    merges.to_csv(os.path.join(csv_dir, 'new_merges.csv'), index=False)

    # Copies next to the baseline -- new, condition-suffixed names only.
    base_csv = os.path.join(baseline_ks4_dir, 'csv')
    for src, dst in [('clusters_to_check.csv', f'clusters_to_check_{name}.csv'),
                     ('overmerge_splits.csv', f'overmerge_splits_{name}.csv')]:
        assert dst not in ('clusters_to_check.csv', 'metrics.csv', 'waveform_metrics.csv')
        pd.read_csv(os.path.join(csv_dir, src)).to_csv(os.path.join(base_csv, dst), index=False)

    stats = label_stats(unit_ids, labels, metrics,
                        os.path.join(cond_dir, 'ks4', 'sorter_output', 'cluster_KSLabel.tsv'))
    known_split = splits.loc[splits['known_overmerged'], 'baseline_cluster_id'].tolist()
    summary = {
        'condition': name,
        'params': json.dumps(CONDITIONS[name]),
        **stats,
        'n_keep_pipeline': int(curation['keep_pipeline'].sum()),
        'n_baseline_split': len(splits),
        'n_split_pieces_cleaner': int(splits['pieces_cleaner'].sum()),
        'n_new_merges': len(merges),
        'known_overmerged_split': ';'.join(map(str, known_split)),
        'ks4_runtime_s': round(ks_elapsed),
        'analysis_runtime_s': round((datetime.now() - t0).total_seconds()),
    }
    with open(done, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    log(name, 'analysis', 'done', f'{stats["n_units"]} units; {len(splits)} baseline clusters split')


def write_summary(base):
    rows = [{'condition': 'baseline (existing sort)', 'params': '{}', **base['stats'],
             'n_keep_pipeline': int((base['labels'] == 'good').sum())}]
    for name in CONDITIONS:
        p = os.path.join(out_root, name, 'condition.json')
        if os.path.isfile(p):
            with open(p) as f:
                rows.append(json.load(f))
    pd.DataFrame(rows).to_csv(SUMMARY_CSV, index=False)
    print(f'Summary -> {SUMMARY_CSV}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('conditions', nargs='*', help=f'subset of {list(CONDITIONS)}')
    parser.add_argument('--check', action='store_true', help='integrity check only')
    args = parser.parse_args()

    if args.check:
        sys.exit(0 if check_integrity() else 1)

    conditions = args.conditions or run_conditions
    unknown = [c for c in conditions if c not in CONDITIONS]
    if unknown:
        parser.error(f'unknown condition(s) {unknown}; choose from {list(CONDITIONS)}')

    import matplotlib
    matplotlib.use('Agg')
    import spikeinterface.full as si
    si.set_global_job_kwargs(n_jobs=12, chunk_duration='1s', progress_bar=True)

    record_integrity_before()
    base = load_baseline()
    rec = get_destriped_recording()

    for name in conditions:
        try:
            process_condition(name, rec, base)
        except Exception as e:
            log(name, 'condition', 'FAILED', f'{type(e).__name__}: {e}')
            traceback.print_exc()
        write_summary(base)

    check_integrity()


if __name__ == '__main__':
    main()
