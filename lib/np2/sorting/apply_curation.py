#!/usr/bin/env python3
"""
Apply a SpikeInterface-GUI curation (merges, splits, removals) to a pipeline
session and re-run Bombcell.

  conda activate spikeinterface
  python apply_curation.py "<...>\\imec0_ks4"

Run after curating with sigui_curate.py and pressing "Save in analyzer", which
writes sorting_analyzer/spikeinterface_gui/curation_data.json.

The curation is applied with SpikeInterface's apply_curation, in hard mode (the
curated units' extensions are recomputed from the recording), with
censor_ms = CENSOR_MS: within each merged unit, a spike less than CENSOR_MS
after the previous one is discarded and the first is kept. Kilosort 4 does
the same after its own merges (remove_duplicates, duplicate_spike_ms = 0.25):
KS sometimes fits a second template to the residual of a spike, and the two
detections can end up in different clusters.

Outputs:
  sorting_analyzer/            the curated analyzer; the previous one is kept as
                               sorting_analyzer_before_<timestamp>/ (it still holds
                               spikeinterface_gui/curation_data.json)
  bombcell/                    unit_labels.csv, all_metrics.csv and figures, re-run
                               with this pipeline's bombcell_thresholds;
                               curation_units.csv lists every new unit (merge or
                               split), its parent units, spikes censored, and the
                               labels set in the GUI; applied_curation.json is the
                               curation that was applied
  cluster_groups.csv, csv/metrics.csv, csv/waveform_metrics.csv
                               rewritten (save_rc2_compatible_files)
  spike_times.npy, spike_clusters.npy, spike_templates.npy, amplitudes.npy,
  spike_positions.npy, spike_detection_templates.npy
                               (root, read by RC2Format) rebuilt from the curated
                               analyzer: every root spike is matched to its curated
                               spike on (sample, parent unit); spikes of removed
                               units and censored spikes are dropped. Checked spike
                               by spike against the curated analyzer before
                               anything is written.

Everything overwritten is first copied to curation_backup_<timestamp>/.
Removed units disappear from all outputs (and so from clusters_to_check.csv);
to bring one back, restore the session from the backup and re-curate.

The session can be curated again afterwards (sigui_curate.py on the curated
analyzer, then this script): each run gets its own timestamped backup.

Then, in MATLAB:
  ctl.create_check_clusters_csv(probe_id)   % new csv/clusters_to_check.csv

phy/ is not updated and is out of date after a curation.

Do NOT re-run spikeGLX_pipeline_np2.py on a session after this (any
start_step): copy_ks4_outputs_to_parent() overwrites the root .npy files with
the uncurated Kilosort output.

Bombcell thresholds come from spikeGLX_pipeline_np2.py as it is NOW -- if they
were changed since the session was sorted, the curated session is labelled with
the new ones.
"""

import os
import gc
import sys
import json
import shutil
import argparse
from datetime import datetime

import numpy as np
import pandas as pd

import spikeGLX_pipeline_np2 as pl   # read only: thresholds + helper functions


# Kilosort 4's duplicate_spike_ms (the value used by this pipeline's sorting)
CENSOR_MS = 0.25

# root .npy files with one entry per spike (kept_spikes.npy is Kilosort's own
# mask over its pre-duplicate-removal spikes, a different length -- untouched)
PER_SPIKE_FILES = ['spike_times.npy', 'spike_clusters.npy', 'spike_templates.npy',
                   'amplitudes.npy', 'spike_positions.npy', 'spike_detection_templates.npy']


def unit_time_pairs(unit_ids, times):
    """(unit, sample) pairs sorted, for comparing two spike sets exactly."""
    unit_ids = np.asarray(unit_ids, dtype=np.int64)
    times = np.asarray(times, dtype=np.int64)
    order = np.lexsort((times, unit_ids))
    return unit_ids[order], times[order]


def analyzer_spikes(analyzer):
    sv = analyzer.sorting.to_spike_vector()
    return (np.asarray(analyzer.unit_ids)[sv['unit_index']].astype(np.int64),
            sv['sample_index'].astype(np.int64))


def same_spikes(units_a, times_a, units_b, times_b):
    ua, ta = unit_time_pairs(units_a, times_a)
    ub, tb = unit_time_pairs(units_b, times_b)
    return np.array_equal(ua, ub) and np.array_equal(ta, tb)


def find_parents(curation, orig_unit_ids, root_times, root_units, cur_units, cur_times):
    """
    {curated unit: (origin, [parent unit ids])}. Parents come from the curation
    itself (merge groups, split units); a new unit id is matched to its
    merge/split entry by spike overlap.
    """
    entries = [('merge', [int(u) for u in m.unit_ids]) for m in curation.merges]
    entries += [('split', [int(s.unit_id)]) for s in curation.splits]
    in_entries = {u for _, units in entries for u in units}
    removed = {int(u) for u in curation.removed}
    entry_times = [np.unique(root_times[np.isin(root_units, units)]) for _, units in entries]

    parents = {}
    for u in np.unique(cur_units):
        u = int(u)
        if u in orig_unit_ids and u not in in_entries and u not in removed:
            parents[u] = ('unchanged', [u])
            continue
        times_u = cur_times[cur_units == u]
        overlap = [np.isin(times_u, t).mean() for t in entry_times]
        best = int(np.argmax(overlap))
        if overlap[best] < 0.99:
            sys.exit(f'ERROR: curated unit {u} matches no merge/split of the curation')
        parents[u] = entries[best]
    return parents


def assign_root_spikes(parents, orig_index, root_times, root_units, cur_units, cur_times):
    """
    New cluster id for every root spike (-1 = dropped: removed unit or censored).
    A root spike goes to the curated unit whose parents include its original unit
    and which has a spike at the same sample. Where two parents of a merged unit
    have a spike at the same sample, the one with the lower unit index is kept,
    as in SpikeInterface's spike vector (sorted by sample, then unit index).
    """
    new_units = np.full(len(root_times), -1, dtype=np.int64)
    for u, (_, units) in parents.items():
        idx = np.flatnonzero(np.isin(root_units, units))
        idx = idx[np.isin(root_times[idx], cur_times[cur_units == u])]
        idx = idx[np.lexsort((orig_index[idx], root_times[idx]))]
        first = np.ones(len(idx), bool)
        first[1:] = root_times[idx][1:] != root_times[idx][:-1]
        idx = idx[first]
        if np.any(new_units[idx] != -1):
            sys.exit(f'ERROR: root spikes assigned to more than one curated unit ({u})')
        new_units[idx] = u
    return new_units


def backup(ks4_dir, backup_dir):
    """Copy every file this script overwrites; written to a temp folder first so a
    crash never leaves a partial backup."""
    files = PER_SPIKE_FILES + ['cluster_groups.csv']
    for sub in ('csv', 'bombcell'):
        files += [os.path.join(sub, f) for f in os.listdir(os.path.join(ks4_dir, sub))
                  if os.path.isfile(os.path.join(ks4_dir, sub, f))]
    tmp_dir = backup_dir + '_incomplete'
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    for f in files:
        src, dst = os.path.join(ks4_dir, f), os.path.join(tmp_dir, f)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        assert os.path.getsize(dst) == os.path.getsize(src), f'backup of {f} is incomplete'
    os.rename(tmp_dir, backup_dir)
    print(f'  {len(files)} files backed up to {backup_dir}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('ks4_dir', help='session imec0_ks4 folder (contains sorting_analyzer/)')
    ks4_dir = parser.parse_args().ks4_dir

    # Agg before anything imports pyplot -- see main() in spikeGLX_pipeline_np2.py
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import spikeinterface.full as si
    import spikeinterface.curation as sc

    si.set_global_job_kwargs(n_jobs=12, chunk_duration='1s', progress_bar=True)

    analyzer_dir = os.path.join(ks4_dir, 'sorting_analyzer')
    curated_dir = os.path.join(ks4_dir, 'sorting_analyzer_curated')
    curation_json = os.path.join(analyzer_dir, 'spikeinterface_gui', 'curation_data.json')
    applied_json = os.path.join(curated_dir, 'applied_curation.json')
    bc_dir = os.path.join(ks4_dir, 'bombcell')
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # ---- the curation ----
    if not os.path.isfile(curation_json):
        sys.exit(f'No curation found ({curation_json}) -- save it in the GUI with "Save in analyzer"')
    curation = sc.load_curation(curation_json)
    if not (curation.merges or curation.splits or curation.removed):
        sys.exit('The curation has no merges, splits or removals -- nothing to apply')
    print(f'Curation: {len(curation.merges)} merges, {len(curation.splits)} splits, '
          f'{len(curation.removed)} removed units')
    for m in curation.merges:
        print(f'  merge  {[int(u) for u in m.unit_ids]}')
    for s in curation.splits:
        print(f'  split  {int(s.unit_id)}')
    if curation.removed:
        print(f'  removed {[int(u) for u in curation.removed]}')

    # ---- checks: root files hold exactly the analyzer's spikes ----
    analyzer = si.load_sorting_analyzer(analyzer_dir)
    assert analyzer.has_recording(), 'recording of the SortingAnalyzer not found'
    root = {f: np.load(os.path.join(ks4_dir, f)) for f in PER_SPIKE_FILES}
    n_spikes = len(root['spike_times.npy'])
    assert all(len(a) == n_spikes for a in root.values()), 'root per-spike files differ in length'
    root_times = root['spike_times.npy'].ravel().astype(np.int64)
    root_units = root['spike_clusters.npy'].ravel().astype(np.int64)
    if not same_spikes(root_units, root_times, *analyzer_spikes(analyzer)):
        sys.exit('ERROR: root spike files do not match sorting_analyzer/ -- not a pipeline session in a '
                 'consistent state (was the pipeline re-run after a curation?)')
    orig_unit_ids = [int(u) for u in analyzer.unit_ids]
    # position of each root spike's unit in analyzer.unit_ids (= SpikeInterface unit index)
    lut = np.full(max(orig_unit_ids) + 1, -1, dtype=np.int64)
    lut[orig_unit_ids] = np.arange(len(orig_unit_ids))
    orig_index = lut[root_units]

    # ---- apply the curation to the SortingAnalyzer ----
    with open(curation_json) as f:
        curation_text = f.read()
    if os.path.exists(curated_dir):
        with open(applied_json) as f:
            if f.read() != curation_text:
                sys.exit(f'ERROR: {curated_dir} was made from a different curation -- delete it and re-run')
        print(f'\nReusing curated SortingAnalyzer from an earlier run: {curated_dir}')
        curated = si.load_sorting_analyzer(curated_dir)
    else:
        print(f'\nApplying curation (apply_curation, hard mode, censor_ms={CENSOR_MS})...')
        curated = sc.apply_curation(analyzer, curation, censor_ms=CENSOR_MS,
                                    merging_mode='hard', verbose=True)
        # apply_curation returns an in-memory analyzer -- save it next to the original
        curated = curated.save_as(format='binary_folder', folder=curated_dir)
        with open(applied_json, 'w') as f:
            f.write(curation_text)

    # Bombcell and save_rc2_compatible_files need both metric extensions --
    # recompute any the curation dropped, with the same arguments as the pipeline
    if not curated.has_extension('template_metrics'):
        print('    template_metrics missing after curation -- recomputing...')
        curated.compute('template_metrics', include_multi_channel_metrics=True)
    if not curated.has_extension('quality_metrics'):
        print('    quality_metrics missing after curation -- recomputing...')
        curated.compute('quality_metrics',
                        metric_names=analyzer.get_extension('quality_metrics').params['metric_names'])
    unit_ids = [int(u) for u in curated.unit_ids]

    # ---- root files: matched spike by spike to the curated analyzer ----
    cur_units, cur_times = analyzer_spikes(curated)
    parents = find_parents(curation, set(orig_unit_ids), root_times, root_units, cur_units, cur_times)
    new_units = assign_root_spikes(parents, orig_index, root_times, root_units, cur_units, cur_times)
    keep = new_units >= 0
    if not same_spikes(new_units[keep], root_times[keep], cur_units, cur_times):
        sys.exit('ERROR: rebuilt root spikes do not match the curated SortingAnalyzer -- nothing written')
    print(f'\nRoot spikes matched to the curated analyzer spike for spike '
          f'({n_spikes} -> {int(keep.sum())} spikes)')

    # ---- Bombcell, same thresholds and options as the pipeline ----
    print('\nBombcell...')
    labels_df = sc.bombcell_label_units(
        curated, thresholds=pl.bombcell_thresholds, split_non_somatic_good_mua=True
    )
    label_col = 'bombcell_label' if 'bombcell_label' in labels_df.columns else labels_df.columns[-1]
    labels = labels_df[label_col].tolist()
    all_metrics_df = curated.get_metrics_extension_data()
    all_metrics_df.index.name = 'unit_id'
    labels_df = labels_df.join(pl.bombcell_failed_thresholds(all_metrics_df, pl.bombcell_thresholds))

    # one row per new unit: origin, parents, censored spikes, GUI labels
    rows = []
    for u, (origin, units) in parents.items():
        if origin == 'unchanged':
            continue
        n_parents = int(np.isin(root_units, units).sum())
        n_now = int((new_units == u).sum())
        row = {'cluster_id': u, 'origin': origin, 'parent_ids': ' '.join(map(str, units)),
               'n_spikes_censored': n_parents - n_now if origin == 'merge' else 0,
               'bombcell_label': labels_df.loc[u, label_col],
               'failed_thresholds': labels_df.loc[u, 'failed_thresholds']}
        for key in curation.label_definitions:
            if key in curated.sorting.get_property_keys():
                row[f'gui_{key}'] = curated.sorting.get_property(key)[unit_ids.index(u)]
        rows.append(row)
    units_df = pd.DataFrame(rows).sort_values('cluster_id') if rows else pd.DataFrame()

    print('\nWriting outputs...')
    backup(ks4_dir, os.path.join(ks4_dir, f'curation_backup_{stamp}'))

    labels_df.to_csv(os.path.join(bc_dir, 'unit_labels.csv'))
    all_metrics_df.to_csv(os.path.join(bc_dir, 'all_metrics.csv'))
    units_df.to_csv(os.path.join(bc_dir, 'curation_units.csv'), index=False)
    shutil.copy2(curation_json, os.path.join(bc_dir, 'applied_curation.json'))
    print(f'    Labels: {labels_df[label_col].value_counts().to_dict()}')

    try:
        pl.plot_bombcell_metric_histograms(
            all_metrics_df, pl.bombcell_thresholds, os.path.join(bc_dir, 'metric_histograms.png'))
        pl.plot_bombcell_upset(curated, labels_df[label_col], pl.bombcell_thresholds, bc_dir)
        pl.plot_bombcell_waveform_classification(
            curated, labels_df[label_col], os.path.join(bc_dir, 'waveform_classification.png'))
    except Exception as e:
        print(f'    Warning: Bombcell summary figures failed ({e})')
    plt.close('all')

    pl.save_rc2_compatible_files(curated, labels, ks4_dir)

    root['spike_clusters.npy'] = new_units.astype(root['spike_clusters.npy'].dtype).reshape(
        root['spike_clusters.npy'].shape)
    for f, arr in root.items():
        np.save(os.path.join(ks4_dir, f), arr[keep])
    print(f'  Wrote {", ".join(PER_SPIKE_FILES)} '
          f'({len(orig_unit_ids)} -> {len(unit_ids)} clusters)')

    # ---- curated analyzer becomes sorting_analyzer/ ----
    before_dir = os.path.join(ks4_dir, f'sorting_analyzer_before_{stamp}')
    del analyzer, curated
    gc.collect()
    try:
        os.rename(analyzer_dir, before_dir)
        os.rename(curated_dir, analyzer_dir)
        print(f'  sorting_analyzer/ -> {os.path.basename(before_dir)}/, '
              f'sorting_analyzer_curated/ -> sorting_analyzer/')
    except OSError as e:
        print(f'  Warning: could not swap the analyzer folders ({e}) -- rename them by hand:\n'
              f'    sorting_analyzer -> {os.path.basename(before_dir)}\n'
              f'    sorting_analyzer_curated -> sorting_analyzer')

    if len(units_df):
        print('\nNew units:')
        print(units_df.to_string(index=False))
    print('\nphy/ is now out of date (not updated by this script).')
    print('Done. Next, in MATLAB: ctl.create_check_clusters_csv(probe_id)')


if __name__ == '__main__':
    main()
