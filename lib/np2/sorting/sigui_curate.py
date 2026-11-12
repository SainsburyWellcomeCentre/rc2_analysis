#!/usr/bin/env python3
"""
Open a pipeline session in SpikeInterface-GUI for manual curation.

  conda activate spikeinterface
  python sigui_curate.py "<...>\\imec0_ks4"

Equivalent to `sigui --curation <ks4_dir>\\sorting_analyzer`, plus:
  - the unit list also shows each unit's Bombcell label (bombcell/unit_labels.csv)
    and rp_contamination
  - and its N_MATCHES most similar units (template_similarity, cosine), as
    match_1 ... match_N = "unit_id (similarity, Bombcell label)", plus
    best_similarity (sortable) -- the GUI's merge search runs over all units
    and cannot be restricted to e.g. good units
  - correlograms / isi_histograms (the GUI's correlogram and ISI views) are
    computed and saved into the analyzer first if missing -- sessions sorted
    before spikeGLX_pipeline_np2.py computed them

In the GUI: merge, split, remove and label units, then press "Save in
analyzer" in the curation panel. This writes
sorting_analyzer/spikeinterface_gui/curation_data.json -- nothing else is
changed until apply_curation.py is run on the session.
"""

import os
import sys
import argparse

import numpy as np
import pandas as pd


N_MATCHES = 5


def top_matches(analyzer, labels, n=N_MATCHES):
    """Unit-list columns with each unit's n most similar units (template_similarity)."""
    sim = np.array(analyzer.get_extension('template_similarity').get_data(), dtype=float)
    np.fill_diagonal(sim, -np.inf)
    unit_ids = np.asarray(analyzer.unit_ids)
    order = np.argsort(-sim, axis=1)[:, :n]
    cols = {'best_similarity': np.round(sim[np.arange(len(unit_ids)), order[:, 0]], 3)}
    for k in range(n):
        cols[f'match_{k + 1}'] = np.array(
            [f'{unit_ids[j]} ({sim[i, j]:.2f}, {labels.get(int(unit_ids[j]), "?")})'
             for i, j in enumerate(order[:, k])], dtype=str)
    return cols


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('ks4_dir', help='session imec0_ks4 folder (contains sorting_analyzer/)')
    ks4_dir = parser.parse_args().ks4_dir

    import spikeinterface.full as si
    from spikeinterface_gui import run_mainwindow

    si.set_global_job_kwargs(n_jobs=12, chunk_duration='1s', progress_bar=True)

    analyzer = si.load_sorting_analyzer(os.path.join(ks4_dir, 'sorting_analyzer'))

    for ext in ('correlograms', 'isi_histograms'):
        if not analyzer.has_extension(ext):
            print(f'Computing {ext} (saved into the analyzer, done once per session)...')
            analyzer.compute(ext)

    extra = {}
    labels = {}
    labels_csv = os.path.join(ks4_dir, 'bombcell', 'unit_labels.csv')
    if os.path.isfile(labels_csv):
        bc = pd.read_csv(labels_csv, index_col=0)['bombcell_label']
        unit_ids = [int(u) for u in analyzer.unit_ids]
        if set(bc.index) == set(unit_ids):
            labels = bc.to_dict()
            extra['bombcell_label'] = np.array([labels[u] for u in unit_ids], dtype=str)
        else:
            print('Warning: bombcell/unit_labels.csv does not match the analyzer units -- label column not shown')
    if analyzer.has_extension('template_similarity'):
        extra.update(top_matches(analyzer, labels))

    # GUI defaults (spikeinterface.widgets.sorting_summary) + rp_contamination; extra
    # columns are appended by the GUI itself
    displayed = ['firing_rate', 'num_spikes', 'x', 'y', 'amplitude_median', 'snr',
                 'rp_violations', 'rp_contamination']

    run_mainwindow(analyzer, mode='desktop', curation=True,
                   displayed_unit_properties=displayed, extra_unit_properties=extra or None)


if __name__ == '__main__':
    sys.exit(main())
