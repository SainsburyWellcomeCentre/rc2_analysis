#!/usr/bin/env python3
"""
SpikeInterface-based electrophysiology pipeline for NP2 probes (SpikeGLX data).
Tested with SpikeInterface 0.104.x.

Preprocessing chain (no external tools -- SpikeInterface end to end):
  [a] bandpass_filter (300-9000 Hz) + phase_shift (ADC sample-delay correction;
      inter_sample_shift is read automatically from the .meta). Saved to disk
      once, as bandpass_only.ap/ -- source for cross-session matching (see
      lib/np2/matching/), which needs the spatial footprint step [b] removes.
      No transient-artifact repair -- see README; if ever needed,
      spikeinterface.preprocessing.detect_and_remove_artifacts.
  [b] detect_bad_channels -> interpolate_bad_channels -> highpass_spatial_filter
             (IBL destriping; handles non-uniform stripe noise better than median-based CAR)
             Applied PER SHANK (channel group) -- see ibl_destripe_by_shank().
             Both detect_bad_channels (coherence+psd) and highpass_spatial_filter
             are spatial and must be computed within a single shank, so on
             multi-shank NP2.0 probes (2- or 4-shank) the recording is split by
             group, each shank destriped, then reaggregated. Shank count is read
             from the data (NP2.0 is the only fixed probe assumption).
             This is the recording actually passed to Kilosort4.

Sorting:
  Kilosort4 via SpikeInterface (do_CAR=False -- destriping already handled above)

Post-processing:
  SortingAnalyzer: waveforms, templates, unit locations, spike amplitudes,
                   template metrics, PCA, quality metrics
  Bombcell: automated unit classification (good / mua / noise / non_soma)
  Phy export: for visualization

Output (rc2_analysis / FileManager.m compatible):
  imec{prb}_ks4/
    spike_times.npy, spike_templates.npy, spike_clusters.npy,
    amplitudes.npy, templates.npy, channel_map.npy, channel_positions.npy, params.py
    cluster_groups.csv          <- Bombcell labels
    csv/
      metrics.csv               <- every quality_metrics column, SpikeInterface's own names
      waveform_metrics.csv      <- every template_metrics column, SpikeInterface's own names
    sorting_analyzer/           <- SortingAnalyzer binary folder (reloadable)
    phy/                        <- Phy visualization output
    bombcell/                   <- Bombcell thresholds and results JSON
  bandpass_only.ap/             <- band-pass + phase-shifted, NOT destriped;
                                    source recording for cross-session matching

TPrime: REMOVED. RC2's sync pulse is wired into the probe's own AP stream
(channel SY0), same clock as the neural data, so there's no cross-clock drift
to correct -- TPrime was never doing anything here. Would matter again only
for multiple simultaneous probes / a separate NI-DAQ stream: re-extract SY
edges with SpikeInterface (threshold + write rising-edge times to .txt) and
feed them to TPrime.exe.

"""

import os
import re
import shutil
import fnmatch
import numpy as np
import pandas as pd
from datetime import datetime


# ============================================================
# Auto-filled per session by MATLAB -- do not edit these here
# ============================================================
#
# NOTE: this file is a TEMPLATE, not run directly by the MATLAB pipeline.
# SortingHelper.m (overwrite_sorting_script) fills in the 5 values below for
# the session being run and writes the result to
# lib/np2/sorting/_generated/spikeGLX_pipeline_session.py, which is what
# actually gets executed. That generated file is overwritten on every run
# and is not tracked in git. Editing these values here only changes what a
# STANDALONE run (no MATLAB) uses -- see sorting/README.txt, "OPTION B".

# logName: log file name (saved in output_dest)
# npx_directory: raw data dir (parent of the SpikeGLX run folder)
# run_specs: [run_name, gate, trigger_string, probe_string]
# output_dest: all pipeline output goes under this directory
# start_step: see run_from_step / run_sorting_from_step in RC2Preprocess.m
logName = 'pipeline_log.csv'
npx_directory = r'D:\data\myrecording'
run_specs = [['myrecording', '0', '0,0', '0']]

# (SortingHelper.m's run_specs replacement consumes up to the next '#' --
# keep a comment here so it doesn't also swallow output_dest below.)
output_dest = r'D:\data\myrecording\output'
start_step = 'preprocess'


# ============================================================
# User input -- edit as needed
# ============================================================
# None of this is required reading for a normal run -- every value below
# already has a justified default (see the comment on each). Only change
# something here if you have a specific reason to (a different probe
# geometry, adjusting Bombcell thresholds for an atypical brain region --
# see the main README, "Adjusting Bombcell thresholds" -- or a different
# machine's tool paths).

# ---- Kilosort 4 settings (all 3 are the documented KS4 defaults) ----
ks_nblocks = 5        # non-rigid drift correction blocks; KS4 docs recommend 5 for long probes
ks_Th_universal = 9   # spike detection threshold, universal (not-yet-learned) templates
ks_Th_learned = 8     # spike detection threshold, templates learned from this recording

# ---- Bombcell thresholds (curation used in Step [6]) ----
# These start from Bombcell's own defaults (spikeinterface.curation.
# bombcell_get_default_thresholds()), written out here instead of called at
# runtime so they are visible and editable without digging through the rest
# of this file. 'greater'/'less' are inclusive pass bounds (None = that side
# is unconstrained); a unit fails a category ("noise"/"mua") if it fails ANY
# ONE metric in it -- see the main README, "Adjusting Bombcell thresholds",
# for what each metric means and when you might want to relax one (e.g.
# cerebellar Purkinje cells). Set any bound to None to disable that side of a
# check.
#
# rp_contamination and snr were relaxed from the SI/Bombcell defaults
# (0.1 -> 0.15, 5 -> 4.5) after checking, on 3 animals (CAA-1124370,
# CAA-1124371, CAA-1123244), that rp_contamination alone accounted for the
# large majority of units failing "good" status (267/857 mua units failed
# ONLY this criterion; amplitude_median failed 0 units in isolation, contra
# Bombcell's own general guidance to tune amplitude first). Units recovered
# by this relaxation were spot-checked in Phy (ISI view) before deciding on
# these values. See lib/np2/sorting/README.txt, "Adjusting Bombcell
# thresholds", for the full analysis and the recovered unit IDs.
bombcell_thresholds = {
    'noise': {
        'num_positive_peaks':        {'greater': None,   'less': 2},
        'num_negative_peaks':        {'greater': None,   'less': 1},
        'peak_to_trough_duration':   {'greater': 0.0001, 'less': 0.00115},
        'waveform_baseline_flatness':{'greater': None,   'less': 0.5},
        'peak_after_to_trough_ratio':{'greater': None,   'less': 0.8},
        'exp_decay':                 {'greater': 0.01,   'less': 0.1},
    },
    'mua': {
        'amplitude_median':  {'greater': 30,  'less': None, 'abs': True},
        'snr':                {'greater': 4.5, 'less': None},
        'amplitude_cutoff':   {'greater': None,'less': 0.2},
        'num_spikes':         {'greater': 300, 'less': None},
        'rp_contamination':   {'greater': None,'less': 0.15},
        'presence_ratio':     {'greater': 0.7, 'less': None},
        'drift_ptp':          {'greater': None,'less': 100},
    },
    'non-somatic': {
        'peak_before_to_trough_ratio':    {'greater': None,   'less': 3},
        'peak_before_width':              {'greater': 0.00015,'less': None},
        'trough_width':                   {'greater': 0.0002, 'less': None},
        'peak_before_to_peak_after_ratio':{'greater': None,   'less': 3},
        'main_peak_to_trough_ratio':      {'greater': None,   'less': 0.8},
    },
}

# ============================================================
# End of user input
# ============================================================


# ---- Utility functions (no ecephys dependency) ----

def parse_probe_str(probe_string):
    """Parse '0', '0,3', or '0:3' -> list of probe index strings."""
    prb_list = []
    for substr in probe_string.split(','):
        if ':' in substr:
            lo, hi = substr.split(':')
            prb_list.extend(str(i) for i in range(int(lo), int(hi) + 1))
        else:
            prb_list.append(substr)
    return prb_list


def check_single_trigger_file(raw_dir):
    """Raise if raw_dir has more than one '_t<N>.ap.bin' file.

    SpikeGLX writes a new _t<N> file each time acquisition restarts within
    the same gate (e.g. after a pause/crash without closing the gate).
    SpikeInterface does not concatenate multiple trigger files -- so surface
    this loudly instead of silently sorting only one of them.
    """
    t_files = [f for f in os.listdir(raw_dir) if fnmatch.fnmatch(f, '*_t*.imec*.ap.bin')]
    if len(t_files) > 1:
        raise RuntimeError(
            f'{raw_dir} has {len(t_files)} trigger files ({t_files}) -- '
            'multi-file concatenation is not supported by this pipeline.'
        )


def copy_ks4_outputs_to_parent(ks4_output_dir, recording_preproc=None):
    """
    KS4 run via SpikeInterface puts its output files in sorter_output/ subfolder.
    Copy them to ks4_output_dir/ so rc2_analysis can find them directly.
    Also ensures spike_clusters.npy exists (created from spike_templates.npy if absent).

    pc_features.npy / pc_feature_ind.npy are excluded: they are large
    (~2 GB) and unused by rc2_analysis (FileManager.ks4_npy is only ever
    called with spike_clusters/spike_templates/spike_times/amplitudes/
    templates/channel_map/channel_positions). They remain available in
    sorter_output/ if ever needed.

    channel_map.npy as written by Kilosort4 is ALWAYS np.arange(n_chan) (see
    spikeinterface/sorters/external/kilosort4.py) -- it never reflects the
    original SpikeGLX channel_id, even in the single-shank case. On multi-shank
    NP2.0 probes this is silently wrong: ibl_destripe_by_shank() splits by shank
    and re-aggregates (si.aggregate_channels), which reorders channels to
    [shank0 channels..., shank1 channels..., ...] instead of SpikeGLX's native
    interleaved order. rc2_analysis (SpikeGLXMetaData.electrode_id_from_channel_id)
    looks up shank_id from the *raw* imroTbl using this channel_map value as if it
    were the true channel_id, so most clusters on shank > 0 get assigned the wrong
    shank_id/depth. If recording_preproc is given (the recording actually passed to
    Kilosort4, before its channel_ids get discarded), overwrite channel_map.npy with
    the true SpikeGLX AP channel numbers in the recording's current channel order.
    """
    sorter_out = os.path.join(ks4_output_dir, 'sorter_output')
    skip_files = {'pc_features.npy', 'pc_feature_ind.npy'}
    if os.path.isdir(sorter_out):
        n = 0
        for fname in os.listdir(sorter_out):
            if fname in skip_files:
                continue
            if fname.endswith('.npy') or fname == 'params.py':
                shutil.copy2(os.path.join(sorter_out, fname),
                             os.path.join(ks4_output_dir, fname))
                n += 1
        print(f'  Copied {n} KS4 output files from sorter_output/ to parent')

    if recording_preproc is not None:
        true_channel_map = np.array(
            [int(re.search(r'AP(\d+)', str(ch)).group(1)) for ch in recording_preproc.channel_ids]
        )
        np.save(os.path.join(ks4_output_dir, 'channel_map.npy'), true_channel_map)
        print('  Fixed channel_map.npy to true SpikeGLX AP channel numbers')

    # spike_clusters.npy: created by Phy on first open; pre-create it here
    sc_path = os.path.join(ks4_output_dir, 'spike_clusters.npy')
    st_path = os.path.join(ks4_output_dir, 'spike_templates.npy')
    if not os.path.exists(sc_path) and os.path.exists(st_path):
        shutil.copy2(st_path, sc_path)
        print('  Created spike_clusters.npy (copy of spike_templates.npy)')


def plot_ks4_motion(ks4_output_dir, out_path):
    """Plot Kilosort4's internal drift-correction estimate (ops.npy: dshift,
    yblk) -- one curve per non-rigid block, showing the estimated shift (um)
    over time at that depth. This is KS4's own correction, not a scatter of
    raw spike positions (see create_driftmap in RC2Preprocess.m for that).

    Use this to judge how much drift a session actually has before deciding
    whether to rerun it with nblocks=0 for cross-session matching -- see the
    main README, "Kilosort4 drift correction and cross-session matching".
    """
    import matplotlib.pyplot as plt

    ops_path = os.path.join(ks4_output_dir, 'ops.npy')
    if not os.path.isfile(ops_path):
        print(f'    Warning: ops.npy not found at {ops_path}, skipping motion plot')
        return

    ops = np.load(ops_path, allow_pickle=True).item()
    dshift = ops.get('dshift')
    yblk = ops.get('yblk')
    if dshift is None or yblk is None:
        print('    Warning: ops.npy has no dshift/yblk (rigid registration?), skipping motion plot')
        return

    fs = ops.get('fs', 30000.0)
    batch_size = ops.get('batch_size', 60000)
    t = np.arange(dshift.shape[0]) * batch_size / fs

    fig, ax = plt.subplots(figsize=(10, 4))
    for block_idx in range(dshift.shape[1]):
        ax.plot(t, dshift[:, block_idx], label=f'{yblk[block_idx]:.0f} um', linewidth=1)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Estimated shift (um)')
    # Clip the y-range to the 0.5-99.5 percentile: KS4 can produce a handful
    # of large one-batch outliers unrelated to real probe motion (e.g. a
    # short/partial batch), which would otherwise dominate the axis scale.
    lo, hi = np.percentile(dshift, [0.5, 99.5])
    pad = 0.1 * max(hi - lo, 1.0)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_title(f"KS4 drift correction (mean_drift={ops.get('mean_drift', float('nan')):.1f} um)")
    ax.legend(title='Block depth', fontsize=8, ncol=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'    Saved motion plot -> {out_path}')


def ibl_destripe_by_shank(recording):
    """IBL destriping (detect + interpolate bad channels, then
    highpass_spatial_filter) applied PER SHANK (channel group).

    This is SpikeInterface's documented "Processing a Recording by Channel
    Group" workflow -- SI has no single multi-shank preprocessing call, only
    the primitives (recording.split_by('group') + si.aggregate_channels), so
    we chain them here. Both spatial steps must be computed within one shank:

      * highpass_spatial_filter is a spatial high-pass ALONG probe depth and
        SI raises 'The recording contains multiple groups!' on a multi-group
        recording.
      * detect_bad_channels (default method 'coherence+psd') uses coherence
        across depth and, per the SI/IBL docs, "must be run on individual
        probes/shanks separately".

    Multi-shank NP2.0 probes (this lab uses 2- and 4-shank NP2.0) expose one
    channel group per shank. The shank count is READ FROM THE RECORDING, so 1-,
    2- and 4-shank NP2.0 probes are all handled with no hardcoded shank count.
    A single-group recording is processed directly (no split/aggregate).

    Returns (preprocessed_recording, bad_channel_ids).
    """
    import spikeinterface.full as si
    import spikeinterface.preprocessing as spre

    def _destripe_one(rec):
        bad_ids, _ = spre.detect_bad_channels(rec)
        rec_i = spre.interpolate_bad_channels(rec, bad_ids)
        return spre.highpass_spatial_filter(rec_i), list(bad_ids)

    groups = np.unique(recording.get_channel_groups())
    if len(groups) == 1:
        return _destripe_one(recording)

    print(f'    multi-shank probe: destriping each of {len(groups)} shank(s) separately')
    filtered, all_bad = [], []
    for r in recording.split_by(property='group').values():   # one recording per shank
        rec_f, bad_ids = _destripe_one(r)
        filtered.append(rec_f)
        all_bad.extend(bad_ids)
    return si.aggregate_channels(filtered), all_bad


def save_rc2_compatible_files(analyzer, labels, ks4_output_dir):
    """
    Write output files expected by rc2_analysis (FileManager.m paths):

      ks4_output_dir/
        cluster_groups.csv           <- cluster_id + group (good/mua/noise from Bombcell)
        csv/
          metrics.csv                <- every quality_metrics column, native SI names
          waveform_metrics.csv       <- every template_metrics column, native SI names,
                                         plus peak_channel/amplitude (not part of that
                                         extension -- derived from the mean templates below)

    Column names are SpikeInterface's own (quality_metrics/template_metrics), NOT
    renamed to match the old ecephys_spike_sorting/Kilosort2 pipeline's metrics.csv
    columns (e.g. old 'isi_viol' is SI's 'isi_violations_ratio', old 'max_drift' is
    SI's 'drift_ptp', old 'duration' is SI's 'peak_to_trough_duration'). That renaming
    used to happen here, but several of those pairs are not the same calculation
    (e.g. 'silhouette_score' was a full pairwise silhouette, SI's default 'silhouette'
    is the simplified/centroid-based one -- see SpikeInterface's silhouette_score.rst).
    Reusing the old name implied an equivalence that often doesn't hold. rc2_analysis
    (RC2Format.m) reads these CSVs directly by column name, so it lists SI's own names.
    """
    csv_dir = os.path.join(ks4_output_dir, 'csv')
    os.makedirs(csv_dir, exist_ok=True)
    unit_ids = list(analyzer.unit_ids)

    # ---- metrics.csv (every quality_metrics column, as computed) ----
    qm_ext = analyzer.get_extension('quality_metrics')
    if qm_ext is not None:
        qm_df = qm_ext.get_data().copy()
        qm_df.index.name = 'unit_id'
        qm_df = qm_df.reset_index().rename(columns={'unit_id': 'cluster_id'})
        qm_df.to_csv(os.path.join(csv_dir, 'metrics.csv'), index=False)
        print(f'  Saved metrics.csv ({len(qm_df)} units, {len(qm_df.columns) - 1} metrics)')

    # ---- waveform_metrics.csv (every template_metrics column + peak_channel/amplitude) ----
    tm_ext = analyzer.get_extension('template_metrics')
    if tm_ext is not None:
        tm_df = tm_ext.get_data().copy()
        tm_df.index.name = 'unit_id'
        tm_df = tm_df.reset_index().rename(columns={'unit_id': 'cluster_id'})

        # peak_channel/amplitude are not template_metrics outputs -- derive them
        # from the mean templates (same computation as before).
        templates_ext = analyzer.get_extension('templates')
        if templates_ext is not None:
            tmpl = templates_ext.get_templates(operator='average')  # (n_units, n_samples, n_ch)
            peak_channels, amplitudes_list = [], []
            for i in range(tmpl.shape[0]):
                amp_per_ch = tmpl[i].max(axis=0) - tmpl[i].min(axis=0)
                best_idx = int(np.argmax(amp_per_ch))
                peak_channels.append(best_idx)
                amplitudes_list.append(float(amp_per_ch[best_idx]))

            tm_df = tm_df.merge(
                pd.DataFrame({'cluster_id': unit_ids, 'peak_channel': peak_channels,
                               'amplitude': amplitudes_list}),
                on='cluster_id', how='left'
            )

        wf_path = os.path.join(csv_dir, 'waveform_metrics.csv')
        tm_df.to_csv(wf_path, index=False)
        print(f'  Saved waveform_metrics.csv ({len(tm_df)} units, {len(tm_df.columns) - 1} metrics)')

    # ---- cluster_groups.csv (Bombcell labels) ----
    group_values = labels if labels is not None else (['unsorted'] * len(unit_ids))
    pd.DataFrame({'cluster_id': unit_ids, 'group': group_values}).to_csv(
        os.path.join(ks4_output_dir, 'cluster_groups.csv'), index=False
    )
    print(f'  Saved cluster_groups.csv ({len(unit_ids)} units)')


# Short display name for each Bombcell/SpikeInterface metric column, used
# consistently across every plot this pipeline produces (metric_histograms,
# the upset plots, waveform_classification) AND in the "User input" section
# at the top of this file -- so a name seen on any plot can be found in the
# same form everywhere else, including the threshold dict a user would edit.
METRIC_DISPLAY_NAMES = {
    'num_positive_peaks':               '# peaks',
    'num_negative_peaks':                '# troughs',
    'waveform_baseline_flatness':        'baseline flatness',
    'peak_to_trough_duration':           'waveform duration',
    'peak_after_to_trough_ratio':        'peak$_2$/trough',
    'exp_decay':                         'spatial decay',
    'amplitude_median':                  'amplitude',
    'snr':                               'SNR',
    'rp_contamination':                  'frac. RPVs',
    'num_spikes':                        '# spikes',
    'presence_ratio':                    'presence ratio',
    'amplitude_cutoff':                  'spikes missing',
    'drift_ptp':                         'maximum drift',
    'peak_before_to_peak_after_ratio':   'peak$_1$/peak$_2$',
    'main_peak_to_trough_ratio':         'peak$_{main}$/trough',
    'peak_before_to_trough_ratio':       'peak$_1$/trough',
    'peak_before_width':                 'peak$_1$ width',
    'trough_width':                      'trough width',
    'drift_std':                         'cum. drift',
    'isolation_distance':                'isolation dist.',
    'l_ratio':                           'L-ratio',
}

# Colour convention used on every plot in this pipeline: which label a metric
# drives the unit TOWARDS when it fails is what determines its colour, not
# just pass/fail -- green always means "passes, no effect on the label";
# red/orange/blue mean "fails, and pushes the unit towards noise/mua/non-soma"
# respectively; grey means "not used in labelling at all, no effect either way".
LABEL_COLORS = {
    'pass':        'tab:green',
    'noise':       'tab:red',
    'mua':         'tab:orange',
    'non-somatic': 'tab:blue',
    'neutral':     '0.4',  # grey -- additional/inspection-only metrics, and non-somatic comparisons
}


def plot_bombcell_metric_histograms(metrics_df, thresholds, out_path):
    """
    Recreates the native Bombcell (MATLAB/Python) quality_metrics_distribution
    layout -- short human-readable axis labels (METRIC_DISPLAY_NAMES, shared
    with the upset plots and waveform_classification.png), a fraction-of-units
    y-axis, a min/max legend per panel, and a green/red/orange/blue bar under
    each x-axis showing which range of that metric is accepted vs. rejected --
    grouped into labelled rows by what each metric is actually USED for
    (noise / mua / non-somatic / not used at all), coloured to match (green =
    passes; red/orange/blue = fails and pushes the unit towards noise/mua/
    non-somatic respectively; grey = no threshold, inspection only), so a
    reader can see at a glance which panels affect which part of the
    good/mua/noise/non-soma decision before touching a threshold.

    21 panels total:
      - 6 noise + 7 mua + 5 non-somatic = 18 metrics that drive the labelling
        (bombcell_label_units, see bombcell_get_default_thresholds).
      - 3 additional metrics (drift_std, isolation_distance, l_ratio) that
        Bombcell computes and displays for manual inspection but that never
        affect good/mua/noise/non-soma (confirmed against the native Bombcell
        classification.py: these three do not appear in its labelling logic
        at all, only in an optional summary table).

    Two of the 18 labelling metrics -- num_positive_peaks/num_negative_peaks
    are shared with 'noise' by name only; every metric otherwise appears in
    exactly one row. peak_before_to_trough_ratio, peak_before_width and
    trough_width are part of Bombcell's own non-somatic decision (see
    bombcell_curation.py: is_non_somatic = (ratio_conditions AND
    width_conditions) OR large_main_peak) but are absent from the native
    Bombcell plotting function itself (get_metric_info_list has no entry for
    them at all) -- included here anyway so every metric that can flip a
    unit's label is visible somewhere.

    SpikeInterface's own sw.plot_metric_histograms uses raw SI column names
    as axis labels (e.g. 'peak_before_width' in seconds, unreadable at 1e-4
    scale), has no colored accept/reject bar, and only plots metrics that
    have a threshold in `thresholds` -- so it never shows the 3 additional
    metrics above. This function reads the same `thresholds` dict (from
    bombcell_get_default_thresholds) plus the 3 additional metrics directly.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # (SI column, short label, unit scale factor, unit suffix, take_abs,
    #  upper percentile clip, integer_valued) -- take_abs mirrors the
    # 'abs': True flag bombcell_get_default_thresholds sets for
    # amplitude_median (amplitude is signed in SI, Bombcell thresholds it
    # unsigned). upper percentile clip guards metrics like isolation_distance
    # that can have a handful of near-infinite outliers (isolated/near-empty
    # clusters) that would otherwise squash the whole histogram into one bin.
    # integer_valued: use one bin per integer instead of a fixed 30 bins --
    # # peaks/# troughs only take small integer values (0, 1, 2, 3...), and
    # 30 evenly-spaced bins over that range slices individual integers into
    # several thin, unreadable bars instead of the wide/clear per-value bars
    # the native Bombcell plot shows.
    #
    # Grouped into labelled sections -- each row of panels is one section,
    # so the reader sees which part of the good/mua/noise/non-soma decision
    # a given panel feeds into (see docstring above).
    # (SI column, unit scale factor, unit suffix, take_abs, upper percentile
    #  clip, integer_valued) -- display name comes from METRIC_DISPLAY_NAMES
    #  (shared with the upset plots and waveform_classification), not repeated
    #  here. take_abs mirrors the 'abs': True flag bombcell_get_default_
    #  thresholds sets for amplitude_median (amplitude is signed in SI,
    #  Bombcell thresholds it unsigned). upper percentile clip guards metrics
    #  like isolation_distance that can have a handful of near-infinite
    #  outliers (isolated/near-empty clusters) that would otherwise squash the
    #  whole histogram into one bin. integer_valued: one bin per integer
    #  instead of a fixed 30 bins -- # peaks/# troughs only take small integer
    #  values (0, 1, 2, 3...), and 30 evenly-spaced bins over that range
    #  slices individual integers into several thin, unreadable bars instead
    #  of the wide/clear per-value bars the native Bombcell plot shows.
    #
    # Grouped into labelled sections -- each row of panels is one section, so
    # the reader sees which part of the good/mua/noise/non-soma decision a
    # given panel feeds into (see docstring above). 'label_color' is the
    # LABEL_COLORS key used for that row's accept/reject bar.
    sections = [
        ('NOISE metrics (fail any one → labelled "noise")', 'noise', [
            ('num_positive_peaks',             1,   '',   False, None, True),
            ('num_negative_peaks',             1,   '',   False, None, True),
            ('waveform_baseline_flatness',     1,   '',   False, None, False),
            ('peak_to_trough_duration',        1e6, 'µs', False, 99,   False),
            ('peak_after_to_trough_ratio',     1,   '',   False, 99,   False),
            ('exp_decay',                      1,   '',   False, 99,   False),
        ]),
        ('MUA metrics (not noise, fail any one → labelled "mua")', 'mua', [
            ('amplitude_median',               1,   'µV', True,  99,   False),
            ('snr',                            1,   '',   False, 99,   False),
            ('rp_contamination',               1,   '',   False, None, False),
            ('num_spikes',                     1,   '',   False, None, False),
            ('presence_ratio',                 1,   '',   False, None, False),
            ('amplitude_cutoff',               100, '%',  False, None, False),
            ('drift_ptp',                      1,   'µm', False, None, False),
        ]),
        ('NON-SOMATIC metrics (combined rule → "non_soma_good"/"non_soma_mua", shown together)', 'non-somatic', [
            ('peak_before_to_peak_after_ratio',1,   '',   False, 99,   False),
            ('main_peak_to_trough_ratio',      1,   '',   False, 99,   False),
            ('peak_before_to_trough_ratio',    1,   '',   False, 99,   False),
            ('peak_before_width',              1e6, 'µs', False, None, False),
            ('trough_width',                   1e6, 'µs', False, None, False),
        ]),
        ('Additional metrics (computed for inspection, not used in labelling)', 'neutral', [
            ('drift_std',                      1,   'µm', False, 99, False),
            # isolation_distance can carry a handful of near-numerically-infinite
            # outliers (division by a near-zero covariance for isolated/sparse
            # clusters) -- up to 1e15 on real data, dwarfing every other unit's
            # value. A 90th-percentile clip (rather than 99th) is needed to keep
            # the histogram readable; the outlier units themselves are unaffected
            # (still in all_metrics.csv / the actual Bombcell threshold check).
            ('isolation_distance',             1,   '',   False, 90, False),
            ('l_ratio',                        1,   '',   False, 95, False),
        ]),
    ]
    # Flatten noise/mua/non-somatic sections into one lookup, same as
    # bombcell_failed_thresholds -- greater/less bounds per SI metric name.
    flat_thresh = {}
    for section in thresholds.values():
        flat_thresh.update(section)

    n_cols = max(len(panels) for _, _, panels in sections)
    n_rows = len(sections)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.2 * n_cols, 2.9 * n_rows))
    axes = np.atleast_2d(axes)

    for row, (section_title, label_color_key, panels) in enumerate(sections):
        fail_color = LABEL_COLORS[label_color_key]
        pass_color = LABEL_COLORS['pass']
        for col_idx in range(n_cols):
            ax = axes[row, col_idx]
            if col_idx >= len(panels):
                ax.axis('off')
                continue

            col, scale, suffix, take_abs, upper_pct, integer_valued = panels[col_idx]
            short_label = METRIC_DISPLAY_NAMES[col]
            if col not in metrics_df.columns:
                ax.set_title(f'{short_label}\n(not computed)')
                ax.axis('off')
                continue

            values = metrics_df[col].to_numpy(dtype=float)
            if take_abs:
                values = np.abs(values)
            values = values * scale
            values = values[np.isfinite(values)]
            if upper_pct is not None and len(values) > 0:
                values = values[values <= np.percentile(values, upper_pct)]
            if len(values) == 0:
                ax.set_title(f'{short_label}\n(no valid data)')
                ax.axis('off')
                continue

            if integer_valued:
                lo, hi = int(np.floor(values.min())), int(np.ceil(values.max())) + 1
                bins = np.arange(lo, hi + 1) - 0.5  # bin edges centered on each integer
            else:
                bins = 30
            counts, bin_edges = np.histogram(values, bins=bins)
            frac = counts / counts.sum() if counts.sum() > 0 else counts
            ax.bar(bin_edges[:-1], frac, width=np.diff(bin_edges), align='edge',
                   color=plt.cm.tab20(col_idx % 20), edgecolor='black', linewidth=0.5)

            bounds = flat_thresh.get(col, {})
            greater = bounds.get('greater', None)
            less = bounds.get('less', None)
            xmin, xmax = float(values.min()), float(values.max())
            xspan = max(xmax - xmin, 1e-12)
            pad = 0.03 * xspan
            xlo, xhi = xmin - pad, xmax + pad

            def _scaled(v, scale=scale):
                return v * scale if v is not None else None

            g = _scaled(greater)
            l = _scaled(less)
            # 3-segment accept/reject bar: green = passes (no effect on the
            # label), fail_color = fails and pushes the unit towards this
            # row's label (red=noise, orange=mua, blue=non-somatic) --
            # neutral grey throughout for the "additional metrics" row, which
            # has no threshold at all (never fails a unit either way).
            y0 = ax.get_ylim()
            bar_y = -0.04 * (y0[1] if y0[1] > 0 else 1)
            if label_color_key == 'neutral':
                ax.plot([xlo, xhi], [bar_y, bar_y], color=LABEL_COLORS['neutral'], lw=4, solid_capstyle='butt')
            elif g is not None and l is not None:
                ax.plot([xlo, g], [bar_y, bar_y], color=fail_color, lw=4, solid_capstyle='butt')
                ax.plot([g, l], [bar_y, bar_y], color=pass_color, lw=4, solid_capstyle='butt')
                ax.plot([l, xhi], [bar_y, bar_y], color=fail_color, lw=4, solid_capstyle='butt')
            elif g is not None:
                ax.plot([xlo, g], [bar_y, bar_y], color=fail_color, lw=4, solid_capstyle='butt')
                ax.plot([g, xhi], [bar_y, bar_y], color=pass_color, lw=4, solid_capstyle='butt')
            elif l is not None:
                ax.plot([xlo, l], [bar_y, bar_y], color=pass_color, lw=4, solid_capstyle='butt')
                ax.plot([l, xhi], [bar_y, bar_y], color=fail_color, lw=4, solid_capstyle='butt')
            else:
                # no threshold defined for this metric even though its row
                # normally has one (shouldn't happen for the 18 labelling
                # metrics, but guards against a future mismatch) -- neutral.
                ax.plot([xlo, xhi], [bar_y, bar_y], color=LABEL_COLORS['neutral'], lw=4, solid_capstyle='butt')

            ax.set_xlim(xlo, xhi)

            # Headroom above the tallest bar so the min/max legend (below)
            # never overlaps a bar, even when the tallest bar is near the top
            # of the axis (e.g. presence ratio, frac. RPVs).
            y0, y1 = ax.get_ylim()
            ax.set_ylim(y0, y1 * 1.22)

            # Small threshold legend, top-right of each panel -- the
            # accept/reject bar shows WHERE the pass range falls but not its
            # exact numeric bound(s), which colour alone can't convey
            # precisely. Shows the PASS-RANGE THRESHOLD(S) (g/l, already
            # computed above and scaled to match the axis), not the min/max
            # of the data -- e.g. '# peaks' shows 'max=2' (the bombcell_
            # thresholds bound), not the largest value actually observed.
            # No threshold at all (the "additional metrics" row) -> no legend.
            legend_parts = []
            if g is not None:
                legend_parts.append(f'min={g:.3g}')
            if l is not None:
                legend_parts.append(f'max={l:.3g}')
            if legend_parts:
                ax.text(0.98, 0.97, '\n'.join(legend_parts), transform=ax.transAxes, ha='right', va='top',
                         fontsize=6.5, color='0.3', linespacing=1.3)

            # Adaptive tick spacing: matplotlib's default locator often picks
            # too few/too coarse ticks when the real data only spans a small
            # fraction of a panel with a long outlier tail (e.g. peak2/trough
            # mostly 0-3 got ticks every 5 up to 10) -- ask for more
            # candidate ticks over the ACTUAL data range so the axis reflects
            # what the histogram really shows instead of a generic default.
            if integer_valued:
                ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=min(10, hi - lo)))
            else:
                ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=8, min_n_ticks=5))

            if col == 'num_spikes':
                # scientific notation (x10^5) instead of raw ticks -- 6-digit
                # spike counts as bare numbers ("100000, 200000...") are hard
                # to read at a glance and don't match this panel's own axis
                # label convention (unit/scale factor shown once, not per tick).
                ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
            else:
                # Trim trailing zeros on every tick ('0.20' -> '0.2') without
                # touching ticks that need their full precision ('0.25' stays
                # '0.25') -- '%g' drops only the zeros that carry no
                # information, on every panel, not just ones with a unit
                # suffix (the unit itself, if any, lives in the axis label,
                # not on the ticks -- see xlabel below).
                ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:g}'))

            if col_idx == 0:
                ax.set_ylabel('frac. units')
            ax.set_xlabel(f'{short_label} ({suffix})' if suffix else short_label)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    # extra vertical spacing between rows so a section title (added below,
    # after the layout is final) has room above its row of panels without
    # overlapping the x-axis labels of the row above it
    fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=3.5)

    # section titles, centred above each row of panels, bold -- added AFTER
    # tight_layout so axes positions (used to place each title) are final
    for row, (section_title, _, _) in enumerate(sections):
        row_axes = [a for a in axes[row] if a.get_subplotspec() is not None]
        if not row_axes:
            continue
        left = row_axes[0].get_position().x0
        right = row_axes[-1].get_position().x1
        top = max(a.get_position().y1 for a in row_axes)
        fig.text((left + right) / 2, top + 0.012, section_title,
                  ha='center', va='bottom', fontsize=10, fontweight='bold')

    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_bombcell_upset(analyzer, unit_labels, thresholds, out_dir):
    """
    One UpSet plot per label category (noise / mua / non-somatic), showing
    which combinations of failed metrics occur together within that category
    -- reimplemented rather than using SpikeInterface's own
    sw.plot_bombcell_labels_upset (spikeinterface.widgets.BombcellUpsetPlotWidget)
    for two reasons:
      1. That widget always keeps non_soma_good/non_soma_mua as two separate
         plots; this pipeline shows a single combined "non-somatic" plot
         instead (the axonal/dendritic waveform-shape failure reasons are the
         same regardless of whether the unit was otherwise good or mua).
      2. That widget labels each row with the raw SpikeInterface column name
         (e.g. 'peak_after_to_trough_ratio'), not the short names
         (METRIC_DISPLAY_NAMES) used on metric_histograms.png -- reusing its
         internal helpers (_get_metrics_for_unit_label, _build_failure_table)
         but renaming columns before handing off to `upsetplot` keeps the two
         plot types easy to cross-reference.

    Writes one file per category that has at least one failing unit:
      noise_units_upset.png, mua_units_upset.png, non_somatic_units_upset.png
    Returns the list of (category, path) pairs actually written.
    """
    import warnings
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from spikeinterface.widgets.bombcell_curation import BombcellUpsetPlotWidget

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=FutureWarning, module='upsetplot')
            from upsetplot import UpSet, from_memberships
    except ImportError:
        return []

    metrics = analyzer.get_metrics_extension_data()
    unit_labels = pd.Series(unit_labels).to_numpy()

    # non-somatic is a virtual category here: real labels are non_soma_good /
    # non_soma_mua (or non_soma if split_non_somatic_good_mua=False) -- fold
    # them into one boolean mask + one shared metric list (the 5 non-somatic
    # thresholds are identical regardless of the unit's good/mua quality).
    categories = {
        'noise': (unit_labels == 'noise', thresholds.get('noise', {})),
        'mua':   (unit_labels == 'mua',   thresholds.get('mua', {})),
        'non-somatic': (
            np.isin(unit_labels, ['non_soma', 'non_soma_good', 'non_soma_mua']),
            thresholds.get('non-somatic', {}),
        ),
    }

    # Reuse SpikeInterface's own failure-table builder (handles NaN/abs/
    # greater-less bounds identically to bombcell_label_units) rather than
    # duplicating that logic here.
    dummy = BombcellUpsetPlotWidget.__new__(BombcellUpsetPlotWidget)
    failure_table = dummy._build_failure_table(metrics, thresholds)
    failure_table = failure_table.rename(columns=METRIC_DISPLAY_NAMES)

    written = []
    for category, (mask, cat_thresholds) in categories.items():
        n_units = int(np.sum(mask))
        if n_units == 0:
            continue

        relevant_cols = [METRIC_DISPLAY_NAMES.get(m, m) for m in cat_thresholds]
        available_cols = [c for c in relevant_cols if c in failure_table.columns]
        if not available_cols:
            continue

        unit_failures = failure_table.loc[mask, available_cols]
        memberships = [
            unit_failures.columns[unit_failures.loc[idx]].tolist()
            for idx in unit_failures.index
        ]
        memberships = [m for m in memberships if m]
        if not memberships:
            continue

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=FutureWarning, module='upsetplot')
            upset_data = from_memberships(memberships)
            if len(upset_data) == 0:
                continue
            fig = plt.figure(figsize=(12, 6))
            UpSet(upset_data, subset_size='count', show_counts=True,
                  sort_by='cardinality', sort_categories_by='cardinality').plot(fig=fig)
        fig.suptitle(f'{category} (n={n_units})', fontsize=14, y=1.02)

        out_name = {'noise': 'noise_units_upset', 'mua': 'mua_units_upset',
                    'non-somatic': 'non_somatic_units_upset'}[category]
        out_path = os.path.join(out_dir, f'{out_name}.png')
        # bbox_inches='tight' so the suptitle (e.g. 'noise (n=48)') is not
        # cropped out of frame.
        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        written.append((category, out_path))

    return written


def plot_bombcell_waveform_classification(analyzer, unit_labels, out_path):
    """
    One panel per label category (good / mua / noise / non-somatic), each
    showing every unit's template waveform (best channel) overlaid, so a
    reader can compare the shapes Bombcell put in each bucket at a glance.

    Reimplemented rather than using SpikeInterface's own sw.plot_unit_labels
    (spikeinterface.widgets.WaveformOverlayByLabelWidget) for two reasons:
      1. That widget crashes on this pipeline's data with
         "AttributeError: 'NoneType' object has no attribute 'set_visible'"
         (spikeinterface/widgets/unit_labels.py:126) -- a grid-sizing bug
         triggered when the number of distinct labels present exactly fills
         the subplot grid it computes, leaving no "extra" axis to hide.
      2. That widget always keeps non_soma_good/non_soma_mua as two separate
         panels; this pipeline shows a single combined "non-somatic" panel
         instead (4 panels total: good / mua / noise / non-somatic), matching
         plot_bombcell_upset's category grouping.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    unit_labels = pd.Series(unit_labels).to_numpy()
    templates_ext = analyzer.get_extension('templates')
    templates = templates_ext.get_data()  # (n_units, n_samples, n_channels)

    categories = {
        'good':  unit_labels == 'good',
        'mua':   unit_labels == 'mua',
        'noise': unit_labels == 'noise',
        'non-somatic': np.isin(unit_labels, ['non_soma', 'non_soma_good', 'non_soma_mua']),
    }

    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5), sharey=True)
    for ax, (category, mask) in zip(axes, categories.items()):
        n_units = int(np.sum(mask))
        if n_units == 0:
            ax.set_title(f'{category} (n=0)')
            ax.text(0.5, 0.5, 'No units', ha='center', va='center', transform=ax.transAxes)
        else:
            alpha = max(0.05, min(0.3, 10 / n_units))
            for unit_idx in np.where(mask)[0]:
                template = templates[unit_idx]
                best_chan = np.argmax(np.max(np.abs(template), axis=0))
                ax.plot(template[:, best_chan], color='black', alpha=alpha, linewidth=0.5)
            ax.set_title(f'{category} (n={n_units})')
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def bombcell_failed_thresholds(metrics_df, thresholds):
    """
    Per unit, list which Bombcell threshold metrics that unit fails, e.g.
    'noise: num_positive_peaks, num_negative_peaks | mua: snr'. Empty string
    for units that pass every threshold in both sections.

    Only covers the 'noise' and 'mua' threshold sections, whose pass/fail
    rule is a simple AND across metrics (a unit fails the section if it
    fails ANY metric in it) -- this matches SpikeInterface's
    threshold_metrics_label_units exactly, verified against
    bombcell_label_units' own noise/mua labels.

    Deliberately excludes the 'non-somatic' section: Bombcell's actual rule
    there is NOT a simple per-metric AND (see bombcell_label_units source) --
    it combines a width OR-check, a ratio AND-check and a peak-ratio check
    into (narrow_width AND large_ratio) OR large_main_peak. Reimplementing
    that here risked silently mislabeling a unit's non-somatic reason, so it
    is left out rather than shown with unverified per-metric detail.
    """
    import numpy as np
    per_unit_reasons = {uid: [] for uid in metrics_df.index}

    def failed_metrics_for(uid, section_thresholds, nan_fails):
        failed_metrics = []
        for metric_name, bounds in section_thresholds.items():
            if metric_name not in metrics_df.columns:
                continue
            min_value = bounds.get('greater', None)
            max_value = bounds.get('less', None)
            if min_value is None and max_value is None:
                continue
            value = metrics_df.at[uid, metric_name]
            if bounds.get('abs', False):
                value = abs(value)
            if np.isnan(value):
                if nan_fails:
                    failed_metrics.append(metric_name)
                continue
            if min_value is not None and value < min_value:
                failed_metrics.append(metric_name)
            elif max_value is not None and value > max_value:
                failed_metrics.append(metric_name)
        return failed_metrics

    # mua is only evaluated on units that already passed noise (matches
    # bombcell_label_units: mua_labels only covers non_noise_indices), and
    # its NaN policy is 'ignore' (a NaN metric neither passes nor fails it)
    # rather than noise's 'fail'.
    noise_thresholds = thresholds.get('noise', {})
    mua_thresholds = thresholds.get('mua', {})
    for uid in metrics_df.index:
        noise_failed = failed_metrics_for(uid, noise_thresholds, nan_fails=True)
        if noise_failed:
            per_unit_reasons[uid].append(f'noise: {", ".join(noise_failed)}')
            continue
        mua_failed = failed_metrics_for(uid, mua_thresholds, nan_fails=False)
        if mua_failed:
            per_unit_reasons[uid].append(f'mua: {", ".join(mua_failed)}')

    return pd.Series(
        {uid: ' | '.join(reasons) for uid, reasons in per_unit_reasons.items()},
        name='failed_thresholds',
    )


# ============================================================
# Main pipeline
# ============================================================

def main():
    # matplotlib backend must be forced to Agg (headless) BEFORE anything
    # that imports pyplot -- spikeinterface.widgets does so as a side effect
    # of its own import, which would otherwise lock in whatever interactive
    # backend matplotlib picks by default on this machine (e.g. TkAgg/QtAgg).
    # That backend can then fail at first render time, deep inside the
    # Bombcell summary-figure calls, rather than at import time -- seen in
    # practice as 'bombcell done' in the log with elapsed=0 and no figures
    # written, because the failure landed inside that step's try/except.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    import spikeinterface.full as si
    import spikeinterface.preprocessing as spre
    import spikeinterface.sorters as ss
    import spikeinterface.exporters as sexp
    import spikeinterface.curation as sc

    # Parallelise SortingAnalyzer extension computation (waveforms,
    # spike_amplitudes, spike_locations, ...) and export_to_phy -- these
    # default to n_jobs=1 (single-threaded) otherwise. n_jobs is kept below
    # the machine's full core count to leave headroom for other work (e.g.
    # MATLAB, Phy) running at the same time.
    si.set_global_job_kwargs(n_jobs=12, chunk_duration='1s', progress_bar=True)

    valid_start_steps = ('preprocess', 'kilosort4', 'postprocess', 'bombcell')
    if start_step not in valid_start_steps:
        raise ValueError(f"start_step must be one of {valid_start_steps}, got {start_step!r}")
    do_preprocess        = start_step == 'preprocess'
    do_kilosort4         = start_step in ('preprocess', 'kilosort4')
    do_sorting_analyzer  = start_step in ('preprocess', 'kilosort4', 'postprocess')
    print(f'start_step = {start_step!r}  '
          f'(Preprocess: {do_preprocess}, Kilosort4: {do_kilosort4}, SortingAnalyzer: {do_sorting_analyzer})')

    os.makedirs(output_dest, exist_ok=True)
    logFullPath = os.path.join(output_dest, logName)
    with open(logFullPath, 'w') as f:
        f.write('session_id,timestamp,step,status,n_units,elapsed_s\n')

    def log_step(sid, step, status, n_units='', elapsed=''):
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(logFullPath, 'a') as f:
            f.write(f'{sid},{ts},{step},{status},{n_units},{elapsed}\n')

    for spec in run_specs:
        run_name, gate, trig_str, probe_str = spec[0], spec[1], spec[2], spec[3]

        print(f'\n{"="*60}')
        print(f'Run: {run_name}  gate: {gate}  triggers: {trig_str}  probes: {probe_str}')
        print(f'{"="*60}')

        prb_list = parse_probe_str(probe_str)

        # ---- Per-probe processing ----
        for prb in prb_list:
            session_id = f'{run_name}_imec{prb}'
            print(f'\n{"--"*30}')
            print(f'Probe {prb}  ({session_id})')
            print(f'{"--"*30}')

            run_str      = f'{run_name}_g{gate}'
            raw_dir      = os.path.join(npx_directory, run_str, f'{run_str}_imec{prb}')
            data_dir     = os.path.join(output_dest, f'preprocessed_{run_str}', f'{run_str}_imec{prb}')
            bandpass_dir = os.path.join(data_dir, 'bandpass_only.ap')
            ks4_out_dir  = os.path.join(data_dir, f'imec{prb}_ks4')

            if not os.path.isdir(raw_dir):
                print(f'ERROR: raw SpikeGLX data not found: {raw_dir}')
                log_step(session_id, 'read_recording', 'error_no_bin')
                continue
            check_single_trigger_file(raw_dir)
            os.makedirs(data_dir, exist_ok=True)

            # ---- Step [1] Preprocessing: bandpass + phase_shift, then destripe ----
            if do_preprocess:
                print(f'\n[1] Reading raw SpikeGLX data: {raw_dir}')
                recording_raw = si.read_spikeglx(
                    raw_dir, stream_id=f'imec{prb}.ap', load_sync_channel=False
                )
                n_ch  = recording_raw.get_num_channels()
                dur_s = recording_raw.get_num_frames() / recording_raw.get_sampling_frequency()
                print(f'    {n_ch} channels, {dur_s:.1f} s @ {recording_raw.get_sampling_frequency():.0f} Hz')

                print('    bandpass_filter + phase_shift...')
                recording_bandpass = spre.phase_shift(
                    spre.bandpass_filter(recording_raw, freq_min=300, freq_max=9000)
                )
                # Saved once: source recording for cross-session matching
                # (lib/np2/matching/), which needs the spatial footprint that
                # destriping below deliberately removes.
                recording_bandpass = recording_bandpass.save(
                    folder=bandpass_dir, overwrite=True
                )

                print('\n[2] SpikeInterface preprocessing (IBL destriping)')
                # Per-shank: detect_bad_channels/highpass_spatial_filter must be
                # computed within a single shank. Shank count read from the
                # recording, so 1-/2-/4-shank NP2.0 all work.
                recording_preproc, bad_ids = ibl_destripe_by_shank(recording_bandpass)
                print(f'    Bad channels ({len(bad_ids)}): {list(bad_ids)}')
            else:
                print(f'\n[1-2] Loading existing bandpass_only recording: {bandpass_dir}')
                if not os.path.isdir(bandpass_dir):
                    print(f'ERROR: no existing bandpass_only recording found at {bandpass_dir}')
                    log_step(session_id, 'preprocess', 'error_no_output')
                    continue
                recording_bandpass = si.load(bandpass_dir)
                recording_preproc, bad_ids = ibl_destripe_by_shank(recording_bandpass)
                print(f'    Bad channels ({len(bad_ids)}): {list(bad_ids)}')

            # ---- Step 4: Kilosort4 ----
            if do_kilosort4:
                print(f'\n[4] Kilosort4')
                t0 = datetime.now()
                sorting = ss.run_sorter(
                    sorter_name='kilosort4',
                    recording=recording_preproc,
                    folder=ks4_out_dir,      # SI 0.104: 'folder' (was 'output_folder' in older SI)
                    do_CAR=False,            # destriping already handled by highpass_spatial_filter
                    nblocks=ks_nblocks,
                    Th_universal=ks_Th_universal,
                    Th_learned=ks_Th_learned,
                    remove_existing_folder=True,
                )
                elapsed_ks = (datetime.now() - t0).total_seconds()
                n_units = len(sorting.unit_ids)
                print(f'    Done in {elapsed_ks:.0f}s  |  {n_units} units found')
                log_step(session_id, 'kilosort4', 'done', n_units=n_units, elapsed=f'{elapsed_ks:.0f}')
            else:
                print(f'\n[4] Loading existing Kilosort4 output: {ks4_out_dir}')
                if not os.path.isdir(os.path.join(ks4_out_dir, 'sorter_output')):
                    print(f'ERROR: no existing Kilosort4 output found at {ks4_out_dir}')
                    log_step(session_id, 'kilosort4', 'error_no_output')
                    continue
                sorting = si.read_sorter_folder(ks4_out_dir)
                sorting.register_recording(recording_preproc)
                n_units = len(sorting.unit_ids)
                print(f'    Loaded {n_units} units')
                log_step(session_id, 'kilosort4', 'loaded_existing', n_units=n_units)

            # Copy KS4 .npy output files to ks4_out_dir root (rc2_analysis reads from there)
            copy_ks4_outputs_to_parent(ks4_out_dir, recording_preproc)
            plot_ks4_motion(ks4_out_dir, os.path.join(ks4_out_dir, 'ks4_motion.png'))

            # ---- Step 5: SortingAnalyzer ----
            print(f'\n[5] SortingAnalyzer')
            t0 = datetime.now()
            analyzer_folder = os.path.join(ks4_out_dir, 'sorting_analyzer')

            if do_sorting_analyzer:
                analyzer = si.create_sorting_analyzer(
                    sorting=sorting,
                    recording=recording_preproc,
                    format='binary_folder',
                    folder=analyzer_folder,
                    overwrite=True,
                )

                print('    random_spikes + waveforms + templates...')
                analyzer.compute('random_spikes', method='uniform', max_spikes_per_unit=500)
                analyzer.compute('waveforms', ms_before=1.5, ms_after=2.5, dtype='float32')
                analyzer.compute('templates', operators=['average', 'std'])

                print('    unit_locations + spike_amplitudes + spike_locations + noise_levels + template_similarity...')
                analyzer.compute('unit_locations', method='monopolar_triangulation')
                analyzer.compute('spike_amplitudes')
                analyzer.compute('spike_locations')   # required for the 'drift' quality metric
                analyzer.compute('noise_levels')      # required for the 'snr' quality metric
                analyzer.compute('template_similarity')

                print('    template_metrics (shape + velocity + spread)...')
                analyzer.compute('template_metrics', include_multi_channel_metrics=True)

                print('    principal_components (for PC-based quality metrics)...')
                analyzer.compute('principal_components', n_components=5, mode='by_channel_local')

                print('    quality_metrics...')
                # IMPORTANT: Bombcell (bombcell_label_units, step 6) aborts if any of
                # its required metrics are missing. Its defaults need num_spikes, snr,
                # amplitude_median, rp_contamination (from 'rp_violation') and
                # drift_ptp (from 'drift') -- these MUST be in this list, and 'snr'
                # needs noise_levels / 'drift' needs spike_locations (computed above).
                analyzer.compute('quality_metrics', metric_names=[
                    'num_spikes',
                    'firing_rate', 'presence_ratio', 'snr',
                    'isi_violation',    # -> isi_violations_ratio
                    'rp_violation',     # -> rp_contamination (Bombcell)
                    'amplitude_cutoff', 'amplitude_median',
                    'drift',            # -> drift_ptp -> max_drift, drift_std -> cumulative_drift
                    'mahalanobis',      # -> isolation_distance + l_ratio
                    'd_prime',
                    'nearest_neighbor', # -> nn_hit_rate + nn_miss_rate
                    'silhouette',       # -> silhouette (renamed to silhouette_score)
                ])

                elapsed_an = (datetime.now() - t0).total_seconds()
                print(f'    SortingAnalyzer done in {elapsed_an:.0f}s')
                log_step(session_id, 'sorting_analyzer', 'done', elapsed=f'{elapsed_an:.0f}')
            else:
                print(f'    Loading existing SortingAnalyzer: {analyzer_folder}')
                if not os.path.isdir(analyzer_folder):
                    print(f'ERROR: no existing SortingAnalyzer found at {analyzer_folder}')
                    log_step(session_id, 'sorting_analyzer', 'error_no_output')
                    continue
                analyzer = si.load_sorting_analyzer(analyzer_folder)
                elapsed_an = (datetime.now() - t0).total_seconds()
                print(f'    Loaded existing SortingAnalyzer in {elapsed_an:.0f}s')
                log_step(session_id, 'sorting_analyzer', 'loaded_existing', elapsed=f'{elapsed_an:.0f}')

            # ---- Step 6: Bombcell automated curation ----
            print('\n[6] Bombcell curation')
            t0 = datetime.now()
            labels = None
            figures_ok = False
            try:
                # bombcell_thresholds: defined in the "User input" section
                # near the top of this file (edit there to change a threshold).
                # split_non_somatic_good_mua=True: keep the good/mua distinction
                # for non-somatic (axonal/dendritic) units instead of collapsing
                # them into a single 'non_soma' label -- otherwise a unit's
                # pre-existing quality (good vs mua) is lost once it is flagged
                # as non-somatic.
                labels_df = sc.bombcell_label_units(
                    analyzer, thresholds=bombcell_thresholds, split_non_somatic_good_mua=True
                )
                # SI returns the labels in the 'bombcell_label' column
                # (good / mua / noise / non_soma_good / non_soma_mua).
                label_col = 'bombcell_label' if 'bombcell_label' in labels_df.columns else labels_df.columns[-1]
                labels = labels_df[label_col].tolist()
                label_summary = dict(zip(*np.unique(labels, return_counts=True)))
                print(f'    Labels: {label_summary}')
                # Full per-unit metrics table (quality_metrics + template_metrics,
                # every column, unfiltered/unrenamed) -- needed to inspect and
                # tune Bombcell thresholds per the Bombcell team's own guidance
                # (their default thresholds are a starting point, not fixed).
                # csv/metrics.csv only keeps a small rc2_analysis-renamed subset
                # and does not include most Bombcell threshold metrics (snr,
                # rp_contamination, num_positive_peaks, waveform_baseline_flatness,
                # ...), so it is not enough for this on its own.
                all_metrics_df = analyzer.get_metrics_extension_data()
                all_metrics_df.index.name = 'unit_id'

                # Per-unit breakdown of which threshold(s) a unit failed, for
                # the same reason -- labels_upset_*.png only shows this
                # aggregated across the population, not per individual unit.
                failed = bombcell_failed_thresholds(all_metrics_df, bombcell_thresholds)
                labels_df = labels_df.join(failed)

                bc_folder = os.path.join(ks4_out_dir, 'bombcell')
                os.makedirs(bc_folder, exist_ok=True)
                labels_df.to_csv(os.path.join(bc_folder, 'unit_labels.csv'))
                all_metrics_df.to_csv(os.path.join(bc_folder, 'all_metrics.csv'))
            except Exception as e:
                print(f'    Warning: Bombcell labelling failed ({e})')
                print('    Units will be labelled "unsorted"')

            # Summary figures (population-level view across all units -- not
            # available from bombcell_label_units alone). Labelling above
            # already succeeded if we get here with labels is not None; kept
            # in its own try/except so a plotting failure (e.g. a metric
            # missing from this SortingAnalyzer) never masks the labels that
            # were already computed and saved, and is logged/reported on its
            # own instead of silently downgrading the whole step to 'done'.
            if labels is not None:
                try:
                    print('    Saving summary figures...')
                    plot_bombcell_metric_histograms(
                        all_metrics_df, bombcell_thresholds,
                        os.path.join(bc_folder, 'metric_histograms.png'),
                    )

                    # noise_units_upset.png / mua_units_upset.png /
                    # non_somatic_units_upset.png (skips any category with 0
                    # failing units) -- see plot_bombcell_upset's docstring
                    # for why this pipeline reimplements SpikeInterface's own
                    # sw.plot_bombcell_labels_upset instead of calling it.
                    plot_bombcell_upset(analyzer, labels_df[label_col], bombcell_thresholds, bc_folder)

                    # See plot_bombcell_waveform_classification's docstring
                    # for why this pipeline reimplements SpikeInterface's own
                    # sw.plot_unit_labels / WaveformOverlayByLabelWidget
                    # instead of calling it directly.
                    plot_bombcell_waveform_classification(
                        analyzer, labels_df[label_col], os.path.join(bc_folder, 'waveform_classification.png')
                    )
                    plt.close('all')

                    figures_ok = True
                except Exception as e:
                    print(f'    Warning: Bombcell summary figures failed ({e})')
                    plt.close('all')

            elapsed_bc = (datetime.now() - t0).total_seconds()
            if labels is None:
                bc_status = 'failed'
            elif not figures_ok:
                bc_status = 'done_no_figures'
            else:
                bc_status = 'done'
            log_step(session_id, 'bombcell', bc_status, elapsed=f'{elapsed_bc:.0f}')

            # ---- Step 7: Phy export ----
            print('\n[7] Phy export')
            t0 = datetime.now()
            phy_folder = os.path.join(ks4_out_dir, 'phy')
            sexp.export_to_phy(
                analyzer,
                output_folder=phy_folder,
                remove_if_exists=True,
                copy_binary=False,
            )
            # Write Bombcell labels to Phy's cluster_group.tsv for visualization in Phy
            if labels is not None:
                with open(os.path.join(phy_folder, 'cluster_group.tsv'), 'w') as f:
                    f.write('cluster_id\tgroup\n')
                    for uid, lbl in zip(analyzer.unit_ids, labels):
                        f.write(f'{uid}\t{lbl}\n')
            elapsed_phy = (datetime.now() - t0).total_seconds()
            print(f'    Done in {elapsed_phy:.0f}s  ->  {phy_folder}')
            log_step(session_id, 'phy_export', 'done', elapsed=f'{elapsed_phy:.0f}')

            # ---- Step 8: Save rc2_analysis-compatible files ----
            print('\n[8] Saving rc2_analysis-compatible files')
            t0 = datetime.now()
            save_rc2_compatible_files(analyzer, labels, ks4_out_dir)
            elapsed_rc2 = (datetime.now() - t0).total_seconds()
            print(f'    Done in {elapsed_rc2:.0f}s')
            log_step(session_id, 'rc2_files', 'done', n_units=n_units, elapsed=f'{elapsed_rc2:.0f}')

            print(f'\nProbe {prb} complete.')
            print(f'  KS4 output  : {ks4_out_dir}')
            print(f'  Phy         : {phy_folder}')
            print(f'  CSV files   : {os.path.join(ks4_out_dir, "csv")}')

    print(f'\nAll runs complete. Log: {logFullPath}')


if __name__ == '__main__':
    main()
