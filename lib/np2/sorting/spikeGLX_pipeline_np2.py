#!/usr/bin/env python3
"""
SpikeInterface-based electrophysiology pipeline for NP2 probes (SpikeGLX data).
Tested with SpikeInterface 0.104.x.

Replaces ecephys_spike_sorting modules with SpikeInterface.

Preprocessing chain:
  CatGT  : bandpass filter (300-9000 Hz Butterworth) + gfix artifact removal
             NOTE: -gbldmx is intentionally removed -- replaced by highpass_spatial_filter
  SI     : detect_bad_channels -> interpolate_bad_channels -> highpass_spatial_filter
             (IBL destriping; handles non-uniform stripe noise better than median-based CAR)
             Applied PER SHANK (channel group) -- see ibl_destripe_by_shank().
             Both detect_bad_channels (coherence+psd) and highpass_spatial_filter
             are spatial and must be computed within a single shank, so on
             multi-shank NP2.0 probes (2- or 4-shank) the recording is split by
             group, each shank destriped, then reaggregated. Shank count is read
             from the data (NP2.0 is the only fixed probe assumption).

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
      metrics.csv               <- quality metrics (firing_rate, isi_viol, ...)
      waveform_metrics.csv      <- template shape metrics (duration, PT_ratio, ...)
      waveform_metrics_fix.csv  <- copy of waveform_metrics.csv
    sorting_analyzer/           <- SortingAnalyzer binary folder (reloadable)
    phy/                        <- Phy visualization output
    bombcell/                   <- Bombcell thresholds and results JSON

TPrime (optional, for behavioral data synchronization):
  Set runTPrime=True and add sync edge extraction flags to catGT_cmd_string.
  Called directly via subprocess -- no ecephys module required.

Column name mapping (SpikeInterface 0.104 -> rc2_analysis):
  Template metrics:
    peak_to_trough_duration  -> duration
    trough_half_width        -> halfwidth
    main_peak_to_trough_ratio-> PT_ratio
    velocity_above, velocity_below, spread, repolarization_slope, recovery_slope: unchanged
  Quality metrics:
    isi_violations_ratio     -> isi_viol
    silhouette               -> silhouette_score
    drift_ptp                -> max_drift
    drift_std                -> cumulative_drift
    isolation_distance, l_ratio, d_prime, nn_hit_rate, nn_miss_rate: unchanged
"""

import os
import re
import sys
import shutil
import subprocess
import fnmatch
import numpy as np
import pandas as pd
from datetime import datetime


# ============================================================
# User input -- Edit this section
# ============================================================
#
# NOTE: this file is a TEMPLATE. It is not run directly by the MATLAB
# pipeline. SortingHelper.m (overwrite_sorting_script) fills in the
# session-specific values below and writes the result to
# lib/np2/sorting/_generated/spikeGLX_pipeline_session.py, which
# is what actually gets executed for a given session. That generated
# file is overwritten on every run and is not tracked in git -- edit
# this template, not the generated copy.

# Log file name (saved in catGT_dest)
logName = 'pipeline_log.csv'

# Raw data directory (parent of the SpikeGLX run folder)
npx_directory = r'D:\data\myrecording'

# run_specs: list of [run_name, gate, trigger_string, probe_string]
#   run_name:       undecorated run name, no g/t specifier (the -run field in CatGT)
#   gate:           gate index as string, e.g. '0'
#   trigger_string: triggers to concatenate, e.g. '0,0' (single file) or 'start,end' (all)
#   probe_string:   probes to process, e.g. '0', '0,3', '0:3'
run_specs = [['myrecording', '0', '0,0', '0']]

# Output destination: all CatGT + KS4 output is written under this directory
catGT_dest = r'D:\data\myrecording\output'

# Where to resume this run from -- for debugging/re-running part of an
# already-sorted session without redoing the expensive earlier steps:
#   'catgt'        : full run from CatGT onwards (default)
#   'kilosort4'     : skip CatGT, read the existing CatGT .ap.bin, destripe,
#                     run Kilosort4 and everything after it
#   'postprocess'   : skip CatGT and Kilosort4, reload the existing KS4
#                     sorter output (re-destriping the recording -- it is
#                     not itself saved to disk, only the sort is), then run
#                     SortingAnalyzer, Bombcell, Phy export and CSV export
# 'kilosort4' and 'postprocess' both require the earlier steps' output to
# already exist under catGT_dest for this run/probe.
start_step = 'catgt'

# ---- CatGT settings ----
run_CatGT = True
catGT_stream_string = '-ap'

# -gbldmx intentionally absent: replaced by SI highpass_spatial_filter (IBL destriping)
# which handles non-uniform stripe noise across probe depth more robustly.
# gfix=0,0.10,0.02 : detect and repair electrical artifacts
#   0    = disable median-threshold for artifact detection (use default)
#   0.10 = exclusion window in seconds around each artifact
#   0.02 = correction (blanking) window in seconds
catGT_cmd_string = '-prb_fld -out_prb_fld -apfilter=butter,12,300,9000 -gfix=0,0.10,0.02 '

# ---- TPrime settings (behavioral data synchronization) ----
# Set runTPrime=True to synchronize timestamps between streams.
# IMPORTANT: also add sync edge extraction to catGT_cmd_string, e.g.:
#   -SY=0,384,6,500    (imec sync line on channel 384, bit 6, threshold 500)
#   -XA=0,1,3,500      (NI analog sync)
# and update toStream_sync_params / niStream_sync_params to match.
runTPrime = False
sync_period = 1.0                          # 1.0 for SYNC wave from imec basestation
toStream_sync_params = 'SY=0,384,6,500'   # copy from catGT_cmd_string, no spaces
niStream_sync_params = 'XA=0,1,3,500'     # set to None if no NI auxiliary data

# ---- Kilosort 4 settings ----
ks_nblocks = 6       # non-rigid drift correction blocks (0=rigid, 6=good for long probes)
ks_Th_universal = 8  # template detection threshold (KS4 default)
ks_Th_learned = 9    # learned template threshold   (KS4 default)

# ---- Tool paths (edit for your computer) ----
catGTPath  = r'C:\Users\Lab\SWC\CatGT-win'
tPrimePath = r'C:\Users\Lab\SWC\TPrime-win'

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


def get_trial_range(gate, prb_folder):
    """Scan prb_folder .bin files and return (min_trial, max_trial) indices."""
    search_str = f'_g{gate}_t'
    min_idx, max_idx = sys.maxsize, 0
    for fname in os.listdir(prb_folder):
        if fnmatch.fnmatch(fname, '*.bin'):
            g_pos = fname.find(search_str)
            if g_pos < 0:
                continue
            t_start = g_pos + len(search_str)
            t_end = fname.find('.', t_start)
            if t_end < 0:
                continue
            try:
                t_idx = int(fname[t_start:t_end])
                min_idx = min(min_idx, t_idx)
                max_idx = max(max_idx, t_idx)
            except ValueError:
                continue
    return min_idx, max_idx


def parse_trig_str(trigger_string, gate, prb_folder):
    """Parse '0,0' or 'start,end' -> (first_trig, last_trig) integers."""
    first_str, last_str = trigger_string.split(',')
    min_idx = max_idx = None
    if 'start' in first_str or 'end' in last_str:
        min_idx, max_idx = get_trial_range(gate, prb_folder)
    first_trig = min_idx if 'start' in first_str else int(first_str)
    last_trig  = max_idx if 'end'   in last_str  else int(last_str)
    return first_trig, last_trig


def parse_catgt_log(log_dir, run_name, gate_string, prb_list):
    """Read CatGT.log and return array of gfix edit rates (edits/sec) per probe."""
    gfix_str = f'{run_name}_{gate_string} Gfix'
    gfix_edits = np.zeros(len(prb_list))
    pfound, gfound = [], []
    log_path = os.path.join(log_dir, 'CatGT.log')
    try:
        with open(log_path, 'r') as f:
            for line in f:
                g = line.find(gfix_str)
                if g > -1:
                    parts = line[g:].split()
                    pfound.append(parts[3])
                    gfound.append(float(parts[5]))
    except FileNotFoundError:
        print(f'  Warning: CatGT.log not found at {log_path}')
    for i, prb in enumerate(prb_list):
        if prb in pfound:
            gfix_edits[i] = gfound[pfound.index(prb)]
    return gfix_edits


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
          metrics.csv                <- quality metrics with rc2_analysis column names
          waveform_metrics.csv       <- template metrics with rc2_analysis column names
          waveform_metrics_fix.csv   <- copy of waveform_metrics.csv

    SpikeInterface 0.104 column names -> rc2_analysis expected names:
      Template metrics:
        peak_to_trough_duration  -> duration
        trough_half_width        -> halfwidth   (half_width metric, trough value)
        main_peak_to_trough_ratio-> PT_ratio    (from waveform_ratios metric)
        velocity_above/below, spread, repolarization/recovery_slope: unchanged
      Quality metrics:
        isi_violations_ratio     -> isi_viol
        silhouette               -> silhouette_score
        drift_ptp                -> max_drift
        drift_std                -> cumulative_drift
        isolation_distance, l_ratio, d_prime, nn_hit_rate, nn_miss_rate: unchanged
    """
    csv_dir = os.path.join(ks4_output_dir, 'csv')
    os.makedirs(csv_dir, exist_ok=True)
    unit_ids = list(analyzer.unit_ids)

    # ---- metrics.csv (quality metrics) ----
    qm_ext = analyzer.get_extension('quality_metrics')
    if qm_ext is not None:
        qm_df = qm_ext.get_data().copy()
        qm_df.index.name = 'unit_id'
        qm_df = qm_df.reset_index().rename(columns={'unit_id': 'cluster_id'})

        qm_df = qm_df.rename(columns={
            'isi_violations_ratio': 'isi_viol',
            'silhouette':           'silhouette_score',
            'drift_ptp':            'max_drift',
            'drift_std':            'cumulative_drift',
        })

        required_qm = [
            'cluster_id', 'firing_rate', 'presence_ratio', 'isi_viol',
            'amplitude_cutoff', 'isolation_distance', 'l_ratio', 'd_prime',
            'nn_hit_rate', 'nn_miss_rate', 'silhouette_score',
            'max_drift', 'cumulative_drift',
        ]
        for col in required_qm:
            if col not in qm_df.columns:
                qm_df[col] = np.nan

        qm_df[required_qm].to_csv(os.path.join(csv_dir, 'metrics.csv'), index=False)
        print(f'  Saved metrics.csv ({len(qm_df)} units)')

    # ---- waveform_metrics.csv (template metrics + SNR + peak_channel + amplitude) ----
    tm_ext = analyzer.get_extension('template_metrics')
    if tm_ext is not None:
        tm_df = tm_ext.get_data().copy()
        tm_df.index.name = 'unit_id'
        tm_df = tm_df.reset_index().rename(columns={'unit_id': 'cluster_id'})

        # Rename SI 0.104 names -> rc2_analysis expected names
        tm_df = tm_df.rename(columns={
            'peak_to_trough_duration':  'duration',
            'trough_half_width':        'halfwidth',
            'main_peak_to_trough_ratio': 'PT_ratio',
        })

        # Add SNR from quality metrics
        if qm_ext is not None:
            qm_data = qm_ext.get_data().copy()
            qm_data.index.name = 'unit_id'
            qm_data = qm_data.reset_index().rename(columns={'unit_id': 'cluster_id'})
            if 'snr' in qm_data.columns:
                tm_df = tm_df.merge(qm_data[['cluster_id', 'snr']], on='cluster_id', how='left')

        # Compute peak_channel and amplitude from mean templates
        templates_ext = analyzer.get_extension('templates')
        if templates_ext is not None:
            tmpl = templates_ext.get_templates(operator='average')  # (n_units, n_samples, n_ch)
            chan_ids = analyzer.channel_ids
            peak_channels, amplitudes_list = [], []
            for i in range(tmpl.shape[0]):
                amp_per_ch = tmpl[i].max(axis=0) - tmpl[i].min(axis=0)
                best_idx = int(np.argmax(amp_per_ch))
                # channel_ids may be strings or ints; use index for simplicity
                peak_channels.append(best_idx)
                amplitudes_list.append(float(amp_per_ch[best_idx]))

            tm_df = tm_df.merge(
                pd.DataFrame({'cluster_id': unit_ids, 'peak_channel': peak_channels,
                               'amplitude': amplitudes_list}),
                on='cluster_id', how='left'
            )

        # Ensure all rc2_analysis-expected waveform columns are present
        required_wf = [
            'cluster_id', 'peak_channel', 'snr', 'duration', 'halfwidth',
            'PT_ratio', 'repolarization_slope', 'recovery_slope', 'amplitude',
            'spread', 'velocity_above', 'velocity_below',
        ]
        for col in required_wf:
            if col not in tm_df.columns:
                tm_df[col] = np.nan

        wf_path     = os.path.join(csv_dir, 'waveform_metrics.csv')
        wf_fix_path = os.path.join(csv_dir, 'waveform_metrics_fix.csv')
        tm_df[required_wf].to_csv(wf_path, index=False)
        shutil.copy2(wf_path, wf_fix_path)
        print(f'  Saved waveform_metrics.csv + waveform_metrics_fix.csv ({len(tm_df)} units)')

    # ---- cluster_groups.csv (Bombcell labels) ----
    group_values = labels if labels is not None else (['unsorted'] * len(unit_ids))
    pd.DataFrame({'cluster_id': unit_ids, 'group': group_values}).to_csv(
        os.path.join(ks4_output_dir, 'cluster_groups.csv'), index=False
    )
    print(f'  Saved cluster_groups.csv ({len(unit_ids)} units)')


def plot_bombcell_metric_histograms(metrics_df, thresholds, out_path):
    """
    Recreates the native Bombcell (MATLAB) quality_metrics_distribution.png
    layout -- 18 metrics, in the native panel order, with short human-readable
    axis labels (matching bc.qm.plotGlobalQualityMetric's defineMetrics), a
    fraction-of-units y-axis, and a red/orange/green bar under each x-axis
    showing which range of that metric is rejected/borderline/accepted.

    SpikeInterface's own sw.plot_metric_histograms uses raw SI column names
    as axis labels (e.g. 'peak_before_width' in seconds, unreadable at 1e-4
    scale) and has no colored accept/reject bar, and only plots metrics that
    have a threshold in `thresholds` -- so it never shows isolation_distance
    or l_ratio, which the native Bombcell histogram does. This function
    reads the same `thresholds` dict (from bombcell_get_default_thresholds)
    plus isolation_distance/l_ratio directly, in the same 18-panel layout as
    the native GUI, for people already familiar with it.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # (SI column, short label, unit scale factor, unit suffix) in the same
    # order as the native plotGlobalQualityMetric.m panel layout
    # (indices_ordered in defineMetrics), skipping the 2 native metrics that
    # are excluded there too (RPV_window_index, %SpikesMissing-symmetric)
    # and the 2 SI has no equivalent for (percentageSpikesMissing_gaussian
    # duplicated as amplitude_cutoff already, mainPeakToTroughRatio's own
    # 'scndPeakToTroughRatio' folded into peak_after_to_trough_ratio).
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
    panels = [
        ('num_positive_peaks',            '# peaks',            1, '', False, None, True),
        ('num_negative_peaks',             '# troughs',          1, '', False, None, True),
        ('waveform_baseline_flatness',     'baseline flatness',  1, '', False, None, False),
        ('peak_to_trough_duration',        'waveform duration',  1e6, ' µs', False, None, False),
        ('peak_after_to_trough_ratio',     'peak$_2$/trough',    1, '', False, None, False),
        ('exp_decay',                      'spatial decay',      1, '', False, None, False),
        ('peak_before_to_peak_after_ratio','peak$_1$/peak$_2$',  1, '', False, 99, False),
        ('main_peak_to_trough_ratio',      'peak$_{main}$/trough', 1, '', False, None, False),
        ('amplitude_median',               'amplitude',          1, ' µV', True, None, False),
        ('snr',                            'SNR',                1, '', False, None, False),
        ('rp_contamination',               'frac. RPVs',         1, '', False, None, False),
        ('num_spikes',                     '# spikes',           1, '', False, None, False),
        ('presence_ratio',                 'presence ratio',     1, '', False, None, False),
        ('amplitude_cutoff',               '% spikes missing',   100, ' %', False, None, False),
        ('drift_ptp',                      'maximum drift',      1, ' µm', False, None, False),
        ('drift_std',                      'cum. drift',         1, ' µm', False, None, False),
        # isolation_distance can carry a handful of near-numerically-infinite
        # outliers (division by a near-zero covariance for isolated/sparse
        # clusters) -- up to 1e15 on real data, dwarfing every other unit's
        # value. A 90th-percentile clip (rather than 99th) is needed to keep
        # the histogram readable; the outlier units themselves are unaffected
        # (still in all_metrics.csv / the actual Bombcell threshold check).
        ('isolation_distance',             'isolation dist.',    1, '', False, 90, False),
        ('l_ratio',                        'L-ratio',            1, '', False, 95, False),
    ]
    # Flatten noise/mua/non-somatic sections into one lookup, same as
    # bombcell_failed_thresholds -- greater/less bounds per SI metric name.
    flat_thresh = {}
    for section in thresholds.values():
        flat_thresh.update(section)

    n = len(panels)
    n_rows = int(np.floor(np.sqrt(n)))
    n_cols = int(np.ceil(n / n_rows))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.2 * n_cols, 2.6 * n_rows))
    axes = np.atleast_1d(axes).ravel()

    for i, (col, short_label, scale, suffix, take_abs, upper_pct, integer_valued) in enumerate(panels):
        ax = axes[i]
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
               color=plt.cm.tab20(i % 20), edgecolor='black', linewidth=0.5)

        bounds = flat_thresh.get(col, {})
        greater = bounds.get('greater', None)
        less = bounds.get('less', None)
        xmin, xmax = float(values.min()), float(values.max())
        xspan = max(xmax - xmin, 1e-12)
        pad = 0.03 * xspan
        xlo, xhi = xmin - pad, xmax + pad

        def _scaled(v):
            return v * scale if v is not None else None

        g = _scaled(greater)
        l = _scaled(less)
        # 3-segment accept/reject bar: red = rejected, green = accepted,
        # orange = the boundary case with only one side constrained.
        y0 = ax.get_ylim()
        bar_y = -0.04 * (y0[1] if y0[1] > 0 else 1)
        if g is not None and l is not None:
            ax.plot([xlo, g], [bar_y, bar_y], color='red', lw=4, solid_capstyle='butt')
            ax.plot([g, l], [bar_y, bar_y], color='green', lw=4, solid_capstyle='butt')
            ax.plot([l, xhi], [bar_y, bar_y], color='red', lw=4, solid_capstyle='butt')
        elif g is not None:
            ax.plot([xlo, g], [bar_y, bar_y], color='red', lw=4, solid_capstyle='butt')
            ax.plot([g, xhi], [bar_y, bar_y], color='green', lw=4, solid_capstyle='butt')
        elif l is not None:
            ax.plot([xlo, l], [bar_y, bar_y], color='green', lw=4, solid_capstyle='butt')
            ax.plot([l, xhi], [bar_y, bar_y], color='red', lw=4, solid_capstyle='butt')
        else:
            ax.plot([xlo, xhi], [bar_y, bar_y], color='orange', lw=4, solid_capstyle='butt')

        ax.set_xlim(xlo, xhi)
        if suffix:
            # Put the unit on the tick labels themselves (not just the axis
            # label) -- '% spikes missing' with a bare 0-1-looking axis
            # (values are genuinely ~0-1% here, not 0-100%) reads as a raw
            # fraction; '0.5 %' on each tick removes the ambiguity.
            # default arg (suffix=suffix) binds THIS iteration's value at
            # definition time -- a bare closure over the loop variable
            # `suffix` would have every panel's formatter see whatever
            # `suffix` happened to be on the LAST loop iteration instead
            # (matplotlib calls the formatter lazily, at draw time, by which
            # point the loop has already finished).
            ax.xaxis.set_major_formatter(
                plt.FuncFormatter(lambda x, _, suffix=suffix: f'{x:g}{suffix}')
            )
        if i % n_cols == 0:
            ax.set_ylabel('frac. units')
        ax.set_xlabel(short_label)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for j in range(n, len(axes)):
        axes[j].axis('off')

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
    import spikeinterface.sorters as ss
    import spikeinterface.exporters as sexp
    import spikeinterface.curation as sc
    import spikeinterface.widgets as sw

    # Parallelise SortingAnalyzer extension computation (waveforms,
    # spike_amplitudes, spike_locations, ...) and export_to_phy -- these
    # default to n_jobs=1 (single-threaded) otherwise. n_jobs is kept below
    # the machine's full core count to leave headroom for other work (e.g.
    # MATLAB, Phy) running at the same time.
    si.set_global_job_kwargs(n_jobs=12, chunk_duration='1s', progress_bar=True)

    valid_start_steps = ('catgt', 'kilosort4', 'postprocess')
    if start_step not in valid_start_steps:
        raise ValueError(f"start_step must be one of {valid_start_steps}, got {start_step!r}")
    do_catgt      = run_CatGT and start_step == 'catgt'
    do_kilosort4  = start_step in ('catgt', 'kilosort4')
    print(f'start_step = {start_step!r}  (CatGT: {do_catgt}, Kilosort4: {do_kilosort4})')

    # Clean stale log files in working directory
    for stale in ('CatGT.log', 'Tprime.log'):
        try:
            os.remove(stale)
        except OSError:
            pass

    os.makedirs(catGT_dest, exist_ok=True)
    logFullPath = os.path.join(catGT_dest, logName)
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

        # Resolve 'start'/'end' trigger wildcards
        run_folder_name = f'{run_name}_g{gate}'
        prb0_fld = os.path.join(npx_directory, run_folder_name,
                                 f'{run_folder_name}_imec{prb_list[0]}')
        first_trig, last_trig = parse_trig_str(trig_str, gate, prb0_fld)
        trig_str_resolved = f'{first_trig},{last_trig}'

        # ---- Step 1: CatGT ----
        if do_catgt:
            print('\n--- CatGT ---')
            catgt_exe = os.path.join(catGTPath, 'CatGT.exe')
            cmd = (
                f'"{catgt_exe}"'
                f' -dir="{npx_directory}"'
                f' -run={run_name}'
                f' -g={gate}'
                f' -t={trig_str_resolved}'
                f' {catGT_stream_string}'
                f' -prb={probe_str}'
                f' {catGT_cmd_string.strip()}'
                f' -dest="{catGT_dest}"'
            )
            print(cmd)
            t0 = datetime.now()
            subprocess.check_call(cmd, shell=True)
            elapsed = (datetime.now() - t0).total_seconds()
            gfix_edits = parse_catgt_log(os.getcwd(), run_name, gate, prb_list)
            for i, prb in enumerate(prb_list):
                print(f'  Probe {prb}: gfix edits/sec = {gfix_edits[i]:.3f}')
            log_step(run_name, 'CatGT', 'done', elapsed=f'{elapsed:.1f}')

            # Keep a copy of CatGT's own log (params, gfix/artifact stats) next
            # to the rest of this run's output -- same location and name as the
            # original ecephys_spike_sorting pipeline -- before it is deleted
            # from the working directory (see cleanup at the top of main()).
            # Kept per-session (rather than discarded) because it is useful
            # for debugging a specific run after the fact.
            catgt_log_src = os.path.join(os.getcwd(), 'CatGT.log')
            if os.path.isfile(catgt_log_src):
                run_str = f'{run_name}_g{gate}'
                catgt_run_dir = os.path.join(catGT_dest, f'catgt_{run_str}')
                os.makedirs(catgt_run_dir, exist_ok=True)
                shutil.copy2(
                    catgt_log_src,
                    os.path.join(catgt_run_dir, f'catgt_{run_str}_prb_{probe_str}_CatGT.log'),
                )
        else:
            print(f'Skipping CatGT (start_step={start_step!r}, run_CatGT={run_CatGT})')

        # ---- Per-probe processing ----
        for prb in prb_list:
            session_id = f'{run_name}_imec{prb}'
            print(f'\n{"--"*30}')
            print(f'Probe {prb}  ({session_id})')
            print(f'{"--"*30}')

            run_str     = f'{run_name}_g{gate}'
            data_dir    = os.path.join(catGT_dest, f'catgt_{run_str}', f'{run_str}_imec{prb}')
            catgt_bin   = os.path.join(data_dir, f'{run_str}_tcat.imec{prb}.ap.bin')
            ks4_out_dir = os.path.join(data_dir, f'imec{prb}_ks4')

            if not os.path.exists(catgt_bin):
                print(f'ERROR: CatGT output not found: {catgt_bin}')
                log_step(session_id, 'read_recording', 'error_no_bin')
                continue

            # ---- Step 2: Read CatGT output via SpikeInterface ----
            print(f'\n[2] Reading CatGT output: {data_dir}')
            recording_raw = si.read_spikeglx(
                data_dir, stream_id=f'imec{prb}.ap', load_sync_channel=False
            )
            n_ch  = recording_raw.get_num_channels()
            dur_s = recording_raw.get_num_frames() / recording_raw.get_sampling_frequency()
            print(f'    {n_ch} channels, {dur_s:.1f} s @ {recording_raw.get_sampling_frequency():.0f} Hz')

            # ---- Step 3: SpikeInterface preprocessing (IBL destriping) ----
            print('\n[3] SpikeInterface preprocessing')

            # Detect + interpolate bad channels and highpass_spatial_filter,
            # all applied PER SHANK (channel group). Both detect_bad_channels
            # (coherence+psd) and highpass_spatial_filter are spatial steps that
            # SI/IBL require to be computed within a single shank. Shank count is
            # read from the recording, so 1-/2-/4-shank NP2.0 all work.
            # Destriping replaces CatGT -gbldmx (non-uniform stripe removal).
            # Always redone even when resuming at 'postprocess': the destriped
            # recording is intentionally never written to disk (it would
            # duplicate the CatGT .ap.bin), and this step is fast (minutes)
            # compared to Kilosort4, so recomputing it on resume is cheap.
            recording_preproc, bad_ids = ibl_destripe_by_shank(recording_raw)
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

            # ---- Step 5: SortingAnalyzer ----
            print(f'\n[5] SortingAnalyzer')
            t0 = datetime.now()
            analyzer_folder = os.path.join(ks4_out_dir, 'sorting_analyzer')
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

            # ---- Step 6: Bombcell automated curation ----
            print('\n[6] Bombcell curation')
            t0 = datetime.now()
            labels = None
            figures_ok = False
            try:
                bombcell_thresholds = sc.bombcell_get_default_thresholds()
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
                    # NOTE: sw.plot_unit_labels is literally an alias for
                    # WaveformOverlayByLabelWidget (see spikeinterface's
                    # widget_list.py: plot_unit_labels = WaveformOverlayByLabelWidget),
                    # i.e. the exact same figure as waveform_classification.png
                    # below -- not called separately here to avoid saving the
                    # same plot twice under two names.
                    plot_bombcell_metric_histograms(
                        all_metrics_df, bombcell_thresholds,
                        os.path.join(bc_folder, 'metric_histograms.png'),
                    )

                    # BombcellUpsetPlotWidget builds one figure per label in
                    # unit_labels_to_plot, in that order, skipping any label
                    # with 0 units -- so figures[] and this filtered name
                    # list stay in lockstep, and each file gets an explicit,
                    # unambiguous name instead of a positional index.
                    present_labels = set(labels_df[label_col].unique())
                    # non_soma_good and non_soma_mua can both be present at
                    # once (split_non_somatic_good_mua=True), each getting
                    # its own figure -- suffix their filenames so one never
                    # overwrites the other, unlike noise/mua which are always
                    # singular per run.
                    upset_labels_to_names = [
                        ('noise', 'noise_units_upset'),
                        ('mua', 'mua_units_upset'),
                        ('non_soma_good', 'non_somatic_units_upset_good'),
                        ('non_soma_mua', 'non_somatic_units_upset_mua'),
                        ('non_soma', 'non_somatic_units_upset'),
                    ]
                    labels_to_plot = [lbl for lbl, _ in upset_labels_to_names if lbl in present_labels]
                    w = sw.plot_bombcell_labels_upset(
                        analyzer, unit_labels=labels_df[label_col], thresholds=bombcell_thresholds,
                        unit_labels_to_plot=labels_to_plot,
                    )
                    names_for_plotted = [
                        name for lbl, name in upset_labels_to_names if lbl in labels_to_plot
                    ]
                    figs = w.figures if hasattr(w, 'figures') else [w.figure]
                    for fig, out_name in zip(figs, names_for_plotted):
                        # bbox_inches='tight' so the per-label suptitle
                        # (e.g. 'noise (n=121)') is not cropped out of frame.
                        fig.savefig(os.path.join(bc_folder, f'{out_name}.png'), bbox_inches='tight')
                    plt.close('all')

                    w = sw.WaveformOverlayByLabelWidget(analyzer, labels_df[label_col].to_numpy())
                    w.figure.savefig(os.path.join(bc_folder, 'waveform_classification.png'))
                    plt.close(w.figure)

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

        # ---- Step 9: TPrime (optional, runs once per run) ----
        if runTPrime:
            print('\n--- TPrime ---')
            run_str = f'{run_name}_g{gate}'
            ref_prb = prb_list[0]
            ref_dir = os.path.join(catGT_dest, f'catgt_{run_str}', f'{run_str}_imec{ref_prb}')
            out_dir = os.path.join(catGT_dest, f'catgt_{run_str}')

            tprime_exe = os.path.join(tPrimePath, 'TPrime.exe')
            cmd = (
                f'"{tprime_exe}"'
                f' -syncperiod={sync_period}'
                f' -tostream="{ref_dir}",{toStream_sync_params}'
                f' -dest="{out_dir}"'
            )
            if niStream_sync_params:
                ni_dir = out_dir
                cmd += f' -fromstream="{ni_dir}",{niStream_sync_params}'
            print(cmd)
            subprocess.check_call(cmd, shell=True)
            log_step(run_name, 'TPrime', 'done')

    print(f'\nAll runs complete. Log: {logFullPath}')


if __name__ == '__main__':
    main()
