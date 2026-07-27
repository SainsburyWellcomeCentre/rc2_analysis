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
# lib/spikeinterface/np2/_generated/spikeGLX_pipeline_session.py, which
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


def copy_ks4_outputs_to_parent(ks4_output_dir):
    """
    KS4 run via SpikeInterface puts its output files in sorter_output/ subfolder.
    Copy them to ks4_output_dir/ so rc2_analysis can find them directly.
    Also ensures spike_clusters.npy exists (created from spike_templates.npy if absent).

    pc_features.npy / pc_feature_ind.npy are excluded: they are large
    (~2 GB) and unused by rc2_analysis (FileManager.ks4_npy is only ever
    called with spike_clusters/spike_templates/spike_times/amplitudes/
    templates/channel_map/channel_positions). They remain available in
    sorter_output/ if ever needed.
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


# ============================================================
# Main pipeline
# ============================================================

def main():
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
    import matplotlib
    matplotlib.use('Agg')  # headless: this script has no display, only saves figures
    import matplotlib.pyplot as plt

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
        if run_CatGT:
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
            print('Skipping CatGT (run_CatGT=False)')

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
            recording_preproc, bad_ids = ibl_destripe_by_shank(recording_raw)
            print(f'    Bad channels ({len(bad_ids)}): {list(bad_ids)}')

            # ---- Step 4: Kilosort4 ----
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

            # Copy KS4 .npy output files to ks4_out_dir root (rc2_analysis reads from there)
            copy_ks4_outputs_to_parent(ks4_out_dir)

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
                # Save Bombcell results to JSON
                bc_folder = os.path.join(ks4_out_dir, 'bombcell')
                os.makedirs(bc_folder, exist_ok=True)
                labels_df.to_csv(os.path.join(bc_folder, 'unit_labels.csv'))

                # Summary figures (population-level view across all units --
                # not available from bombcell_label_units alone):
                #   - unit_labels: units on the probe, coloured by label
                #   - metric_histograms: distribution of each metric with thresholds
                #   - labels_upset: which metric(s) caused each noise/mua rejection
                print('    Saving summary figures...')
                w = sw.plot_unit_labels(analyzer, labels_df[label_col])
                w.figure.suptitle('Bombcell labels')
                w.figure.savefig(os.path.join(bc_folder, 'unit_labels.png'))
                plt.close(w.figure)

                w = sw.plot_metric_histograms(analyzer, bombcell_thresholds, figsize=(15, 10))
                w.figure.savefig(os.path.join(bc_folder, 'metric_histograms.png'))
                plt.close(w.figure)

                w = sw.plot_bombcell_labels_upset(
                    analyzer, unit_labels=labels_df[label_col], thresholds=bombcell_thresholds,
                    unit_labels_to_plot=['noise', 'mua'],
                )
                for i, fig in enumerate(plt.get_fignums()):
                    plt.figure(fig).savefig(os.path.join(bc_folder, f'labels_upset_{i}.png'))
                plt.close('all')
            except Exception as e:
                print(f'    Warning: Bombcell failed ({e})')
                print('    Units will be labelled "unsorted"')
            elapsed_bc = (datetime.now() - t0).total_seconds()
            log_step(session_id, 'bombcell',
                     'done' if labels is not None else 'failed', elapsed=f'{elapsed_bc:.0f}')

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
