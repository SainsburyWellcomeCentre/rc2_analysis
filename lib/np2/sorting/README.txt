Sorting pipeline for Neuropixels 2.0 (4-shank) SpikeGLX data
==============================================================

Per recording: CatGT + SpikeInterface (destriping) + Kilosort 4 + Bombcell.

  spikeGLX_pipeline_np2.py   -- sort one recording (CatGT + SI + KS4 + Bombcell)

This pipeline does NOT use the old ecephys_spike_sorting / Kilosort2 fork.

For environment setup and installation, see the main repo README:
  ../../../README.md

For cross-session unit tracking (DeepUnitMatch / UnitMatchPy), see:
  ../matching/README.txt


======================================================================
 RUNNING
======================================================================

Run this once per recording. Two equivalent ways:

  OPTION A -- from MATLAB (full RC2 preprocessing; recommended in this lab)
    ctl = RC2Preprocess();
    ctl.preprocess_step_1(probe_id)                       % full stage-1 pipeline
    % sorting only:
    ctl.run_sorting_from_step(probe_id)
    % re-sort already-CatGT'd data (skip CatGT):
    ctl.run_sorting_from_step(probe_id, 'run_catgt', false)

  OPTION B -- Python standalone (no MATLAB)
    conda activate spikeinterface
    cd C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\np2\sorting
    # edit the "User input" section of spikeGLX_pipeline_np2.py:
    #   npx_directory, run_specs, catGT_dest, run_CatGT, catGTPath
    python spikeGLX_pipeline_np2.py

  Inspect a result in Phy (optional):
    phy template-gui <...>/imec0_ks4/phy/params.py


======================================================================
 OUTPUT (compatible with rc2_analysis / FileManager.m)
======================================================================

catGT_dest/catgt_{run}_g{gate}/{run}_g{gate}_imec{prb}/imec{prb}_ks4/
  spike_times.npy / spike_clusters.npy / spike_templates.npy
  amplitudes.npy / templates.npy / channel_map.npy / channel_positions.npy
  params.py
  cluster_groups.csv       <- Bombcell labels (cluster_id, group)
  cluster_group.tsv        <- Bombcell labels in Phy/UnitMatch format
  sorting_analyzer/        <- SI SortingAnalyzer (reloadable)
  phy/                     <- Phy visualization files
  bombcell/                <- Bombcell unit_labels.csv
  RawWaveforms/            <- (created by the matching stage) per-unit UnitMatch waveforms
  csv/
    metrics.csv            <- quality metrics (firing_rate, isi_viol, ...)
    waveform_metrics.csv   <- template metrics (duration, PT_ratio, snr, ...)
    waveform_metrics_fix.csv


======================================================================
 WHAT CHANGED vs. THE ORIGINAL (ecephys_spike_sorting / KS2) PIPELINE
======================================================================

REMOVED (ecephys_spike_sorting modules):
  - catGT_helper        -> CatGT called directly via subprocess
  - ks4_helper          -> KS4 run via spikeinterface.sorters.run_sorter()
  - kilosort_postprocessing -> SI handles redundant spike removal
  - noise_templates     -> replaced by Bombcell curation
  - mean_waveforms      -> SI SortingAnalyzer computes waveforms + template metrics
  - quality_metrics     -> SI SortingAnalyzer computes quality metrics

ADDED:
  - detect_bad_channels + interpolate_bad_channels (SI preprocessing)
  - highpass_spatial_filter (IBL destriping) REPLACES CatGT -gbldmx
  - SortingAnalyzer: waveforms, templates, template metrics, PCA, quality metrics
  - Bombcell curation (good/mua/noise/non_soma labels)

KEPT:
  - CatGT band-pass (300-9000 Hz) + gfix
  - TPrime available (off by default)
  - Phy via export_to_phy()
  - KS4 settings: do_CAR=False, nblocks=6, Th_universal=8, Th_learned=9

WHY destriping instead of -gbldmx:
  -gbldmx does median subtraction per 32-channel ADC block (hardware-specific).
  highpass_spatial_filter (IBL destriping) is a spatial high-pass along probe
  depth that handles non-uniform stripes better (IBL paper Fig. 4).
  NOTE: destriping is used for SORTING only. UnitMatch RawWaveforms are taken
  from the CatGT bin (band-pass only) so the cross-channel spatial footprint --
  which UnitMatch matches on -- is preserved.


======================================================================
 FIXES & COMPATIBILITY NOTES  (2026-07, SpikeInterface 0.104.x)
======================================================================

These were found and fixed while bringing the pipeline up on real 2-shank NP2.0
data. Read this if a fresh install misbehaves.

1. MULTI-SHANK DESTRIPING (per-shank, not whole-probe)
   detect_bad_channels (coherence+psd) and highpass_spatial_filter are SPATIAL
   steps that must be computed WITHIN one shank. On a multi-shank NP2.0 probe,
   highpass_spatial_filter raises "The recording contains multiple groups!" and
   the coherence bad-channel detection is invalid across shanks.
   FIX: ibl_destripe_by_shank() in spikeGLX_pipeline_np2.py splits the recording
   by channel group (= shank), runs detect -> interpolate -> highpass_spatial_
   filter on each shank, then si.aggregate_channels() reassembles the probe.
   This is SpikeInterface's documented "Processing a Recording by Channel Group"
   workflow (SI has no single multi-shank call, only split_by + aggregate).
   The shank count is READ FROM THE RECORDING, so 1-/2-/4-shank NP2.0 all work
   with no hardcoded number. (NP2.0 is the only permanent probe assumption.)

2. run_sorter() KWARG: 'folder' NOT 'output_folder'
   SpikeInterface 0.104 renamed the run_sorter output argument from
   'output_folder' to 'folder'. Passing 'output_folder' makes SI treat it as a
   Kilosort parameter and abort with "Invalid parameters: ['output_folder']".
   FIX: the KS4 call uses folder=ks4_out_dir.

3. BLAS / LAPACK MUST BE CONSISTENT
   A conda MKL-numpy + too-new mkl (e.g. 2026) + pip OpenBLAS-scipy mix makes
   np.linalg.solve hard-crash (Windows fatal exception 0xc06d007f), which kills
   every SI step silently. Keep numpy+scipy from one source; verify with:
       python -c "import numpy as np; print(np.linalg.solve([[3,1],[1,2]],[9,8]))"
   See the main repo README (Step 6) for the full explanation and fix options.

4. BOMBCELL NEEDS NO SEPARATE PACKAGE
   spikeinterface.curation.bombcell_label_units is SI's own reimplementation of
   the Bombcell labelling logic (from quality_metrics + template_metrics). The
   standalone `bombcell` PyPI package is NOT imported. Do not `pip install
   bombcell` for this pipeline.

4b. BOMBCELL REQUIRES A SPECIFIC SET OF QUALITY METRICS
   bombcell_label_units aborts ("Metric(s) [...] not present in the quality
   metrics DataFrame") if its default thresholds reference metrics that were not
   computed. It needs, in particular: num_spikes, snr, amplitude_median,
   rp_contamination (from the 'rp_violation' metric) and drift_ptp (from the
   'drift' metric). Therefore the SortingAnalyzer step computes:
     - noise_levels        (extension)  -> enables 'snr'
     - spike_locations      (extension)  -> enables 'drift'
     - quality_metrics with metric_names including num_spikes, snr,
       amplitude_median, rp_violation and drift (see spikeGLX_pipeline_np2.py).
   If Bombcell still fails, every unit is labelled 'unsorted' and the MATLAB
   automated curation would then keep zero clusters -- so keep this list in sync
   with what bombcell_label_units expects.

5. DISK: a run needs ~2x the recording size, temporarily
   KS4-via-SI writes the destriped recording to a temporary binary
   (imec{prb}_ks4/sorter_output/recording.dat, ~= the raw recording size)
   before sorting. Ensure the OUTPUT drive has ~2x the recording free. The
   temp .dat is removed/overwritten on the next run.

6. MULTI-SHANK channel_map.npy MUST BE FIXED UP AFTER aggregate_channels
   Kilosort4 always writes channel_map.npy as np.arange(n_chan) -- it never
   reflects the true SpikeGLX channel_id, even on single-shank probes. On
   multi-shank NP2.0 probes this matters: ibl_destripe_by_shank() splits by
   shank and reaggregates (si.aggregate_channels), which reorders channels to
   [shank0 channels..., shank1 channels..., ...] instead of SpikeGLX's native
   interleaved order. rc2_analysis (SpikeGLXMetaData.electrode_id_from_channel_id)
   looks up shank_id from the raw imroTbl using channel_map's value as if it
   were the true channel_id, so most clusters on shank > 0 would get the wrong
   shank_id/depth if this were left uncorrected.
   FIX: copy_ks4_outputs_to_parent() overwrites channel_map.npy with the true
   SpikeGLX AP channel numbers, read from recording_preproc.channel_ids (the
   recording actually passed to Kilosort4, before SI discards that mapping).


======================================================================
 NOTES
======================================================================

GPU / CPU:
  KS4 uses the GPU if torch+CUDA is installed.
  python -c "import torch; print(torch.cuda.is_available())"
  On CPU it still runs, ~10x slower.

Reload a finished SortingAnalyzer without re-sorting:
  import spikeinterface.full as si
  ana = si.load_sorting_analyzer(r'<...>/imec0_ks4/sorting_analyzer')
  print(ana.get_extension('quality_metrics').get_data().head())
