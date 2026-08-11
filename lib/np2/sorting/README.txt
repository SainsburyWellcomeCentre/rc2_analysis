Sorting pipeline for Neuropixels 2.0 (4-shank) SpikeGLX data
==============================================================

Per recording: SpikeInterface end to end (bandpass+phase_shift, destriping,
Kilosort 4, Bombcell) -- no external tools.

  spikeGLX_pipeline_np2.py   -- sort one recording

This pipeline does NOT use the old ecephys_spike_sorting / Kilosort2 fork,
nor any external SpikeGLX post-processing tool (see "WHAT CHANGED" below for why).

For environment setup and installation, see the main repo README:
  ../../../README.md

For cross-session unit tracking (UnitMatchPy), see:
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
    % resume from a later step (see 'help RC2Preprocess.run_sorting_from_step'):
    ctl.run_sorting_from_step(probe_id, 'start_step', 'kilosort4')

  OPTION B -- Python standalone (no MATLAB)
    conda activate spikeinterface
    cd C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\np2\sorting
    # edit the "User input" section of spikeGLX_pipeline_np2.py:
    #   npx_directory, run_specs, output_dest
    python spikeGLX_pipeline_np2.py

  Inspect a result in Phy (optional):
    phy template-gui <...>/imec0_ks4/phy/params.py


======================================================================
 OUTPUT (compatible with rc2_analysis / FileManager.m)
======================================================================

output_dest/preprocessed_{run}_g{gate}/{run}_g{gate}_imec{prb}/
  bandpass_only.ap/        <- bandpass + phase_shift, NOT destriped; source
                               recording for cross-session matching
  imec{prb}_ks4/
    spike_times.npy / spike_clusters.npy / spike_templates.npy
    amplitudes.npy / templates.npy / channel_map.npy / channel_positions.npy
    params.py
    cluster_groups.csv       <- Bombcell labels (cluster_id, group)
    cluster_group.tsv        <- Bombcell labels in Phy/UnitMatch format
    sorting_analyzer/        <- SI SortingAnalyzer (reloadable)
    phy/                     <- Phy visualization files
    bombcell/                <- Bombcell unit_labels.csv
    ks4_motion.png           <- Kilosort4's own drift-correction estimate
                                 (see "KILOSORT4 DRIFT CORRECTION" below)
    csv/
      metrics.csv            <- every SortingAnalyzer quality_metrics column,
                                 SpikeInterface's own names (num_spikes,
                                 firing_rate, isi_violations_ratio, ...) --
                                 NOT the old ecephys_spike_sorting/Kilosort2
                                 names (see "WHAT CHANGED" below)
      waveform_metrics.csv   <- every SortingAnalyzer template_metrics
                                 column, SpikeInterface's own names
                                 (peak_to_trough_duration, ...), plus
                                 peak_channel/amplitude


======================================================================
 WHAT CHANGED vs. THE ORIGINAL (ecephys_spike_sorting / KS2) PIPELINE
======================================================================

REMOVED (ecephys_spike_sorting modules):
  - ks4_helper          -> KS4 run via spikeinterface.sorters.run_sorter()
  - kilosort_postprocessing -> Kilosort4 does its own within-cluster duplicate
                          spike removal + template merging internally (see
                          kilosort.postprocessing/template_matching)
  - noise_templates     -> replaced by Bombcell curation
  - mean_waveforms      -> SI SortingAnalyzer computes waveforms + template metrics
  - quality_metrics     -> SI SortingAnalyzer computes quality metrics

CHANGED (2026-08, metrics.csv/waveform_metrics.csv column names):
  These CSVs used to rename SpikeInterface's quality_metrics/template_metrics
  columns to match the old ecephys_spike_sorting/Kilosort2 pipeline's names
  (e.g. isi_violations_ratio -> isi_viol, drift_ptp -> max_drift, silhouette
  -> silhouette_score, peak_to_trough_duration -> duration). Several of those
  pairs are NOT the same calculation -- e.g. silhouette_score was a full
  pairwise silhouette, SI's default 'silhouette' is the simplified/centroid-
  based one (see SpikeInterface's silhouette_score.rst); drift_ptp is based
  on actual spike-location estimates, the old max_drift on PC1-weighted
  channel depth. Reusing the old name implied an equivalence that often
  doesn't hold, so these CSVs now use SpikeInterface's own column names, and
  every quality_metrics/template_metrics column is written (not just the
  subset the old pipeline had) -- including the metrics Bombcell's own
  thresholds are actually computed from (amplitude_median, snr,
  rp_contamination, num_positive_peaks, waveform_baseline_flatness, ...),
  previously absent from metrics.csv entirely. RC2Format.m (rc2_analysis)
  was updated to match -- see its `qm`/`waveform_metrics` variable lists.

REMOVED (2026-08, no external SpikeGLX post-processing tool is used anymore):
  - the band-pass + tshift step -> bandpass_filter + phase_shift (SpikeInterface).
                          No replacement for the old transient-artifact repair
                          step -- see main README.
  - stream synchronisation -> this lab's RC2 sync pulse is wired into the
                          probe's own AP stream (same clock as the neural
                          data), so there was never any cross-clock drift to
                          correct. Revisit only if a future experiment records
                          multiple probes at once or a separate NI-DAQ stream.

ADDED:
  - detect_bad_channels + interpolate_bad_channels (SI preprocessing)
  - highpass_spatial_filter (IBL destriping) REPLACES the old per-32-channel
    median-subtraction spatial step
  - SortingAnalyzer: waveforms, templates, template metrics, PCA, quality metrics
  - Bombcell curation (good/mua/noise/non_soma labels)

KEPT:
  - Phy via export_to_phy()
  - KS4 settings: do_CAR=False, nblocks=5, Th_universal=9, Th_learned=8

WHY destriping instead of per-32-channel median subtraction:
  Per-channel-block median subtraction is hardware-specific (32-channel ADC
  blocks). highpass_spatial_filter (IBL destriping) is a spatial high-pass
  along probe depth that handles non-uniform stripes better (IBL paper Fig. 4).
  NOTE: destriping is used for SORTING only. Cross-session matching (UnitMatch)
  reads the bandpass_only recording instead, so the cross-channel spatial
  footprint it matches on is preserved.


======================================================================
 KILOSORT4 DRIFT CORRECTION AND CROSS-SESSION MATCHING
======================================================================

Kilosort4 corrects for probe drift WITHIN a session (nblocks/do_correction,
see the "Kilosort 4 settings" comment in spikeGLX_pipeline_np2.py). UnitMatch
(cross-session matching, ../matching/) does its own, separate drift
correction BETWEEN sessions (e.g. the probe having shifted a few um
overnight). The two are not redundant -- per UnitMatch's own author (Chris
Halcrow, personal communication, 2026-08): "Unit Match does its own drift
correction, which is meant to find an overall shift between two sessions...
this is a separate thing."

Whether KS4's within-session correction needs to be turned off (nblocks=0)
before matching is NOT something this pipeline decides automatically -- per
the same source: "check and see how much drift you do have. If it's small,
try running with drift turned off and have a look at the results... there
probably isn't too much difference." This is a per-session, empirical call
for the researcher to make when preparing a matching analysis, not a
pipeline default.

To check: every sorting run saves ks4_motion.png (imec{prb}_ks4/), plotting
Kilosort4's own drift estimate over time, one curve per non-rigid block
(ops.npy: dshift, yblk) -- not the same as the driftmap PDF from
create_driftmap (RC2Preprocess.m), which scatters raw spike positions
instead of KS4's correction estimate. Look at ks4_motion.png for the
sessions you plan to match; if the drift is small, no action needed.

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
 ADJUSTING BOMBCELL THRESHOLDS
======================================================================

bombcell_thresholds (top of spikeGLX_pipeline_np2.py) starts from Bombcell's
own defaults. Bombcell's own guidance is to tune thresholds per-dataset
rather than trust the defaults, in particular amplitude_median and snr
("proxy metrics... we recommend tuning for your specific dataset" -- see
Bombcell wiki, "Detailed overview of quality metrics").

On this lab's data (checked on 3 animals: CAA-1124370, CAA-1124371,
CAA-1123244; protocols passive_protocol_motion_clouds / sparse_noise), the
empirical picture did not match that general guidance:

  Of 857 "mua" units across the 3 animals, checking which single threshold
  each one failed (bombcell/unit_labels.csv, 'failed_thresholds' column):
    rp_contamination alone: 267 units (would become "good" if only this
                             threshold were relaxed)
    snr alone:                 9 units
    amplitude_median alone:    0 units -- never the sole blocker on any
                             of the 857 units checked, despite being
                             Bombcell's top general recommendation to tune

  rp_contamination is capped at 1.0 by construction (SpikeInterface's
  Llobet-based formula returns a complex number past a certain violation
  rate and clips it to 1.0 -- see SpikeInterface docs, "Inter-spike-interval
  (ISI) violations"). 69-81% of units across the 3 animals sat exactly at
  this cap. This is NOT a low-spike-count artefact: capped units had a
  HIGHER median spike count than non-capped ones (Mann-Whitney p=4.6e-6) --
  consistent with genuine, correlated contamination (e.g. two real units
  insufficiently separated by Kilosort) rather than a small-sample fluke.

  Note: SpikeInterface's rp_contamination uses the Llobet et al. 2022
  formula. Bombcell's own default (hillOrLlobetMethod=True) uses the
  Hill et al. 2011 formula instead -- SpikeInterface's bombcell_label_units
  does not expose a way to switch to Hill for this particular threshold
  (it does compute the Hill-based number under a different name,
  isi_violations_ratio, but bombcell_label_units's rp_contamination
  threshold is wired to the Llobet-based column). The two are highly
  correlated (Spearman r=0.76 on CAA-1124370) but disagree on ~3% of units,
  so this is a real source of divergence from Bombcell's documented
  defaults, not just a relabelling.

Given the above, snr and rp_contamination -- not amplitude_median -- were
identified as the practical levers for this lab's data. Candidate relaxed
thresholds (rp_contamination 0.1->0.15, snr 5->4.5) were checked for their
effect before adopting them:

  Units that would flip from "mua" to "good" (failing only some subset of
  {rp_contamination, snr}, and passing both at the relaxed thresholds):
    CAA-1124370: 45 -> 56 good (+11)
    CAA-1124371: 37 -> 40 good (+3)
    CAA-1123244: 31 -> 37 good (+6)
    Total:      113 -> 133 good (+20, +18%)

  A sample of the newly-passing units was spot-checked in Phy (ISI view,
  waveform view) before adopting these values -- e.g. on CAA-1124370:
  cluster_id 100, 45, 270, 422, 17, 447, 171, 221, 77, 200, 212.

Current values (bombcell_thresholds in spikeGLX_pipeline_np2.py):
  rp_contamination: less=0.15   (was 0.1)
  snr:               greater=4.5 (was 5)
  all other thresholds: unchanged Bombcell/SpikeInterface defaults

This is a per-lab empirical choice, not a universal recommendation --
re-derive it (same method: bombcell/unit_labels.csv's 'failed_thresholds'
column, cross-checked in Phy) if the data profile changes substantially
(different brain region, probe, or protocol).


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
