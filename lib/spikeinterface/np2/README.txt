SpikeInterface-based pipeline for Neuropixels 2.0 (4-shank) SpikeGLX data
========================================================================

Two stages:
  1. SORTING        : per recording -> CatGT + SpikeInterface + Kilosort4 + Bombcell
  2. UNIT TRACKING  : across recordings of one animal -> DeepUnitMatch (optional)

Scripts in this directory:

  spikeGLX_pipeline_np2.py   -- STAGE 1: sort one recording (CatGT + SI + KS4 + Bombcell)
  extract_raw_waveforms.py   -- STAGE 2 prep: write RawWaveforms/ for (Deep)UnitMatch
  run_deep_unit_match.py     -- STAGE 2: cross-session tracking with DeepUnitMatch (recommended)
  run_unit_match.py          -- STAGE 2: cross-session tracking with classic UnitMatchPy (fallback)

This pipeline does NOT use the old ecephys_spike_sorting / Kilosort2 fork.


======================================================================
 ENVIRONMENT -- what you need, and which tool needs it
======================================================================

Everything runs in ONE conda environment (named `spikeinterface` here),
Python 3.10, on Windows with an NVIDIA GPU. Below is each tool, what it is
used for, and how it is installed.

  TOOL              USED FOR                         INSTALLED AS
  ----------------- -------------------------------- ----------------------------
  conda (Miniconda) the python environment           installer (see step 1)
  SpikeInterface    orchestration, preprocessing     pip install spikeinterface[full]
                    (destriping), SortingAnalyzer,
                    quality metrics, Phy export
  Kilosort 4        spike sorting                    pip install kilosort
  PyTorch (CUDA)    GPU for KS4 AND DeepUnitMatch     pip install torch (CUDA build)
  Bombcell          automated curation labels        pip install bombcell
                    (good/mua/noise) -- LABEL ONLY
  CatGT             band-pass + gfix on the AP data   external .exe (download)
  TPrime            stream sync (OFF by default;      external .exe (download)
                    trigger is read from the SY
                    channel instead -- see note)
  Phy               manual cluster inspection (opt.)  pip install phy / conda
  UnitMatchPy       cross-session matching + I/O      pip install UnitMatchPy
                    (also pulls mtscomp, mat73)
  DeepUnitMatch     deep-NN matching + trained model  git clone of UnitMatch repo
                    (model is NP2.0 4-shank only)

Notes:
  * Bombcell runs via spikeinterface.curation.bombcell_label_units = LABELS ONLY.
    It does NOT write RawWaveforms. Those are produced by extract_raw_waveforms.py
    (called automatically by the matchers).
  * UnitMatchPy + DeepUnitMatch are imported from the cloned UnitMatch repo (so we
    get the DeepUnitMatch package + pretrained model), but `pip install UnitMatchPy`
    is still the easy way to pull the runtime deps (mtscomp, mat73, scikit-learn...).
  * The same GPU is used by Kilosort4 and by DeepUnitMatch.


======================================================================
 INSTALLATION ON A NEW COMPUTER (do once)
======================================================================

1. INSTALL CONDA (Miniconda or Anaconda)
   https://docs.conda.io/en/latest/miniconda.html

2. CREATE + ACTIVATE THE ENVIRONMENT
   conda create -n spikeinterface python=3.10 -y
   conda activate spikeinterface

3. SPIKEINTERFACE (sorting orchestration + preprocessing + metrics)
   pip install spikeinterface[full]

   >> CRITICAL -- keep the BLAS/LAPACK stack CONSISTENT (Windows) <<
   `pip install spikeinterface[full]` pulls numpy + scipy as pip wheels
   (which bundle their own OpenBLAS). If numpy is LATER replaced by a conda
   numpy (MKL-backed) -- or conda pulls a newer `mkl` than that numpy was
   built for -- you get a MIXED, INCOMPATIBLE linear-algebra stack. The
   symptom is a HARD CRASH (no Python traceback) on the first LAPACK call:

       Windows fatal exception: code 0xc06d007f   (in numpy.linalg.solve)

   Everything that filters/sorts (detect_bad_channels, highpass_spatial_filter,
   Kilosort4, quality metrics, UnitMatch) hits LAPACK, so the whole pipeline
   dies silently. `matmul` still works (that's BLAS, not LAPACK), which makes
   it look like numpy is fine -- it is not.

   PREVENTION (do ONE, and do not mix numpy from conda with scipy from pip):
     - EASIEST: after step 3, install numpy+scipy from ONE source and let it
       own BLAS. All-pip (OpenBLAS, self-contained) is simplest:
           pip install --force-reinstall numpy scipy
       or all-conda (MKL):
           conda install -n spikeinterface numpy scipy    # then DON'T pip-reinstall them
     - If you use conda MKL numpy, PIN mkl to what numpy was built against
       (a too-new mkl is the usual culprit):
           conda install -n spikeinterface "mkl<2025"
     - Or force OpenBLAS instead of MKL for the whole env:
           conda install -n spikeinterface "libblas=*=*openblas" "libcblas=*=*openblas" "liblapack=*=*openblas"

   ALWAYS verify LAPACK before trusting the env (see VERIFY block below):
       python -c "import numpy as np; print(np.linalg.solve([[3,1],[1,2]],[9,8]))"
   If that line crashes instead of printing [2. 3.], STOP and fix BLAS first.

4. KILOSORT 4 (+ GPU PyTorch -- strongly recommended)
   pip install kilosort
   pip install torch --index-url https://download.pytorch.org/whl/cu118
   (match cu118 to your CUDA; check with: nvidia-smi)
   NOTE: torch ships its own OpenMP/BLAS; installing it can perturb the
   numpy BLAS stack, so re-run the LAPACK verify line above afterwards.

5. BOMBCELL (curation labels) -- NO SEPARATE INSTALL NEEDED
   The pipeline uses spikeinterface.curation.bombcell_label_units, which is
   SpikeInterface's OWN reimplementation of the Bombcell labelling logic
   (good / mua / noise / non_soma). It derives labels purely from the
   SortingAnalyzer's quality_metrics + template_metrics -- it does NOT import
   the standalone `bombcell` PyPI package. So `pip install bombcell` is NOT
   required and is intentionally omitted. (It ships with spikeinterface[full].)

6. CATGT (external tool: band-pass + gfix)
   Download CatGT-win: https://billkarsh.github.io/SpikeGLX/#catgt
   Extract to e.g. C:\Users\Lab\SWC\CatGT-win
   Set catGTPath in spikeGLX_pipeline_np2.py to that folder.

7. TPRIME (external tool, OPTIONAL -- only if syncing a separate NI stream)
   Download TPrime-win: https://billkarsh.github.io/SpikeGLX/#tprime
   Extract to e.g. C:\Users\Lab\SWC\TPrime-win ; set tPrimePath.
   Left OFF by default (runTPrime=False): the stimulus trigger is recorded on
   the probe's own SY (sync) channel -- the 385th channel -- so events already
   share the spike timebase and no cross-stream alignment is needed.

8. PHY (OPTIONAL, manual inspection)
   pip install phy

9. (DEEP)UNITMATCH (OPTIONAL -- only for chronic, multi-session tracking)
   pip install UnitMatchPy            (matching + I/O; pulls mtscomp, mat73, ...)
   pip install torch                  (only if not already installed in step 4)
   git clone https://github.com/EnnyvanBeest/UnitMatch    (e.g. into C:\Users\Lab\SWC\UnitMatch)

   The clone provides the DeepUnitMatch package and the pretrained model
   (UnitMatchPy/DeepUnitMatch/utils/model, trained for Npix 2.0 4-shank).
   The default clone location is set in the scripts as DEFAULT_UNITMATCH_REPO;
   override per run with --unitmatch-repo.

   If you did NOT `pip install UnitMatchPy`, install its deps directly:
     pip install mtscomp mat73 scikit-learn joblib h5py

VERIFY THE INSTALL (each line should print OK / a version, no error):
   conda activate spikeinterface
   python -c "import spikeinterface; print('SI', spikeinterface.__version__)"
   python -c "import kilosort; print('KS4 OK')"
   python -c "import torch; print('CUDA', torch.cuda.is_available())"
   python -c "from spikeinterface.curation import bombcell_label_units; print('Bombcell OK')"
   # DeepUnitMatch (set the path to your clone):
   python -c "import sys; sys.path.insert(0, r'C:\Users\Lab\SWC\UnitMatch\UnitMatchPy'); \
              from DeepUnitMatch.testing import test; test.load_trained_model(device='cpu'); print('DeepUnitMatch OK')"


======================================================================
 STAGE 1 -- SORT EACH RECORDING
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
    cd C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\spikeinterface\np2
    # edit the "User input" section of spikeGLX_pipeline_np2.py:
    #   npx_directory, run_specs, catGT_dest, run_CatGT, catGTPath
    python spikeGLX_pipeline_np2.py

  Inspect a result in Phy (optional):
    phy template-gui <...>/imec0_ks4/phy/params.py


======================================================================
 STAGE 2 -- TRACK UNITS ACROSS RECORDINGS (DeepUnitMatch, optional)
======================================================================

Do this only after ALL recordings of one animal are sorted (Stage 1).
This step is OPTIONAL and INDEPENDENT -- it never runs during sorting.

You do NOT edit any paths. Point the script at ONE folder that contains all
the recordings you want to match; it finds every imec*_ks4 session under it:

    conda activate spikeinterface
    cd C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\spikeinterface\np2
    python run_deep_unit_match.py  D:\path\to\recordings_root

What happens:
  1. discovers every sorted session (imec*_ks4) under the folder and PRINTS the
     list in match order (sorted by path = chronological if folders are named by
     date/sequence -- always check the printed list);
  2. extracts each session's RawWaveforms/ from its CatGT bin + KS4 sorting
     (automatic; skipped if already present);
  3. runs DeepUnitMatch (model -> similarity -> drift-correct + Bayes -> matches);
  4. saves results to <recordings_root>\unit_match_deep\.

Options:
    --save-dir DIR        output folder (default <recordings_root>\unit_match_deep)
    --threshold 0.5       match-probability threshold
    --dist-thresh 50      max drift-corrected centroid distance (um)
    --include-mua         match mua units too (default: good units only)
    --re-extract          force-rebuild RawWaveforms/
    --unitmatch-repo DIR  clone of EnnyvanBeest/UnitMatch (default already set)

Classic UnitMatchPy instead of the deep model (same interface):
    python run_unit_match.py  D:\path\to\recordings_root

Rebuild RawWaveforms only (rarely needed; the matchers do it automatically):
    python extract_raw_waveforms.py  D:\path\to\recordings_root

Outputs (in the save dir):
    MatchTable.csv          unit pairs with match probability
    UniqueIDConversion.*    cluster IDs with cross-session unique IDs
    MatchingOverview.png    similarity / probability / final-match matrices

Requirements recap for this stage: each session must already have KS4 output and
Bombcell labels (cluster_group.tsv); RawWaveforms are generated for you. The
trigger / sorting itself is NOT re-run -- matching uses the spike times and
clusters from Stage 1 plus the CatGT voltage.


======================================================================
 OUTPUT OF STAGE 1 (compatible with rc2_analysis / FileManager.m)
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
  RawWaveforms/            <- (created by Stage 2) per-unit UnitMatch waveforms
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
  - DeepUnitMatch cross-session unit tracking (this directory)

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

3. BLAS / LAPACK MUST BE CONSISTENT (see INSTALL step 3)
   A conda MKL-numpy + too-new mkl (e.g. 2026) + pip OpenBLAS-scipy mix makes
   np.linalg.solve hard-crash (Windows fatal exception 0xc06d007f), which kills
   every SI step silently. Keep numpy+scipy from one source; verify with the
   np.linalg.solve one-liner. Full explanation and fix in INSTALL step 3.

4. BOMBCELL NEEDS NO SEPARATE PACKAGE
   spikeinterface.curation.bombcell_label_units is SI's own reimplementation of
   the Bombcell labelling logic (from quality_metrics + template_metrics). The
   standalone `bombcell` PyPI package is NOT imported. Do not `pip install
   bombcell` for this pipeline. (See INSTALL step 5.)

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


======================================================================
 NOTES
======================================================================

GPU / CPU:
  KS4 and DeepUnitMatch use the GPU if torch+CUDA is installed.
  python -c "import torch; print(torch.cuda.is_available())"
  On CPU they still run, ~10x slower.

Reload a finished SortingAnalyzer without re-sorting:
  import spikeinterface.full as si
  ana = si.load_sorting_analyzer(r'<...>/imec0_ks4/sorting_analyzer')
  print(ana.get_extension('quality_metrics').get_data().head())
