# rc2_analysis

Preprocessing and analysis pipeline for electrophysiological data acquired with the
rollercoaster (RC2) setup on **Neuropixels 2.0 (4-shank)** probes recorded with **SpikeGLX**.

The pipeline runs entirely on a modern, actively supported software stack:

| Stage | Tool |
|---|---|
| AP band-pass + artifact repair (`gfix`) | [CatGT](https://billkarsh.github.io/SpikeGLX/#catgt) |
| Orchestration, destriping, waveforms, quality metrics, Phy export | [SpikeInterface](https://github.com/SpikeInterface/spikeinterface) |
| Spike sorting | [Kilosort 4](https://github.com/MouseLand/Kilosort) |
| Automated curation (good / mua / noise / non-soma labels) | [Bombcell](https://github.com/Julie-Fabre/bombcell) *(via SpikeInterface's `bombcell_label_units`)* |
| Manual cluster inspection (optional) | [Phy](https://github.com/cortex-lab/phy) |
| Cross-session unit tracking (optional) | [DeepUnitMatch / UnitMatch](https://github.com/EnnyvanBeest/UnitMatch) |
| Stream synchronisation (optional, off by default) | [TPrime](https://billkarsh.github.io/SpikeGLX/#tprime) |

All spike sorting now happens in **Python** — MATLAB only drives the pipeline and does the
downstream formatting, quality control and analysis.

> [!IMPORTANT]
> ### Reference and reproducibility
>
> This codebase accompanies the publication:
> **Velez-Fort, Cossell, Porta, Clopath, Margrie (202x), *Title*, Journal.**
>
> The pipeline has evolved since the analyses in the paper. The published results were
> produced with an **earlier version of this pipeline** (a Kilosort + `ecephys_spike_sorting`
> workflow). That version is preserved in the git history — check out an earlier commit / tag
> if you need to reproduce the original figures exactly.
>
> The guide below documents the **current, maintained pipeline** (SpikeInterface + Kilosort 4 +
> Bombcell). If you have questions, please open an issue.

> [!CAUTION]
> This codebase is provided as a reference and is not actively maintained. Installation can be
> involved; follow the steps in order.

---

# Installation Guide

This guide walks you through installing everything needed to run
[rc2_analysis](https://github.com/SainsburyWellcomeCentre/rc2_analysis) on a fresh Windows
machine. Follow the steps in order.

## Requirements | Version compatibility

The table below shows the combination we tested and recommend. **The correct CUDA version is
determined by your GPU** — check your GPU (Step 2) before installing anything else.

| Component | Tested version |
|---|---|
| OS | Windows 10 |
| MATLAB | R2025b (R2021a or newer should work) |
| Python (conda env) | 3.10 |
| CUDA Toolkit | 11.8 / 12.x (match your GPU + PyTorch build) |
| SpikeInterface | 0.104.x |
| Kilosort | 4 |
| PyTorch | CUDA build matching your CUDA |
| CatGT / TPrime | latest |

> [!NOTE]
> Compared with the old pipeline, MATLAB **no longer** compiles Kilosort MEX files. You
> therefore **do not need Visual Studio, the MATLAB Kilosort repo, `mexGPUall`, or the MATLAB
> Engine for Python**. The GPU (CUDA) is now used by **Python** — Kilosort 4 and DeepUnitMatch
> run on PyTorch.

---

## Step 1 — Set up the working directory

Create a working directory and clone `rc2_analysis` into it:

```bash
mkdir C:\SWC
cd C:\SWC
git clone https://github.com/SainsburyWellcomeCentre/rc2_analysis
```

The guide uses `C:\SWC\` as a reference — adapt paths to your machine. A typical layout:

```
SWC/
├── rc2_analysis/                 # this repo (contains lib/spikeinterface/np2 pipeline)
├── original_pipeline/
│   ├── npy-matlab/               # read .npy files in MATLAB
│   ├── spikes/                   # cortex-lab spikes (driftmap plotting)
│   ├── CatGT-win/                # external SpikeGLX tool
│   └── TPrime-win/               # external SpikeGLX tool (optional)
├── UnitMatch/                    # optional: cross-session tracking (Step 9)
└── data/
    ├── raw_data/
    ├── processed_data/
    │   └── formatted_data/
    ├── figures/
    └── temp/
```

---

## Step 2 — Check your GPU and install CUDA Toolkit

Kilosort 4 and DeepUnitMatch run on the GPU via PyTorch, and need CUDA.

1. Find your GPU: right-click the desktop → Display settings → Advanced display → note the GPU.
2. Check the CUDA versions it supports on the [CUDA GPU table](https://developer.nvidia.com/cuda-gpus).
3. Install a matching CUDA Toolkit from [developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads).

You will install a PyTorch build that matches this CUDA version in Step 6. Confirm your driver
with `nvidia-smi`.

> The pipeline still runs on CPU if no GPU/CUDA is available, roughly 10× slower.

---

## Step 3 — Install MATLAB

Download and install MATLAB from [mathworks.com](https://mathworks.com). Use the version in the
Requirements table (R2025b) if possible; R2021a or newer should also work.

**Required toolboxes** (these are the only ones the analysis code actually uses):

- Signal Processing Toolbox — filtering and `bandpower` (HF power profile)
- Statistics and Machine Learning Toolbox — `prctile` and related stats

To list installed toolboxes, run `ver` in the MATLAB command window.

> Curve Fitting, Image Processing and Parallel Computing toolboxes are **no longer required** —
> they were only needed by the old MATLAB-based Kilosort, which has been replaced by Kilosort 4
> in Python.

---

## Step 4 — Clone the MATLAB helper repositories

Two small MATLAB repositories are needed so MATLAB can read Kilosort `.npy` output and plot
driftmaps. In [Git Bash](https://git-scm.com/downloads), from `C:\SWC\original_pipeline\`:

```bash
git clone https://github.com/kwikteam/npy-matlab
git clone https://github.com/cortex-lab/spikes
```

`npy-matlab` reads `.npy` files in MATLAB; `spikes` (cortex-lab) provides `ksDriftmap` /
`plotDriftmap` used to build driftmaps. No compilation is needed.

---

## Step 5 — Download CatGT and TPrime

These are pre-compiled SpikeGLX tools — just download and extract (no installation).
From the [SpikeGLX download page](https://billkarsh.github.io/SpikeGLX/):

- **CatGT** → e.g. `C:\SWC\original_pipeline\CatGT-win`  *(required)*
- **TPrime** → e.g. `C:\SWC\original_pipeline\TPrime-win`  *(optional — only if syncing a
  separate NIDAQ stream; off by default, see note in Step 8)*

---

## Step 6 — Create the Python (SpikeInterface) environment

The whole sorting pipeline runs in **one conda environment**, Python 3.10, on Windows with an
NVIDIA GPU. If you don't have conda, install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.

```bash
conda create -n spikeinterface python=3.10 -y
conda activate spikeinterface

# 1. SpikeInterface (orchestration, preprocessing, metrics, Phy export, Bombcell labels)
pip install spikeinterface[full]

# 2. Kilosort 4
pip install kilosort

# 3. PyTorch with CUDA (match cuXXX to your CUDA from Step 2; check with nvidia-smi)
pip install torch --index-url https://download.pytorch.org/whl/cu118

# 4. Phy for manual inspection (optional)
pip install phy
```

> [!IMPORTANT]
> **Keep the BLAS/LAPACK stack consistent (Windows).** `pip install spikeinterface[full]`
> pulls numpy + scipy wheels that bundle their own OpenBLAS. If a **conda** numpy (MKL) is later
> mixed in — or conda pulls a too-new `mkl` — you get a hard crash *with no Python traceback*
> on the first LAPACK call (`Windows fatal exception: code 0xc06d007f` in `numpy.linalg.solve`).
> Every filtering/sorting step hits LAPACK, so the pipeline dies silently.
>
> Simplest prevention — after installing, own BLAS from one source:
> ```bash
> pip install --force-reinstall numpy scipy
> ```
> Then **verify** (should print `[2. 3.]`, not crash):
> ```bash
> python -c "import numpy as np; print(np.linalg.solve([[3,1],[1,2]],[9,8]))"
> ```
> Installing `torch` can also perturb the numpy BLAS stack — re-run this check afterwards.

> [!NOTE]
> **Bombcell needs no separate install.** The pipeline uses
> `spikeinterface.curation.bombcell_label_units`, SpikeInterface's own reimplementation of the
> Bombcell labelling logic. It ships with `spikeinterface[full]`; do **not** `pip install bombcell`.
> It produces *labels only* (good / mua / noise / non-soma) — it does not write raw waveforms.

**Verify the environment** (each line should print a version / OK, no error):

```bash
conda activate spikeinterface
python -c "import spikeinterface; print('SI', spikeinterface.__version__)"
python -c "import kilosort; print('KS4 OK')"
python -c "import torch; print('CUDA', torch.cuda.is_available())"
python -c "from spikeinterface.curation import bombcell_label_units; print('Bombcell OK')"
```

---

## Step 7 — Point the pipeline at CatGT/TPrime and your Python

The sorting script `lib/spikeinterface/np2/spikeGLX_pipeline_np2.py` has a small **Tool paths**
section near the top. Set it to your CatGT (and, if used, TPrime) folders:

```python
catGTPath  = r'C:\SWC\original_pipeline\CatGT-win'
tPrimePath = r'C:\SWC\original_pipeline\TPrime-win'
```

MATLAB launches this script with the Python executable from your conda env — set that path in
`path_config.m` (`si_np2_python_exe`, see [Configuration](#configuration) below), e.g.
`C:\Users\<you>\miniconda3\envs\spikeinterface\python.exe`.

Everything else in the script's *User input* section (recording directory, run specs, output
destination, CatGT on/off) is filled in automatically per session when you run from MATLAB. Edit
it by hand only if you run the Python script standalone (see
[the pipeline README](lib/spikeinterface/np2/README.txt)).

---

## Step 8 — (Optional) Cross-session unit tracking: DeepUnitMatch

Only needed for chronic / multi-session recordings where you want to track the same units across
recordings. Skip this if you only sort single recordings.

```bash
conda activate spikeinterface
pip install UnitMatchPy          # matching + I/O; pulls mtscomp, mat73, scikit-learn, ...
git clone https://github.com/EnnyvanBeest/UnitMatch    C:\SWC\UnitMatch
```

The clone provides the **DeepUnitMatch** package and its pretrained model
(`UnitMatchPy/DeepUnitMatch/utils/model`), which is **trained for Neuropixels 2.0 4-shank
probes only** — the probe this lab uses. The same GPU (PyTorch) is used for matching.

See [Cross-session unit tracking](#optional-cross-session-unit-tracking) below and
[lib/spikeinterface/np2/README.txt](lib/spikeinterface/np2/README.txt) for how to run it.

---

# Configuration

Before the first run, edit `path_config.m` to point at your data and tools, run `setup_paths.m`
to load the MATLAB paths, and fill in the experiment list CSV.

## Configuration file (`path_config.m`)

Key entries:

`git_work_tree_dir`
: Directory containing the `.git` folder for rc2_analysis. Used so the commit SHA is saved
  alongside generated data/figures.

`experiment_list_csv`
: Full path to the `.csv` describing each experiment (see below).

`formatted_data_dir`
: Where the formatted `.mat` files are written (one per probe: `<formatted_data_dir>\<probe_id>.mat`).

`raw_probe_dir`
: Root of the raw probe data, with one subfolder per animal. For **Neuropixels 2.0 (4-shank)**
  SpikeGLX runs the AP binary is at:
  `<raw_probe_dir>\<animal_id>\<probe_id>_g0\<probe_id>_g0_imec0\<probe_id>_g0_t0.imec0.ap.bin`
  where `<probe_id> = <animal_id>_<session_suffix1>_<session_suffix2>_...`

`raw_camera_dir`
: Root of the camera `.avi` files (`<raw_camera_dir>\<animal_id>_<session_suffix>\camera0.avi`).

`raw_rc2_dir`
: Root of the RC2 NIDAQ `.bin` files
  (`<raw_rc2_dir>\<animal_id>\<animal_id>\<animal_id>_<session_suffix>_001.bin`).

`processed_probe_fast_dir` / `processed_probe_slow_dir`
: Where preprocessing output goes. The `fast` location should be an SSD (CatGT + Kilosort 4 read
  the large `.ap.bin` files); data can later be moved to `slow` long-term storage.

`processed_camera_fast_dir` / `processed_camera_slow_dir`
: As above, for camera data.

`figure_dir`
: Where figures are written.

`npy_matlab_dir`
: Local clone of https://github.com/kwikteam/npy-matlab (Step 4).

`spikes_dir`
: Local clone of https://github.com/cortex-lab/spikes (Step 4).

`si_np2_scripts_dir`
: Directory holding the SpikeInterface pipeline scripts (this repo's `lib/spikeinterface/np2`).

`si_np2_template`
: Path to the pipeline template `spikeGLX_pipeline_np2.py`. MATLAB copies this per session,
  injects the session paths, and runs the copy.

`si_np2_python_exe`
: Python executable of the `spikeinterface` conda env (Step 6) used to run the pipeline.

`runningmouse_python_exe` / `runningmouse_main_script`
: Python executable and entry script for the camera motion-energy processing (`runningmouse`).

## Setup

From the `rc2_analysis` directory, run `setup_paths.m` to add the required directories to the
MATLAB path:

```matlab
>> setup_paths
```

## Experiment list file

Each RC2 session is one row in the CSV at `experiment_list_csv`. Columns:

- **animal_id** — ID of the animal.
- **probe_id** — probe recording ID, `<animal_id>_<session_suffix1>_<session_suffix2>_...`
- **session_id** — RC2 session ID, `<animal_id>_<session_suffix>_001`
- **protocol** — protocol used (e.g. `mismatch_nov20`, `sparse_noise`).
- **experiment_group** — group a recording belongs to (used to select probe IDs during analysis).
- **probe_type** — `24` for Neuropixels 2.0 (4-shank).
- **git_commit** — SHA of the acquisition commit (reference only).
- **discard** — whether to discard the experiment.

---

# Preprocessing and Analysis Workflow

This section covers running the pipeline after installation: sorting, curation, quality control,
formatting and analysis.

### 0. Fill in the experiment list

Add the session details to `experiment_list_csv` (set in `path_config.m`).

### 1. Create the controller object

```matlab
>> ctl = RC2Preprocess();
```

### 2. Run stage 1 (sorting + preprocessing)

```matlab
>> ctl.preprocess_step_1(<probe_id>);
```

where `<probe_id>` is e.g. `'CAA-1115688_rec1_rec2'`. This runs seven steps in order:

1. `move_raw_to_local` — copy raw data from storage to the fast (SSD) drive.
2. `patch_meta_NP2013` — patch the SpikeGLX meta file where needed.
3. `si_sorting` — run the **SpikeInterface** pipeline: CatGT band-pass + `gfix` → IBL destriping
   (per shank) → **Kilosort 4** → SortingAnalyzer (waveforms, template & quality metrics) →
   **Bombcell** curation labels → Phy export.
4. `create_check_clusters_csv` — write the cluster table to review.
5. `create_trigger_file` — extract the probe sync/trigger channel.
6. `create_driftmap` — build the driftmap.
7. `process_camera_data` — compute camera motion energy.

To resume part-way (e.g. after fixing something), use `run_from_step`, which always runs from the
chosen step **to the end**:

```matlab
>> ctl.run_from_step(<probe_id>, 'si_sorting');
```

To run **only** the sorting step (optionally skipping CatGT, or enabling TPrime):

```matlab
>> ctl.run_sorting_from_step(<probe_id>);                          % full sort
>> ctl.run_sorting_from_step(<probe_id>, 'run_catgt', false);      % re-sort already-CatGT'd data
```

The Kilosort 4 / SpikeInterface output for each probe is written to a directory of the form:

`<processed_probe_fast_dir>\<animal_id>\output\catgt_<probe_id>_g0\<probe_id>_g0_imec0\imec0_ks4`

referred to below as `<ks4_dir>`. It contains the standard Phy/Kilosort `.npy` files plus
`cluster_groups.csv` (Bombcell labels), a `csv/` folder with `metrics.csv` and
`waveform_metrics.csv`, a reloadable `sorting_analyzer/`, and `phy/` and `bombcell/` folders.
See [lib/spikeinterface/np2/README.txt](lib/spikeinterface/np2/README.txt) for the full output
description.

> [!NOTE]
> A sorting run temporarily needs about **2×** the recording size on the output drive
> (SpikeInterface writes the destriped recording to a temporary binary before sorting).

### 3. Check the clusters

Bombcell already labels every unit (good / mua / noise / non-soma) in `cluster_groups.csv`. To
inspect clusters manually:

a. Open the results in **Phy**:

```bash
conda activate spikeinterface
cd <ks4_dir>\phy
phy template-gui params.py
```

b. Or plot cluster information in MATLAB:

```matlab
>> cluster_info = ctl.cluster_info(<probe_id>);
>> cluster_info.plot(<cluster_id>);   % <cluster_id> is an integer, not a string
```

If you curate manually, record your judgements and then create the selected-clusters file:

```matlab
>> ctl.create_selected_clusters_txt(<probe_id>);
```

This writes `selected_clusters.txt` in `<ks4_dir>`. If any reviewer discards a cluster, it is
discarded.

### 3b. Adjusting Bombcell thresholds for specific brain regions

Bombcell's default thresholds (below) were calibrated on generic cortical/hippocampal
recordings. Some regions have neurons with atypical waveform shapes that can be wrongly
flagged, so thresholds may need adjusting depending on where you recorded. All thresholds are
set in [lib/spikeinterface/np2/spikeGLX_pipeline_np2.py](lib/spikeinterface/np2/spikeGLX_pipeline_np2.py),
in the `[6] Bombcell curation` step, via the `thresholds=` argument to `bombcell_label_units`
(edit the template, not the generated per-session script — see the note at the top of the
template's "User input" section).

**Known case — cerebellum (Purkinje cells).** Purkinje cells fire two waveform types from the
same neuron: simple spikes (normal bi/triphasic shape) and complex spikes (one large initial
spike followed by several smaller "spikelets", ~600 Hz, variable in number). Complex spikes
routinely have 2+ positive peaks and will fail Bombcell's default `num_positive_peaks` noise
threshold, even though they are a real, scientifically important signal (cerebellar error
coding), not noise. Kilosort may even split simple and complex spikes from the same cell into
two clusters because the waveforms differ so much. If cerebellar clusters you inspect in Phy
show this large-spike-plus-spikelets shape and low (~1 Hz) firing rate, they are likely
mislabelled complex spikes, not noise — consider relaxing `num_positive_peaks` (see below) for
those sessions, or reviewing them manually rather than trusting the automatic "noise" label.

**All adjustable thresholds** (from `spikeinterface.curation.bombcell_get_default_thresholds()`).
`"greater"`/`"less"` are inclusive PASS bounds (a unit passes a metric if
`greater <= value <= less`, either bound optional); a unit fails the category ("noise" or
"mua") if it fails *any one* of that category's metrics:

| Category | Metric | Default | Meaning |
|---|---|---|---|
| noise | `num_positive_peaks` | `less: 2` | passes with 0-1 positive peaks; 2+ → noise |
| noise | `num_negative_peaks` | `less: 1` | passes with 0-1 negative peaks (troughs); 2+ → noise |
| noise | `peak_to_trough_duration` | `greater: 0.0001, less: 0.00115` (s) | spike duration must be in this window |
| noise | `waveform_baseline_flatness` | `less: 0.5` | baseline before/after the spike must be reasonably flat |
| noise | `peak_after_to_trough_ratio` | `less: 0.8` | rebound peak must not be too large relative to the trough |
| noise | `exp_decay` | `greater: 0.01, less: 0.1` | waveform decay time constant must be in this window |
| mua | `amplitude_median` | `greater: 30` (µV) | minimum median spike amplitude |
| mua | `snr` | `greater: 5` | minimum signal-to-noise ratio |
| mua | `amplitude_cutoff` | `less: 0.2` | estimated fraction of missed (sub-threshold) spikes must be low |
| mua | `num_spikes` | `greater: 300` | minimum spike count over the recording |
| mua | `rp_contamination` | `less: 0.1` | refractory-period violation rate must be low |
| mua | `presence_ratio` | `greater: 0.7` | unit must be present through most of the recording |
| mua | `drift_ptp` | `less: 100` (µm) | peak-to-peak drift must be limited |
| non-somatic | `peak_before_to_trough_ratio` | `less: 3` | ratio used to detect axonal/dendritic waveform shape |
| non-somatic | `peak_before_width` | `greater: 0.00015` (s) | pre-trough peak width |
| non-somatic | `trough_width` | `greater: 0.0002` (s) | trough width |
| non-somatic | `peak_before_to_peak_after_ratio` | `less: 3` | ratio used to detect axonal/dendritic waveform shape |
| non-somatic | `main_peak_to_trough_ratio` | `less: 0.8` | large positive peak relative to trough flags non-somatic origin |

Set any bound to `None` to disable that side of a check. Example, relaxing the peak-count
noise check for a cerebellar session:

```python
from spikeinterface.curation import bombcell_get_default_thresholds

thresholds = bombcell_get_default_thresholds()
thresholds['noise']['num_positive_peaks']['less'] = 4  # allow complex-spike-like waveforms

labels_df = sc.bombcell_label_units(
    analyzer, thresholds=thresholds, split_non_somatic_good_mua=True,
)
```

`split_non_somatic_good_mua=True` (used in this pipeline) additionally keeps the good/mua
distinction for non-somatic units — labels are `non_soma_good` / `non_soma_mua` instead of a
single `non_soma`, so a unit's pre-existing quality is not lost just because it was flagged as
axonal/dendritic. `RestrictClusters.m` (`curation_table`/`curation_mua_table`) surfaces this as
an `is_non_somatic` column and excludes non-somatic units from `keep` by default — hand-edit
`keep` to include a specific non-somatic unit (e.g. for axonal signal analyses).

### 4. Check the trigger

RC2 sessions that were accidentally started (then stopped) appear as extra triggers on the sync
channel but have no protocol. Remove them with:

```matlab
>> ctl.correct_trigger_file(<probe_id>);
```

This opens a small GUI where you drag across the triggers to remove and save the result.

### 5. Check LFP power profile and anatomy

Run per shank (provide the shank ID even for a single shank — it would be `0`):

```matlab
>> hf_power = ctl.hf_power(<probe_id>, <shank_id>);
>> hf_power.run();
>> ctl.save_hf_power(hf_power);
```

Output goes to `<ks4_dir>\tracks\`:

- `offset_<shank_id>.txt` — offset between anatomical and electrophysiological layer 5 (0 if no
  anatomy present).
- `hf_power_<shank_id>.pdf` — HF power depth profile (per-batch, averaged, and multi-unit
  histogram), with the electrophysiological peak marked. Region boundaries are shown if anatomy
  is present.
- `hf_power_<shank_id>.mat` — parameters and per-channel power (reloadable).

*Anatomy* — optional. If anatomy from `brainreg-segment` is available, place it in
`<ks4_dir>\tracks\track_<shank_id>.csv` and re-run the steps above. When tracing the probe track:
always start at the pia, and choose a number of sampling points equal to the insertion depth in
µm (e.g. 1750 points for a 1750 µm insertion).

*Selecting batches* — most recordings need some manual exploration:

```matlab
>> hf_power.plot_raw_batches();                    % see which batches have a clear L5 peak
>> hf_power.batches_to_use = [2, 4, 6, 8:10];      % restrict to these batches
>> hf_power.search_above = 600;                    % restrict L5 peak search range (µm above tip)
>> hf_power.search_below = 1000;
>> hf_power.interactive_plot();                    % interactively shift anatomy vs. power
```

### 6. Format the data

```matlab
>> ctl.format(<probe_id>);
```

Creates `<formatted_data_dir>\<probe_id>.mat` — a single structure with the whole recording,
including trial splitting, synchronisation to the probe, and allocation of clusters to brain
regions (from the anatomy `track_<shank_id>.csv` and any offset from step 5).

Summary CSVs are also written:

- `<formatted_data_dir>\csvs\trial_matched_offsets\<probe_id>.csv` — per replay-trial alignment
  offset (samples).
- `<formatted_data_dir>\csvs\stationary_vs_motion_fr\<probe_id>.csv` — per-cluster firing rate in
  motion vs. stationary periods.

### 7. Load the formatted data

```matlab
>> experiment = ctl.load_formatted_data(<probe_id>);
```

Returns a `FormattedData` object with access to the full recording.

> Work in progress: if you use "protocol always vis", add the `<session_id>` in
> `lib/data_cleansing.m`.

---

## Optional: cross-session unit tracking

After **all** recordings of one animal are sorted (stage 1), you can track units across
recordings with **DeepUnitMatch**. This is optional and independent — it never runs during
sorting. Point the script at one folder containing the recordings to match; it discovers every
`imec*_ks4` session under it (no path editing):

```bash
conda activate spikeinterface
cd C:\SWC\rc2_analysis\lib\spikeinterface\np2
python run_deep_unit_match.py  D:\path\to\recordings_root
```

It extracts each session's raw waveforms (from the CatGT bin, automatically), runs DeepUnitMatch,
and writes results to `<recordings_root>\unit_match_deep\`. Use `run_unit_match.py` for the
classic UnitMatchPy model instead. Full options and outputs are documented in
[lib/spikeinterface/np2/README.txt](lib/spikeinterface/np2/README.txt).

---

## Figures

Common plotting scripts under `rc2_analysis\scripts\plot\`:

1. Trial structure — `trial_structure.m`
2. Overlay aligned trials — `overlay_aligned_trials.m`
3. Rasters around motion onset / mismatch — `rasters.m`
4. Heatmaps + population average trace — `heatmaps.m`
5. Unity plots (`unity_plots`) — mismatch baseline vs. response; stationary vs. motion; motion vs.
   motion; for population and single clusters.
6. MI vs. depth plots (`unity_plots`) — same comparisons as above.
7. Tuning curves.
