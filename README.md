

Cohérence du pipeline — vue d'ensemble
Côté scientifique, la chaîne se tient et chaque choix a une justification vérifiée (pas supposée) :

Filtrage/déphasage/destriping : équivalents directs et sans perte à ce que faisait l'ancienne pipeline, alignés sur la méthode IBL de référence
gfix sacrifié sciemment, avec un raisonnement mesuré (taux d'artefacts variable mais faible en moyenne sur les runs réels, Bombcell rattrape les pires cas) plutôt qu'un oubli
Le point le plus solide scientifiquement : la bascule vers un SortingAnalyzer sparse pour le matching n'est pas un compromis de facilité — le papier UnitMatch lui-même confirme que ça n'a jamais eu besoin d'être dense (rayon de matching ≤150µm), donc c'est une correction d'une fausse prémisse de l'ancien script, pas juste un alignement sur l'API
Côté développement logiciel, l'architecture reste cohérente avec les patterns déjà en place :

RC2Preprocess continue d'orchestrer, SortingHelper continue de faire le pont MATLAB→Python, le script Python reste un template rempli par regex — rien de nouveau architecturalement, juste moins d'étapes
Le point le plus fragile qui subsiste : le nom de dossier de sortie est dupliqué littéralement à deux endroits indépendants (spikeGLX_pipeline_np2.py et FileManager.generate_imec0_ks4), sans source commune — exactement le même risque de désynchronisation qui existait déjà avec catgt_... avant. Rien ne garantit qu'un futur renommage des deux côtés reste synchronisé, à part la discipline de grep.
Ce qui reste non vérifié en conditions réelles (le point le plus important avant de committer) : aucun run complet n'a encore été exécuté avec ce nouveau code. Le point le plus à risque est le nouveau bloc bandpass_filter → phase_shift → save() → reload → destripe — jamais testé pour de vrai, contrairement au reste qui reprend des patterns déjà validés dans la session (le .save()/si.load() a été vérifié API par API, mais pas exécuté sur un vrai run).

Recommandation concrète avant de committer : lancer ctl.run_sorting_from_step sur une session réelle (par exemple en pointant npx_directory vers CAA-1124370 brut) pour confirmer que le bloc de preprocessing tourne de bout en bout, avant de considérer que le refactor est fonctionnellement validé — jusqu'ici tout a été vérifié par lecture de code et compilation syntaxique, pas par exécution.

# rc2_analysis

Preprocessing and analysis pipeline for electrophysiological data acquired with the
rollercoaster (RC2) setup on **Neuropixels 2.0 (4-shank)** probes recorded with **SpikeGLX**.

The pipeline runs entirely on a modern, actively supported software stack:

| Stage | Tool |
|---|---|
| AP band-pass + phase shift, destriping, orchestration, waveforms, quality metrics, Phy export | [SpikeInterface](https://github.com/SpikeInterface/spikeinterface) |
| Spike sorting | [Kilosort 4](https://github.com/MouseLand/Kilosort) |
| Automated curation (good / mua / noise / non-soma labels) | [Bombcell](https://github.com/Julie-Fabre/bombcell) *(via SpikeInterface's `bombcell_label_units`)* |
| Manual cluster inspection (optional) | [Phy](https://github.com/cortex-lab/phy) |
| Cross-session unit tracking (optional) | [UnitMatch](https://github.com/EnnyvanBeest/UnitMatch) *(via SpikeInterface's `SortingAnalyzer` integration)* |

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

> [!NOTE]
> Compared with the old pipeline, MATLAB **no longer** compiles Kilosort MEX files. You
> therefore **do not need Visual Studio, the MATLAB Kilosort repo, `mexGPUall`, or the MATLAB
> Engine for Python**. The GPU (CUDA) is now used by **Python** — Kilosort 4 runs on PyTorch.

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
├── rc2_analysis/                 # this repo (contains lib/np2 pipeline: sorting/ + matching/)
├── original_pipeline/
│   ├── npy-matlab/               # read .npy files in MATLAB
│   └── spikes/                   # cortex-lab spikes (driftmap plotting)
├── UnitMatch/                    # optional: cross-session tracking (Step 7)
└── data/
    ├── raw_data/
    ├── processed_data/
    │   └── formatted_data/
    ├── figures/
    └── temp/
```

---

## Step 2 — Check your GPU and install CUDA Toolkit

Kilosort 4 runs on the GPU via PyTorch, and needs CUDA.

1. Find your GPU: right-click the desktop → Display settings → Advanced display → note the GPU.
2. Check the CUDA versions it supports on the [CUDA GPU table](https://developer.nvidia.com/cuda-gpus).
3. Install a matching CUDA Toolkit from [developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads).

You will install a PyTorch build that matches this CUDA version in Step 5. Confirm your driver
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

## Step 5 — Create the Python (SpikeInterface) environment

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

## Step 6 — Point the pipeline at your Python executable

MATLAB launches the sorting script with the Python executable from your conda env — set that
path in `path_config.m` (`si_np2_python_exe`, see [Configuration](#configuration) below), e.g.
`C:\Users\<you>\miniconda3\envs\spikeinterface\python.exe`.

Everything else in `lib/np2/sorting/spikeGLX_pipeline_np2.py`'s *User input* section (recording
directory, run specs, output destination) is filled in automatically per session when you run
from MATLAB. Edit it by hand only if you run the Python script standalone (see
[the pipeline README](lib/np2/sorting/README.txt)).

---

## Step 7 — (Optional) Cross-session unit tracking: UnitMatch

Only needed for chronic / multi-session recordings where you want to track the same units across
recordings. Skip this if you only sort single recordings.

```bash
conda activate spikeinterface
pip install -U UnitMatchPy       # matching + I/O; pulls mtscomp, mat73, scikit-learn, ...
git clone https://github.com/EnnyvanBeest/UnitMatch    C:\SWC\UnitMatch
```

The clone is only used for its UnitMatchPy source (imported directly, so a `git pull` there
picks up upstream fixes without a reinstall).

See [Cross-session unit tracking](#optional-cross-session-unit-tracking) below and
[lib/np2/matching/README.txt](lib/np2/matching/README.txt) for how to run it.

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
: Where preprocessing output goes. The `fast` location should be an SSD (SpikeInterface and
  Kilosort 4 read the large `.ap.bin` files); data can later be moved to `slow` long-term storage.

`processed_camera_fast_dir` / `processed_camera_slow_dir`
: As above, for camera data.

`figure_dir`
: Where figures are written.

`npy_matlab_dir`
: Local clone of https://github.com/kwikteam/npy-matlab (Step 4).

`spikes_dir`
: Local clone of https://github.com/cortex-lab/spikes (Step 4).

`si_np2_scripts_dir`
: Directory holding the sorting pipeline scripts (this repo's `lib/np2/sorting`).

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
3. `si_sorting` — run the **SpikeInterface** pipeline: bandpass filter + phase shift → IBL
   destriping (per shank) → **Kilosort 4** → SortingAnalyzer (waveforms, template & quality
   metrics) → **Bombcell** curation labels → Phy export.
4. `create_check_clusters_csv` — write the cluster table to review.
5. `create_trigger_file` — extract the probe sync/trigger channel.
6. `create_driftmap` — build the driftmap.
7. `process_camera_data` — compute camera motion energy.

To resume part-way (e.g. after fixing something), use `run_from_step`, which always runs from the
chosen step **to the end**:

```matlab
>> ctl.run_from_step(<probe_id>, 'si_sorting');
```

To run **only** the sorting step, optionally resuming from a later sub-step (see
`help RC2Preprocess.run_sorting_from_step`):

```matlab
>> ctl.run_sorting_from_step(<probe_id>);                                    % full sort
>> ctl.run_sorting_from_step(<probe_id>, 'start_step', 'kilosort4');         % re-sort, skip preprocessing
```

The Kilosort 4 / SpikeInterface output for each probe is written to a directory of the form:

`<processed_probe_fast_dir>\<animal_id>\output\preprocessed_<probe_id>_g0\<probe_id>_g0_imec0\imec0_ks4`

referred to below as `<ks4_dir>`. It contains the standard Phy/Kilosort `.npy` files plus
`cluster_groups.csv` (Bombcell labels), a `csv/` folder with `metrics.csv` and
`waveform_metrics.csv`, a reloadable `sorting_analyzer/`, and `phy/` and `bombcell/` folders.
Next to it, `bandpass_only.ap/` holds the band-pass + phase-shifted (not destriped) recording
used for cross-session matching. See [lib/np2/sorting/README.txt](lib/np2/sorting/README.txt)
for the full output description.

> [!NOTE]
> A sorting run temporarily needs about **3×** the recording size on the output drive:
> `bandpass_only.ap/` (kept) plus the destriped recording SpikeInterface writes to a temporary
> binary before sorting (removed after).

### 3. Check the clusters

Bombcell already labels every unit (good / mua / noise / non-soma) in `cluster_groups.csv`.
Curation is fully automated by default — no manual step is required to get a usable
`selected_clusters.txt`:

```matlab
>> ctl.create_check_clusters_csv(<probe_id>);
```

This writes `clusters_to_check.csv` in `<ks4_dir>` (one row per cluster: `cluster_id`,
`bombcell_group`, `is_non_somatic`, `keep_pipeline`, `keep`), and immediately generates
`selected_clusters.txt` from the automated `keep` decision (Bombcell `good`, somatic).

**Manual override (optional).** Inspect clusters if you want to double check the automated call:

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

Then hand-edit the **`keep`** column of `clusters_to_check.csv` (`1` to force a cluster in, `0` to
force it out) and re-run:

```matlab
>> ctl.create_selected_clusters_txt(<probe_id>);
```

to regenerate `selected_clusters.txt` from your edits.

> [!NOTE]
> `keep_pipeline` sits next to `keep` and always holds the original, untouched automated
> decision — it is only ever used as a reference to check what you changed later
> (`sum(tbl.keep ~= tbl.keep_pipeline)` in MATLAB), never edit it by hand. Only the `keep` column
> is read by `create_selected_clusters_txt`; `clusters_to_check.csv` is fully regenerated (both
> columns reset to the automated decision) every time `create_check_clusters_csv` runs, so make a
> copy of your hand-edited file before re-running it if you want to keep your edits.

### 3b. Adjusting Bombcell thresholds for specific brain regions

Bombcell's default thresholds (below) were calibrated on generic cortical/hippocampal
recordings. Some regions have neurons with atypical waveform shapes that can be wrongly
flagged, so thresholds may need adjusting depending on where you recorded. All thresholds are
in the `bombcell_thresholds` dictionary in the **"User input"** section near the top of
[lib/np2/sorting/spikeGLX_pipeline_np2.py](lib/np2/sorting/spikeGLX_pipeline_np2.py) — edit the
template there, not the generated per-session script (see the note at the top of the file).

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
"mua") if it fails *any one* of that category's metrics. **Plot label** is how each metric is
labelled in `bombcell/metric_histograms.png` — use it to find the matching panel when adjusting
a threshold below:

| Category | Metric | Plot label | Default | Meaning |
|---|---|---|---|---|
| noise | `num_positive_peaks` | `# peaks` | `less: 2` | passes with 0-1 positive peaks; 2+ → noise |
| noise | `num_negative_peaks` | `# troughs` | `less: 1` | passes with 0-1 negative peaks (troughs); 2+ → noise |
| noise | `peak_to_trough_duration` | `waveform duration` | `greater: 0.0001, less: 0.00115` (s) | spike duration must be in this window |
| noise | `waveform_baseline_flatness` | `baseline flatness` | `less: 0.5` | baseline before/after the spike must be reasonably flat |
| noise | `peak_after_to_trough_ratio` | `peak₂/trough` | `less: 0.8` | rebound peak must not be too large relative to the trough |
| noise | `exp_decay` | `spatial decay` | `greater: 0.01, less: 0.1` | waveform decay time constant must be in this window |
| mua | `amplitude_median` | `amplitude` | `greater: 30` (µV) | minimum median spike amplitude |
| mua | `snr` | `SNR` | `greater: 5` | minimum signal-to-noise ratio |
| mua | `amplitude_cutoff` | `spikes missing (%)` | `less: 0.2` (i.e. 20%) | estimated fraction of missed (sub-threshold) spikes must be low |
| mua | `num_spikes` | `# spikes` | `greater: 300` | minimum spike count over the recording |
| mua | `rp_contamination` | `frac. RPVs` | `less: 0.1` | refractory-period violation rate must be low |
| mua | `presence_ratio` | `presence ratio` | `greater: 0.7` | unit must be present through most of the recording |
| mua | `drift_ptp` | `maximum drift` | `less: 100` (µm) | peak-to-peak drift must be limited |
| non-somatic | `peak_before_to_trough_ratio` | `peak₁/trough` | `less: 3` | ratio used to detect axonal/dendritic waveform shape |
| non-somatic | `peak_before_width` | `peak₁ width` | `greater: 0.00015` (s) | pre-trough peak width |
| non-somatic | `trough_width` | `trough width` | `greater: 0.0002` (s) | trough width |
| non-somatic | `peak_before_to_peak_after_ratio` | `peak₁/peak₂` | `less: 3` | ratio used to detect axonal/dendritic waveform shape |
| non-somatic | `main_peak_to_trough_ratio` | `peak_main/trough` | `less: 0.8` | large positive peak relative to trough flags non-somatic origin |

**Not adjustable — for inspection only.** Bombcell also computes and plots these 3 metrics, but
(confirmed against Bombcell's own `classification.py`) none of them ever affects the
good/mua/noise/non-soma label — there is no threshold to set for them:

| Metric | Plot label | Meaning |
|---|---|---|
| `drift_std` | `cum. drift` | cumulative drift over the whole recording (vs. `drift_ptp`'s single worst excursion) |
| `isolation_distance` | `isolation dist.` | how well-separated this unit's spikes are from neighbouring units in PCA feature space |
| `l_ratio` | `L-ratio` | estimated fraction of borderline spikes that may belong to a neighbouring unit |

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
axonal/dendritic. `ClusterSelectionTable.m` (`curation_table`/`curation_mua_table`) surfaces this
as an `is_non_somatic` column and excludes non-somatic units from `keep` by default — hand-edit
`keep` to include a specific non-somatic unit (e.g. for axonal signal analyses).

After changing thresholds, re-run just Bombcell (skipping preprocessing, Kilosort4 and
SortingAnalyzer) with `ctl.run_from_step(probe_id, 'bombcell')` — see
[Tips for debugging](#tips-for-debugging).

### 3c. No Bombcell GUI, and what UnitRefine would add

**No interactive GUI.** The standalone `bombcell` Python/MATLAB package ships its own interactive
curation GUI (`unit_quality_gui.py`, native package only), but this pipeline uses
`spikeinterface.curation.bombcell_label_units` — SpikeInterface's own reimplementation of
Bombcell's labelling logic, not the native package (see [Fixes & compatibility
notes](lib/np2/sorting/README.txt), point 4) — which exposes no GUI, only the programmatic
`thresholds=` interface described above. `metric_histograms.png`, the upset plots and
`waveform_classification.png` in `bombcell/` are this pipeline's substitute for browsing results
interactively.

**UnitRefine, and why this pipeline doesn't use it (yet).** SpikeInterface's own guide
([Auto-label units](https://spikeinterface.readthedocs.io/en/latest/how_to/auto_label_units.html))
recommends running **both** Bombcell and UnitRefine, as complementary — not competing — methods:
Bombcell applies fixed, documented thresholds to individual metrics (what this pipeline does);
UnitRefine instead uses a classifier pre-trained on hand-curated data, which can catch
patterns a fixed threshold on any single metric would miss, at the cost of being less
transparent about *why* a unit was labelled a given way. This pipeline does not currently run
UnitRefine — Bombcell alone is the curation step.

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

## Tips for debugging

The commands above are meant to be enough for a normal run — the pipeline is designed to run
start to finish in one call. The following are for when something goes wrong partway through and
you need to re-run only part of it, or see what actually failed.

**See the Python traceback on failure.** By default the sorting script's console window closes as
soon as it finishes, whether it succeeded or failed, so a Python error is easy to miss. Set:

```matlab
>> ctl.leave_window_open_on_error = true;
```

before calling `preprocess_step_1`/`run_from_step`/`run_sorting_from_step` and the window stays
open on completion (success or failure) so you can read the traceback. Close it manually when
done.

**Resume from a specific step**, instead of re-running everything from scratch:

```matlab
>> ctl.run_from_step(<probe_id>, <start_step>);
```

`help RC2Preprocess` lists the valid `start_step` names, in execution order, and how to resume
partway through sorting specifically -- e.g. re-run SortingAnalyzer + Bombcell without redoing
Kilosort4 via `run_from_step(probe_id, 'postprocess')`, or, if you only changed a Bombcell
threshold, skip SortingAnalyzer too via `run_from_step(probe_id, 'bombcell')` (its metrics do not
depend on Bombcell's thresholds).

**Run only the sorting step**, without chaining to the later stage-1 steps (trigger, driftmap,
camera data):

```matlab
>> ctl.run_sorting_from_step(<probe_id>);
```

See `help RC2Preprocess.run_sorting_from_step` for its `start_step` options.

---

## Optional: cross-session unit tracking

For chronic recordings, track the same units across several sessions of one animal with
**UnitMatch**. This is optional and independent from sorting — any session not yet sorted is
sorted first:

```matlab
>> ctl.match_sessions({probe_id_1, probe_id_2, ...});   % chronological order, N >= 2
```

Each session's SortingAnalyzer for matching is built on its `bandpass_only.ap/` recording (not
the destriped one used for sorting), since destriping removes the cross-channel spatial footprint
UnitMatch matches on. Results are saved next to the first session by default (`MatchTable.csv`,
`MatchingOverview.png`, ...). Full options and outputs are documented in
[lib/np2/matching/README.txt](lib/np2/matching/README.txt).

> [!TIP]
> Before matching, check each session's `ks4_motion.png` (in its `imec0_ks4/`) to see how much
> probe drift Kilosort4 corrected for within that session — UnitMatch does its own, separate
> drift correction *between* sessions, so this is only relevant if you're deciding whether to
> re-sort with `nblocks=0`. See
> [lib/np2/sorting/README.txt, "KILOSORT4 DRIFT CORRECTION AND CROSS-SESSION MATCHING"](lib/np2/sorting/README.txt).

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
