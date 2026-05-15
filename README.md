# rc2_analysis
Preprocessing and analysis pipeline for electrophysiological data acquired with the rollercoaster setup. It runs CatGT (preprocessing), Kilosort 2.5 (spike sorting), and ecephys postprocessing modules on Neuropixels 2.0 4-shank recordings acquired with SpikeGLX.

> [!IMPORTANT]
> ### Reference and reproducibility
>
> This codebase is provided as a reference for the publication:
> **Velez-Fort, Cossell, Porta, Clopath, Margrie (202x), *Title* Journal.**
>
> The pipeline has evolved since the analyses presented in the paper.  
> To ensure **full reproducibility of the original results**, the following 
> legacy configuration may be required:
>
> - MATLAB R2021a (Windows 10)
> - Kilosort2: commit `2a399268d6e1710f482aed5924ba90d52718452a`
>
> These versions correspond to earlier iterations of the pipeline used in 
> the publication and may differ from the current recommended setup.
>
> The installation guide below reflects a more recent and maintained 
> configuration (e.g. Kilosort 2.5) and has been thoroughly documented. 
> If you have any questions, please open an issue.


> [!CAUTION]
> This codebase will not be actively maintained. Installation may be complex and could potentially be unsuccessful.
# Installation Guide

This guide walks you through installing everything needed to run the 
[rc2_analysis](https://github.com/SainsburyWellcomeCentre/rc2_analysis) spike sorting pipeline on a new Windows machine. Follow the steps in order.

> [!NOTE]
> **Neuropixels Phase 3A compatibility:**
>
> If you are working with **Neuropixels Phase 3A data**, you must use the 
> following version of the SpikeGLX fork of ecephys_spike_sorting:
>
> - ecephys_spike_sorting (SpikeGLX fork):  
>   commit `55ff943892577cee497e28dd6a201f4f101b0777`
>
> This version contains the appropriate support for Phase 3A probes and is 
> required for correct pipeline execution.


---



## Requirements | Version compatibility

The table below shows the combination we tested and recommend. 
**The correct CUDA version is determined by your GPU, which then constrains the compatible versions of Visual Studio and MATLAB — check your GPU (Step 2) before installing anything else.**

| Component | Tested version | 
|---|---|
| MATLAB |  R2025b |
| Python | 3.9 |
| CUDA Toolkit | 12.2 |
| Visual Studio | Community 2019 (16.11.54) |
| CatGT / TPrime / C_Waves | 4.8 / latest |
| Kilosort | 2.5 |


---

## Step 1 — Set up the working directory

Before starting, create your working directory and clone ```rc2_analysis``` into it:


```bash
mkdir C:\SWC
cd C:\SWC
git clone https://github.com/SainsburyWellcomeCentre/rc2_analysis
cd rc2_analysis
git checkout 64bf4517a2135f1591d52394e9877475c9649b90
```
>We used commit 64bf4517 (26 March 2026) — you can use a different commit, but compatibility with the rest of this guide is not guaranteed.

All other repositories and tools will be installed alongside it. The guide uses C:\SWC\ as a reference — adapt to your machine.


```
SWC/
├── rc2_analysis/
├── original_pipeline/
│   ├── Kilosort/
│   ├── npy-matlab/
│   ├── spikes/
│   ├── ecephys_spike_sorting/
│   ├── CatGT-win/
│   ├── TPrime-win/
│   └── C_Waves-win/
└── data/
    ├── raw_data/
    ├── processed_data/
    │   └── formatted_data/
    ├── figures/
    └── temp/
        └── kilosort/
```


Each step in this guide specifies where to clone or extract its component.



---


## Step 2 — Check your GPU and install CUDA Toolkit

Kilosort runs on the GPU and needs CUDA to compile its MEX files.

1. Find your GPU model: right-click the desktop → Display settings → 
Advanced display → note the GPU name
2. Check which CUDA version your GPU supports on the 
[CUDA compatibility table](https://developer.nvidia.com/cuda-gpus)
3. Download the matching CUDA Toolkit from 
[developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads)

We tested with **CUDA Toolkit 12.2**, which is backward compatible with 
older hardware.



---
## Step 3 — Install MATLAB

Download and install MATLAB from [mathworks.com](https://mathworks.com). 
Ideally, use the same version as listed in the Requirements table (R2025b), 
but any version from R2021a onwards should work. Make sure your GPU is 
supported — check the [GPU Computing Requirements](https://uk.mathworks.com/help/parallel-computing/gpu-computing-requirements.html) 
page for details on supported compute capabilities.


**Make sure to include these toolboxes during installation:**

- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox
- Image Processing Toolbox
- Parallel Computing Toolbox

To check which toolboxes are already installed, run this in the MATLAB 
command window:
```matlab
ver
```

Once MATLAB is installed, verify your GPU is recognised — open the MATLAB 
command window and run:

```matlab
gpuDevice()
```
This should return your GPU name and memory without errors. If you get a driver error, update your NVIDIA drivers via Device Manager in Windows.


---
## Step 4 — Install Visual Studio 2019

Kilosort needs a C++ compiler to build its MEX files. On Windows, this 
means **Visual Studio Community 2019** specifically — newer versions may 
not be compatible with your CUDA version.

Download it here: 
[visualstudio.microsoft.com/vs/older-downloads](https://visualstudio.microsoft.com/vs/older-downloads)

During installation, select the **"Desktop development with C++"** workload.

![Capture d’écran 2026-05-05 à 09.38.31](https://hackmd.io/_uploads/HJEIT7w0-e.png)

---

## Step 5 — Clone the required MATLAB repositories

You need three additional MATLAB repositories at **specific commits**. 
We recommend downloading and using [Git Bash](https://git-scm.com/downloads) for these commands.

Open Git Bash, navigate to `C:\SWC\original_pipeline\` and run:

```bash
git clone https://github.com/MouseLand/Kilosort
cd Kilosort
git checkout v2.5.2
cd ..

git clone https://github.com/kwikteam/npy-matlab
cd npy-matlab
git checkout b7b0a4ef6ba26d98a8c54e651d5444083c88311c
cd ..

git clone https://github.com/cortex-lab/spikes
cd spikes
git checkout 8dc617c6dd5b279341506b7411004772ba05f4cc
cd ..
```

> **Why specific commits?** The pipeline was written against these exact versions. Newer versions may have breaking changes.

Then compile Kilosort's MEX files: open MATLAB, navigate to`C:\SWC\original_pipeline\Kilosort\` and run:

```matlab
mexGPUall
```
>If compilation succeeded, you should see no error messages in the MATLAB command window and new `.mexw64` files will have been created in `C:\SWC\original_pipeline\Kilosort\CUDA\`.


---



## Step 6 — Set up a Python environment via Miniconda

The pipeline calls Kilosort and the ecephys tools from Python. The correct 
Python version must satisfy two constraints:

- Compatible with your MATLAB version — check the [MathWorks Python compatibility table](https://mathworks.com/support/requirements/python-compatibility.html)
- Compatible with the ecephys pipeline — tested on Python 3.8.10 and 3.9 only

For MATLAB R2025b, **Python 3.9** is the only version that satisfies both constraints.

If you don't have Python installed, download and install Miniconda from [docs.conda.io](https://docs.conda.io/en/latest/miniconda.html). Then open an Anaconda Prompt and create a dedicated environment:

```bash
conda create --name rc2_env python=3.9
conda activate rc2_env
```

---

## Step 7 — Install ecephys_spike_sorting


Still in the Anaconda Prompt, make sure to activate the environment, and navigate to `C:\SWC\original_pipeline\`:

```bash
conda activate rc2_env
cd C:\SWC\original_pipeline
```

Then clone the correct fork at the commit tested for NP 2.0 4-shank probes:

```bash
git clone https://github.com/jenniferColonell/ecephys_spike_sorting
cd ecephys_spike_sorting
git checkout 60c40251ed568fb036b4364e615a261d3afc4800
```

Then install it with the exact dependency versions required — conflicts will occur with newer versions:

```bash
pip install -e .
pip uninstall setuptools -y && pip install setuptools==59.8.0
pip uninstall argschema -y && pip install argschema==1.17.5
pip uninstall marshmallow -y && pip install marshmallow==2.19.2
pip install h5py
pip install phylib
```

See the [ecephys README](https://github.com/jenniferColonell/ecephys_spike_sorting) for more details on the pipeline modules.

---
## Step 8 — Download CatGT, TPrime, and C_Waves

These are pre-compiled tools that require no installation — just download and extract. 
Download the latest and available Windows versions of CatGT, TPrime and C_Waves from the [SpikeGLX download page](https://billkarsh.github.io/SpikeGLX/) and extract each to a convenient location, for example:


`C:\SWC\original_pipeline\CatGT-win`
`C:\SWC\original_pipeline\TPrime-win`
`C:\SWC\original_pipeline\C_Waves-win`



---

## Step 9 — Install the MATLAB Engine for Python

This allows the pipeline to launch MATLAB from Python. Open an **Anaconda Prompt as administrator** (right-click → Run as administrator).
![Capture d’écran 2026-05-05 à 14.08.40](https://hackmd.io/_uploads/Byji3PP0be.png)


Then run the following commands:

```bash
conda activate rc2_env
cd "C:\Program Files\MATLAB\R2025b\extern\engines\python"
python setup.py install
```

Replace `R2025b` with your actual MATLAB version folder name.

> **Note:** this step must be repeated if you ever rebuild your conda 
> environment.

Close the administrator prompt when done. Use a normal Anaconda Prompt 
from here on.

---





## Step 10 — Configure create_input_json.py

Open`C:\SWC\original_pipeline\ecephys_spike_sorting\scripts\create_input_json.py` and update the hardcoded paths at the top of the file to match your machine:

```python
ecephys_directory = r'C:\SWC\original_pipeline\ecephys_spike_sorting\ecephys_spike_sorting'
kilosort_repository = r'C:\SWC\original_pipeline\Kilosort'
KS2ver = '2.5'   # must match the repository above
npy_matlab_repository = r'C:\SWC\original_pipeline\npy-matlab'
catGTPath = r'C:\SWC\original_pipeline\CatGT-win'
tPrime_path = r'C:\SWC\original_pipeline\TPrime-win'
cWaves_path = r'C:\SWC\original_pipeline\C_Waves-win'
kilosort_output_tmp = r'D:\kilosort_datatemp'
```
> [!IMPORTANT]
> **CatGT command string must be updated manually.** The version of CatGT 
> (4.8 or later) uses a different command syntax than 
> what is hardcoded in `spikeGLX_pipeline_np2.py`. You must update the 
> `catGT_cmd_string` from:
> ```
> -prb_fld -out_prb_fld -aphipass=300 -aplopass=9000 -gbldmx -gfix=0,0.10,0.02 -SY=0,384,0,0
> ```
> To:
> ```
> -g=0 -t=0,0 -prb=0 -ap -gbldmx -prb_fld -out_prb_fld -apfilter=butter,12,300,9000 -gfix=0.4,0.10,0.02
> ```
> Edit `spikeGLX_pipeline_np2.py` directly — not the generated copy.

> [!IMPORTANT]
> `rc2_analysis` will overwrite the run script at each 
> preprocessing run by copying `spikeGLX_pipeline_np2.py` and inserting 
> session-specific parameters. If you want to change the CatGT command 
> string or any Kilosort settings permanently, edit 
> `spikeGLX_pipeline_np2.py` or `create_input_json.py` — **not** the 
> generated copy.

---
## Step 11 — Install Phy for manual cluster curation

We recommend to install Phy in its **own separate environment** to avoid dependency 
conflicts. Copy this section in GitBash:

```bash
conda create --name phy2 python=3.9
conda activate phy2
pip install phylib==2.6.0
pip install phy==2.0b5
pip install joblib==1.1.1
pip install numpy==1.23.5
```

If Phy crashes on startup, delete the `.phy` folder in your home 
directory (`C:\Users\<you>\`) and try again.

>Note: Phy can also be run on a Mac (including M-chip machines) — follow the platform-specific installation instructions in the [Phy repository](https://github.com/cortex-lab/phy).


---

# Configuration 


Before running the pipeline for the first time, configure `path_config.m` to point to your data directories, run `setup_paths.m` to load the MATLAB paths, and populate the experiment list CSV with your session details.


---
## Configuration file

You must point to the following directories in the `path_config.m` file.

git_work_tree_dir
: The directory which contains the .git folder tracking rc2_analysis. This is important when saving data/figures so that we also save the commit SHA1 used when generating the data/figures.

experiment_list_csv
: Full path to `.csv` containing a table of experiment details (full description below).

formatted_data_dir
: Directory containing the eventual formatted `.mat` files. After preprocessing of the data, we combine it into a formatted data file. Paths to the formatted data files are of form `<formatted_data_dir>\<probe_id>.mat`

raw_probe_dir
: Directory containing the raw probe data. This directory contains subdirectories named with the animal IDs. Paths to the raw data depend on which probe was used for recording. (in the following the <probe_id> has the form `<animal_id>_<session_suffix1>_<session_suffix2>_...`)

: For an older version of SpikeGLX used with NP Phase 3A probes:
  `<raw_probe_dir>\<animal_id>\<probe_id>_g0_t0.imec.ap.bin`
  
: For a newer version of SpikeGLX with NP 2.0, 4 shank probes:
  `<raw_probe_dir>\<animal_id>\<probe_id>_g0\<probe_id>_g0_imec0\<probe_id>_g0_t0.imec0.ap.bin`

raw_camera_dir
: Directory containing the camera .avi files. This directory contains subdirectories with the session IDs (i.e. of the form `<animal_id>_<session_suffix>`). Paths to the raw data are of the form: `<raw_camera_dir>\<animal_id>_<session_suffix>\camera0.avi`

raw_rc2_dir
: Directory containing the RC2 .bin files. This directory contains subdirectories with the animal IDs, which themselves contain raw data aquired on the NIDAQ of RC2. Paths to the raw data files are of form: `<raw_rc2_dir>\<animal_id>\<animal_id>\<animal_id>_<session_suffix>_001.bin`

processed_probe_fast_dir
: Directory in which we will save the initial pre-processing of the raw probe data. Pre-processing of the data requires access to the large .ap.bin files therefore they benefit from being on a 'fast' storage device (SSD). Raw data is initially moved here to carry out the preprocessing steps before being moved (manually) to long-term storage.

processed_probe_slow_dir
: Directory in which we will save the pre-processed data for more long term storage. After formatting of the data, we need to access this data less, so it can go on slower, long-term storage devices.

processed_camera_fast_dir/processed_camera_slow_dir
: Similar to above but for the processing of the camera data.

figure_dir
: Directory in which figures will be stored.

npy_matlab_dir
:   Directory with local clone of https://github.com/kwikteam/npy-matlab
 		
spikes_dir
: Directory with local clone of https://github.com/cortex-lab/spikes

ecephys_scripts_dir
: Directory where we will move a modified copy of the template script <ecephys_template>. Used for Neuropixels Phase 3A processing.

ecephys_template
: Full path to a template file which will be used to run ecephys_spike_sorting. During preprocessing this file is modified and moved to the <ecephys_scripts_dir> and used. Used for Neuropixels Phase 3A processing. 

ecephys_np2_python_exe
: The python executable used to start ecephys_spike_sorting (used as we run it in a python virtual environment). Used for Neuropixels Phase 3A processing.

ecephys_np2_scripts_dir
: Directory where we will move a modified copy of the template script <ecephys_template>. Used for Neuropixels 2.0 processing.

ecephys_np2_template
: Full path to a template file which will be used to run ecephys_spike_sorting. During preprocessing this file is modified and moved to the <ecephys_scripts_dir> and used. Used for Neuropixels 2.0 processing. 

ecephys_python_exe
: The python executable used to start ecephys_spike_sorting (used as we run it in a python virtual environment). Used for Neuropixels 2.0 processing.

runningmouse_python_exe
: Python executable used to start runningmouse

runningmouse_main_script
: Main start script for starting runningmouse



## Setup

Navigate to the `rc2_analysis` directory (e.g. `C:\Users\lee\Documents\mvelez\rc2_analysis`).
Run `setup_paths.m`, which adds required directories to the MATLAB path:
```
>> setup_paths
```

## Experiment list file

A few important details needs to be manually entered in a `.csv`, the full path to which is listed in `path_config.m` in `<experiment_list_csv>`.

Each RC2 session has a row in the `.csv`. The details that need to be inserted are:

- *animal_i*d: 		ID of the animal used for the session.
- *probe_id*: 		ID of the probe recording linked to the RC2 session (of form 

   `<probe_id> = <animal_id>_<session_suffix1>_<session_suffix2>_...`
- *session_id*: 		ID of the RC2 session (of form `<session_id> = <animal_id>_<session_suffix>_001`)
- *protocol*:		The protocol used for the session (e.g. mismatch_nov20 or sparse_noise)
- *experiment_group*:	Each probe recording is part of an 'experiment group' which is used when determining which probe IDs to load during analysis
- *probe_type*:		Type of probe used ('3A' or '24')
- *git_commit*: 		SHA1 of the git commit used for the RC2 acquisition (only for reference)
- *discard*: 		Whether we are discarding the experiment


---
# Preprocessing and Analysis Workflow 

This section describes how to run the pipeline after installation, including:
- preprocessing of raw data
- spike sorting and curation
- quality control steps
- data formatting and analysis

---


## Preprocessing data

>[!IMPORTANT]**Update the metadata files**
>
>Before running the pipeline, the `.meta` files from SpikeGLX need four 
>manual corrections for NP 2.0 probes. Open each `.meta` file in a text 
>editor (e.g. Notepad) and change:
>
>| Field | Old value | New value |
>|---|---|---|
>| `imDatPrb_pn` | `NP2013` | `NP2010` |
>| `imDatPrb_type` | `2013` | `24` |
>| `imroTbl` | `(2013,...)` | `(24,...)` |
>| `snsGeomMap` | `(NP2013, 4, ...)` | `(4,...)` — >remove `NP2013,` |
>
> This needs to be done for every recording (one `.meta` file per animal per session).

---



### 0. Experiment list file
Enter the experiment details into `<experiment_list_csv>` (set in  `path_config.m`, e.g. for me it is, `D:\mvelez\experiment_list.csv`).


### 1. Create controller object
```
>> ctl = RC2Preprocess();
```


### 2. Run the first step of preprocessing:
```
>> ctl.preprocess_step_1(<probe_id>);
```

where <probe_id> is of the form 'CAA-1115688_rec1_rec2'.

This will:
- move files from winstor to local drive (fast SSD) (from `<raw_probe_dir>` to `<processed_probe_fast_dir>`)
- run ecephys_spike_sorting (incl. kilosort2)
- move csvs to a single location
- extract the probe sync channel (TTL) to separate file
- create a driftmap
- process the camera data
    
Recently, since upgrading to Windows 10, the Windows console has occassionally been getting stuck after some of the modules in ecephys_spike_sorting. It can be resumed by pressing Ctrl-D at the terminal. But of course this is annoying as we cannot just leave it to run.
 
The output of `ecephys_spike_sorting`/`kilosort2` is put in a directory with the form:

`<processed_probe_fast_dir>\<animal_id>\output\catgt_<probe_id>_g0\<probe_id>_g0_imec0\imec0_ks2`

which will be referred to below as `<kilosort_dir>`.


### 3. Checking the clusters

a. Open file containing 'good' clusters in spreadsheet viewer (`clusters_janelia.csv`). Path is of form:
`<kilosort_dir>\csv\clusters_janelia.csv`

b. Open results in Phy (e.g. at Windows Command Prompt)
```
 > cd Documents\phy2
 >.venv\Scripts\activate
 > D:
 > cd D:\mvelez\mateoData_probe\janelia_pipeline\CAA-1115689\output\catgt_CAA-1115689_rec1_rec2_g0\CAA-1115689_rec1_rec2_g0_imec0\imec0_ks2
 > phy template-gui params.py
```

c. Plot cluster information in MATLAB (e.g. at MATLAB prompt):
```
>> cluster_info = ctl.cluster_info(<probe_id>);
>> cluster_info.plot(<cluster_id>); 
```
<cluster_id> is an integer, not a string.

Historically, we have created a new file `clusters_janelia.xlsx`, and added our judgements as separate columns (with headers 'mateo' and 'lee').
After going through the clusters labelling each cluster with:  g - green (clearly good), w - white (passable), b - brown (discard)

After saving this new table as 'clusters_janelia.xlsx' we then run:
```
>> ctl.create_selected_clusters_txt(<probe_id>);
```

Which creates a file `selected_clusters.txt` in the main kilosort directory. If at least one experimenter discards a cluster (labelled 'b') the cluster is discarded.

###4. Manually check the trigger

RC2 sessions which were accidentally started (and therefore promptly stopped), appear as a sequence of triggers on the probe sync channel. However, they are not associated with a protocol. Therefore, we must remove them from protocol related analysis. To do this we can run:

```
>> ctl.correct_trigger_file(<probe_id>);
```

This opens a small gui in which you can click and drag across the triggers you want to remove and save the results.


###5. Check LFP power profile and anatomy

The following should be run for each shank separately:
```
>> hf_power = ctl.hf_power(<probe_id>, <shank_id>);
```
which returns a HighFrequencyPowerProfile object for the probe and shank ID. (NOTE: even if there is a single shank you have to provide the ID, which would be 0)

To compute the high-frequency power profile:
```
>> hf_power.run();
```
To save the results, pass the `hf_power` object back to the controller.
```
>> ctl.save_hf_power(hf_power);
```

The output is put in the directory `<kilosort_dir>\tracks\`.
The following information is saved:

- `<kilosort_dir>\tracks\offset_<shank_id>.txt`

   Stores the offset between anatomical and electrophysiological layer 5.
       If no anatomy is present (see below) when the above code is run, a zero will be saved there. 

- `<kilosort_dir>\tracks\hf_power_<shank_id>.pdf`
   
   Figure with high frequency power depth profile for all batches (left), averaged power profile (middle) and histogram of multi-unit count along the probe (right).
   Electrophysiological peak is indicated with red line.
   
   If anatomy is also present (see below) the brain region boundaries are shown unshifted (left) and shifted (middle and right) according to the offset difference between anatomy and electrophysiology. 
- `<kilosort_dir>\tracks\hf_power_<shank_id>.mat`
   
   Parameters used to create the HF power information as well as the power on each channel (can be loaded at a later date to browse the results).

*Anatomy*

The above code can be run without any anatomy information.
However, if the anatomy from `brainreg-segment` is ready it should be placed in the tracks directory `<kilosort_dir>\tracks`
in a file with name `track_<shank_id>.csv`. NOTE: At this point, the anatomy data acquired in `brainreg-segment`has to be taken according to the following rules: 1) Always start to trace the probe track from the pia; 2) Make sure the number of sampling points chosen to fit the traced probe track corresponds to 1um - e.g. if the probe was introduced 1750um deep into the brain, make sure to choose 1750 sampling points when fitting the traced probe track.
The above steps can be run again, when the anatomy information is present, and the region boundaries will be plotted and the offset
between anatomical layer 5 and electrophysiological layer 5 will be computed and saved to a file of the form `<kilosort_dir>\tracks\offset_<shank_id>.txt`.

*Selecting batches*

There are more options for exploring the high-frequency power profile using the HighFrequencyPowerProfile class. 
Indeed, most recordings will need some manual exploration of the data to get a good idea of the L5 peak and other features.
This usually includes:

  - choosing which batches to use for the computation
  - choosing a depth range in which to look for the peak in cortex
   
To choose batches run:

```
>> hf_power.plot_raw_batches();
```

This plots the power profile computed in each batch and gives an idea of which batches have clear peaks in cortex.

To restrict analysis using particular batches set the `batches_to_use` property:
```
>> hf_power.batches_to_use = [2, 4, 6, 8:10];  % this will use batches 2, 4, 6, 8, 9, 10  
```

To restrict the region to use to look for the peak of L5 set the `search_above` and `search_below` properties.

```
>> hf_power.search_above = 600;
>> hf_power.search_below = 1000;
```
Will restrict search to between 600 and 1000 um above the probe tip.

A further useful function is the interactive viewer, which is started:
```
>> hf_power.interactive_plot();
```

This brings up a figure which allows you to shift the boundaries of the anatomy relative to the power profile to explore other possible good relationships between anatomy and electrophysiology.

### 6. Format the data

Once the above steps are complete we create a formatted `.mat` file which will contain a single data structure with all the information of the probe recording.
```
>> ctl.format(<probe_id>);
```
Saved to `<formatted_data_dir>` with name `<probe_id>.mat`.

The format includes:

- splitting of the experimental conditions contained in each protocol into 'trials'
- synchronization to the probe recording with trials
- allocation of clusters to brain regions according to the anatomy `track_<shank_id>.csv` file and any offset from step 5.

In addition, `.csv` files are created in separate location with further summary data. These are:

`<formatted_data_dir>\csvs\trial_matched_offsets\<probe_id>.csv`
Contains a table with the offset (in samples) for each replay trial, in order to align it with its original trial.

`<formatted_data_dir>\csvs\stationary_vs_motion_fr\<probe_id>.csv`
Contains a table with the firing rate of each cluster in each trial during the motion and stationary periods.


###7. Load the formatted data

To work with the data run:
```
>> experiment = ctl.load_formatted_data(<probe_id>);
```
   
This creates a `FormattedData` object in the workspace with access to all data in the recording.

Work in progress: If you are using "protocol always vis", add the <session_id> in lib/data_cleansing.m file.


## Figures

There are several scripts which can be used to plot commonly used figures:

Useful for all experiments:

1. Trial structure

    `rc2_analysis\scripts\plot\trial_structure.m`
2. Overlay aligned trials
    
    `rc2_analysis\scripts\plot\overlay_aligned_trials.m`
3. Rasters around motion onset or mismatch
    
    `rc2_analysis\scripts\plot\rasters.m`

4. Heatmaps + population average trace around motion onset or mismatch

    `rc_analysis\scripts\plot\heatmaps.m`

5. Unity plots
    
    `rc2_analysis\scripts\plot\unity_plots`

    a. mismatch, baseline vs. response

    b. Stationary vs. motion (where motion is protocol dependent)

    c. Motion vs. motion

    d. Each of above for population and single clusters

6. MI vs. depth plots

    `rc2_analysis\scripts\plot\unity_plots`

    a. mismatch, baseline vs. response

    b. Stationary vs. motion (where motion is protocol dependent)

    c. Motion vs. motion

7. Tuning curves

    

	
	
	
	