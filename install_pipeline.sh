#!/usr/bin/env bash
#
# install_pipeline.sh — fully automate the rc2_analysis install (see README.md / INSTALL.md).
#
# Run from Git Bash, from inside an already-cloned rc2_analysis checkout, AS ADMINISTRATOR
# if you want CUDA Toolkit installed automatically (Step 2 needs elevation on Windows).
#
# What this DOES automate:
#   - detect GPU / CUDA, download + silently install a matching CUDA Toolkit  [Step 2]
#   - clone the MATLAB helper repos (npy-matlab, spikes)                     [Step 4]
#   - create the `spikeinterface` conda env and pip-install spikeinterface,
#     kilosort, torch (CUDA build auto-detected), phy, upsetplot            [Step 5]
#   - pin numpy/scipy to a single BLAS source                               [Step 5, BLAS note]
#   - verify the environment                                                [Step 5, Verify]
#   - clone + set up `runningmouse` in its OWN conda env                    [Step 6b, optional]
#   - install UnitMatchPy + clone UnitMatch                                 [Step 7, optional]
#   - interactively ask where to put data/outputs, then WRITE path_config.m
#     with every path filled in                                             [Configuration]
#
# What this DOES NOT automate (still manual — see INSTALL.md):
#   - installing MATLAB + toolboxes (Step 3) -- MathWorks requires a license/GUI install,
#     no unattended installer is provided by this script
#   - running setup_paths.m from MATLAB, and filling in experiment_list_csv
#
# Usage:
#   ./install_pipeline.sh [options]
#
# Options:
#   --swc-dir <path>          Root working directory (prompted interactively if not given)
#   --env-name <name>         Conda env name for the main pipeline (default: spikeinterface)
#   --torch-cuda <tag>        Force a PyTorch CUDA tag (cu118, cu121, cu124, cpu) instead of auto-detect
#   --skip-cuda               Skip CUDA Toolkit download/install (use if already installed)
#   --skip-clone              Skip cloning npy-matlab / spikes (already have them)
#   --skip-runningmouse        Skip runningmouse install (camera motion-energy processing)
#   --skip-unitmatch           Skip UnitMatchPy / UnitMatch install (cross-session tracking)
#   --non-interactive          Don't prompt for paths -- use --swc-dir (or its default) for everything
#   -h, --help                 Show this help and exit
#
# Example:
#   ./install_pipeline.sh --swc-dir /c/SWC
#
# Note on CUDA auto-install: downloading + silently installing the NVIDIA CUDA Toolkit
# requires Administrator privileges and can take a while (multi-GB download); the machine
# may need a reboot before the driver is fully active. If you'd rather do this by hand,
# pass --skip-cuda and follow INSTALL.md Step 2.

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults / argument parsing
# ---------------------------------------------------------------------------

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SWC_DIR_DEFAULT="$(cd "${REPO_DIR}/.." && pwd)"
SWC_DIR=""
ENV_NAME="spikeinterface"
PY_VERSION="3.10"
TORCH_CUDA=""            # auto-detected if empty
SKIP_CUDA=0
SKIP_CLONE=0
SKIP_RUNNINGMOUSE=0
SKIP_UNITMATCH=0
NON_INTERACTIVE=0

usage() {
    grep '^#' "${BASH_SOURCE[0]}" | sed '1d;s/^# \{0,1\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --swc-dir)            SWC_DIR="$2"; shift 2 ;;
        --env-name)           ENV_NAME="$2"; shift 2 ;;
        --torch-cuda)         TORCH_CUDA="$2"; shift 2 ;;
        --skip-cuda)          SKIP_CUDA=1; shift ;;
        --skip-clone)         SKIP_CLONE=1; shift ;;
        --skip-runningmouse)  SKIP_RUNNINGMOUSE=1; shift ;;
        --skip-unitmatch)     SKIP_UNITMATCH=1; shift ;;
        --non-interactive)    NON_INTERACTIVE=1; shift ;;
        -h|--help)            usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$1"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$1" >&2; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$1" >&2; exit 1; }

# to_posix <path> -> converts a Windows-style path (C:\foo\bar or C:/foo/bar)
# typed at a prompt into the POSIX form Git Bash needs (/c/foo/bar). A bare
# backslash path left as-is gets silently mangled by bash (backslash is an
# escape character), producing garbage nested folders like "C:foobar" --
# normalise it instead of trusting the user typed /c/... correctly.
to_posix() {
    local p="$1"
    if [[ "$p" =~ ^([A-Za-z]):[\\/](.*)$ ]]; then
        local drive="${BASH_REMATCH[1],,}" rest="${BASH_REMATCH[2]}"
        rest="${rest//\\//}"
        echo "/${drive}/${rest}"
    else
        echo "$p"
    fi
}

# prompt "question" "default" -> echoes the answer (default if empty/non-interactive),
# normalising a Windows-style path if one was typed
prompt() {
    local question="$1" default="$2" answer
    if [[ "${NON_INTERACTIVE}" -eq 1 ]]; then
        echo "${default}"
        return
    fi
    read -r -p "${question} [${default}]: " answer >&2 || true
    to_posix "${answer:-${default}}"
}

# ---------------------------------------------------------------------------
# Ask where everything should live
# ---------------------------------------------------------------------------

log "rc2_analysis pipeline install — where should things go?"

if [[ -z "${SWC_DIR}" ]]; then
    SWC_DIR="$(prompt "Root working directory (SWC-style layout)" "${SWC_DIR_DEFAULT}")"
fi
mkdir -p "${SWC_DIR}"
SWC_DIR="$(cd "${SWC_DIR}" && pwd)"

# Raw data (probe, camera, RC2, motion-clouds) is commonly on the lab's shared
# Ceph storage (Z:\mvelez\...), but that's a per-machine network mapping, not
# guaranteed to exist or be mounted the same way everywhere -- ask, with the
# common lab default suggested, rather than hardcoding it.
RAW_PROBE_DIR="$(prompt "Raw probe data directory" "Z:\\mvelez\\mateoData_probe")"
RAW_CAMERA_DIR="$(prompt "Raw camera data directory" "Z:\\mvelez\\mateoData_cameras")"
RAW_RC2_DIR="$(prompt "Raw RC2 NIDAQ data directory" "Z:\\mvelez\\mateoData_rc2")"
MOTION_CLOUDS_ROOT="$(prompt "Motion Clouds root directory (only used if you run that protocol)" "Z:\\mvelez\\mateoData_mc")"
PROCESSED_PROBE_DIR="$(prompt "Processed probe output directory (fast/SSD ideally)" "${SWC_DIR}/data/processed_data/probe")"
PROCESSED_CAMERA_DIR="$(prompt "Processed camera output directory" "${SWC_DIR}/data/processed_data/cameras")"
FORMATTED_DATA_DIR="$(prompt "Formatted data output directory" "${SWC_DIR}/data/processed_data/formatted_data")"
FIGURE_DIR="$(prompt "Figures output directory" "${SWC_DIR}/data/figures")"
EXPERIMENT_LIST_CSV="$(prompt "Path to experiment_list csv (existing file if you have one, otherwise a new path to create)" "${SWC_DIR}/data/experiment_list.csv")"

for d in "${PROCESSED_PROBE_DIR}" "${PROCESSED_CAMERA_DIR}" "${FORMATTED_DATA_DIR}" \
         "${FIGURE_DIR}" "$(dirname "${EXPERIMENT_LIST_CSV}")"; do
    mkdir -p "${d}"
done

if [[ -f "${EXPERIMENT_LIST_CSV}" ]]; then
    warn "Using existing experiment_list csv at ${EXPERIMENT_LIST_CSV} -- left untouched."
elif [[ ! -f "${EXPERIMENT_LIST_CSV}" ]]; then
    echo "animal_id,date,probe_id,session_id,protocol,experiment_group,np_probe_type,git_commit,discard" \
        > "${EXPERIMENT_LIST_CSV}"
    log "Created empty experiment_list csv with the correct header: ${EXPERIMENT_LIST_CSV}"
fi

echo "  repo dir       : ${REPO_DIR}"
echo "  SWC dir        : ${SWC_DIR}"
echo "  conda env      : ${ENV_NAME} (python ${PY_VERSION})"

# ---------------------------------------------------------------------------
# GPU / CUDA detection (also used to pick the PyTorch build)
# ---------------------------------------------------------------------------

log "Detecting GPU / CUDA"

HAS_GPU=0
DETECTED_CUDA_MAJOR_MINOR=""

if command -v nvidia-smi >/dev/null 2>&1; then
    if nvidia-smi >/dev/null 2>&1; then
        HAS_GPU=1
        # nvidia-smi reports the *driver's max supported* CUDA version, e.g. "CUDA Version: 12.4"
        DETECTED_CUDA_MAJOR_MINOR="$(nvidia-smi | grep -oE 'CUDA Version: [0-9]+\.[0-9]+' | grep -oE '[0-9]+\.[0-9]+' | head -1 || true)"
        echo "  GPU detected. Driver supports up to CUDA ${DETECTED_CUDA_MAJOR_MINOR:-unknown}"
    fi
fi

if [[ "${HAS_GPU}" -eq 0 ]]; then
    warn "No NVIDIA GPU detected via nvidia-smi. Falling back to CPU-only PyTorch (~10x slower, see README)."
fi

if [[ -z "${TORCH_CUDA}" ]]; then
    if [[ "${HAS_GPU}" -eq 1 && -n "${DETECTED_CUDA_MAJOR_MINOR}" ]]; then
        MAJOR="${DETECTED_CUDA_MAJOR_MINOR%%.*}"
        MINOR="${DETECTED_CUDA_MAJOR_MINOR##*.}"
        if   [[ "${MAJOR}" -ge 13 || ( "${MAJOR}" -eq 12 && "${MINOR}" -ge 4 ) ]]; then
            TORCH_CUDA="cu124"; CUDA_INSTALLER_VERSION="12.4"
        elif [[ "${MAJOR}" -eq 12 ]]; then
            TORCH_CUDA="cu121"; CUDA_INSTALLER_VERSION="12.1"
        elif [[ "${MAJOR}" -eq 11 ]]; then
            TORCH_CUDA="cu118"; CUDA_INSTALLER_VERSION="11.8"
        else
            warn "Unrecognised CUDA major version '${MAJOR}', defaulting to cu118"
            TORCH_CUDA="cu118"; CUDA_INSTALLER_VERSION="11.8"
        fi
    else
        TORCH_CUDA="cpu"
    fi
fi

echo "  Selected PyTorch build: ${TORCH_CUDA}"
echo "  (override any time with --torch-cuda <cu118|cu121|cu124|cpu>)"

# ---------------------------------------------------------------------------
# CUDA Toolkit — download + silent install
# ---------------------------------------------------------------------------

if [[ "${SKIP_CUDA}" -eq 1 ]]; then
    log "Skipping CUDA Toolkit install (--skip-cuda)"
elif [[ "${TORCH_CUDA}" == "cpu" ]]; then
    log "No GPU detected -- skipping CUDA Toolkit install (CPU-only run)"
else
    log "Installing CUDA Toolkit ${CUDA_INSTALLER_VERSION}"
    warn "This step needs Administrator privileges and downloads several GB. A reboot may be" \
         2>/dev/null || true
    warn "required before the driver/toolkit is fully active. Pass --skip-cuda to do this by hand."

    if command -v nvcc >/dev/null 2>&1; then
        EXISTING_CUDA="$(nvcc --version 2>/dev/null | grep -oE 'release [0-9]+\.[0-9]+' | grep -oE '[0-9]+\.[0-9]+')"
        warn "CUDA Toolkit already installed (nvcc reports ${EXISTING_CUDA:-unknown}) -- skipping" \
             "download/install. Pass --torch-cuda to match the PyTorch build to it if needed."
    else
        case "${CUDA_INSTALLER_VERSION}" in
            12.4) CUDA_URL="https://developer.download.nvidia.com/compute/cuda/12.4.1/local_installers/cuda_12.4.1_551.78_windows.exe" ;;
            12.1) CUDA_URL="https://developer.download.nvidia.com/compute/cuda/12.1.1/local_installers/cuda_12.1.1_531.14_windows.exe" ;;
            11.8) CUDA_URL="https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda_11.8.0_522.06_windows.exe" ;;
            *)    CUDA_URL="" ;;
        esac

        if [[ -z "${CUDA_URL}" ]]; then
            warn "No known installer URL for CUDA ${CUDA_INSTALLER_VERSION} -- install manually from" \
                 "https://developer.nvidia.com/cuda-downloads (see INSTALL.md Step 2)"
        else
            CUDA_INSTALLER="${TEMP:-/tmp}/cuda_installer.exe"
            log "Downloading CUDA ${CUDA_INSTALLER_VERSION} installer (this can take a while)..."
            if command -v curl >/dev/null 2>&1; then
                curl -L -o "${CUDA_INSTALLER}" "${CUDA_URL}" \
                    || warn "CUDA download failed -- install manually (INSTALL.md Step 2)"
            else
                warn "curl not found -- cannot auto-download CUDA. Install manually (INSTALL.md Step 2)."
            fi

            if [[ -f "${CUDA_INSTALLER}" ]]; then
                log "Running CUDA installer silently (-s) -- requires Administrator; this can take 10-20 min"
                # -s = silent install of every component; NVIDIA installer must run elevated.
                if "${CUDA_INSTALLER}" -s; then
                    warn "CUDA Toolkit installed. A reboot is recommended before running the sorting pipeline."
                else
                    warn "CUDA silent install failed or needs elevation -- rerun this script as" \
                         "Administrator (or --skip-cuda and install manually, see INSTALL.md Step 2)."
                fi
                rm -f "${CUDA_INSTALLER}"
            fi
        fi
    fi
fi

# ---------------------------------------------------------------------------
# MATLAB helper repositories (npy-matlab, spikes)
# ---------------------------------------------------------------------------

# All external repos the pipeline depends on (npy-matlab, spikes, runningmouse,
# UnitMatch) live side by side here -- not part of rc2_analysis itself, but
# required to run it.
HELPERS_DIR="${SWC_DIR}/helpers"
mkdir -p "${HELPERS_DIR}"

if [[ "${SKIP_CLONE}" -eq 0 ]]; then
    log "Cloning MATLAB helper repositories"

    if [[ -d "${HELPERS_DIR}/npy-matlab" ]]; then
        warn "npy-matlab already exists, skipping clone"
    else
        git clone https://github.com/kwikteam/npy-matlab "${HELPERS_DIR}/npy-matlab"
    fi

    if [[ -d "${HELPERS_DIR}/spikes" ]]; then
        warn "spikes already exists, skipping clone"
    else
        git clone https://github.com/cortex-lab/spikes "${HELPERS_DIR}/spikes"
    fi
else
    log "Skipping helper repo cloning (--skip-clone)"
fi

# ---------------------------------------------------------------------------
# Python (SpikeInterface) conda environment
# ---------------------------------------------------------------------------

log "Creating conda environment '${ENV_NAME}'"

# Try PATH first; if not found, fall back to common install locations rather
# than giving up -- a fresh Git Bash session often doesn't have conda's shell
# hook loaded yet even when conda itself is installed and even after adding
# `source .../conda.sh` to ~/.bashrc (that only affects new interactive
# shells, not necessarily this script's invocation).
CONDA_SH=""
if command -v conda >/dev/null 2>&1; then
    CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
else
    for candidate in \
        "${USERPROFILE:-$HOME}/miniconda3" \
        "${USERPROFILE:-$HOME}/anaconda3" \
        "/c/Users/${USER:-$USERNAME}/miniconda3" \
        "/c/Users/${USER:-$USERNAME}/anaconda3" \
        "/c/ProgramData/miniconda3" \
        "/c/ProgramData/anaconda3"
    do
        if [[ -f "${candidate}/etc/profile.d/conda.sh" ]]; then
            CONDA_SH="${candidate}/etc/profile.d/conda.sh"
            break
        fi
    done
fi

[[ -n "${CONDA_SH}" && -f "${CONDA_SH}" ]] || die "conda not found on PATH or in common install locations. Install Miniconda first: https://docs.conda.io/en/latest/miniconda.html (or activate it manually, then rerun this script)"

# shellcheck disable=SC1091
source "${CONDA_SH}"

if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    warn "conda env '${ENV_NAME}' already exists, reusing it"
else
    conda create -n "${ENV_NAME}" python="${PY_VERSION}" -y
fi

conda activate "${ENV_NAME}"

log "Installing SpikeInterface, Kilosort 4, PyTorch (${TORCH_CUDA}), Phy, upsetplot"

pip install "spikeinterface[full]"
pip install kilosort

if [[ "${TORCH_CUDA}" == "cpu" ]]; then
    pip install torch
else
    pip install torch --index-url "https://download.pytorch.org/whl/${TORCH_CUDA}"
fi

pip install phy
pip install upsetplot

log "Pinning numpy/scipy to a single BLAS source (see README BLAS/LAPACK note)"
pip install --force-reinstall numpy scipy

log "Verifying numpy LAPACK works (should print [2. 3.], not crash)"
python -c "import numpy as np; print(np.linalg.solve([[3,1],[1,2]],[9,8]))" \
    || die "numpy LAPACK check failed — see README's BLAS/LAPACK note (Step 5)"

log "Verifying the environment"

python -c "import spikeinterface; print('SI', spikeinterface.__version__)"
python -c "import kilosort; print('KS4 OK')"
python -c "import torch; print('CUDA available:', torch.cuda.is_available())"
python -c "from spikeinterface.curation import bombcell_label_units; print('Bombcell OK')"
python -c "import upsetplot; print('upsetplot OK')"

SPIKEINTERFACE_PYTHON_EXE="$(conda run -n "${ENV_NAME}" python -c "import sys; print(sys.executable)")"

# ---------------------------------------------------------------------------
# runningmouse (camera motion-energy processing) — its own env, per upstream
# repo having no pinned dependencies and pulling its own ffmpeg via `av`
# ---------------------------------------------------------------------------

RUNNINGMOUSE_ENV="runningmouse"
RUNNINGMOUSE_PYTHON_EXE=""
RUNNINGMOUSE_MAIN_SCRIPT=""

if [[ "${SKIP_RUNNINGMOUSE}" -eq 1 ]]; then
    log "Skipping runningmouse install (--skip-runningmouse)"
else
    log "Installing runningmouse (camera motion-energy processing) in its own env '${RUNNINGMOUSE_ENV}'"

    if conda env list | grep -qE "^${RUNNINGMOUSE_ENV}[[:space:]]"; then
        warn "conda env '${RUNNINGMOUSE_ENV}' already exists, reusing it"
    else
        conda create -n "${RUNNINGMOUSE_ENV}" python="${PY_VERSION}" -y
    fi

    conda activate "${RUNNINGMOUSE_ENV}"
    # av via pip (not conda) so it ships its own static ffmpeg and doesn't touch
    # this env's BLAS/numpy resolution; runningmouse itself has no version pins.
    pip install numpy scipy pillow matplotlib av
    conda deactivate

    RUNNINGMOUSE_DIR="${HELPERS_DIR}/runningmouse"
    if [[ -d "${RUNNINGMOUSE_DIR}" ]]; then
        warn "runningmouse already exists at ${RUNNINGMOUSE_DIR}, skipping clone"
    else
        git clone -b master https://github.com/neuroinformatics-unit/runningmouse "${RUNNINGMOUSE_DIR}"
    fi

    RUNNINGMOUSE_PYTHON_EXE="$(conda run -n "${RUNNINGMOUSE_ENV}" python -c "import sys; print(sys.executable)")"
    RUNNINGMOUSE_MAIN_SCRIPT="${RUNNINGMOUSE_DIR}/difference_video/main.py"

    log "Verifying runningmouse imports"
    ( cd "${RUNNINGMOUSE_DIR}" && conda run -n "${RUNNINGMOUSE_ENV}" python -c "import tools.plot, tools.tools; print('runningmouse OK')" ) \
        || warn "runningmouse import check failed -- check ${RUNNINGMOUSE_DIR} manually"
fi

# ---------------------------------------------------------------------------
# UnitMatchPy + UnitMatch — required by lib/np2/matching/match_sessions.py
# (hard `import UnitMatchPy...` in that script, not merely optional), so this
# installs by default; pass --skip-unitmatch if you will never run match_sessions.
# ---------------------------------------------------------------------------

UNITMATCH_DIR=""

if [[ "${SKIP_UNITMATCH}" -eq 1 ]]; then
    log "Skipping UnitMatchPy / UnitMatch install (--skip-unitmatch)"
    warn "lib/np2/matching/match_sessions.py hard-imports UnitMatchPy -- it will fail to run" \
         "until you install it (see INSTALL.md Step 7)."
else
    log "Installing UnitMatchPy (required by match_sessions.py) + cloning UnitMatch"
    conda activate "${ENV_NAME}"
    pip install -U UnitMatchPy

    UNITMATCH_DIR="${HELPERS_DIR}/UnitMatch"
    if [[ -d "${UNITMATCH_DIR}" ]]; then
        warn "UnitMatch already exists at ${UNITMATCH_DIR}, skipping clone"
    else
        git clone https://github.com/EnnyvanBeest/UnitMatch "${UNITMATCH_DIR}"
    fi
fi

# ---------------------------------------------------------------------------
# Write path_config.m with every path now known
# ---------------------------------------------------------------------------

log "Writing path_config.m"

PATH_CONFIG="${REPO_DIR}/path_config.m"

# to_win <posix_path> -> Windows-style backslash path for MATLAB string literals
to_win() {
    local p="$1"
    if command -v cygpath >/dev/null 2>&1; then
        cygpath -w "$p"
    else
        # best-effort: /c/foo/bar -> C:\foo\bar
        echo "$p" | sed -E 's#^/([a-zA-Z])/#\1:/#' | sed 's#/#\\#g'
    fi
}

RAW_PROBE_DIR_WIN="$(to_win "${RAW_PROBE_DIR}")"
RAW_CAMERA_DIR_WIN="$(to_win "${RAW_CAMERA_DIR}")"
RAW_RC2_DIR_WIN="$(to_win "${RAW_RC2_DIR}")"
MOTION_CLOUDS_ROOT_WIN="$(to_win "${MOTION_CLOUDS_ROOT}")"

REPO_DIR_WIN="$(to_win "${REPO_DIR}")"
EXPERIMENT_LIST_CSV_WIN="$(to_win "${EXPERIMENT_LIST_CSV}")"
FORMATTED_DATA_DIR_WIN="$(to_win "${FORMATTED_DATA_DIR}")"
PROCESSED_PROBE_DIR_WIN="$(to_win "${PROCESSED_PROBE_DIR}")"
PROCESSED_CAMERA_DIR_WIN="$(to_win "${PROCESSED_CAMERA_DIR}")"
FIGURE_DIR_WIN="$(to_win "${FIGURE_DIR}")"
NPY_MATLAB_DIR_WIN="$(to_win "${HELPERS_DIR}/npy-matlab")"
SPIKES_DIR_WIN="$(to_win "${HELPERS_DIR}/spikes")"
SI_NP2_SCRIPTS_DIR_WIN="$(to_win "${REPO_DIR}/lib/np2/sorting")"
SI_NP2_TEMPLATE_WIN="$(to_win "${REPO_DIR}/lib/np2/sorting/spikeGLX_pipeline_np2.py")"
SI_PYTHON_EXE_WIN="$(to_win "${SPIKEINTERFACE_PYTHON_EXE}")"

if [[ -n "${RUNNINGMOUSE_PYTHON_EXE}" ]]; then
    RUNNINGMOUSE_PYTHON_EXE_WIN="$(to_win "${RUNNINGMOUSE_PYTHON_EXE}")"
    RUNNINGMOUSE_MAIN_SCRIPT_WIN="$(to_win "${RUNNINGMOUSE_MAIN_SCRIPT}")"
else
    RUNNINGMOUSE_PYTHON_EXE_WIN='<not installed -- run without --skip-runningmouse, or set manually>'
    RUNNINGMOUSE_MAIN_SCRIPT_WIN='<not installed -- run without --skip-runningmouse, or set manually>'
fi

if [[ -n "${UNITMATCH_DIR}" ]]; then
    UNITMATCH_DIR_WIN="$(to_win "${UNITMATCH_DIR}")"
else
    UNITMATCH_DIR_WIN='<not installed -- run without --skip-unitmatch, or set manually>'
fi

cat > "${PATH_CONFIG}" <<EOF
function config = path_config()
% PATH_CONFIG Configuration information on the system
%
%   CONFIG = path_config()
%   returns a list of paths on the system allowing the user access data and code
%   See INSTALL.md / README for a description of the entries.
%
%   Generated by install_pipeline.sh on $(date '+%Y-%m-%d %H:%M:%S').

% the path containing the .git for rc2_analysis
config.git_work_tree_dir        = '${REPO_DIR_WIN}';

config.experiment_list_csv      = '${EXPERIMENT_LIST_CSV_WIN}';
config.formatted_data_dir       = '${FORMATTED_DATA_DIR_WIN}';

config.raw_probe_dir            = '${RAW_PROBE_DIR_WIN}';
config.raw_camera_dir           = '${RAW_CAMERA_DIR_WIN}';
config.raw_rc2_dir              = '${RAW_RC2_DIR_WIN}';

config.processed_probe_fast_dir = '${PROCESSED_PROBE_DIR_WIN}';
config.processed_probe_slow_dir = '${PROCESSED_PROBE_DIR_WIN}';

config.processed_camera_fast_dir = '${PROCESSED_CAMERA_DIR_WIN}';
config.processed_camera_slow_dir = '${PROCESSED_CAMERA_DIR_WIN}';

config.figure_dir               = '${FIGURE_DIR_WIN}';

config.npy_matlab_dir           = '${NPY_MATLAB_DIR_WIN}';
config.spikes_dir               = '${SPIKES_DIR_WIN}';

% SpikeInterface-based NP2 sorting pipeline.
config.si_np2_scripts_dir  = '${SI_NP2_SCRIPTS_DIR_WIN}';
config.si_np2_template     = '${SI_NP2_TEMPLATE_WIN}';
config.si_np2_python_exe   = '${SI_PYTHON_EXE_WIN}';

config.runningmouse_python_exe  = '${RUNNINGMOUSE_PYTHON_EXE_WIN}';
config.runningmouse_main_script = '${RUNNINGMOUSE_MAIN_SCRIPT_WIN}';

% Clone of github.com/EnnyvanBeest/UnitMatch, used by RC2Preprocess.match_sessions
% (cross-session unit tracking). Passed explicitly to match_sessions.py so it
% never falls back to that script's own hardcoded DEFAULT_UNITMATCH_REPO.
config.unitmatch_repo_dir       = '${UNITMATCH_DIR_WIN}';

config.motion_clouds_root       = '${MOTION_CLOUDS_ROOT_WIN}';
EOF

log "path_config.m written to ${PATH_CONFIG}"

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

log "Done."
cat <<EOF

Automated:
  - CUDA Toolkit ${CUDA_INSTALLER_VERSION:-"(skipped)"}
  - MATLAB helper repos cloned under: ${HELPERS_DIR}/
  - conda env '${ENV_NAME}' with SpikeInterface, Kilosort 4, PyTorch (${TORCH_CUDA}), Phy, upsetplot
  - runningmouse: $([ -n "${RUNNINGMOUSE_PYTHON_EXE}" ] && echo "installed in env '${RUNNINGMOUSE_ENV}'" || echo "skipped")
  - UnitMatchPy / UnitMatch: $([ -n "${UNITMATCH_DIR}" ] && echo "installed, cloned to ${UNITMATCH_DIR}" || echo "skipped")
  - path_config.m written with every path above filled in

Still required manually (see INSTALL.md):
  1. If CUDA was just installed, reboot before running the sorting pipeline.
  2. Install MATLAB + required toolboxes: Signal Processing, Statistics and Machine
     Learning, Image Processing, Parallel Computing, Optimization (Step 3).
  3. Run setup_paths.m from MATLAB in ${REPO_DIR}.
  4. Fill in the experiment_list csv (header already written for you):
       ${EXPERIMENT_LIST_CSV}

EOF
