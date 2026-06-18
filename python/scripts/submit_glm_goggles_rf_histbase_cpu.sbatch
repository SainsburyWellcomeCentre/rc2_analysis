#!/bin/bash
#SBATCH --job-name=glm_histbase_cpu
#SBATCH --partition=gpu_lowp
#SBATCH --account=swc-ac
#SBATCH --cpus-per-task=96
#SBATCH --mem=360G
#SBATCH --time=4:00:00
#SBATCH --output=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.out
#SBATCH --error=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.err
#
# Goggles rf_local GLM, History-in-baseline — MASSIVELY PARALLEL CPU variant
# (maximise-compute rule: never run the serial path). The IRLS CPU backend has
# no GPU JIT/transfer overhead (which dominates the small per-cluster fits — see
# feedback_maximise_compute), so we process-parallelise the per-cluster fits
# across ~96 loky workers on one fat gpu_lowp node (192 cores / 760 G), with
# BLAS/OMP pinned to 1 thread per worker (no oversubscription + removes the IRLS
# cv-bps thread-nondeterminism). Both probes run in ONE job so the driver's
# aggregate() + root figs/ + diagnostics/ fire once, after all fits finish.
# Outputs land on ceph (~/local_data/motion_clouds symlink); scp the CSVs/PDFs
# back to the laptop afterwards. CPU-only (no --gres=gpu): the node's GPUs stay
# free for GPU jobs.
#
# Prereqs (as for submit_glm_goggles_rf_histbase.sbatch): repo at $REPO on
# branch glm-improvements; venv at $REPO/python/.venv; ceph data staged.

set -euo pipefail

REPO="$HOME/Source/github/SainsburyWellcomeCentre/rc2_analysis"
mkdir -p /ceph/margrie/laura/motion_clouds/_glm_logs

source "$REPO/python/.venv/bin/activate"
cd "$REPO/python"

# One BLAS/OMP thread per loky worker (workers x 1 thread == cores requested).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

N_JOBS="${SLURM_CPUS_PER_TASK:-96}"

python scripts/run_glm_goggles_rf_sfor_20ms.py \
    --history-in-baseline \
    --backend irls \
    --n-jobs "$N_JOBS"
