#!/bin/bash
#SBATCH --job-name=glm_histbase_cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
#SBATCH --time=2:00:00
#SBATCH --output=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.out
#SBATCH --error=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.err
#
# Goggles rf_local GLM, History-in-baseline — PARALLEL CPU run (maximise-compute:
# never the serial path). The IRLS CPU backend has no GPU JIT/transfer overhead
# (which dominates the small per-cluster fits — feedback_maximise_compute), so we
# process-parallelise the per-cluster fits across the node's cores with loky,
# BLAS/OMP pinned to 1 thread per worker (no oversubscription + removes the IRLS
# cv-bps thread-nondeterminism). Both probes in ONE job so aggregate() + root
# figs/ + diagnostics/ fire once, after all fits finish.
#
# One node on the cpu partition (20 cores / 180 G). NOTE on granularity: the
# parallelism is over CLUSTERS within a node (loky); fanning across many nodes
# would re-load the 1.4 GB probe .mat per node, so the data load would dominate.
# For a ~212-fit cohort, one node at full core-parallelism is the sweet spot
# (a 2-node per-probe split is the most that helps; beyond that is wasteful).
# CPU partition, no GPU. Outputs land on ceph; scp the CSVs/PDFs back after.
#
# Prereqs: repo at $REPO on branch glm-improvements; venv at $REPO/python/.venv;
# ceph data staged (as for submit_glm_goggles_rf_histbase.sbatch).

set -euo pipefail

REPO="$HOME/Source/github/SainsburyWellcomeCentre/rc2_analysis"
mkdir -p /ceph/margrie/laura/motion_clouds/_glm_logs

source "$REPO/python/.venv/bin/activate"
cd "$REPO/python"

# One BLAS/OMP thread per loky worker (workers x 1 thread == cores requested).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

N_JOBS="${SLURM_CPUS_PER_TASK:-20}"

python scripts/run_glm_goggles_rf_sfor_20ms.py \
    --history-in-baseline \
    --backend irls \
    --n-jobs "$N_JOBS"
