#!/bin/bash
#SBATCH --job-name=glm_histbase_cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
# mem/time right-sized 2026-06-23: 180G locked the job out of the idle 60G nodes
# (peak RSS ~18G even in both-probes mode); 2h timed out a both-probes run. 56G fits
# an idle 60G node and schedules instantly; 6h covers the serial both-probes path.
# Override per-run with --mem/--time as needed (and prefer probe-split: PROBE=<id>).
#SBATCH --mem=56G
#SBATCH --time=06:00:00
#SBATCH --output=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.out
#SBATCH --error=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.err
#
# Goggles rf_local GLM, History-in-baseline — PARALLEL CPU run (maximise-compute:
# never the serial path). IRLS CPU backend (no GPU JIT, which dominates the small
# per-cluster fits), per-cluster fits process-parallelised with loky, BLAS/OMP
# pinned to 1 thread/worker (no oversubscription + removes the IRLS cv-bps
# thread-nondeterminism).
#
# MODES via env vars:
#   N_SHARDS=<M>          -> SHARD ARRAY task: needs a job array of 0..(2*M-1);
#                           probe = idx/M, shard = idx%M -> fit clusters[shard::M]
#                           of that probe into _runs/<probe>/_shards/ (no aggregate).
#                           Many small tasks pack into a fragmented cpu partition.
#   PROBE=<id>            -> fit ONLY that probe (writes _runs/<probe>/, no aggregate)
#   RUN_MODE=aggregate    -> skip fitting; merge shards + aggregate() + diagnostics
#   RUN_MODE=diagnostics  -> only (re)compute variance_partition + diagnostic figs
#                           (parallel over the node's cores; a GLM-refit workload)
#   (none)                -> full run: both probes + aggregate (single node)
# Submit the SHARD ARRAY (M shards/probe, few cores each) + dependent aggregate:
#   J=$(sbatch --parsable --array=0-$((2*M-1)) --cpus-per-task=6 \
#              --export=ALL,N_SHARDS=$M this.sbatch)
#   sbatch --dependency=afterok:$J --cpus-per-task=8 --mem=64G \
#          --export=ALL,RUN_MODE=aggregate this.sbatch
#
# Prereqs: repo at $REPO on branch glm-improvements; venv at $REPO/python/.venv;
# ceph data staged.

set -euo pipefail

REPO="$HOME/Source/github/SainsburyWellcomeCentre/rc2_analysis"
mkdir -p /ceph/margrie/laura/motion_clouds/_glm_logs

source "$REPO/python/.venv/bin/activate"
cd "$REPO/python"

# One BLAS/OMP thread per loky worker (workers x 1 thread == cores requested).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

N_JOBS="${SLURM_CPUS_PER_TASK:-20}"
DRIVER="scripts/run_glm_goggles_rf_sfor_20ms.py"

PROBES_ARR=(CAA-1124370_rec1_rec2_rec3 CAA-1124371_rec1_rec2_rec3)

if [ -n "${N_SHARDS:-}" ]; then
    IDX="${SLURM_ARRAY_TASK_ID:?N_SHARDS set but this is not a job array}"
    PR="${PROBES_ARR[$((IDX / N_SHARDS))]}"
    SID="$((IDX % N_SHARDS))"
    python "$DRIVER" --history-in-baseline --backend irls --n-jobs "$N_JOBS" \
        --probe "$PR" --n-shards "$N_SHARDS" --shard-id "$SID"
elif [ "${RUN_MODE:-}" = "diagnostics" ]; then
    # variance_partition refits the GLM ~40x/cluster (LOO-unique x CV) -> a
    # cluster-class job, not a plotting one. Parallelise over the node's cores.
    # Run via -m so the loky worker (_partition_one) pickles cleanly.
    python -m scripts.diagnostics_rf_sfor --run histbase --n-jobs "$N_JOBS"
elif [ "${RUN_MODE:-}" = "aggregate" ]; then
    python "$DRIVER" --history-in-baseline --backend irls --aggregate-only
elif [ -n "${PROBE:-}" ]; then
    python "$DRIVER" --history-in-baseline --backend irls --n-jobs "$N_JOBS" --probe "$PROBE"
else
    python "$DRIVER" --history-in-baseline --backend irls --n-jobs "$N_JOBS"
fi
