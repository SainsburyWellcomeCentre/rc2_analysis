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
# never the serial path). IRLS CPU backend (no GPU JIT, which dominates the small
# per-cluster fits), per-cluster fits process-parallelised with loky, BLAS/OMP
# pinned to 1 thread/worker (no oversubscription + removes the IRLS cv-bps
# thread-nondeterminism).
#
# THREE MODES via env vars (PER-PROBE SPLIT = the way to use 2 nodes):
#   PROBE=<id>            -> fit ONLY that probe (writes _runs/<probe>/, no aggregate)
#   RUN_MODE=aggregate    -> skip fitting; aggregate() + diagnostics over both probes
#   (neither)             -> full run: both probes serial + aggregate (single node)
# Submit the split as:
#   A=$(sbatch --parsable --export=ALL,PROBE=<p1> this.sbatch)
#   B=$(sbatch --parsable --export=ALL,PROBE=<p2> this.sbatch)
#   sbatch --dependency=afterok:$A:$B --export=ALL,RUN_MODE=aggregate this.sbatch
# Both probe jobs run in parallel on 2 cpu nodes; aggregate fires once both finish.
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

if [ "${RUN_MODE:-}" = "aggregate" ]; then
    python "$DRIVER" --history-in-baseline --backend irls --aggregate-only
elif [ -n "${PROBE:-}" ]; then
    python "$DRIVER" --history-in-baseline --backend irls --n-jobs "$N_JOBS" --probe "$PROBE"
else
    python "$DRIVER" --history-in-baseline --backend irls --n-jobs "$N_JOBS"
fi
