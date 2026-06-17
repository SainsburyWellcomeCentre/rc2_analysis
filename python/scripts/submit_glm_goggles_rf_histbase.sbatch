#!/bin/bash
#SBATCH --job-name=glm_histbase
#SBATCH --partition=gpu_lowp
#SBATCH --account=swc-ac
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.out
#SBATCH --error=/ceph/margrie/laura/motion_clouds/_glm_logs/%j.err
#
# Goggles rf_local GLM, History-in-baseline, on the JAX/nemos GPU backend.
# One job, BOTH probes (serial): the driver only runs aggregate() + the root
# figs/ + diagnostics/ once, after both probes finish, so a per-probe array
# would fire aggregate on incomplete data. Outputs land under
# ~/local_data/motion_clouds (symlinked to /ceph/margrie/laura/motion_clouds),
# i.e. on ceph; the CSVs + PDFs are scp'd back to the laptop afterwards.
#
# Prereqs (set up once by the cluster-setup step, see CHANGELOG):
#   - rc2_analysis cloned at $REPO below (branch glm-improvements)
#   - a venv at $REPO/python/.venv with rc2 deps + nemos + jax[cuda12]
#   - ~/local_data/motion_clouds -> /ceph/margrie/laura/motion_clouds symlink
#   - sequence/.mat, image_folders/.mat, RF parquet cohort staged on ceph

set -euo pipefail

REPO="$HOME/Source/github/SainsburyWellcomeCentre/rc2_analysis"
mkdir -p /ceph/margrie/laura/motion_clouds/_glm_logs

source "$REPO/python/.venv/bin/activate"
cd "$REPO/python"

# Both probes, History as an always-on baseline nuisance, GPU JAX engine.
python scripts/run_glm_goggles_rf_sfor_20ms.py \
    --history-in-baseline \
    --backend nemos \
    --device gpu
