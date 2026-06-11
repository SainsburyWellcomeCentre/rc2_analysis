#!/bin/bash
#SBATCH --job-name=gabor_mc
#SBATCH --partition=gpu_lowp
#SBATCH --account=swc-ac
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --array=0-35
#SBATCH --output=/ceph/margrie/laura/goggle_gabor/logs/%A_%a.out
#SBATCH --error=/ceph/margrie/laura/goggle_gabor/logs/%A_%a.err

set -euo pipefail
source ~/goggle_gabor/.venv/bin/activate
python ~/goggle_gabor/gabor_extract_gpu.py --cloud-index "$SLURM_ARRAY_TASK_ID"
