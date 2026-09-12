#!/usr/bin/env bash
#SBATCH --job-name=chiral4form
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=chiral4form-%A_%a.log
set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
: "${SLURM_ARRAY_TASK_ID:?Use a job array}"
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONPATH="$PWD/src"
.venv/bin/python -m chiral4form worker \
  --prime-index "$SLURM_ARRAY_TASK_ID" \
  --config "${CONFIG:-configs/degree8.json}" \
  --output "${OUTPUT:-runs/cluster-degree8}"
