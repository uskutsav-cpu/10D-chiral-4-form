#!/bin/bash
#SBATCH --job-name=sd5-deg12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=0-63
#SBATCH --output=logs/deg12-%A_%a.out
#SBATCH --error=logs/deg12-%A_%a.err
set -euo pipefail

# Degree twelve is NOT claimed in the manuscript. This job exists so that the
# extension can be run reproducibly, not because any result depends on it.

REPO_COMMIT="${REPO_COMMIT:?set REPO_COMMIT}"
REPO_URL="${REPO_URL:?set REPO_URL}"
ARCHIVE_DIR="${ARCHIVE_DIR:?set ARCHIVE_DIR}"

SHARD="${SLURM_ARRAY_TASK_ID}"
NSHARDS="${SLURM_ARRAY_TASK_COUNT}"

# Output is immutable: a shard that already produced a result refuses to run
# again, so a resubmitted array cannot silently overwrite completed work.
OUTDIR="${SLURM_SUBMIT_DIR}/deg12/shard-${SHARD}"
if [ -f "${OUTDIR}/COMPLETE" ]; then
    echo "shard ${SHARD} already complete; refusing to overwrite"
    exit 0
fi
mkdir -p "${OUTDIR}" logs

WORK="${SLURM_SUBMIT_DIR}/build-${SLURM_ARRAY_JOB_ID}"
if [ ! -d "${WORK}/repo" ]; then
    mkdir -p "${WORK}"
    git clone --quiet "${REPO_URL}" "${WORK}/repo"
    ( cd "${WORK}/repo" && git checkout --quiet "${REPO_COMMIT}" \
      && python -m venv .venv \
      && .venv/bin/pip install --quiet -r requirements-lock.txt \
      && .venv/bin/pip install --quiet opt_einsum networkx joblib python-igraph )
fi
cd "${WORK}/repo"

# Deterministic per-shard seed: derived from the shard index alone, so a rerun
# of shard k explores exactly the same candidates.
SEED=$(( 20260731 + SHARD * 7919 ))

export OMP_NUM_THREADS=1        # joblib supplies the parallelism

.venv/bin/python spinor_trace_bridge/scripts/run_exact_jacobian.py \
    --primes 32749,32719 \
    --seeds "${SEED}" \
    --degrees 4,6,8,10,12 \
    --patience 60 \
    --max-graphs 4000 \
    --out "${OUTDIR}/exact_jacobian.json"

touch "${OUTDIR}/COMPLETE"
echo "shard ${SHARD} of ${NSHARDS} complete -> ${OUTDIR}"
