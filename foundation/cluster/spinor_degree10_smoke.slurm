#!/bin/bash
#SBATCH --job-name=sd5-smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/smoke-%j.out
#SBATCH --error=logs/smoke-%j.err
set -euo pipefail

# ---- pinned inputs; edit these three, nothing else -------------------------
REPO_COMMIT="${REPO_COMMIT:?set REPO_COMMIT to the exact commit to reproduce}"
REPO_URL="${REPO_URL:?set REPO_URL}"
ARCHIVE_DIR="${ARCHIVE_DIR:?set ARCHIVE_DIR to your own copy of the spinor archive}"

WORK="${SLURM_SUBMIT_DIR}/run-${SLURM_JOB_ID}"
mkdir -p "${WORK}" logs
cd "${WORK}"

git clone --quiet "${REPO_URL}" repo
cd repo
git checkout --quiet "${REPO_COMMIT}"

python -m venv .venv
.venv/bin/pip install --quiet -r requirements-lock.txt
.venv/bin/pip install --quiet opt_einsum networkx joblib python-igraph

# 1. the environment must reproduce what a laptop already established
.venv/bin/python -m pytest -q
( cd spinor_trace_bridge && ../.venv/bin/python -m pytest -q )

# 2. verify the third-party archive against its manifest before using it
.venv/bin/python scripts/verify_archive.py --archive "${ARCHIVE_DIR}" \
    --manifest release_candidate/spinor-archive/MANIFEST.json

# 3. the actual smoke calculation: low degrees with the stopping rule disabled
.venv/bin/python spinor_trace_bridge/scripts/run_comparison.py \
    --primes 32749 --degrees 4,6,8,10 \
    --out "${WORK}/smoke_comparison.json"

echo "SMOKE OK. Do not submit the production job unless the ranks above read"
echo "1, 2, 7, 14 and the span comparison matches the committed certificate."
