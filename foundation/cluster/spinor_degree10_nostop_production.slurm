#!/bin/bash
#SBATCH --job-name=sd5-deg10-nostop
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/deg10-nostop-%j.out
#SBATCH --error=logs/deg10-nostop-%j.err
set -euo pipefail

# The one long computation this project attempted and did not finish.
#
# Every archived scan stopped at the Hilbert target, a number supplied from
# outside the computation, so none of them says anything about the ansatz. With
# --no-hilbert-stop the terminal status becomes informative: degrees 4, 6 and 8
# reached candidate_exhaustion on a laptop. Degree 10 was killed after roughly
# four hours, at rank 12 of 14, still inside its first batch of 64 candidates.
#
# What this job would settle, in either direction:
#   * terminal status candidate_exhaustion at rank 14 -- the port-graph ansatz
#     is sufficient at degree 10, and the degree-8 shortfall does not recur;
#   * candidate_exhaustion at rank < 14 -- it is not, and the missing directions
#     are located, as at degree 8 where structured tensor words supplied them.
#
# Either outcome is a result. A time-limited stop is not, and the recorder marks
# it as such rather than reporting the rank reached.

# ---- pinned inputs; edit these three, nothing else -------------------------
REPO_COMMIT="${REPO_COMMIT:?set REPO_COMMIT to the exact commit to reproduce}"
REPO_URL="${REPO_URL:?set REPO_URL}"
ARCHIVE_DIR="${ARCHIVE_DIR:?set ARCHIVE_DIR to your own copy of the spinor archive}"

# Do not run this before spinor_degree10_smoke.slurm has passed on this cluster.
WORK="${SLURM_SUBMIT_DIR}/deg10-nostop-${SLURM_JOB_ID}"
OUT="${SLURM_SUBMIT_DIR}/deg10-nostop-out"       # immutable: never reused
mkdir -p "${WORK}" logs
if [ -e "${OUT}" ]; then
    echo "refusing to overwrite ${OUT}; move it aside or change OUT" >&2
    exit 1
fi
cd "${WORK}"

git clone --quiet "${REPO_URL}" repo
cd repo
git checkout --quiet "${REPO_COMMIT}"

python -m venv .venv
.venv/bin/pip install --quiet -r requirements-lock.txt
.venv/bin/pip install --quiet opt_einsum networkx joblib python-igraph

# The archive is third-party and is not redistributed. Verify the copy the
# operator supplied before computing anything from it.
.venv/bin/python scripts/verify_archive.py --archive "${ARCHIVE_DIR}" \
    --manifest release_candidate/spinor-archive/MANIFEST.json

cd "${ARCHIVE_DIR}"

# --optimize random-greedy is not a preference. opt_einsum's `auto` switches to
# exact dynamic programming around twenty operands and a degree-10 amputated
# contraction has nineteen; the path SEARCH then runs away before any array is
# allocated. That is what silently stalled two earlier runs.
#
# --dense-I trades memory for a better contraction order on the topologies where
# the factorised plan forces a bad one. At 64G that trade is available; on the
# laptop it was not.
"${WORK}/repo/.venv/bin/python" -m sd5_invariants.cli scan \
    --degrees 10 \
    --no-hilbert-stop \
    --candidate-mode exhaustive \
    --optimize random-greedy \
    --dense-I \
    --batch-size 64 \
    --n-jobs "${SLURM_CPUS_PER_TASK}" \
    --joblib-backend threading \
    --seed 20260801 \
    --outdir "${OUT}"

# Turn the run into a certificate, including if it stopped early. The recorder
# distinguishes candidate_exhaustion from a resource stop and refuses to report
# a rank as though it were terminal when it is not.
cd "${WORK}/repo"
.venv/bin/python spinor_trace_bridge/scripts/record_no_stop_run.py \
    --rundir "${OUT}" \
    --out-json verification/spinor_degree10_no_stop.json \
    --out-md verification/SPINOR_DEGREE10_NO_STOP.md

echo
echo "Read stopped_by in verification/spinor_degree10_no_stop.json before"
echo "quoting any rank. candidate_exhaustion is a result; time_limit is not."
