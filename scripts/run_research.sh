#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
PYTHON="${PYTHON:-python3}"
if [[ -x "$ROOT/.venv/bin/python" ]]; then PYTHON="$ROOT/.venv/bin/python"; fi
export PYTHON
case "${1:-smoke}" in
  smoke)
    make verify smoke PYTHON="$PYTHON"
    ;;
  core)
    make verify smoke degree8 completion8 map-proof symbolic PYTHON="$PYTHON"
    "$PYTHON" -m chiral4form lift runs/degree8/fields_prime30011.json \
      runs/degree8/fields_prime30013.json runs/degree8/fields_prime30029.json \
      --bound 500 --output runs/degree8-rational-candidate.json
    "$PYTHON" -m chiral4form complete runs/degree8-generalized/fields_prime30011.json \
      --candidates I6_2 I8_2 I8_3 I8_4 I8_5 I8_6 --target-rank 10 --depth 3 \
      --output runs/completion8-finite-search.json
    ;;
  degree12)
    "$PYTHON" -m chiral4form source --fetch --checkout .cache/upstream --output data/imported
    "$PYTHON" -m chiral4form plan --source .cache/upstream --config configs/degree12.json
    "$PYTHON" -m chiral4form run --source .cache/upstream --config configs/degree12.json --output runs/degree12
    ;;
  *) printf 'Usage: bash scripts/run_research.sh {smoke|core|degree12}\n' >&2; exit 2 ;;
esac
