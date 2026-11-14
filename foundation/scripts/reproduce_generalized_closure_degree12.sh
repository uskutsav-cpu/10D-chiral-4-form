#!/bin/sh
# Reproduce the Phase 1 generalized-flow deficit result.
# Expects the interacting certificates to be present; see
# docs/generalized_flow_complete_through_degree12.md.
set -e
PY=${PY:-.venv/bin/python}

echo "== degree 10 =="
$PY scripts/find_missing_flow_directions.py --degree 10 \
    --out results/generalized_flow/degree10_missing_directions.json

echo "== degree 12 =="
$PY scripts/find_missing_flow_directions.py --degree 12 --single-prime-scan \
    --out results/generalized_flow/degree12_missing_directions.json

echo "== targeted tests =="
$PY -m pytest tests/test_generalized_flow_deficits.py -q

echo
echo "expected: degree 10 -> I10_6, I10_7, I10_12"
echo "          degree 12 -> I12_59, I12_60, I12_61, I12_62"
