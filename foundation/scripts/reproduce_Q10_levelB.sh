#!/usr/bin/env bash
# Reproduce the degree-10 published-candidate result.
#
#   QUICK  (default)  verify committed artifacts and cheap invariants
#   FULL              additionally re-evaluate the projection from scratch
#
# QUICK is designed to run in a fresh clone with no cached checkpoints and to
# finish in about a minute. It verifies claims against committed artifacts and
# recomputes only what is cheap and exact; it does NOT re-evaluate the atlas,
# which is what FULL is for.
#
# Usage:  scripts/reproduce_Q10_levelB.sh [QUICK|FULL]
set -euo pipefail

MODE="${1:-QUICK}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PY="${PYTHON:-}"
if [ -z "$PY" ]; then
  if [ -x .venv/bin/python ]; then PY=.venv/bin/python; else PY=python3; fi
fi

echo "=== reproduce_Q10_levelB.sh [$MODE] ==="
echo "repo    : $ROOT"
echo "commit  : $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo')"
echo "python  : $($PY --version 2>&1)"
echo

# Checkpoints must never be written to an eventually-consistent filesystem.
# A fresh clone may sit anywhere, so pin the checkpoint root to local temp.
export SDINV_CKPT_ROOT="${SDINV_CKPT_ROOT:-${TMPDIR:-/tmp}/sdinv_ckpt_repro}"
echo "ckpt    : $SDINV_CKPT_ROOT"
echo

fail() { echo "FAIL: $*" >&2; exit 1; }

# --- 1. all twelve candidates are registered and evaluable -----------------
echo "--- 1. registration of all twelve equation-(4.24) candidates"
$PY - <<'PYEOF' || fail "registration check"
import sys
sys.path.insert(0, "src")
from sdinv.published_degree10_invariants import (
    AMBIGUITY_VARIANTS, NOT_IMPLEMENTED, PUBLISHED_DEGREE10)
expected = {f"P10_{i:02d}" for i in range(1, 13)}
assert set(PUBLISHED_DEGREE10) | set(NOT_IMPLEMENTED) == expected, "coverage"
assert not (set(PUBLISHED_DEGREE10) & set(NOT_IMPLEMENTED)), "disjoint"
print(f"  implemented {len(PUBLISHED_DEGREE10)}/12, "
      f"blocked {len(NOT_IMPLEMENTED)}, variants {len(AMBIGUITY_VARIANTS)}")
for name, spec in sorted(PUBLISHED_DEGREE10.items()):
    if spec.get("ambiguity"):
        print(f"    {name}: reading recorded under {spec['ambiguity']}")
PYEOF

# --- 2. every candidate is a Lorentz scalar --------------------------------
# This is the check that caught the non-scalar P10_07. A rotation cannot see a
# metric misplacement; a boost can.
echo "--- 2. boost invariance and homogeneity (the P10_07 check)"
$PY -m pytest -q \
  tests/test_published_degree10_invariants.py::test_every_published_candidate_is_boost_invariant \
  tests/test_published_degree10_invariants.py::test_p10_07_alpha_edge_placement_is_free \
  tests/test_published_degree10_invariants.py::test_the_original_p10_07_placement_really_was_broken \
  tests/test_published_degree10_invariants.py::test_homogeneity_degree_ten \
  || fail "boost / homogeneity"

# --- 3. the positive control: the projector can return rank 3 --------------
# Without this, a projector stuck at zero would produce the same table as a
# genuine null result.
echo "--- 3. positive control: Level-A representatives must span Q10"
$PY scripts/positive_control_degree10_quotient.py || fail "positive control"

# --- 4. committed projection artifact --------------------------------------
echo "--- 4. committed projection vectors and Q10 rank"
$PY - <<'PYEOF' || fail "artifact check"
import json, pathlib, sys
p = pathlib.Path("results/intrinsic_candidates/published_degree10_map.json")
if not p.exists():
    print("  no artifact committed yet"); sys.exit(0)
d = json.loads(p.read_text())
rank = d["Q10_rank_from_implemented_published"]
print(f"  schema {d.get('schema')}  candidates {len(d.get('implemented', []))}"
      f"  primes {sorted(d['per_prime'])}")
print(f"  consistent={d['consistent']}  Q10 rank={rank}")
assert d["consistent"] is True, "rank disagrees between primes"
for prime, rec in sorted(d["per_prime"].items()):
    assert rec["dim_Q10"] == 3, f"dim Q10 != 3 at {prime}"
    bad = [k for k, v in rec["projections"].items() if v["status"] != "solved"]
    assert not bad, f"unsolved at {prime}: {bad} -- rank would be uninformative"
print("  every projection solved; a zero is a statement about the candidate")
PYEOF

# --- 4b. Level-B basis, minimality, holdout agreement ----------------------
echo "--- 4b. Level-B basis, removal minimality, holdout agreement"
$PY - <<'PYEOF' || fail "Level-B basis check"
import json, pathlib, sys
root = pathlib.Path(".")
b = root / "results/intrinsic_candidates/intrinsic_Q10_levelB_basis.json"
m = root / "results/intrinsic_candidates/published_degree10_map.json"
c = root / "results/intrinsic_candidates/degree10_positive_control.json"
if not b.exists():
    print("  no basis artifact committed yet"); sys.exit(0)
basis = json.loads(b.read_text())
amap = json.loads(m.read_text())
ctrl = json.loads(c.read_text())

triple = basis["preferred_basis"]
assert triple and len(triple) == 3, "basis is not three elements"
print(f"  preferred basis: {triple}")
assert basis["triples_agree_across_primes"] is True, \
    "the independent-triple set disagrees between primes"

# removal minimality: forced members must appear in every independent triple
for prime, rec in sorted(basis["per_prime"].items()):
    assert rec["rank"] == 3, f"Q10 rank {rec['rank']} at {prime}, expected 3"
    forced = rec["forced_members"]
    for f in forced:
        assert all(f in t for t in rec["independent_triples"]), \
            f"{f} recorded as forced but is missing from some triple"
    print(f"  prime {prime}: rank 3/3, forced {forced}, "
          f"{len(rec['independent_triples'])} independent triples")

# fit / holdout agreement on the rank
roles = {p: amap["per_prime"][p]["role"] for p in amap["per_prime"]}
ranks = {r: {amap["per_prime"][p]["Q10_rank_from_published"]
             for p in roles if roles[p] == r} for r in set(roles.values())}
print(f"  ranks by role: {ranks}")
assert all(s == {3} for s in ranks.values()), "fit/holdout rank disagreement"

# Level-A in Level-B must round-trip exactly
labels = ["Q10_A", "Q10_B", "Q10_C"]
for prime, ent in sorted(basis.get("levelA_in_levelB", {}).items()):
    p = int(prime)
    rows = [amap["per_prime"][prime]["projections"][n]["quotient_vector"]
            for n in triple]
    for lbl, coef in sorted(ent.items()):
        target = ctrl["per_prime"][prime]["vectors"][labels.index(lbl)]
        recon = [sum(coef[j] * rows[j][i] for j in range(3)) % p
                 for i in range(3)]
        assert recon == [t % p for t in target], \
            f"{lbl} does not reconstruct at {prime}"
    print(f"  prime {prime}: Level-A -> Level-B round-trips for all three")
print("  ambiguity-robust: " + str(basis["ambiguity_robust"]))
PYEOF

# --- 5. mutation cover ------------------------------------------------------
echo "--- 5. bracket / index / normalisation mutation tests"
$PY -m pytest -q \
  tests/test_published_degree10_invariants.py::test_p10_05_and_p10_08_black_brackets_are_not_vacuous \
  tests/test_published_degree10_invariants.py::test_p10_05_and_p10_08_index_mutation_is_detected \
  tests/test_published_degree10_invariants.py::test_p10_04_red_stage_is_not_vacuous_and_readings_are_distinguished \
  tests/test_index_symmetry_ops.py \
  || fail "mutation tests"

# --- 6. modular overflow guards --------------------------------------------
echo "--- 6. modular overflow guards"
$PY -m pytest -q tests/test_modular_overflow_guards.py || fail "overflow guards"

if [ "$MODE" = "FULL" ]; then
  echo "--- 7. FULL: re-evaluate the projection (hours; checkpointed)"
  echo "    restrict with SDINV_PRIMES=32749,32717 for a shorter pass"
  $PY scripts/project_published_degree10_ckpt.py || fail "projection"
fi

echo
echo "=== $MODE reproduction PASSED ==="
