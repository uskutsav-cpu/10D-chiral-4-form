"""Closed families of the pure stress flow, pinned exactly.

The static free-stress span and the set a FLOW can reach are different
objects, and conflating them is the main interpretive trap in this project.
This module pins both, so a future change that silently merges them fails.

Reached dimensions, identical over every available prime:

    degree   full   static span   free-seed closure
      4        1         1               1
      6        2         1               1
      8        7         2               3
     10       14         2              11
     12       72         4              67

Two independent facts live in that table:

1. The closure is much LARGER than the static span at degrees 8, 10, 12.
   A no-go argument built on the static span alone would be wrong.
2. At degree 6 the closure is still 1, not 2. The missing direction is
   exactly the intrinsic quotient class K6. See test_k6_flow_role.py.
"""

import json
import os
from pathlib import Path
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from stress_flow_closure import closure  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_DIR = ROOT / "results" / "stress_flow" / "certificates"
PRIMES = (32749, 32719, 32717, 32693)

FULL_DIMENSION = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}
STATIC_STRESS_SPAN = {4: 1, 6: 1, 8: 2, 10: 2, 12: 4}
FREE_SEED_CLOSURE = {4: 1, 6: 1, 8: 3, 10: 11, 12: 67}
K6_SEEDED_CLOSURE = {4: 1, 6: 2, 8: 3, 10: 11, 12: 68}


def _certificates():
    out = []
    for prime in PRIMES:
        path = CERTIFICATE_DIR / f"interacting_degree12_{prime}.json"
        if path.exists():
            with path.open() as stream:
                out.append((prime, json.load(stream)))
    assert len(out) >= 3, "need at least three interacting certificates"
    return out


def test_free_seed_closure_is_exact_and_prime_independent():
    for prime, certificate in _certificates():
        dims, _ = closure(certificate, {4: ["I4_1"]}, prime)
        assert dims == FREE_SEED_CLOSURE, f"prime {prime}: {dims}"


def test_k6_seed_adds_exactly_the_sextic_and_one_degree12_direction():
    for prime, certificate in _certificates():
        dims, _ = closure(
            certificate, {4: ["I4_1"], 6: ["I6_2"]}, prime)
        assert dims == K6_SEEDED_CLOSURE, f"prime {prime}: {dims}"
        assert dims[6] - FREE_SEED_CLOSURE[6] == 1
        assert dims[12] - FREE_SEED_CLOSURE[12] == 1


def test_closure_strictly_exceeds_the_static_span():
    """Guards against restating the static result as a dynamical no-go."""
    for degree in (8, 10, 12):
        assert FREE_SEED_CLOSURE[degree] > STATIC_STRESS_SPAN[degree], (
            "the flow closure must not be equated with the static span")


def test_closure_is_a_proper_subspace_at_every_degree_above_four():
    """There IS a genuine obstruction -- the closure is not everything."""
    for degree in (6, 8, 10, 12):
        assert FREE_SEED_CLOSURE[degree] < FULL_DIMENSION[degree]


def test_seeding_j6_changes_nothing():
    """J6 already lies in the free-seed closure, so seeding it is inert."""
    for prime, certificate in _certificates():
        dims, _ = closure(
            certificate, {4: ["I4_1"], 6: ["I6_1"]}, prime)
        assert dims == FREE_SEED_CLOSURE, f"prime {prime}: {dims}"
