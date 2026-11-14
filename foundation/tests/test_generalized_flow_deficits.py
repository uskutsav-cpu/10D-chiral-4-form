"""The degree-10 and degree-12 generalized-flow deficits, pinned exactly.

Phase 1 result. Starting from the free seed with K6 and the four degree-8
directions adjoined, the closure is (1,2,7,11,68) against full (1,2,7,14,72).
The scan over every basis direction at each degree finds exactly the deficit
many raisers:

    degree 10:  I10_6, I10_7, I10_12                    (deficit 3)
    degree 12:  I12_59, I12_60, I12_61, I12_62          (deficit 4)

These are SEED directions, not flow generators -- see docs/CLAIM_LEDGER.md
C-MIN-03 -- and they are graph labels, not intrinsic tensors, which is why the
Phase 1 gate is not yet met.
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
GENERALIZED = ROOT / "results" / "generalized_flow"
PRIMES = (32749, 32719, 32717, 32693)

FULL = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
DEGREE10_MISSING = ["I10_6", "I10_7", "I10_12"]
DEGREE12_MISSING = ["I12_59", "I12_60", "I12_61", "I12_62"]


def _certificates():
    out = []
    for prime in PRIMES:
        path = CERTIFICATE_DIR / f"interacting_degree12_{prime}.json"
        if path.exists():
            with path.open() as stream:
                out.append((prime, json.load(stream)))
    assert len(out) >= 3
    return out


def _seed(deg10=(), deg12=()):
    seed = {k: list(v) for k, v in SEED8.items()}
    if deg10:
        seed[10] = list(deg10)
    if deg12:
        seed[12] = list(deg12)
    return seed


def test_degree10_deficit_is_closed_by_exactly_these_three():
    for prime, certificate in _certificates():
        base, _ = closure(certificate, _seed(), prime)
        assert base[10] == 11, f"prime {prime}: base degree-10 closure moved"
        dims, _ = closure(certificate, _seed(DEGREE10_MISSING), prime)
        assert dims[10] == FULL[10], f"prime {prime}: degree 10 not closed"


def test_degree10_each_direction_is_necessary():
    for prime, certificate in _certificates():
        for omitted in DEGREE10_MISSING:
            kept = [d for d in DEGREE10_MISSING if d != omitted]
            dims, _ = closure(certificate, _seed(kept), prime)
            assert dims[10] < FULL[10], (
                f"prime {prime}: {omitted} was supposed to be necessary")


def test_degree12_deficit_is_closed_by_exactly_these_four():
    for prime, certificate in _certificates():
        base, _ = closure(certificate, _seed(DEGREE10_MISSING), prime)
        assert base[12] == 68, f"prime {prime}: base degree-12 closure moved"
        dims, _ = closure(
            certificate, _seed(DEGREE10_MISSING, DEGREE12_MISSING), prime)
        assert dims[12] == FULL[12], f"prime {prime}: degree 12 not closed"


def test_degree12_each_direction_is_necessary():
    for prime, certificate in _certificates():
        for omitted in DEGREE12_MISSING:
            kept = [d for d in DEGREE12_MISSING if d != omitted]
            dims, _ = closure(
                certificate, _seed(DEGREE10_MISSING, kept), prime)
            assert dims[12] < FULL[12], (
                f"prime {prime}: {omitted} was supposed to be necessary")


def test_the_two_functional_dependencies_are_among_the_unreachable():
    """CJ-01 evidence, pinned so it cannot silently drift.

    This is a CONJECTURE (2 of 2, converse false), not a theorem. The test
    exists to detect a change in the underlying data, not to assert a law.
    """
    with (ROOT / "results" / "10d_order12.json").open() as stream:
        atlas = json.load(stream)
    dependencies = atlas["discovery"]["functional_dependencies"]
    assert dependencies == ["I12_61", "I12_62"]
    assert set(dependencies) <= set(DEGREE12_MISSING)
    # the converse must remain false, or CJ-01 changes character entirely
    assert set(DEGREE12_MISSING) - set(dependencies) == {"I12_59", "I12_60"}


def test_committed_artifacts_match_the_pinned_sets():
    for degree, expected in ((10, DEGREE10_MISSING), (12, DEGREE12_MISSING)):
        path = GENERALIZED / f"degree{degree}_missing_directions.json"
        if not path.exists():
            continue
        with path.open() as stream:
            payload = json.load(stream)
        assert payload["missing_directions"] == expected
        assert payload["degree_closed"] is True
        assert "SEED directions" in payload["caveat"], (
            "the seed-vs-generator caveat must stay in the artifact")
