"""Minimality of the generalized-flow completion, proved by removal.

Through degree 8 the free-seed stress flow needs exactly five extra
directions to close: K6 and the four degree-8 complement directions
I8_3..I8_6. Minimality is asserted the only way it can be honestly asserted
-- by showing that dropping any one of them reopens a gap.

See docs/minimal_generalized_flow.md.
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

FULL = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}
REQUIRED_DEGREE8 = ("I8_3", "I8_4", "I8_5", "I8_6")
INERT_DEGREE8 = ("I8_1", "I8_2", "I4_1^2")


def _certificates():
    out = []
    for prime in PRIMES:
        path = CERTIFICATE_DIR / f"interacting_degree12_{prime}.json"
        if path.exists():
            with path.open() as stream:
                out.append((prime, json.load(stream)))
    assert len(out) >= 3
    return out


def _seed(extra_six=False, degree8=()):
    seed = {4: ["I4_1"]}
    if extra_six:
        seed[6] = ["I6_2"]
    if degree8:
        seed[8] = list(degree8)
    return seed


def test_full_five_directions_close_degrees_six_and_eight():
    for prime, certificate in _certificates():
        dims, _ = closure(
            certificate, _seed(True, REQUIRED_DEGREE8), prime)
        assert dims[6] == FULL[6], f"prime {prime}: degree 6 open"
        assert dims[8] == FULL[8], f"prime {prime}: degree 8 open"


def test_dropping_any_degree8_direction_reopens_degree_eight():
    for prime, certificate in _certificates():
        for omitted in REQUIRED_DEGREE8:
            kept = tuple(d for d in REQUIRED_DEGREE8 if d != omitted)
            dims, _ = closure(certificate, _seed(True, kept), prime)
            assert dims[8] < FULL[8], (
                f"prime {prime}: {omitted} was supposed to be necessary but "
                f"degree 8 closed without it")


def test_dropping_k6_reopens_degree_six_and_degree_twelve():
    for prime, certificate in _certificates():
        with_k6, _ = closure(
            certificate, _seed(True, REQUIRED_DEGREE8), prime)
        without, _ = closure(
            certificate, _seed(False, REQUIRED_DEGREE8), prime)
        assert without[6] < with_k6[6], "K6 not needed at degree 6"
        assert without[12] < with_k6[12], "K6 not needed at degree 12"


def test_inert_degree8_directions_are_not_part_of_a_minimal_set():
    """I8_1, I8_2, I4_1^2 already lie in the free-seed closure."""
    for prime, certificate in _certificates():
        base, _ = closure(certificate, _seed(True), prime)
        for name in INERT_DEGREE8:
            dims, _ = closure(certificate, _seed(True, (name,)), prime)
            assert dims == base, (
                f"prime {prime}: {name} was expected to be inert")


def test_degrees_ten_and_twelve_remain_open():
    """The completion through degree 8 does NOT close the higher degrees."""
    for prime, certificate in _certificates():
        dims, _ = closure(
            certificate, _seed(True, REQUIRED_DEGREE8), prime)
        assert dims[10] < FULL[10]
        assert dims[12] < FULL[12]
