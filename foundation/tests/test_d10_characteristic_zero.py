"""Regression test: the exact characteristic-zero D10/Q10 certificate."""
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
D10 = ROOT / "results/stress_flow/D10_characteristic_zero.json"
Q10 = ROOT / "results/stress_flow/Q10_characteristic_zero.json"


def load(p):
    if not p.exists():
        pytest.skip(f"{p.name} absent")
    return json.loads(p.read_text())


def test_lift_was_validated_against_a_held_out_prime():
    """The lift is only trustworthy because a prime was withheld from it."""
    d = load(D10)
    lift = d["lift"]
    assert lift["holdout_prime"] not in lift["fitting_primes"]
    assert lift["failed"] == [], lift["failed"]
    assert lift["holdout_mismatches"] == [], lift["holdout_mismatches"]


def test_minor_is_nonzero_over_the_integers():
    d = load(D10)
    m = d["lower_bound_certificate"]
    assert m is not None and m["nonzero"] is True
    assert int(m["integer_determinant"]) != 0
    assert m["size"] == d["D10_dim_over_Q"]


def test_exact_dimensions_are_consistent():
    d, q = load(D10), load(Q10)
    assert d["dims_over_Q"]["10"] == d["D10_dim_over_Q"]
    assert q["A10_dim_over_Q"] - q["D10_dim_over_Q"] == q["Q10_dim_over_Q"]


def test_equality_wording_requires_a_settled_certificate():
    """`= 3` may only be claimed when the certificate says settled."""
    d, q = load(D10), load(Q10)
    if not d["settled"]:
        assert q["status"] != "exact"
        assert "<=" in q["permitted_wording"]
    else:
        assert q["status"] == "exact"
