"""The assembled interacting flow equations, with rational coefficients.

`results/stress_flow/interacting_flow_equations.json` is the rationally
reconstructed coefficient system for

    dV/dlambda = f(T, lambda)

through five-form degree 12, validated on an independent holdout prime.

The key quantity is the NEW FORCING span: the directions the flow can
genuinely create, computed with `Tr(tau)` excluded because that generator
only rescales interaction coordinates already present (it is the homogeneity
operator, `Tr(tau)[V_d] = 10*(d-2)*V_d`).

The obstruction is the quotient of the full space by the new-forcing span:

    degree    4    6    8   10   12
    full      1    2    7   14   72
    forcing   1    1    3    5   21
    quotient  0    1    4    9   51

Degree 4 closes. Degree 6 leaves exactly one direction unforced, and that
direction is the intrinsic sextic quotient class K6 -- which is the same
conclusion reached by the completely independent closure computation in
tests/test_stress_flow_closure.py. Two different code paths, same answer.
"""

import json
import os
from pathlib import Path
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

ROOT = Path(__file__).resolve().parents[1]
EQUATIONS = ROOT / "results" / "stress_flow" / "interacting_flow_equations.json"

FULL = {"4": 1, "6": 2, "8": 7, "10": 14, "12": 72}
NEW_FORCING = {"4": 1, "6": 1, "8": 3, "10": 5, "12": 21}
QUOTIENT = {"4": 0, "6": 1, "8": 4, "10": 9, "12": 51}

FIT_PRIMES = [32749, 32719, 32693, 32771, 32713]
HOLDOUT_PRIME = 32717


def _load():
    assert EQUATIONS.exists(), (
        "run scripts/assemble_interacting_stress_adapted.py first")
    with EQUATIONS.open() as stream:
        return json.load(stream)


def test_flow_equations_are_rationally_reconstructed_and_validated():
    d = _load()
    validation = d["exact_validation"]
    assert validation["all_modular_and_rational_holdouts_passed"] is True
    assert validation["fit_primes"] == FIT_PRIMES
    assert validation["independent_validation_prime"] == HOLDOUT_PRIME
    assert HOLDOUT_PRIME not in validation["fit_primes"], (
        "the validation prime must be independent of the fit set")
    assert validation["degree12_sha256"] == (
        "9a784dc56a2bc8186a4abb59e9177051be05640e95b4d7d24a4538bb45335113")


def test_new_forcing_and_obstruction_dimensions():
    d = _load()
    assert d["basis_dimensions"] == FULL
    assert d["new_forcing_dimension_by_degree"] == NEW_FORCING
    assert d["new_forcing_quotient_dimension_by_degree"] == QUOTIENT
    for degree in FULL:
        assert NEW_FORCING[degree] + QUOTIENT[degree] == FULL[degree]


def test_degree_four_closes_but_nothing_above_it_does():
    d = _load()
    quotient = d["new_forcing_quotient_dimension_by_degree"]
    assert quotient["4"] == 0, "degree 4 should be fully forced"
    for degree in ("6", "8", "10", "12"):
        assert quotient[degree] > 0, (
            f"degree {degree} unexpectedly closed; the obstruction would "
            f"then be absent and the classification claim wrong")


def test_sextic_obstruction_is_exactly_one_dimensional():
    """The K6 statement, read off the assembled equations."""
    d = _load()
    assert d["new_forcing_dimension_by_degree"]["6"] == 1
    assert d["new_forcing_quotient_dimension_by_degree"]["6"] == 1


def test_tr_tau_is_excluded_from_forcing_for_the_stated_reason():
    d = _load()
    assert "Tr(tau)" in d["new_forcing_definition"]
    assert "already present" in d["new_forcing_definition"]
    assert d["flow_generator_count"] == 18
    assert d["target_count"] == 192


def test_coefficient_count_is_not_confused_with_the_invariant_count():
    """96 truncated coefficients vs 81 cumulative generic functional rank."""
    d = _load()
    assert d["coefficient_direction_count"] == sum(FULL.values()) == 96
    assert d["generic_functional_rank"] == 81
    assert "81" in d["coordinate_count_explanation"]
