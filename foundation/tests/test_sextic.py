"""Exact intrinsic identification of the two sextic directions."""

from fractions import Fraction
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.exactmap import sample_selfdual_five_form  # noqa: E402
from sdinv.invariant_registry import load_verified_registry  # noqa: E402
from sdinv.sextic import (  # noqa: E402
    INTRINSIC_TO_REGISTRY,
    REGISTRY_TO_INTRINSIC,
    paper_i6_1,
    paper_i6_2,
    paper_i6_2_direct,
    paper_i6_2_graph_expansion,
    registry_to_intrinsic_coefficients,
    sextic_quotient_coordinate,
)


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PRIMES = (32749, 32719, 32693)
SAMPLE_SEEDS = (20260921, 20260922, 20260923, 20260924)


def _residue(fraction, prime):
    fraction = Fraction(fraction)
    return (
        fraction.numerator
        * pow(fraction.denominator, -1, prime)
        % prime
    )


def test_intrinsic_sextic_change_of_basis_is_exact_and_invertible():
    product = tuple(tuple(
        sum(
            INTRINSIC_TO_REGISTRY[row][middle]
            * REGISTRY_TO_INTRINSIC[middle][column]
            for middle in range(2)
        )
        for column in range(2)
    ) for row in range(2))
    assert product == ((Fraction(1), Fraction(0)),
                       (Fraction(0), Fraction(1)))
    assert (
        INTRINSIC_TO_REGISTRY[0][0]
        * INTRINSIC_TO_REGISTRY[1][1]
        - INTRINSIC_TO_REGISTRY[0][1]
        * INTRINSIC_TO_REGISTRY[1][0]
    ) == Fraction(32, 125)
    assert registry_to_intrinsic_coefficients((7, 11)) == (
        Fraction(7) * Fraction(3, 32) + Fraction(11, 288),
        Fraction(1375, 3),
    )
    assert sextic_quotient_coordinate((7, 11)) == Fraction(1375, 3)


def test_1050_cubic_has_a_reproducible_contraction_graph_expansion():
    expansion = paper_i6_2_graph_expansion()
    assert len(expansion) == 8
    assert sum(
        (item["coefficient"] for item in expansion), Fraction()
    ) == Fraction(11, 250)
    assert len({item["graph"] for item in expansion}) == len(expansion)


@pytest.mark.parametrize("prime", PRIMES)
def test_intrinsic_sextics_match_registry_on_four_samples(prime):
    registry = load_verified_registry(ROOT)
    for seed in SAMPLE_SEEDS:
        five_form = sample_selfdual_five_form(seed, prime)
        graph_i61, graph_i62 = registry.evaluate_degree(
            6, five_form, prime)
        intrinsic_i61 = paper_i6_1(five_form, prime)
        intrinsic_i62 = paper_i6_2(five_form, prime)
        assert intrinsic_i61 == (
            _residue(Fraction(32, 3), prime) * graph_i61
        ) % prime
        assert intrinsic_i62 == (
            _residue(Fraction(-1, 1125), prime) * graph_i61
            + _residue(Fraction(3, 125), prime) * graph_i62
        ) % prime


@pytest.mark.parametrize("prime", PRIMES)
def test_1050_cubic_graph_expansion_matches_direct_tensor_formula(prime):
    five_form = sample_selfdual_five_form(20260921, prime)
    assert paper_i6_2(five_form, prime) == paper_i6_2_direct(
        five_form, prime)

