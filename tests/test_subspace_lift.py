from fractions import Fraction

from chiral4form.finite_field import rref
from chiral4form.modular_reconstruction import reduce_fraction
from chiral4form.subspace_lift import canonical_row_spaces, lift_stable_row_space


def reduce_matrix(matrix, p):
    return [[reduce_fraction(Fraction(x), p) for x in row] for row in matrix]


def test_lift_simple_row_space_with_holdout():
    qq = [
        [1, 0, Fraction(1, 3), Fraction(-5, 7)],
        [0, 1, Fraction(2, 5), Fraction(11, 13)],
    ]
    primes = [1009, 1013, 1019, 1021]
    rows = {p: reduce_matrix(qq, p) for p in primes}
    spaces, pivots = canonical_row_spaces(rows, ncols=4, expected_rank=2)
    out = lift_stable_row_space(spaces, pivots, ncols=4, holdout_primes=[1021])
    assert out["status"] == "lifted_and_verified"
    assert out["qq_verified"] is True
    assert out["pivots"] == [0, 1]
