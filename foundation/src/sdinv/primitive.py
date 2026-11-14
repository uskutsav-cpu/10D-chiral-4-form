"""Homogeneous polynomial evidence, separate from functional Jacobian rank.

The Hilbert coefficients bound vector-space dimensions. Multiplicative
factor exponents are not automatically counts of generators in degrees
where generators and relations can coexist.
"""

import numpy as np

from .modp import RankSieve

HILBERT_DIMENSIONS = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72, 14: 247, 16: 1364}
HILBERT_FACTOR_EXPONENTS = {4: 1, 6: 2, 8: 6, 10: 12, 12: 62, 14: 221}


def product_monomials(generators, degree):
    """Enumerate distinct lower-generator monomials of a prescribed degree.

    Products may satisfy relations; callers must rank their evaluations.
    The returned list is not a claim of product-space dimension.
    """
    items = sorted((name, d) for name, d in generators.items() if 0 < d < degree)
    result = []

    def visit(start, remaining, factors):
        if remaining == 0:
            if len(factors) > 1:
                result.append(tuple(factors))
            return
        for k in range(start, len(items)):
            name, d = items[k]
            if d <= remaining:
                visit(k, remaining - d, factors + [name])

    visit(0, degree, [])
    return result


class PolynomialSieve:
    """Rank value vectors, or gradients concatenated at several points.

    A pivot proves polynomial linear independence. A failed pivot at finitely
    many points alone does not prove an identity. Gradient pivots also prove
    polynomial independence; characteristic must exceed the positive degree.
    """

    def __init__(self, degree, ncols, prime, method="values"):
        if degree <= 0 or degree >= prime:
            raise ValueError("require 0 < degree < prime")
        if method not in {"values", "stacked_gradients"}:
            raise ValueError("unknown polynomial evidence method")
        self.degree, self.method = degree, method
        self.sieve = RankSieve(ncols, prime)
        self.product_ids, self.primitive_ids = [], []

    def add(self, invariant_id, vector, *, product=False):
        if product and self.primitive_ids:
            raise ValueError("seed the entire product span before candidates")
        row = np.asarray(vector)
        if row.shape != (self.sieve.ncols,):
            raise ValueError("polynomial evidence vector has wrong dimension")
        if self.sieve.add(row):
            (self.product_ids if product else self.primitive_ids).append(invariant_id)
            return True
        return False

    @property
    def rank(self):
        return self.sieve.rank

    @property
    def complete(self):
        return self.rank == HILBERT_DIMENSIONS.get(self.degree)
