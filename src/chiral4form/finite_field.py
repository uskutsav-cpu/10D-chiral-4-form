"""Small exact linear-algebra helpers over prime fields.

The research-scale pipelines can swap in faster backends.  These routines are
kept intentionally small and auditable for certificates, unit tests, and
cross-checks.
"""

from __future__ import annotations

from typing import Iterable, Sequence

Matrix = list[list[int]]


def _check_prime_field_modulus(p: int) -> None:
    if p <= 2:
        raise ValueError("p must be an odd prime used as a field modulus")


def _rectangular(matrix: Sequence[Sequence[int]]) -> tuple[int, int]:
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    if any(len(row) != cols for row in matrix):
        raise ValueError("matrix must be rectangular")
    return rows, cols


def rref(matrix: Sequence[Sequence[int]], p: int) -> tuple[Matrix, list[int]]:
    """Return reduced row-echelon form and pivot columns over F_p."""

    _check_prime_field_modulus(p)
    rows, cols = _rectangular(matrix)
    a = [[int(x) % p for x in row] for row in matrix]
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [(x * inv) % p for x in a[r]]
        for i in range(rows):
            if i == r:
                continue
            factor = a[i][c] % p
            if factor:
                a[i] = [(x - factor * y) % p for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return a, pivots


def matrix_rank(matrix: Sequence[Sequence[int]], p: int) -> int:
    """Exact matrix rank over F_p."""

    return len(rref(matrix, p)[1])


def transpose(matrix: Sequence[Sequence[int]]) -> Matrix:
    rows, cols = _rectangular(matrix)
    if rows == 0:
        return []
    return [[int(matrix[i][j]) for i in range(rows)] for j in range(cols)]


def nullspace(matrix: Sequence[Sequence[int]], p: int) -> Matrix:
    """Basis vectors for the right nullspace of a matrix over F_p."""

    rr, pivots = rref(matrix, p)
    _, cols = _rectangular(matrix)
    pivot_set = set(pivots)
    free = [c for c in range(cols) if c not in pivot_set]
    out: Matrix = []
    for free_col in free:
        v = [0] * cols
        v[free_col] = 1
        for row, pivot_col in enumerate(pivots):
            v[pivot_col] = (-rr[row][free_col]) % p
        out.append(v)
    return out


def annihilator_basis(row_generators: Sequence[Sequence[int]], p: int) -> Matrix:
    """Return linear functionals annihilating a row-generated subspace.

    A functional `ell` is represented by a column-coordinate vector satisfying
    `R @ ell = 0`, so this is the right nullspace of the generator matrix.
    """

    return nullspace(row_generators, p)


def matvec(matrix: Sequence[Sequence[int]], vector: Sequence[int], p: int) -> list[int]:
    _, cols = _rectangular(matrix)
    if len(vector) != cols:
        raise ValueError("dimension mismatch")
    return [sum(int(a) * int(b) for a, b in zip(row, vector)) % p for row in matrix]


def span_contains(
    row_generators: Sequence[Sequence[int]], vector: Sequence[int], p: int
) -> bool:
    """Exact membership test for a row span."""

    rows, cols = _rectangular(row_generators)
    if rows == 0:
        return all(int(x) % p == 0 for x in vector)
    if len(vector) != cols:
        raise ValueError("dimension mismatch")
    before = matrix_rank(row_generators, p)
    after = matrix_rank([*map(list, row_generators), list(vector)], p)
    return before == after


def normalize_projective(vector: Iterable[int], p: int) -> list[int]:
    """Canonical projective normalization: first nonzero entry is one."""

    _check_prime_field_modulus(p)
    v = [int(x) % p for x in vector]
    first = next((x for x in v if x), None)
    if first is None:
        return v
    inv = pow(first, -1, p)
    return [(x * inv) % p for x in v]
