"""Auditable exact prime-field algebra; no tolerance and no machine-int overflow."""
from __future__ import annotations
from functools import lru_cache
from numbers import Integral
from typing import Iterable, Sequence

Matrix = list[list[int]]

def integer(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("exact integer required (floats and booleans are rejected)")
    return int(value)

@lru_cache(maxsize=128)
def _check_prime_field_modulus(p: int) -> None:
    p = integer(p)
    if not 2 < p < 2**64:
        raise ValueError("use an odd prime smaller than 2**64")
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if p == q:
            return
        if p % q == 0:
            raise ValueError(f"composite modulus: {p}")
    if p % 2 == 0:
        raise ValueError(f"composite modulus: {p}")
    d, s = p - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    # Deterministic Miller--Rabin bases for unsigned 64-bit integers.
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % p == 0:
            continue
        x = pow(a, d, p)
        if x in (1, p - 1):
            continue
        for _ in range(s - 1):
            x = x * x % p
            if x == p - 1:
                break
        else:
            raise ValueError(f"composite modulus: {p}")


def _rectangular(matrix: Sequence[Sequence[int]], ncols: int | None = None) -> tuple[int,int]:
    rows = len(matrix)
    cols = len(matrix[0]) if rows else (0 if ncols is None else integer(ncols))
    if cols < 0 or (ncols is not None and ncols != cols):
        raise ValueError("column count mismatch")
    if any(len(row) != cols for row in matrix):
        raise ValueError("matrix must be rectangular")
    return rows, cols


def rref(matrix: Sequence[Sequence[int]], p: int, *, ncols: int | None = None) -> tuple[Matrix,list[int]]:
    _check_prime_field_modulus(p)
    rows, cols = _rectangular(matrix, ncols)
    a = [[integer(x) % p for x in row] for row in matrix]
    pivots = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [x * inv % p for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                factor = a[i][c]
                a[i] = [(x - factor*y) % p for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return a, pivots


def matrix_rank(matrix: Sequence[Sequence[int]], p: int, *, ncols: int | None = None) -> int:
    return len(rref(matrix, p, ncols=ncols)[1])


def transpose(matrix: Sequence[Sequence[int]]) -> Matrix:
    rows, cols = _rectangular(matrix)
    return [[integer(matrix[i][j]) for i in range(rows)] for j in range(cols)]


def nullspace(matrix: Sequence[Sequence[int]], p: int, *, ncols: int | None = None) -> Matrix:
    rr, pivots = rref(matrix, p, ncols=ncols)
    _, cols = _rectangular(matrix, ncols)
    out = []
    for free in (c for c in range(cols) if c not in pivots):
        v = [0]*cols
        v[free] = 1
        for row, pivot in enumerate(pivots):
            v[pivot] = -rr[row][free] % p
        out.append(v)
    return out


def annihilator_basis(row_generators: Sequence[Sequence[int]], p: int, *, ncols: int | None = None) -> Matrix:
    return nullspace(row_generators, p, ncols=ncols)


def matvec(matrix: Sequence[Sequence[int]], vector: Sequence[int], p: int) -> list[int]:
    _check_prime_field_modulus(p)
    _, cols = _rectangular(matrix, len(vector))
    v = [integer(x) for x in vector]
    return [sum(integer(a)*b for a,b in zip(row,v)) % p for row in matrix]


def matmul(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]], p: int) -> Matrix:
    _check_prime_field_modulus(p)
    _, n = _rectangular(a)
    m, k = _rectangular(b)
    if n != m:
        raise ValueError("matrix dimensions do not align")
    return [[sum(integer(row[j])*integer(b[j][c]) for j in range(n)) % p
             for c in range(k)] for row in a]


def span_contains(rows: Sequence[Sequence[int]], vector: Sequence[int], p: int) -> bool:
    _check_prime_field_modulus(p)
    _rectangular(rows, len(vector))
    v = [integer(x) for x in vector]
    return matrix_rank([*map(list,rows),v],p) == matrix_rank(rows,p)


def normalize_projective(vector: Iterable[int], p: int) -> list[int]:
    _check_prime_field_modulus(p)
    v = [integer(x)%p for x in vector]
    first = next((x for x in v if x), None)
    return v if first is None else [x*pow(first,-1,p)%p for x in v]


def independent_rows(rows: Sequence[Sequence[int]], p: int) -> list[int]:
    if not rows:
        return []
    return rref(transpose(rows),p)[1]


def determinant(matrix: Sequence[Sequence[int]], p: int) -> int:
    _check_prime_field_modulus(p)
    n, m = _rectangular(matrix)
    if n != m:
        raise ValueError("determinant requires a square matrix")
    a = [[integer(x)%p for x in row] for row in matrix]
    det = 1
    for c in range(n):
        k = next((i for i in range(c,n) if a[i][c]),None)
        if k is None:
            return 0
        if k != c:
            a[c],a[k] = a[k],a[c]
            det = -det
        pivot = a[c][c]
        det = det*pivot%p
        inv = pow(pivot,-1,p)
        for i in range(c+1,n):
            f = a[i][c]*inv%p
            for j in range(c+1,n):
                a[i][j] = (a[i][j]-f*a[c][j])%p
    return det%p


def solve_many(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]], p: int) -> Matrix:
    """Unique exact solution of A X = B; reject inconsistent/underdetermined fits."""
    m,n = _rectangular(a)
    mb,k = _rectangular(b)
    if m != mb or not m or not n:
        raise ValueError("empty or incompatible linear system")
    rr,piv = rref([list(x)+list(y) for x,y in zip(a,b)],p)
    if any(c >= n for c in piv):
        raise ValueError("inconsistent exact linear system")
    if len(piv) < n:
        raise ValueError("underdetermined exact linear system")
    return [row[n:n+k] for row in rr[:n]]
