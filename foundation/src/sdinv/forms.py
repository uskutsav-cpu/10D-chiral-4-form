"""Antisymmetric p-forms in d dimensions, over F_p.

A form is stored two ways:
  - "sparse": a vector of length C(d,p) indexed by SORTED index tuples.
    This is the honest count of independent components (20 for a 3-form in
    6D, 252 for a 5-form in 10D).
  - "dense": a full (d,)*p numpy array, needed for einsum contraction.

The Hodge dual is computed on the sparse representation via index
complement -- O(C(d,p)) work. Never build the d^d epsilon tensor.
"""

import itertools
import numpy as np

from .modp import P, inv


def basis_tuples(d, p):
    """All sorted index tuples -- one per independent component."""
    return list(itertools.combinations(range(d), p))


def metric_signs(d, lorentzian=True):
    """Diagonal metric. Lorentzian = diag(-1, +1, ..., +1)."""
    s = np.ones(d, dtype=np.int64)
    if lorentzian:
        s[0] = -1
    return s


def perm_sign(perm):
    """Sign of a permutation given as a sequence."""
    perm = list(perm)
    n = len(perm)
    seen = [False] * n
    sign = 1
    for i in range(n):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            ln += 1
        if ln % 2 == 0:
            sign = -sign
    return sign


def to_dense(vec, d, p_deg, mod=P):
    """Expand sparse component vector to a dense antisymmetric array."""
    T = np.zeros((d,) * p_deg, dtype=np.int64)
    for k, idx in enumerate(basis_tuples(d, p_deg)):
        v = int(vec[k]) % mod
        if v == 0:
            continue
        for perm in itertools.permutations(range(p_deg)):
            sgn = perm_sign(perm)
            T[tuple(idx[perm[i]] for i in range(p_deg))] = (
                v * sgn) % mod
    return T


def random_form(d, p_deg, rng, mod=P):
    """Random sparse component vector."""
    n = len(basis_tuples(d, p_deg))
    return rng.integers(1, mod, size=n, dtype=np.int64)


def hodge_matrix(d, p_deg, lorentzian=True, mod=P):
    """Matrix of the Hodge star on sparse components.

    (*F)_{I} = (1/(d-p)!) eps_{I J} F^{J}, which on sorted multi-indices is
    just +/- the component on the complementary index set, times the metric
    signs from raising J.
    """
    src = basis_tuples(d, p_deg)
    dst = basis_tuples(d, d - p_deg)
    dst_pos = {t: i for i, t in enumerate(dst)}
    s = metric_signs(d, lorentzian)
    M = np.zeros((len(dst), len(src)), dtype=np.int64)
    for j, I in enumerate(src):
        Ic = tuple(sorted(set(range(d)) - set(I)))
        # sign of the permutation (Ic, I) relative to (0,1,...,d-1)
        perm = list(Ic) + list(I)
        sgn = perm_sign([perm.index(k) for k in range(d)])
        raise_sign = int(np.prod(s[list(I)]))
        M[dst_pos[Ic], j] = (sgn * raise_sign) % mod
    return M


def selfdual_projector(d, p_deg, lorentzian=True, mod=P):
    """P = (1 + *) / 2, valid when *^2 = +1 (true for 5-forms in 10D
    Lorentzian: (-1)^{p(d-p)} * (-1)^s = (-1)^25 * (-1)^1 = +1)."""
    assert d == 2 * p_deg, "star maps p-forms to p-forms only when d = 2p"
    H = hodge_matrix(d, p_deg, lorentzian, mod)
    n = H.shape[0]
    I = np.eye(n, dtype=np.int64)
    half = inv(2, mod)
    return ((I + H) * half) % mod


def check_star_squared(d, p_deg, lorentzian=True, mod=P):
    """Sanity check: returns the scalar c with *^2 = c * identity."""
    H = hodge_matrix(d, p_deg, lorentzian, mod)
    H2 = (H @ H) % mod
    n = H.shape[0]
    c = H2[0, 0]
    assert np.array_equal(H2, (c * np.eye(n, dtype=np.int64)) % mod), \
        "star^2 is not a multiple of the identity"
    return int(c) if c < mod // 2 else int(c) - mod
