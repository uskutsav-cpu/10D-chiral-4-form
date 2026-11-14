"""The six displayed tensor expressions (4.12)--(4.16) of arXiv:2509.14350v2.

Normalized symmetrization is applied after the five-index antisymmetrization
already present in N1050. We implement the displayed expressions literally;
no additional trace subtraction is silently inserted. Lower-degree products
must be included when comparing different primitive representatives.
"""

import itertools

import numpy as np

from .modp import inv, mod_einsum
from .stress import (five_form_moment, matrix_trace_power, composite_n1050,
                     composite_n4125, _raise_axes, _antisymmetrize_axes)

SOURCE = "https://arxiv.org/html/2509.14350v2#S4.SS1.SSS3"
BASIS = [{"id": f"T8_{i}", "equation": eq} for i, eq in
         enumerate(["4.12", "4.12", "4.13", "4.14", "4.15", "4.16"], 1)]
ALTERNATIVE_BASIS = [{"id": f"H8_{i}", "equation": eq} for i, eq in
                     enumerate(["4.18", "4.19 (definition); 4.20 (reduction)", "4.21", "4.22", "4.23", "4.12 (second)"], 1)]
SELECTED_BASIS = BASIS[:5] + [ALTERNATIVE_BASIS[0]]


def _contract_four(tensor, subscripts, prime):
    """Explicit one-metric-per-edge contraction of four covariant tensors."""
    seen, operands = set(), []
    for term in subscripts.split(","):
        raised = [i for i, index in enumerate(term) if index in seen]
        seen.update(term)
        operands.append(_raise_axes(tensor, raised, prime))
    return int(mod_einsum(subscripts + "->", operands, prime))


def evaluate_selected_basis(form, prime):
    """A verified selection across the paper's two displayed lists.

    T8_1,...,T8_5,H8_1 are used explicitly, rather than assuming either
    displayed list supplies a primitive complement under a literal reading.
    """
    first = evaluate_basis(form, prime)
    red = _antisymmetrize_axes(composite_n1050(form, prime), (3, 4, 5), prime)
    hatted1 = -_contract_four(red, "abcdrs,abcduv,urefgh,sefvgh", prime) % prime
    return first[:5] + [hatted1]


def evaluate_alternative_basis(form, prime):
    """Five hatted N1050 expressions plus the second invariant in (4.12)."""
    n = composite_n1050(form, prime)
    red = _antisymmetrize_axes(n, (3, 4, 5), prime)
    words = ["abcdrs,abcduv,urefgh,sefvgh", "abcdrs,abcduv,refugh,segfhv",
             "abcrst,abcuvw,rsefgu,tvfegw", "abcrst,abcuvw,uresfg,tvfweg",
             "urabcd,sabvcd,suefgh,vefrgh"]
    result = [_contract_four(red, word, prime) for word in words]
    result[0] = -result[0] % prime
    _, mixed = five_form_moment(form, prime)
    n4125 = composite_n4125(form, prime)
    result.append(int(mod_einsum("ad,be,cf,defabc->", [mixed, mixed, mixed,
                      _raise_axes(n4125, (3, 4, 5), prime)], prime)))
    return result


def evaluate_basis(form, prime):
    lower, mixed = five_form_moment(form, prime)
    upper = _raise_axes(lower, (0, 1), prime)
    n1050 = composite_n1050(form, prime)
    n4125 = composite_n4125(form, prime)
    q = mod_einsum("abcdmn,abcdrl->mnrl", [n1050, _raise_axes(n1050, (0, 1, 2, 3), prime)], prime)
    sym = np.zeros_like(q)
    for perm in itertools.permutations(range(4)):
        sym = (sym + q.transpose(perm)) % prime
    sym = sym * inv(24, prime) % prime
    q_mixed = _raise_axes(q, (2, 3), prime)
    return [
        matrix_trace_power(mixed, 4, prime),
        int(mod_einsum("ad,be,cf,defabc->", [mixed, mixed, mixed, _raise_axes(n4125, (3, 4, 5), prime)], prime)),
        int(mod_einsum("mn,mnrl,rl->", [upper, q, upper], prime)),
        int(mod_einsum("mnrl,mnrl->", [sym, _raise_axes(q, (0, 1, 2, 3), prime)], prime)),
        int((mod_einsum("am,bn,mnab->", [mixed, mixed, q_mixed], prime)
             - mod_einsum("bm,an,mnab->", [mixed, mixed, q_mixed], prime)) * inv(2, prime) % prime),
        int(mod_einsum("am,bn,mnpqrs,pqrsab->", [mixed, mixed, n1050, _raise_axes(n4125, range(6), prime)], prime)),
    ]
