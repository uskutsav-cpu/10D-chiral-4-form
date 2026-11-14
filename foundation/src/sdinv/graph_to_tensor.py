"""Deterministic translation of a contraction graph into an explicit
Einstein-index tensor formula.

A verified contraction graph already *is* a coordinate-independent Lorentz
scalar; this module makes that explicit rather than leaving it as a label.

Translation rule, unambiguous by construction:

  * every vertex v carries one self-dual five-form  F^(v)_{mu1..mu5};
  * every edge (i,j) of multiplicity m contributes m contracted index PAIRS,
    one slot on vertex i and one on vertex j;
  * every contracted pair carries exactly one inverse metric eta^{..};
  * each vertex has exactly five slots, matching valence 5.

So for a graph on n vertices the scalar is

    I  =  [ prod_v  F_{ s(v,1) ... s(v,5) } ]  x  [ prod_pairs  eta^{ s_a s_b } ]

with n*5/2 metric factors. Dummy names are assigned canonically (edge order,
then multiplicity index) so the output is reproducible.

Sign convention. Contractions of odd-rank forms can pick up an orientation
sign under graph isomorphism, so a translation may differ from the
repository's optimised evaluator by an overall sign. That sign is *measured*
and recorded, never silently absorbed.
"""

import re

import numpy as np

from .forms import metric_signs, to_dense
from .modp import P, mod_einsum

FORM_RANK = 5
DIM = 10

_EDGE = re.compile(r"(\d+)-?(\d+)\^(\d+)")


def parse_graph_label(label):
    """'n10[03^1,05^2,...]' or 'n12[0-4^2,...]' -> symmetric multiplicity matrix.

    Both label dialects appear in the committed artifacts: degree 10 uses
    concatenated single digits, degree 12 uses explicit dashes because vertex
    indices reach two digits.
    """
    head, body = label.split("[", 1)
    n = int(head[1:])
    body = body.rstrip("]")
    matrix = np.zeros((n, n), dtype=np.int64)
    for token in body.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            pair, mult = token.split("^")
            i, j = pair.split("-")
        else:
            pair, mult = token.split("^")
            if n <= 10:
                i, j = pair[0], pair[1]
            else:
                raise ValueError(f"ambiguous vertex parse in {token!r}")
        i, j, m = int(i), int(j), int(mult)
        matrix[i, j] += m
        matrix[j, i] += m
    return matrix


def index_specification(matrix):
    """Canonical dummy-index assignment.

    Returns (slots, pairs) where slots[v] is the ordered list of dummy names
    on vertex v, and pairs is the list of contracted (name, name) couples --
    one inverse metric each.
    """
    n = matrix.shape[0]
    slots = [[] for _ in range(n)]
    pairs = []
    counter = 0
    for i in range(n):
        for j in range(i + 1, n):
            for _ in range(int(matrix[i, j])):
                a = f"i{counter}"
                b = f"j{counter}"
                counter += 1
                slots[i].append(a)
                slots[j].append(b)
                pairs.append((a, b))
    for v, s in enumerate(slots):
        if len(s) != FORM_RANK:
            raise ValueError(
                f"vertex {v} has valence {len(s)}, expected {FORM_RANK}")
    return slots, pairs


def to_latex(matrix, form_symbol="F"):
    """Einstein-notation LaTeX for the scalar."""
    slots, pairs = index_specification(matrix)
    factors = []
    for v, s in enumerate(slots):
        idx = "".join("\\" + name if False else f"{{\\mu_{{{name[1:]}}}}}"
                      for name in s)
        factors.append(f"{form_symbol}_{{{' '.join(s)}}}")
    metrics = " ".join(f"\\eta^{{{a}{b}}}" for a, b in pairs)
    return metrics + " \\, " + " ".join(factors)


def contraction_specification(label):
    """Full machine-readable spec for one graph."""
    matrix = parse_graph_label(label)
    slots, pairs = index_specification(matrix)
    n = matrix.shape[0]
    return {
        "label": label,
        "vertices": n,
        "form_rank": FORM_RANK,
        "field_degree": n,
        "multiplicity_matrix": matrix.tolist(),
        "slots": {str(v): slots[v] for v in range(n)},
        "metric_pairs": [list(p) for p in pairs],
        "metric_factor_count": len(pairs),
        "latex": to_latex(matrix),
        "rule": (
            "I = prod_v F_{slots[v]} x prod_pairs eta^{pair}; one inverse "
            "metric per contracted index pair; vertex valence 5"),
    }


def dense_evaluate(label, five_form_dense, lorentzian=True, mod=P):
    """Evaluate the scalar directly from the explicit index specification.

    Deliberately independent of the repository's graph slot-planner: the
    einsum subscripts are built from the dummy-index assignment above, and the
    metric appears as explicit eta factors rather than being folded into the
    tensors. This is the implementation-independent check on the translation.
    """
    matrix = parse_graph_label(label)
    slots, pairs = index_specification(matrix)
    n = matrix.shape[0]

    # One letter per contracted PAIR, not per slot. Writing eta as an explicit
    # operand would need two letters per pair (60 at degree 12) and blow the
    # 52-letter budget, besides adding 30 operands. Instead the metric is
    # applied by raising the index on exactly one side of each pair, which is
    # the same contraction.
    alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if len(pairs) > len(alphabet):
        raise ValueError(
            f"{len(pairs)} contracted pairs exceeds the einsum letter budget")
    letter = {}
    raise_on = {}
    for k, (a, b) in enumerate(pairs):
        letter[a] = letter[b] = alphabet[k]
        raise_on[b] = True          # raise the second slot of every pair

    signs = metric_signs(DIM, lorentzian).astype(np.int64) % mod

    subs, ops = [], []
    for v in range(n):
        tensor = five_form_dense
        for axis, name in enumerate(slots[v]):
            if raise_on.get(name):
                shape = [1] * tensor.ndim
                shape[axis] = DIM
                tensor = (tensor * signs.reshape(shape)) % mod
        subs.append("".join(letter[s] for s in slots[v]))
        ops.append(tensor)

    return int(mod_einsum(",".join(subs) + "->", ops, mod) % mod)
