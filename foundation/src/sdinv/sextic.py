"""Intrinsic sixth-order invariants of a self-dual D=10 five-form."""

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import itertools

import numpy as np

from .contract import value
from .forms import perm_sign
from .graphs import canonical, graph_from_label, graph_label
from .modp import mod_einsum
from .stress import (
    _raise_axes,
    composite_n1050,
    five_form_moment,
    matrix_trace_power,
)


INTRINSIC_TO_REGISTRY = (
    (Fraction(32, 3), Fraction(-1, 1125)),
    (Fraction(0), Fraction(3, 125)),
)
REGISTRY_TO_INTRINSIC = (
    (Fraction(3, 32), Fraction(1, 288)),
    (Fraction(0), Fraction(125, 3)),
)


def _n1050_shuffle_terms(first_five, last, internal_prefix):
    """Expand ``Lambda^mn_[abc Lambda_de]fmn`` into ten shuffle terms."""
    first_five = tuple(first_five)
    internal = (f"{internal_prefix}0", f"{internal_prefix}1")
    terms = []
    for selected_positions in itertools.combinations(range(5), 3):
        selected_positions = tuple(selected_positions)
        remaining_positions = tuple(
            position for position in range(5)
            if position not in selected_positions
        )
        shuffle = selected_positions + remaining_positions
        coefficient = Fraction(perm_sign(shuffle), 10)
        selected = tuple(first_five[position]
                         for position in selected_positions)
        remaining = tuple(first_five[position]
                          for position in remaining_positions)
        terms.append((
            coefficient,
            (
                internal + selected,
                remaining + (last,) + internal,
            ),
        ))
    return tuple(terms)


def _term_graph(ordered_vertices):
    """Return a graph and the orientation sign of ordered form slots."""
    occurrences = defaultdict(list)
    for vertex, slots in enumerate(ordered_vertices):
        if len(slots) != 5 or len(set(slots)) != 5:
            raise ValueError("each five-form vertex needs five distinct slots")
        for label in slots:
            occurrences[label].append(vertex)
    if any(len(vertices) != 2 for vertices in occurrences.values()):
        raise ValueError("every contraction label must occur exactly twice")

    matrix = np.zeros((6, 6), dtype=np.int64)
    pair_labels = defaultdict(list)
    for label, vertices in occurrences.items():
        left, right = sorted(vertices)
        if left == right:
            raise ValueError("antisymmetric self-contractions vanish")
        matrix[left, right] += 1
        matrix[right, left] += 1
        pair_labels[left, right].append(label)

    evaluator_slots = [[] for _ in range(6)]
    for left in range(6):
        for right in range(left + 1, 6):
            for label in sorted(pair_labels[left, right]):
                evaluator_slots[left].append(label)
                evaluator_slots[right].append(label)

    orientation = 1
    for evaluator, requested in zip(evaluator_slots, ordered_vertices):
        positions = {label: index for index, label in enumerate(evaluator)}
        orientation *= perm_sign([positions[label] for label in requested])
    return matrix, orientation


def _evaluator_slots(matrix):
    """The edge order used by the graph contraction evaluator."""
    slots = [[] for _ in range(matrix.shape[0])]
    for left in range(matrix.shape[0]):
        for right in range(left + 1, matrix.shape[0]):
            for occurrence in range(int(matrix[left, right])):
                edge = (left, right, occurrence)
                slots[left].append(edge)
                slots[right].append(edge)
    return slots


def _isomorphism_sign(source, target):
    """Return the orientation sign relating isomorphic labeled graphs."""
    if source.shape != target.shape:
        raise ValueError("graph orders differ")
    order = source.shape[0]
    target_slots = _evaluator_slots(target)
    for permutation in itertools.permutations(range(order)):
        if not np.array_equal(
            source[np.ix_(permutation, permutation)], target
        ):
            continue
        inverse = {
            source_vertex: target_vertex
            for target_vertex, source_vertex in enumerate(permutation)
        }
        orientation = 1
        source_slots = _evaluator_slots(source)
        for target_vertex, source_vertex in enumerate(permutation):
            relabeled = []
            for left, right, occurrence in source_slots[source_vertex]:
                new_left, new_right = sorted((inverse[left], inverse[right]))
                relabeled.append((new_left, new_right, occurrence))
            positions = {
                edge: index
                for index, edge in enumerate(target_slots[target_vertex])
            }
            orientation *= perm_sign([
                positions[edge] for edge in relabeled
            ])
        return orientation
    raise ValueError("graphs have the same certificate but no isomorphism")


@lru_cache(maxsize=1)
def paper_i6_2_graph_expansion():
    """Expand paper equation (2.14) into exact contraction-graph scalars.

    The convention is

    ``N1050_[abc,de]f = Lambda^mn_[abc Lambda_de]fmn``

    with normalized five-index antisymmetrization.  The three N tensors in
    equation (2.14) are contracted as

    ``N_abcdef N^abc_ghi N^defgih``.
    """
    outer = _n1050_shuffle_terms(
        ("r1", "r2", "r3", "r4", "r5"), "k", "u")
    second = _n1050_shuffle_terms(
        ("r1", "r2", "r3", "a1", "a2"), "a3", "v")
    third = _n1050_shuffle_terms(
        ("r4", "r5", "k", "a1", "a3"), "a2", "w")

    coefficients = defaultdict(Fraction)
    representatives = {}
    for (left_coefficient, left_vertices), (
        middle_coefficient, middle_vertices
    ), (right_coefficient, right_vertices) in itertools.product(
        outer, second, third
    ):
        matrix, orientation = _term_graph(
            left_vertices + middle_vertices + right_vertices)
        certificate = canonical(matrix)
        representative = representatives.setdefault(certificate, matrix)
        relabel_sign = _isomorphism_sign(matrix, representative)
        coefficients[certificate] += (
            left_coefficient
            * middle_coefficient
            * right_coefficient
            * orientation
            * relabel_sign
        )

    result = []
    for certificate, coefficient in coefficients.items():
        if not coefficient:
            continue
        result.append({
            "graph": graph_label(representatives[certificate]),
            "coefficient": coefficient,
        })
    return tuple(sorted(result, key=lambda item: item["graph"]))


def paper_i6_1(five_form, mod, backend="optimized"):
    """Return the paper's first sextic invariant ``Tr(M^3)``."""
    _, mixed = five_form_moment(five_form, mod, backend)
    return matrix_trace_power(mixed, 3, mod)


def paper_i6_2(five_form, mod):
    """Evaluate the intrinsic 1050-cubic invariant in equation (2.14)."""
    result = 0
    for item in paper_i6_2_graph_expansion():
        coefficient = item["coefficient"]
        denominator_inverse = pow(
            int(coefficient.denominator), -1, mod)
        residue = (
            int(coefficient.numerator) * denominator_inverse
        ) % mod
        result = (
            result
            + residue
            * value(
                graph_from_label(item["graph"]),
                five_form,
                10,
                5,
                True,
                mod,
            )
        ) % mod
    return int(result)


def paper_i6_2_direct(five_form, mod, backend="optimized"):
    """Independent direct six-index evaluation of paper equation (2.14)."""
    n1050 = composite_n1050(five_form, mod, backend)
    first_three_upper = _raise_axes(n1050, (0, 1, 2), mod)
    all_upper = _raise_axes(n1050, range(6), mod)
    return int(mod_einsum(
        "abcdef,abcghi,defgih->",
        [n1050, first_three_upper, all_upper],
        mod,
    ))


def registry_to_intrinsic_coefficients(coefficients):
    """Map ``(I6_1,I6_2)`` coefficients to ``(Tr(M^3),K_1050)``."""
    if len(coefficients) != 2:
        raise ValueError("the sextic coefficient vector must have length 2")
    return tuple(
        sum(
            REGISTRY_TO_INTRINSIC[row][column]
            * Fraction(coefficients[column])
            for column in range(2)
        )
        for row in range(2)
    )


def sextic_quotient_coordinate(coefficients):
    """Coefficient of the intrinsic 1050-cubic class modulo ``Tr(M^3)``."""
    return registry_to_intrinsic_coefficients(coefficients)[1]
