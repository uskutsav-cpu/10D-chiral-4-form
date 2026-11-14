"""Exact derivatives on the constrained self-dual five-form field space.

This module deliberately sits outside the immutable degree-12 validation
engine.  It reuses that engine's reverse contraction primitive without
changing any file covered by the committed atlas certificate.
"""

import itertools
from math import factorial

import numpy as np

from .contract import (
    _contraction_plan_profile,
    _network_value_gradient,
    _optimal_contraction_tree,
    _signed,
    _slot_plan,
)
from .forms import (
    basis_tuples,
    hodge_matrix,
    metric_signs,
    perm_sign,
    to_dense,
)
from .modp import P, inv, mod_einsum


def value_and_dense_gradient(matrix, five_form, d, valence,
                             lorentzian=True, mod=P,
                             max_memory_bytes=None):
    """Return a graph value and its ambient lower-input dense gradient."""
    slots, tails = _slot_plan(matrix, valence)
    signs = metric_signs(d, lorentzian) % mod
    terms = ["".join(labels) for labels in slots]
    profile = _contraction_plan_profile(
        terms, [(d,) * valence for _ in terms])
    if (
        max_memory_bytes is not None
        and profile["estimated_peak_bytes"] > int(max_memory_bytes)
    ):
        raise MemoryError(
            "exact contraction plan rejected before signed operands: "
            f"estimated {profile['estimated_peak_bytes']} bytes exceeds "
            f"limit {int(max_memory_bytes)} bytes")
    operands = [
        _signed(five_form, tails[vertex], signs, mod)
        for vertex in range(matrix.shape[0])
    ]
    scalar, operand_gradients = _network_value_gradient(
        terms,
        operands,
        mod,
        plan_profile=profile,
        max_memory_bytes=max_memory_bytes,
    )
    ambient = np.zeros_like(five_form, dtype=np.int64)
    for vertex, gradient in enumerate(operand_gradients):
        ambient = (
            ambient + _signed(gradient, tails[vertex], signs, mod)
        ) % mod
    return int(scalar) % mod, ambient


def covariant_antisymmetric_gradient(ambient, d, valence,
                                     lorentzian=True, mod=P):
    """Antisymmetrize an ambient gradient and lower all of its indices."""
    tensor = np.asarray(ambient, dtype=np.int64) % mod
    expected = (int(d),) * int(valence)
    if tensor.shape != expected:
        raise ValueError(
            f"gradient must have shape {expected}, got {tensor.shape}")
    result = np.zeros_like(tensor, dtype=np.int64)
    for permutation in itertools.permutations(range(int(valence))):
        result = (
            result
            + perm_sign(permutation) * tensor.transpose(permutation)
        ) % mod
    result = result * inv(factorial(int(valence)), mod) % mod
    signs = metric_signs(d, lorentzian) % mod
    for axis in range(int(valence)):
        shape = [1] * int(valence)
        shape[axis] = int(d)
        result = result * signs.reshape(shape) % mod
    return np.asarray(result, dtype=np.int64)


def anti_selfdual_projector(d, valence, lorentzian=True, mod=P):
    """Return the exact middle-form projector ``(1-*)/2``."""
    if d != 2 * valence:
        raise ValueError("anti-self-dual projection requires d=2*valence")
    hodge = hodge_matrix(d, valence, lorentzian, mod)
    identity = np.eye(hodge.shape[0], dtype=np.int64)
    return (identity - hodge) * inv(2, mod) % mod


def invariant_value_and_derivative(matrix, five_form, d=10, valence=5,
                                   lorentzian=True, mod=P,
                                   max_memory_bytes=None):
    """Return ``(I, dI/dLambda^upper)`` on the self-dual field space.

    The paper's constrained variation makes the derivative anti-self-dual.
    Its dense covariant components are normalized so that a homogeneous
    invariant of degree ``n`` obeys
    ``Lambda_upper . derivative_lower = n*I``.
    """
    scalar, ambient = value_and_dense_gradient(
        matrix,
        five_form,
        d,
        valence,
        lorentzian,
        mod,
        max_memory_bytes,
    )
    covariant = covariant_antisymmetric_gradient(
        ambient, d, valence, lorentzian, mod)
    tuples = basis_tuples(d, valence)
    compact = np.asarray(
        [covariant[index] for index in tuples],
        dtype=np.int64,
    )
    compact = (
        anti_selfdual_projector(d, valence, lorentzian, mod) @ compact
    ) % mod
    return scalar, to_dense(compact, d, valence, mod)


def _network_value_gradient_hvp(terms, operands, tangents, mod=P):
    """Forward-over-reverse contraction returning leaf Hessian products."""
    if len(terms) != len(operands) or len(terms) != len(tangents):
        raise ValueError("terms, operands, and tangents must align")
    nodes = [
        {
            "term": term,
            "value": np.asarray(operand, dtype=np.int64) % mod,
            "tangent": np.asarray(tangent, dtype=np.int64) % mod,
            "left": None,
            "right": None,
        }
        for term, operand, tangent in zip(terms, operands, tangents)
    ]
    split, boundary, _ = _optimal_contraction_tree(
        terms, [np.shape(operand) for operand in operands])
    label_bits = {
        label: 1 << index
        for index, label in enumerate(dict.fromkeys("".join(terms)))
    }

    def forward(mask):
        if mask.bit_count() == 1:
            return (mask & -mask).bit_length() - 1
        left_mask, right_mask = split[mask]
        left, right = forward(left_mask), forward(right_mask)
        left_term = nodes[left]["term"]
        right_term = nodes[right]["term"]
        keep_labels = boundary[mask]
        keep = "".join(
            label for label in dict.fromkeys(left_term + right_term)
            if keep_labels & label_bits[label]
        )
        expression = f"{left_term},{right_term}->{keep}"
        value = mod_einsum(
            expression,
            [nodes[left]["value"], nodes[right]["value"]],
            mod,
        )
        tangent = (
            mod_einsum(
                expression,
                [nodes[left]["tangent"], nodes[right]["value"]],
                mod,
            )
            + mod_einsum(
                expression,
                [nodes[left]["value"], nodes[right]["tangent"]],
                mod,
            )
        ) % mod
        parent = len(nodes)
        nodes.append({
            "term": keep,
            "value": value,
            "tangent": tangent,
            "left": left,
            "right": right,
        })
        return parent

    root = forward((1 << len(operands)) - 1)
    if nodes[root]["term"]:
        raise ValueError("network is not fully contracted")
    adjoints = {root: np.array(1, dtype=np.int64)}
    tangent_adjoints = {root: np.array(0, dtype=np.int64)}
    for parent in range(root, len(operands) - 1, -1):
        node = nodes[parent]
        left, right = node["left"], node["right"]
        parent_term = node["term"]
        left_term = nodes[left]["term"]
        right_term = nodes[right]["term"]
        adjoint = adjoints[parent]
        tangent_adjoint = tangent_adjoints[parent]

        left_expression = (
            f"{parent_term},{right_term}->{left_term}")
        adjoints[left] = mod_einsum(
            left_expression,
            [adjoint, nodes[right]["value"]],
            mod,
        )
        tangent_adjoints[left] = (
            mod_einsum(
                left_expression,
                [tangent_adjoint, nodes[right]["value"]],
                mod,
            )
            + mod_einsum(
                left_expression,
                [adjoint, nodes[right]["tangent"]],
                mod,
            )
        ) % mod

        right_expression = (
            f"{parent_term},{left_term}->{right_term}")
        adjoints[right] = mod_einsum(
            right_expression,
            [adjoint, nodes[left]["value"]],
            mod,
        )
        tangent_adjoints[right] = (
            mod_einsum(
                right_expression,
                [tangent_adjoint, nodes[left]["value"]],
                mod,
            )
            + mod_einsum(
                right_expression,
                [adjoint, nodes[left]["tangent"]],
                mod,
            )
        ) % mod
    return (
        int(np.asarray(nodes[root]["value"]).item()) % mod,
        int(np.asarray(nodes[root]["tangent"]).item()) % mod,
        [adjoints[index] for index in range(len(operands))],
        [tangent_adjoints[index] for index in range(len(operands))],
    )


def invariant_value_derivative_hvp(matrix, five_form, direction, d=10,
                                   valence=5, lorentzian=True, mod=P):
    """Return value, directional value, derivative, and derivative HVP.

    ``direction`` must be a dense self-dual tangent five-form.  The last
    output is the directional derivative of
    ``dI/dLambda^upper`` and therefore has the same covariant anti-self-dual
    normalization as :func:`invariant_value_and_derivative`.
    """
    five_form = np.asarray(five_form, dtype=np.int64) % mod
    direction = np.asarray(direction, dtype=np.int64) % mod
    expected = (int(d),) * int(valence)
    if five_form.shape != expected or direction.shape != expected:
        raise ValueError(f"five_form and direction must have shape {expected}")
    slots, tails = _slot_plan(matrix, valence)
    signs = metric_signs(d, lorentzian) % mod
    terms = ["".join(labels) for labels in slots]
    operands = [
        _signed(five_form, tails[vertex], signs, mod)
        for vertex in range(matrix.shape[0])
    ]
    tangents = [
        _signed(direction, tails[vertex], signs, mod)
        for vertex in range(matrix.shape[0])
    ]
    scalar, scalar_tangent, gradients, hvps = (
        _network_value_gradient_hvp(
            terms, operands, tangents, mod))
    ambient_gradient = np.zeros_like(five_form, dtype=np.int64)
    ambient_hvp = np.zeros_like(five_form, dtype=np.int64)
    for vertex in range(matrix.shape[0]):
        ambient_gradient = (
            ambient_gradient
            + _signed(gradients[vertex], tails[vertex], signs, mod)
        ) % mod
        ambient_hvp = (
            ambient_hvp
            + _signed(hvps[vertex], tails[vertex], signs, mod)
        ) % mod

    projector = anti_selfdual_projector(
        d, valence, lorentzian, mod)
    tuples = basis_tuples(d, valence)

    def project(ambient):
        covariant = covariant_antisymmetric_gradient(
            ambient, d, valence, lorentzian, mod)
        compact = np.asarray(
            [covariant[index] for index in tuples],
            dtype=np.int64,
        )
        compact = projector @ compact % mod
        return to_dense(compact, d, valence, mod)

    return scalar, scalar_tangent, project(ambient_gradient), project(
        ambient_hvp)
