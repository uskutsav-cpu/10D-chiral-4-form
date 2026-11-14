"""Independent exact checks of the interacting formal trace expansion."""

import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.exactmap import sample_selfdual_five_form
from sdinv.contract import build_compact_derivative_basis
from sdinv.formal_flow import (
    TRACE_GENERATORS,
    evaluate_matrix_polynomial,
    flow_generator_polynomials,
    interaction_data,
    normalized_stress_polynomial,
    trace_polynomials,
)
from sdinv.forms import selfdual_projector, to_dense
from sdinv.graphs import graph_from_label
from sdinv.interaction import (
    invariant_value_and_derivative,
    invariant_value_derivative_hvp,
)
from sdinv.invariant_registry import load_verified_registry
from sdinv.modp import inv
from sdinv.stress import (
    five_form_moment,
    interacting_stress,
    matrix_trace_power,
    paper_i4_i8_i12,
    stress_mixed,
)


ROOT = Path(__file__).resolve().parents[1]


def test_normalized_sparse_tensor_equals_direct_interacting_stress():
    prime = 32749
    five_form = sample_selfdual_five_form(20261001, prime)
    registry = load_verified_registry(ROOT)
    values, derivatives = interaction_data(five_form, registry, prime)
    polynomial, _ = normalized_stress_polynomial(
        five_form, registry, prime, values, derivatives)
    coefficients = {"I4_1": 7, "I6_1": 11}
    expanded = evaluate_matrix_polynomial(
        polynomial, coefficients, prime)

    derivative = (
        7 * derivatives["I4_1"] + 11 * derivatives["I6_1"]
    ) % prime
    potential = (
        7 * values["I4_1"] + 11 * values["I6_1"]
    ) % prime
    euler = (
        4 * 7 * values["I4_1"] + 6 * 11 * values["I6_1"]
    ) % prime
    direct = 48 * stress_mixed(interacting_stress(
        five_form, derivative, potential, euler, prime), prime) % prime
    assert np.array_equal(expanded, direct)


def test_trace_polynomial_free_terms_and_published_i4_specialization():
    prime = 32749
    five_form = sample_selfdual_five_form(20261002, prime)
    registry = load_verified_registry(ROOT)
    values, derivatives = interaction_data(five_form, registry, prime)
    traces, _ = trace_polynomials(
        five_form, registry, prime, values, derivatives)
    _, m_mixed = five_form_moment(five_form, prime)
    for power in range(2, 7):
        assert traces[power][(2 * power, ())] == matrix_trace_power(
            m_mixed, power, prime)

    i4 = values["I4_1"]
    paper_i8 = paper_i4_i8_i12(five_form, prime)["I8"]
    # The verified graph basis has I4_1=paper_I4/2.
    expected_tau2_c4_squared = (
        8 * paper_i8 + 88 * inv(7, prime) * i4 * i4
    ) % prime
    assert traces[2][(8, ("I4_1", "I4_1"))] == (
        expected_tau2_c4_squared)
    assert traces[3][(8, ("I4_1",))] == 12 * i4 * i4 % prime

    for sextic in ("I6_1", "I6_2"):
        assert traces[3][(10, (sextic,))] == (
            24 * i4 * values[sextic] % prime)


def test_complete_trace_generator_registry_through_degree12():
    assert len(TRACE_GENERATORS) == 18
    assert len({generator.id for generator in TRACE_GENERATORS}) == 18
    assert {generator.trace_powers for generator in TRACE_GENERATORS} == {
        (1,), (2,), (3,),
        (4,), (1, 1), (1, 2), (2, 2),
        (5,), (1, 3), (2, 3),
        (6,), (1, 4), (2, 4), (3, 3),
        (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2),
    }

    prime = 32749
    five_form = sample_selfdual_five_form(20261003, prime)
    registry = load_verified_registry(ROOT)
    generators = flow_generator_polynomials(
        five_form, registry, prime)
    assert set(generators) == {
        generator.id for generator in TRACE_GENERATORS}
    assert all(
        degree <= 12
        for polynomial in generators.values()
        for degree, _ in polynomial
    )


def test_graph_hessian_vector_product_matches_exact_interpolation():
    prime = 32749
    five_form = sample_selfdual_five_form(20261007, prime)
    registry = load_verified_registry(ROOT)
    item = registry.item("I6_2")
    matrix = graph_from_label(item.graph)
    compact_basis = build_compact_derivative_basis(
        10,
        5,
        selfdual_projector(10, 5, True, prime),
        prime,
        independent=True,
    )
    direction = to_dense(
        compact_basis.directions[:, 17], 10, 5, prime)
    value, tangent, derivative, hvp = invariant_value_derivative_hvp(
        matrix, five_form, direction, 10, 5, True, prime)

    # Exact derivative weights for polynomials of degree at most six.
    points = tuple(range(7))
    vandermonde = np.asarray([
        [pow(point, power, prime) for point in points]
        for power in range(7)
    ], dtype=np.int64)
    target = np.zeros(7, dtype=np.int64)
    target[1] = 1
    from sdinv.exactmap import solve_full_column_rank
    weights = solve_full_column_rank(vandermonde, target, prime)

    sampled_values = []
    sampled_derivatives = []
    for point in points:
        shifted = (five_form + point * direction) % prime
        sample_value, sample_derivative = invariant_value_and_derivative(
            matrix, shifted, 10, 5, True, prime)
        sampled_values.append(sample_value)
        sampled_derivatives.append(sample_derivative)
    interpolated_tangent = sum(
        int(weight) * sample
        for weight, sample in zip(weights, sampled_values)
    ) % prime
    interpolated_hvp = np.zeros_like(hvp)
    for weight, sample in zip(weights, sampled_derivatives):
        interpolated_hvp = (
            interpolated_hvp + int(weight) * sample
        ) % prime

    assert value == sampled_values[0]
    assert np.array_equal(derivative, sampled_derivatives[0])
    assert tangent == interpolated_tangent
    assert np.array_equal(hvp, interpolated_hvp)
