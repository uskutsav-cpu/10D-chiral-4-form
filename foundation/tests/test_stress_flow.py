"""Regression gates for the exact stress-flow physics layer."""

from fractions import Fraction
import os
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.exactmap import (  # noqa: E402
    crt,
    rank_mod,
    rational_reconstruct,
    reconstruct_vector,
    solve_full_column_rank,
)
from sdinv.forms import (  # noqa: E402
    basis_tuples,
    hodge_matrix,
    metric_signs,
    random_form,
    selfdual_projector,
    to_dense,
)
from sdinv.invariant_registry import (  # noqa: E402
    InvariantItem,
    degree12_placeholder_items,
    degree12_product_items,
    load_verified_registry,
)
from sdinv.graphs import graph_to_record  # noqa: E402
from sdinv.modp import inv  # noqa: E402
from sdinv.stress import (  # noqa: E402
    _antisymmetrize_axes,
    _raise_axes,
    composite_n,
    composite_n4125,
    five_form_moment,
    interacting_stress,
    interaction_gradient_i4,
    interacting_trace_formula,
    matrix_trace_power,
    modmax_stress,
    modmax_stress_square_formula,
    n4125_mm,
    paper_i4_i8_i12,
    paper_i8_i12_expanded,
    raw_n_mm,
    stress_correction_r,
    stress_mixed,
    stress_traces,
    stress_v_i4,
    stress_v_i4_raw,
    symmetric_inner,
)


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESULT = os.path.join(ROOT, "results", "stress_flow_exact_low_degree.json")
THIRD_PRIME = 32693


def _sample(seed, prime):
    projector = selfdual_projector(10, 5, True, prime)
    compact = (
        projector
        @ random_form(10, 5, np.random.default_rng(seed), prime)
    ) % prime
    return to_dense(compact, 10, 5, prime)


def _transform_covariant_tensor(tensor, transformation, prime):
    result = np.asarray(tensor, dtype=np.int64) % prime
    for axis in range(result.ndim):
        result = np.tensordot(
            transformation, result, axes=(1, axis))
        result = np.moveaxis(result, 0, axis) % prime
    return result


@pytest.mark.parametrize("prime", [32749, 32719, THIRD_PRIME])
def test_paper_m_n_and_relevant_4125_contraction_conventions(prime):
    five_form = _sample(20260901, prime)
    m_lower, m_mixed = five_form_moment(five_form, prime)
    n_lower = composite_n(five_form, prime)
    r_lower = stress_correction_r(five_form, prime)
    projected = n4125_mm(five_form, prime)

    assert np.array_equal(m_lower, m_lower.T)
    assert matrix_trace_power(m_mixed, 1, prime) == 0
    assert np.array_equal(
        n_lower, n_lower.transpose(3, 4, 5, 0, 1, 2))
    for tensor in (r_lower, projected):
        assert np.array_equal(tensor, tensor.T)
        assert matrix_trace_power(stress_mixed(tensor, prime), 1, prime) == 0


def test_simple_reference_m_and_n_match_optimized_backend():
    five_form = _sample(20260902, 32749)
    optimized_m = five_form_moment(five_form, 32749, "optimized")
    reference_m = five_form_moment(five_form, 32749, "reference")
    assert all(np.array_equal(left, right)
               for left, right in zip(optimized_m, reference_m))
    assert np.array_equal(
        composite_n(five_form, 32749, "optimized"),
        composite_n(five_form, 32749, "reference"),
    )


def test_full_4125_projector_symmetries_traces_and_relevant_contraction():
    prime = 32749
    five_form = _sample(20260908, prime)
    projected = composite_n4125(five_form, prime)
    assert np.array_equal(
        projected, projected.transpose(3, 4, 5, 0, 1, 2))
    signs = metric_signs(10, True) % prime
    cross_trace = np.einsum(
        "a,abcade->bcde", signs, projected) % prime
    assert not np.any(cross_trace)
    assert not np.any(
        _antisymmetrize_axes(projected, (0, 1, 2, 3), prime))
    _, m_mixed = five_form_moment(five_form, prime)
    assert np.array_equal(
        raw_n_mm(projected, m_mixed, prime),
        n4125_mm(five_form, prime),
    )


@pytest.mark.parametrize("prime", [32749, 32719, THIRD_PRIME])
def test_equations_3_3_and_3_4_are_exactly_equivalent(prime):
    five_form = _sample(20260903, prime)
    assert np.array_equal(
        stress_v_i4(five_form, v=137, v_i=211, mod=prime),
        stress_v_i4_raw(five_form, v=137, v_i=211, mod=prime),
    )


@pytest.mark.parametrize("prime", [32749, 32719, THIRD_PRIME])
def test_general_interacting_stress_reduces_to_v_i4_formula(prime):
    """Equation (2.33) must reproduce the independent equation (3.3)."""
    five_form = _sample(20260903, prime)
    v, v_i = 137, 211
    i4 = matrix_trace_power(
        five_form_moment(five_form, prime)[1], 2, prime)
    derivative = interaction_gradient_i4(five_form, v_i, prime)
    general = interacting_stress(
        five_form,
        derivative,
        v=v,
        euler=4 * i4 * v_i,
        mod=prime,
    )
    assert np.array_equal(
        general,
        stress_v_i4_raw(five_form, v=v, v_i=v_i, mod=prime),
    )
    derivative_vector = np.asarray([
        derivative[index] for index in basis_tuples(10, 5)
    ])
    assert np.array_equal(
        hodge_matrix(10, 5, True, prime) @ derivative_vector % prime,
        -derivative_vector % prime,
    )
    assert matrix_trace_power(stress_mixed(general, prime), 1, prime) == (
        interacting_trace_formula(v, 4 * i4 * v_i, prime)
    )
    assert np.array_equal(
        interaction_gradient_i4(
            five_form, v_i, prime, backend="reference"),
        derivative,
    )
    assert np.array_equal(
        interacting_stress(
            five_form,
            derivative,
            v=v,
            euler=4 * i4 * v_i,
            mod=prime,
            backend="reference",
        ),
        general,
    )


@pytest.mark.parametrize("prime", [32749, 32719, THIRD_PRIME])
def test_registry_gradient_has_paper_normalization(prime):
    """The graph reverse derivative calibrates to the explicit I4 formula."""
    five_form = _sample(20260909, prime)
    registry = load_verified_registry(ROOT)
    graph_i4, graph_gradient = registry.evaluate_item_with_gradient(
        "I4_1", five_form, prime)
    paper_gradient = interaction_gradient_i4(five_form, 1, prime)
    assert np.array_equal(2 * graph_gradient % prime, paper_gradient)
    assert 2 * graph_i4 % prime == matrix_trace_power(
        five_form_moment(five_form, prime)[1], 2, prime)


@pytest.mark.parametrize("prime", [32749, 32719, THIRD_PRIME])
def test_registry_gradients_are_anti_selfdual_and_obey_euler(prime):
    five_form = _sample(20260910, prime)
    five_form_upper = _raise_axes(five_form, range(5), prime)
    registry = load_verified_registry(ROOT)
    hodge = hodge_matrix(10, 5, True, prime)
    for item_id in ("I4_1", "I6_1", "I6_2", "I4_1^2"):
        item = registry.item(item_id)
        scalar, derivative = registry.evaluate_item_with_gradient(
            item_id, five_form, prime)
        derivative_vector = np.asarray([
            derivative[index] for index in basis_tuples(10, 5)
        ])
        assert np.array_equal(
            hodge @ derivative_vector % prime,
            -derivative_vector % prime,
        )
        assert int(np.sum(
            five_form_upper * derivative, dtype=np.int64
        ) % prime) == item.degree * scalar % prime


@pytest.mark.parametrize("prime", [32749, 32719, THIRD_PRIME])
def test_published_i8_i12_and_modmax_stress_square(prime):
    five_form = _sample(20260904, prime)
    compact = paper_i4_i8_i12(five_form, prime)
    expanded = paper_i8_i12_expanded(five_form, prime)
    assert compact["I8"] == expanded["I8"]
    assert compact["I12"] == expanded["I12"]

    stress = modmax_stress(five_form, b=173, mod=prime)
    formula = modmax_stress_square_formula(
        five_form, b=173, mod=prime)
    assert np.array_equal(stress, stress.T)
    assert matrix_trace_power(stress_mixed(stress, prime), 1, prime) == 0
    assert symmetric_inner(stress, stress, prime) == (
        formula["stress_square"])
    assert formula["stress_square"] == (
        formula["root_term"]
        + formula["delta_I8"]
        + formula["delta_I12"]
    ) % prime


def test_equation_3_16_normalization_identity():
    """The paper uses the 4!-rescaled INZ density in its flow equation."""
    # b=-1/2*sqrt(7/6)*tanh(gamma/2).  After subtracting Delta,
    # T^2=I4*sech(gamma/2)^4/(4*(4!)^2).  Both sides of (3.16)
    # consequently have coefficient sqrt(7/6)*sech^2/4 times sqrt(I4).
    rhs_without_sqrt_7_over_6 = (
        Fraction(24, 2) * Fraction(1, 2 * 24))
    minus_db_dgamma_without_sqrt_7_over_6 = Fraction(1, 4)
    assert rhs_without_sqrt_7_over_6 == (
        minus_db_dgamma_without_sqrt_7_over_6)


def test_lorentz_boost_covariance_of_m_r_and_modmax_t():
    prime = 32749
    parameter = 11
    parameter_inverse = inv(parameter, prime)
    half = inv(2, prime)
    cosine = (parameter + parameter_inverse) * half % prime
    sine = (parameter_inverse - parameter) * half % prime
    boost = np.eye(10, dtype=np.int64)
    boost[0, 0] = cosine
    boost[0, 1] = sine
    boost[1, 0] = sine
    boost[1, 1] = cosine
    eta = np.diag(metric_signs(10, True)).astype(np.int64) % prime
    assert np.array_equal(boost.T @ eta @ boost % prime, eta)

    five_form = _sample(20260905, prime)
    transformed = _transform_covariant_tensor(
        five_form, boost, prime)
    for builder in (
        lambda tensor: five_form_moment(tensor, prime)[0],
        lambda tensor: stress_correction_r(tensor, prime),
        lambda tensor: modmax_stress(tensor, 17, prime),
        lambda tensor: interacting_stress(
            tensor,
            interaction_gradient_i4(tensor, 23, prime),
            v=23 * matrix_trace_power(
                five_form_moment(tensor, prime)[1], 2, prime),
            euler=4 * 23 * matrix_trace_power(
                five_form_moment(tensor, prime)[1], 2, prime),
            mod=prime,
        ),
    ):
        original = builder(five_form)
        expected = boost @ original @ boost.T % prime
        assert np.array_equal(builder(transformed), expected)


def test_stress_trace_cayley_hamilton_interface():
    five_form = _sample(20260906, 32749)
    stress = modmax_stress(five_form, 19, 32749)
    traces = stress_traces(stress, 10, 32749)
    assert tuple(traces) == tuple(range(1, 11))
    with pytest.raises(ValueError, match="CH limit"):
        stress_traces(stress, 11, 32749)


def test_v_equals_c_i4_trace_expansion_through_degree10():
    """Independently multiply the homogeneous T2+T4+T6 components."""
    prime, coupling = 32749, 29
    five_form = _sample(20260907, prime)
    m_lower, m_mixed = five_form_moment(five_form, prime)
    r_lower = stress_correction_r(five_form, prime)
    r_mixed = stress_mixed(r_lower, prime)
    i4 = matrix_trace_power(m_mixed, 2, prime)
    tr_m3 = matrix_trace_power(m_mixed, 3, prime)
    i8 = paper_i4_i8_i12(five_form, prime)["I8"]
    m2r = int(np.trace(
        (m_mixed @ m_mixed % prime) @ r_mixed % prime) % prime)

    def scaled(matrix, *factors):
        result = matrix.copy() % prime
        for factor in factors:
            result = result * (int(factor) % prime) % prime
        return result

    homogeneous = {
        2: scaled(m_mixed, inv(48, prime)),
        4: scaled(
            np.eye(10, dtype=np.int64),
            coupling,
            i4,
            inv(24, prime),
        ),
        6: (
            scaled(
                m_mixed,
                coupling * coupling,
                -2 * inv(7, prime),
                i4,
            )
            + scaled(
                r_mixed,
                coupling * coupling,
                inv(3, prime),
            )
        ) % prime,
    }

    def trace_polynomial(power):
        polynomial = {0: np.eye(10, dtype=np.int64)}
        for _ in range(power):
            product = {}
            for left_degree, left in polynomial.items():
                for right_degree, right in homogeneous.items():
                    degree = left_degree + right_degree
                    product[degree] = (
                        product.get(
                            degree, np.zeros((10, 10), dtype=np.int64))
                        + left @ right
                    ) % prime
            polynomial = product
        return {
            degree: int(np.trace(matrix) % prime)
            for degree, matrix in polynomial.items()
        }

    tr2 = trace_polynomial(2)
    assert tr2[4] == i4 * inv(48 ** 2, prime) % prime
    assert tr2[8] == (
        coupling * coupling
        * (i8 * inv(72, prime)
           + 11 * i4 * i4 * inv(2016, prime))
    ) % prime

    tr3 = trace_polynomial(3)
    assert tr3[8] == (
        coupling * i4 * i4 * inv(8 * 48 ** 2, prime)
    ) % prime
    assert tr3[10] == (
        coupling * coupling * inv(48 ** 2, prime)
        * (m2r - 6 * inv(7, prime) * i4 * tr_m3)
    ) % prime

    tr4 = trace_polynomial(4)
    assert tr4[10] == (
        coupling * i4 * tr_m3 * inv(6 * 48 ** 3, prime)
    ) % prime
    tr5 = trace_polynomial(5)
    assert tr5[10] == (
        matrix_trace_power(m_mixed, 5, prime)
        * inv(48 ** 5, prime)
    ) % prime


def test_verified_registry_and_degree12_extension_contract():
    registry = load_verified_registry(ROOT)
    assert {degree: len(registry.basis(degree))
            for degree in registry.degrees} == {
        4: 1,
        6: 2,
        8: 7,
        10: 14,
    }
    products = degree12_product_items()
    placeholders = degree12_placeholder_items()
    assert len(products) == 10
    assert len(placeholders) == 62
    assert sum(item.kind == "product" for item in products) == 10
    with pytest.raises(ValueError, match="62 primitives"):
        registry.with_degree12_primitives(())
    unsupported = tuple(
        InvariantItem(
            id=f"unsupported_{index}",
            degree=12,
            kind="graph",
            graph=f"n12[unsupported-{index}]",
            source="fixture",
        )
        for index in range(62)
    )
    with pytest.raises(ValueError, match="graph_record"):
        registry.with_degree12_primitives(unsupported)
    degree12_graph = np.zeros((12, 12), dtype=np.int64)
    for vertex in range(12):
        for offset in (1, 2, 6):
            neighbor = (vertex + offset) % 12
            degree12_graph[vertex, neighbor] = 1
            degree12_graph[neighbor, vertex] = 1
    assert np.all(degree12_graph.sum(axis=1) == 5)
    record = graph_to_record(degree12_graph)
    concrete = tuple(
        InvariantItem(
            id=item.id,
            degree=12,
            kind="graph_record",
            graph_record=record,
            source="fixture",
        )
        for item in placeholders
    )
    extended = registry.with_degree12_primitives(concrete)
    assert len(extended.basis(12)) == 72


def test_small_exact_map_linear_algebra_and_rational_reconstruction():
    primes = (32749, 32719, THIRD_PRIME)
    matrix = np.array([[1, 2], [3, 5], [8, 13]], dtype=np.int64)
    rational = (Fraction(-17, 9), Fraction(23, 7))
    vectors = []
    for prime in primes:
        expected = np.array([
            value.numerator * inv(value.denominator, prime) % prime
            for value in rational
        ], dtype=np.int64)
        target = matrix @ expected % prime
        assert rank_mod(matrix, prime) == 2
        vectors.append(solve_full_column_rank(matrix, target, prime))
    assert tuple(reconstruct_vector(vectors, primes)) == rational
    residue, modulus = crt(
        [(-17 * inv(9, prime)) % prime for prime in primes],
        primes,
    )
    assert rational_reconstruct(residue, modulus) == Fraction(-17, 9)


def test_pipeline_rejects_insufficient_crt_and_held_out_configurations():
    script = os.path.join(ROOT, "scripts", "stress_flow_pipeline.py")
    seeds_15 = [str(seed) for seed in range(20260901, 20260916)]
    for unsupported_primes in (
        ("2", "32749", "32719", "32693", "32717"),
        ("1000000007", "1000000009"),
    ):
        unsupported = subprocess.run(
            [
                sys.executable,
                script,
                "--primes", *unsupported_primes,
                "--sample-seeds", *seeds_15,
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert unsupported.returncode == 2
        assert "supported exact-arithmetic range" in unsupported.stderr

    too_small_crt = subprocess.run(
        [
            sys.executable,
            script,
            "--primes", "32749", "32719", "32693",
            "--sample-seeds", *seeds_15,
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert too_small_crt.returncode == 2
    assert "combined CRT modulus is too small" in too_small_crt.stderr

    seeds_14 = [str(seed) for seed in range(20260901, 20260915)]
    no_held_out = subprocess.run(
        [
            sys.executable,
            script,
            "--primes", "32749", "32719", "32693", "32717",
            "--sample-seeds", *seeds_14,
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert no_held_out.returncode == 2
    assert "15 distinct samples" in no_held_out.stderr


def test_saved_exact_maps_and_nonmembership_certificates():
    import json

    with open(RESULT) as stream:
        result = json.load(stream)
    assert result["schema"] == 1
    assert result["paper"]["source_sha256"] == (
        "eb889f3dfd60280a95a907343e39ade60d62df6f56b4cd53c33cd445283ffb4e")
    assert result["paper"]["pdf_sha256"] == (
        "c4f812035fb8d8d07aa8b0dba1a2e55de3e909a77f75397ee17f3fe21f4f5e90")
    assert "Pinned arXiv v2" in result["paper"]["hash_policy"]
    assert len(result["exact_arithmetic"]["primes"]) == 5
    assert not result["exact_arithmetic"]["floating_point_used_for_claims"]
    assert result["stress_trace_registry"]["cayley_hamilton_limit"] == 10
    assert [
        row["leading_five_form_degree"]
        for row in result["stress_trace_registry"]["traces"]
    ] == list(range(4, 21, 2))
    assert result["stress_subalgebra"]["dimension_table"] == [
        {
            "degree": 4,
            "five_form_value_dimension": 1,
            "stress_generated_dimension": 1,
            "complement_dimension": 0,
        },
        {
            "degree": 6,
            "five_form_value_dimension": 2,
            "stress_generated_dimension": 1,
            "complement_dimension": 1,
        },
        {
            "degree": 8,
            "five_form_value_dimension": 7,
            "stress_generated_dimension": 2,
            "complement_dimension": 5,
        },
        {
            "degree": 10,
            "five_form_value_dimension": 14,
            "stress_generated_dimension": 2,
            "complement_dimension": 12,
        },
    ]
    maps = result["degree_maps"]
    assert [entry["text"] for entry in maps["4"]["targets"]["tr_M2"]] == [
        "2"]
    assert [entry["text"] for entry in maps["6"]["targets"]["tr_M3"]] == [
        "32/3", "0"]
    assert [
        entry["text"] for entry in maps["8"]["targets"]["paper_I8"]
    ] == ["2", "-12", "0", "0", "0", "0", "-18/7"]

    def fractions(entries):
        return [
            Fraction(entry["numerator"], entry["denominator"])
            for entry in entries
        ]

    def add(left, right):
        return [a + b for a, b in zip(left, right)]

    def scale(row, factor):
        return [factor * entry for entry in row]

    flow = result["flow_closure_pilot"]
    trace_rows = flow["V_equals_c4_I4_trace_expansion"]
    i8 = fractions(maps["8"]["targets"]["paper_I8"])
    i4_squared = fractions(maps["8"]["targets"]["tr_M2^2"])
    assert fractions(
        trace_rows["TrT2"]["degree8_coefficient_of_c4^2"]["coordinates"]
    ) == add(
        scale(i8, Fraction(1, 72)),
        scale(i4_squared, Fraction(11, 2016)),
    )
    assert fractions(
        trace_rows["TrT3"]["degree8_coefficient_of_c4"]["coordinates"]
    ) == scale(i4_squared, Fraction(1, 8 * 48 ** 2))

    m2r = fractions(maps["10"]["targets"]["tr_M2R"])
    m2m3 = fractions(maps["10"]["targets"]["tr_M2*tr_M3"])
    assert fractions(
        trace_rows["TrT3"][
            "degree10_coefficient_of_c4^2"]["coordinates"]
    ) == scale(
        add(m2r, scale(m2m3, Fraction(-6, 7))),
        Fraction(1, 48 ** 2),
    )
    assert fractions(
        trace_rows["TrT4"][
            "degree10_coefficient_of_c4"]["coordinates"]
    ) == scale(m2m3, Fraction(1, 6 * 48 ** 3))

    general_v6 = flow["general_V6_linear_terms"]
    for suffix in ("1", "2"):
        row = fractions(
            general_v6[
                f"TrT3_degree10_coefficient_of_c6_{suffix}"][
                    "coordinates"])
        nonzero = [entry for entry in row if entry]
        assert nonzero == [Fraction(1, 2 * 48 ** 2)]

    for degree_map in maps.values():
        for prime_text, run in degree_map["per_prime"].items():
            prime = int(prime_text)
            assert run["held_out_sample_seeds"]
            assert run["held_out_validation_passed"]
            basis_values = np.asarray(
                run["sample_values"]["basis"], dtype=np.int64)
            for target, rational_row in degree_map["targets"].items():
                reconstructed = np.asarray([
                    entry["numerator"]
                    * inv(entry["denominator"], prime) % prime
                    for entry in rational_row
                ], dtype=np.int64)
                modular_solution = np.asarray(
                    run["solutions"][target], dtype=np.int64)
                assert np.array_equal(
                    reconstructed % prime, modular_solution % prime)
                target_values = np.asarray(
                    run["sample_values"]["targets"][target],
                    dtype=np.int64,
                )
                assert np.array_equal(
                    basis_values @ modular_solution % prime,
                    target_values % prime,
                )

    for degree in ("4", "6", "8", "10"):
        report = result["stress_subalgebra"]["degrees"][degree]
        expected_stress_rank = report["stress_generated_dimension"]
        for evidence in report[
                "modular_nonmembership_certificates"].values():
            assert evidence["stress_value_rank"] == expected_stress_rank
            assert all(
                rank == expected_stress_rank + 1
                for rank in evidence[
                    "stress_plus_obstruction_ranks"].values()
            )
    for runs in result["stage1_modmax_reproduction"].values():
        for sample in runs:
            assert sample["expanded_I8_matches"]
            assert sample["expanded_I12_matches"]
            assert sample["equation_3_3_equals_3_4"]
            assert sample["stress_symmetric"]
            assert sample["stress_trace"] == 0
            assert sample["stress_square_direct"] == (
                sample["stress_square_formula"])
