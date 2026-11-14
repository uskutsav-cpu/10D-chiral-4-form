"""Regression and adversarial checks on the classification deliverables."""

import copy
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from sdinv.certificate import (context, determinant_mod, digest, evaluate_item,
                               integral_directions, verify_witness)
from sdinv.forms import hodge_matrix
from sdinv.graphs import graph_from_label, graph_label, graph_to_record
from sdinv.latex import graph_to_latex
from sdinv.modp import RankSieve
from sdinv.primitive import PolynomialSieve, product_monomials
from sdinv.serialize import dumps, loads, to_record


def _regular_graph(n):
    matrix = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for d in (-2, -1, 1, 2, n // 2):
            matrix[i, (i + d) % n] = 1
    return matrix


@pytest.mark.parametrize("n", [10, 12, 14])
def test_graph_interchange_preserves_order_and_old_formats(n):
    graph = _regular_graph(n)
    for encoding in (dumps(graph), graph_label(graph), graph_to_record(graph)):
        assert np.array_equal(loads(encoding), graph)
    assert loads(dumps(graph))[0, n - 1] == 1


@pytest.mark.parametrize("record", [
    {"n": 12, "edges": [[0, 11, 1], [0, 11, 2]]},
    {"n": 12, "edges": [[0, 12, 1]]},
    {"n": 12, "edges": [[0, 0, 1]]},
    {"n": 12, "edges": [[0, 11, 1.5]]},
    {"n": 12, "edges": [[0, 11, True]]},
    {"n": -1, "edges": []},
    {"order": 4, "upper_triangle": [0, 1, 1, 1, 1, -1]},
    "n12[011^1]",
])
def test_graph_interchange_rejects_ambiguous_or_lossy_data(record):
    with pytest.raises(ValueError):
        loads(record)


def test_product_and_functional_ranks_are_different():
    # x^2 and y^2 are independent homogeneous polynomials on two variables.
    # Their gradients add no functional directions after x and y.
    polynomial = PolynomialSieve(2, 3, 32749)
    assert polynomial.add("x^2", [1, 4, 9])
    assert polynomial.add("y^2", [4, 9, 25])
    functional = RankSieve(2, 32749)
    assert functional.add([1, 0]) and functional.add([0, 1])
    assert not functional.add([2, 0]) and not functional.add([0, 4])
    assert polynomial.rank == functional.rank == 2


def test_degree12_product_inventory_is_ten_monomials():
    generators = {"I4_1": 4, "I6_1": 6, "I6_2": 6}
    generators.update({f"I8_{k}": 8 for k in range(1, 7)})
    generators.update({f"I10_{k}": 10 for k in range(1, 13)})
    products = product_monomials(generators, 12)
    assert len(products) == len(set(products)) == 10
    assert ("I4_1",) * 3 in products
    assert ("I6_1", "I6_2") in products


def test_integral_coordinates_are_selfdual_and_have_no_denominators():
    directions, indices = integral_directions()
    assert set(np.unique(directions)) == {-1, 0, 1}
    for prime in (32749, 32719):
        assert np.array_equal(hodge_matrix(10, 5, True, prime) @ directions % prime, directions % prime)
        _, basis, compact = context(prime, list(range(126)))
        assert compact[indices].tolist() == list(range(126))
        assert basis.ncols == 126


def test_determinant_handles_row_swaps_and_singular_matrices():
    assert determinant_mod([[0, 2], [3, 4]], 32749) == 32743
    assert determinant_mod([[2, 4], [3, 6]], 32749) == 0
    with pytest.raises(ValueError, match="prime"):
        determinant_mod([[1]], 32745)


def test_integral_jacobian_matches_finite_difference_directional_polynomial():
    # Quartic directional derivative recovered independently by exact
    # interpolation from five values (not a floating difference quotient).
    from sdinv.contract import planned_value
    from sdinv.exactmap import solve_full_column_rank
    item = json.loads((ROOT / "results/10d_order8.json").read_text())["generators"][0]
    graph = graph_from_label(item["graph"])
    record = {"id": item["id"], "degree": 4, "graph": to_record(graph)}
    prime = 32749
    a = np.random.default_rng(178).integers(0, prime, 126)
    direction = np.random.default_rng(179).integers(0, prime, 126)
    evaluated = evaluate_item(record, prime, a)
    values = []
    for t in range(5):
        form, _, _ = context(prime, (a + t * direction) % prime)
        values.append(planned_value(graph, form, 10, 5, True, prime))
    vandermonde = [[pow(t, k, prime) for k in range(5)] for t in range(5)]
    coefficients = solve_full_column_rank(vandermonde, values, prime)
    assert int(np.array(evaluated["jacobian_row"]) @ direction % prime) == coefficients[1]


def test_latex_has_one_upper_and_one_lower_occurrence_per_edge():
    graph = _regular_graph(12)
    formula = graph_to_latex(graph)
    for i, j, mult in to_record(graph)["edges"]:
        for k in range(1, mult + 1):
            assert formula.count(rf"\mu_{{{i},{j},{k}}}") == 2
    assert formula.count("F") == 12


def test_saved_certificate_and_tamper_rejection():
    from verify_rank81 import verify
    basis_path, certificate_path = ROOT / "results/rank81_basis.json", ROOT / "results/rank81_certificate.json"
    if not certificate_path.exists():
        pytest.skip("run scripts/search_rank81.py to produce release certificate")
    manifest, certificate = json.loads(basis_path.read_text()), json.loads(certificate_path.read_text())
    if not certificate.get("complete"):
        pytest.skip("certificate generation is still in progress")
    assert len(verify(manifest, certificate)) >= 4
    witness = certificate["witnesses"][0]
    for key, replacement in [("determinant_mod_p", 0), ("pivot_columns", [0] * 81),
                              ("invariant_ids", list(reversed(witness["invariant_ids"]))),
                              ("rank", 80)]:
        corrupted = copy.deepcopy(witness)
        corrupted[key] = replacement
        with pytest.raises(ValueError):
            verify_witness(manifest, corrupted)
    altered = copy.deepcopy(manifest)
    altered["invariants"][0]["graph"]["edges"][0][2] += 1
    with pytest.raises(ValueError):
        verify_witness(altered, witness)


def test_published_degree10_map_does_not_invent_an_invertible_matrix():
    path = ROOT / "results/order10_change_of_basis.json"
    if not path.exists():
        pytest.skip("run scripts/map_literature_basis.py --degree 10")
    result = json.loads(path.read_text())
    assert result["published_span_rank"] == 12
    assert result["union_rank"] == 13
    assert result["primitive_quotient_rank"] == 11
    assert len(result["matrix_12x14"]) == 12
    assert all(len(row) == 14 for row in result["matrix_12x14"])


def test_degree_runner_resume_and_bounded_workers_reproduce(tmp_path):
    from run_degree import run
    checkpoint = tmp_path / "degree4.json"
    first = run(4, 32749, [178, 179], checkpoint=checkpoint)
    resumed = run(4, 32749, [178, 179], checkpoint=checkpoint)
    parallel = run(4, 32749, [178, 179], workers=2)
    assert first == resumed == parallel
    assert first["polynomial_complete"] and first["functional_rank"] == 1
    with pytest.raises(ValueError, match="identity"):
        run(4, 32749, [180, 181], checkpoint=checkpoint)


def test_octic_map_checks_products_and_current_tensor_evaluator():
    from fractions import Fraction
    from map_literature_basis import reduce_matrix, inverse
    from sdinv.published_degree8_invariants import evaluate_selected_basis
    path = ROOT / "results/order8_change_of_basis.json"
    if not path.exists():
        pytest.skip("generate the octic literature map")
    result = json.loads(path.read_text())
    matrix = [[Fraction(x) for x in row] for row in result["graph_to_literature"]]
    assert inverse(matrix) == [[Fraction(x) for x in row] for row in result["literature_to_graph"]]
    assert len(result["primitive_quotient_matrix_6x6"]) == 6
    assert not result["exact_6x6_without_products"]
    prime = result["holdout_prime"]
    for sample in result["prime_witnesses"][str(prime)]["samples"]:
        assert np.array_equal(reduce_matrix(matrix, prime) @ np.array(sample["tensors"]) % prime, sample["graphs"])
    sample = result["prime_witnesses"][str(prime)]["samples"][0]
    form, _, _ = context(prime, sample["coordinates"])
    assert evaluate_selected_basis(form, prime) == sample["tensors"][:6]
