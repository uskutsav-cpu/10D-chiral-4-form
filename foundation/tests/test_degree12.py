"""Fast structural gates for the committed degree-12 certificate."""

import importlib.util
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.catalog import canonical_graph_id
from sdinv.graphs import graph_from_label, validate_graph
from sdinv.invariant_registry import (
    load_verified_registry_through_degree12,
)


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESULT = os.path.join(ROOT, "results", "10d_order12.json")
BENCHMARKS = os.path.join(ROOT, "results", "degree12_benchmarks.json")
PIPELINE_SPEC = importlib.util.spec_from_file_location(
    "degree12_pipeline",
    os.path.join(ROOT, "scripts", "degree12_pipeline.py"),
)
degree12_pipeline = importlib.util.module_from_spec(PIPELINE_SPEC)
PIPELINE_SPEC.loader.exec_module(degree12_pipeline)


def test_degree12_product_inventory_and_leibniz_rows():
    mod = 101
    values = {
        "I4_1": 3,
        "I6_1": 5,
        "I6_2": 7,
        **{f"I8_{k}": 10 + k for k in range(1, 7)},
    }
    rows = {
        name: np.asarray([index + 1, 2 * index + 3], dtype=np.int64)
        for index, name in enumerate(values)
    }
    products = degree12_pipeline._product_rows(values, rows, mod)
    assert [name for name, _, _ in products] == [
        "I4_1^3",
        "I6_1^2",
        "I6_1*I6_2",
        "I6_2^2",
        *[f"I4_1*I8_{k}" for k in range(1, 7)],
    ]
    assert np.array_equal(
        products[0][2], 3 * values["I4_1"] ** 2 * rows["I4_1"] % mod)
    assert np.array_equal(
        products[2][2],
        (
            values["I6_1"] * rows["I6_2"]
            + values["I6_2"] * rows["I6_1"]
        ) % mod,
    )


def test_semantic_fingerprint_ignores_only_runtime_measurements():
    first = {
        "rank": 72,
        "seconds": 1.0,
        "nested": {"evaluation_seconds": 2.0, "ids": ["a", "b"]},
    }
    second = {
        "rank": 72,
        "seconds": 99.0,
        "nested": {"evaluation_seconds": 88.0, "ids": ["a", "b"]},
    }
    assert (
        degree12_pipeline.semantic_result_sha256(first)
        == degree12_pipeline.semantic_result_sha256(second)
    )
    second["rank"] = 71
    assert (
        degree12_pipeline.semantic_result_sha256(first)
        != degree12_pipeline.semantic_result_sha256(second)
    )


def test_committed_degree12_certificate_is_internally_consistent():
    with open(RESULT) as stream:
        result = json.load(stream)
    generators = result["generators"]
    assert len(generators) == 62
    assert len(result["products"]) == 10
    assert len(result["degree12_basis"]) == 72
    assert len(result["primitive_candidates_through_degree12"]) == 83
    assert result["discovery"]["degree12_polynomial_rank"] == 72
    assert result["discovery"]["cumulative_functional_rank"] == 81
    assert result["discovery"]["functional_dependencies"] == [
        "I12_61",
        "I12_62",
    ]
    assert result["validation_engine_sha256"] == (
        degree12_pipeline._engine_sha256())
    with open(BENCHMARKS) as stream:
        benchmarks = json.load(stream)
    with open(RESULT, "rb") as stream:
        assert benchmarks["result_sha256"] == hashlib.sha256(
            stream.read()).hexdigest()
    assert benchmarks["result_semantic_sha256"] == (
        degree12_pipeline.semantic_result_sha256(result))

    graph_ids = set()
    for item in generators:
        M = graph_from_label(item["graph"])
        validate_graph(M, 5, 4, True)
        assert canonical_graph_id(M) == item["graph_id"]
        assert item["graph_id"] not in graph_ids
        graph_ids.add(item["graph_id"])

    assert set(result["validation"]) == {"32749", "32719", "32693"}
    for run in result["validation"].values():
        assert len(run["samples"]) == 4
        assert run["degree12_polynomial_space"]["rank"] == 72
        assert len(run["degree12_polynomial_space"]["primitive_pivots"]) == 62
        assert len(run["lorentz_boost"]["matched"]) == 62
        for sample in run["samples"]:
            assert sample["cumulative_rank"] == 81
            assert sample["dependency_ids"] == ["I12_61", "I12_62"]


def test_committed_degree12_basis_loads_in_exact_artifact_order():
    with open(RESULT) as stream:
        result = json.load(stream)
    registry = load_verified_registry_through_degree12(ROOT)
    assert registry.metadata["verified_through_degree"] == 12
    assert [item.id for item in registry.basis(12)] == [
        item["id"] for item in result["degree12_basis"]
    ]
    assert len(registry.basis(12)) == 72
