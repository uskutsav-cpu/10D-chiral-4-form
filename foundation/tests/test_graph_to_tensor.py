"""Level-A intrinsic representatives: graph -> explicit Einstein indices.

A verified contraction graph already defines a coordinate-independent Lorentz
scalar. This module pins the deterministic translation into explicit index
form and the independent dense evaluator that validates it.

The dense evaluator is a genuinely different code path from the repository's
slot-planner: it builds einsum subscripts from the dummy-index assignment and
applies the metric by raising one slot of each contracted pair.

Level A is explicit and intrinsic. It is NOT compact or canonical -- these are
10- and 12-fold contractions with 25 and 30 metric factors respectively.
"""

import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.modp import P  # noqa: E402
from sdinv.forms import selfdual_projector, to_dense, random_form  # noqa: E402
from sdinv.contract import value  # noqa: E402
from sdinv.graph_to_tensor import (  # noqa: E402
    contraction_specification, dense_evaluate, index_specification,
    parse_graph_label)

ROOT = Path(__file__).resolve().parents[1]
EXPORT = ROOT / "results" / "intrinsic_candidates" / "explicit_F_contractions.json"

DEG10_LABEL = ("n10[04^2,06^1,08^1,09^1,15^1,16^1,18^1,19^2,25^2,27^2,28^1,"
               "36^3,37^1,39^1,47^1,48^2,57^1,59^1]")
DEG12_LABEL = ("n12[0-4^2,0-6^2,0-8^1,1-5^2,1-7^1,1-8^1,1-9^1,2-7^2,2-9^1,"
               "2-10^1,2-11^1,3-9^2,3-10^2,3-11^1,4-6^1,4-9^1,4-10^1,5-7^1,"
               "5-8^1,5-11^1,6-10^1,6-11^1,7-8^1,8-11^1]")


def _sample(prime, seed):
    projector = selfdual_projector(10, 5, True, prime)
    raw = random_form(10, 5, np.random.default_rng(seed), prime)
    return to_dense((projector @ raw) % prime, 10, 5, prime)


def test_parser_handles_both_label_dialects():
    """Degree 10 concatenates digits; degree 12 needs dashes past vertex 9."""
    a = parse_graph_label(DEG10_LABEL)
    b = parse_graph_label(DEG12_LABEL)
    assert a.shape == (10, 10) and b.shape == (12, 12)
    assert set(a.sum(1).tolist()) == {5}
    assert set(b.sum(1).tolist()) == {5}


def test_index_assignment_gives_five_slots_per_vertex():
    for label, pairs_expected in ((DEG10_LABEL, 25), (DEG12_LABEL, 30)):
        matrix = parse_graph_label(label)
        slots, pairs = index_specification(matrix)
        assert all(len(s) == 5 for s in slots)
        assert len(pairs) == pairs_expected
        # every dummy name appears exactly twice across all slots
        flat = [n for s in slots for n in s]
        assert len(flat) == 2 * pairs_expected
        assert len(set(flat)) == 2 * pairs_expected


def test_dense_evaluator_matches_graph_evaluator():
    """The independent check on the translation."""
    for label in (DEG10_LABEL, DEG12_LABEL):
        matrix = parse_graph_label(label)
        for prime in (32749, 32719):
            form = _sample(prime, 20260729)
            graph = value(matrix, form, 10, 5, True, prime) % prime
            dense = dense_evaluate(label, form, True, prime)
            assert dense == graph, (
                f"{label[:12]} prime {prime}: dense {dense} != graph {graph}")


def test_scalar_is_homogeneous_of_the_right_degree():
    for label, degree in ((DEG10_LABEL, 10), (DEG12_LABEL, 12)):
        prime = 32749
        form = _sample(prime, 20260730)
        scaled = (form * 3) % prime
        base = dense_evaluate(label, form, True, prime)
        got = dense_evaluate(label, scaled, True, prime)
        assert got == (pow(3, degree, prime) * base) % prime


def test_specification_is_complete_and_unambiguous():
    spec = contraction_specification(DEG12_LABEL)
    assert spec["vertices"] == 12
    assert spec["metric_factor_count"] == 30
    assert len(spec["slots"]) == 12
    assert all(len(v) == 5 for v in spec["slots"].values())
    assert "latex" in spec and spec["latex"]


def test_exported_level_a_records_are_validated_and_honestly_labelled():
    if not EXPORT.exists():
        return
    with EXPORT.open() as stream:
        payload = json.load(stream)
    assert payload["all_validated"] is True
    assert len(payload["classes"]) == 7
    for name, record in payload["classes"].items():
        assert record["level"] == "A"
        assert record["validation"]["all_agree"] is True
        assert record["validation"]["homogeneity_degree_check"] is True
        # Level B and C must not be silently claimed
        assert "NOT DERIVED" in record["level_B_status"]
        assert "NOT DERIVED" in record["level_C_status"]
        assert "not be described as compact" in record["caveat"]
