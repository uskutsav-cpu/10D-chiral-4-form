"""The three published degree-12 invariants of equation (4.25).

Source: Cederwall, Hutomo, Kuzenko, Lechner, Sorokin, "Some remarks on
invariants", J. Phys. A 59 (2026) 065203; cross-checked against
arXiv:2509.14350v2, equation (4.25), journal PDF page 18.

Transcribed from the rendered source, not from memory.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.modp import P, ALT_P  # noqa: E402
from sdinv.forms import selfdual_projector, to_dense, random_form  # noqa: E402
from sdinv.published_degree12_invariants import (  # noqa: E402
    PUBLISHED_DEGREE12, evaluate_all)


def _sample(prime, seed):
    projector = selfdual_projector(10, 5, True, prime)
    raw = random_form(10, 5, np.random.default_rng(seed), prime)
    return to_dense((projector @ raw) % prime, 10, 5, prime)


def test_all_three_equation_4_25_structures_are_present():
    assert set(PUBLISHED_DEGREE12) == {"P12_01", "P12_02", "P12_03"}
    for spec in PUBLISHED_DEGREE12.values():
        assert spec["field_degree"] == 12
        assert sum(spec["blocks"].values()) == 6, (
            "six quadratic blocks at five-form degree 12")


def test_homogeneity_degree_twelve():
    for prime in (P, ALT_P):
        form = _sample(prime, 20260729)
        base = evaluate_all(form, prime)
        for c in (2, 3):
            scaled = evaluate_all((form * c) % prime, prime)
            for name, value in base.items():
                assert scaled[name] == (pow(c, 12, prime) * value) % prime, (
                    f"{name} is not homogeneous of degree 12 at prime {prime}")


def test_values_are_nonzero_and_prime_consistent_in_structure():
    """A structure that vanished identically would be uninformative."""
    for prime in (P, ALT_P):
        for seed in (20260729, 20260730):
            values = evaluate_all(_sample(prime, seed), prime)
            assert len(values) == 3
            assert any(v % prime for v in values.values()), (
                "all three vanished, which would make the map vacuous")


def test_backends_agree_where_a_reference_path_exists():
    prime = P
    form = _sample(prime, 20260731)
    optimized = evaluate_all(form, prime, backend="optimized")
    reference = evaluate_all(form, prime, backend="reference")
    assert optimized == reference


def test_certified_q12_projection_is_zero_and_non_vacuous():
    """Both primes must agree, and every structure must actually SOLVE.

    A structure that failed to solve would also contribute rank 0, for an
    entirely uninformative reason. The status field separates the cases.
    """
    import json
    from pathlib import Path
    path = (Path(__file__).resolve().parents[1] / "results"
            / "intrinsic_candidates" / "published_degree12_map.json")
    if not path.exists():
        return
    with path.open() as stream:
        payload = json.load(stream)

    assert payload["consistent"] is True
    assert payload["Q12_rank_from_published"] == 0
    assert len(payload["per_prime"]) >= 2, "need at least two primes"

    for prime, record in payload["per_prime"].items():
        assert record["dim_Q12"] == 4
        assert record["Q12_rank_from_published"] == 0
        assert set(record["projections"]) == {"P12_01", "P12_02", "P12_03"}
        for name, proj in record["projections"].items():
            assert proj["status"] == "solved", (
                f"prime {prime}: {name} did not solve; a rank of 0 would then "
                f"say nothing about the quotient")
            assert len(proj["quotient_vector"]) == 4
            assert all(v == 0 for v in proj["quotient_vector"])
            assert proj["nonzero"] is False


def test_certification_does_not_overclaim():
    """The doc must carry the forbidden-inference list."""
    from pathlib import Path
    doc = (Path(__file__).resolve().parents[1] / "docs"
           / "PUBLISHED_DEGREE12_MAP.md")
    if not doc.exists():
        return
    # normalise whitespace: the permitted wording is a blockquote that wraps
    # mid-phrase, so a raw substring search would fail on the line break
    body = " ".join(doc.read_text().replace(">", " ").split())
    assert "lie inside the computed reachable closure D12" in body
    for forbidden in ("complete M/N degree-12 space has zero image",
                      "no compact tensor representation",
                      "are redundant",
                      "no generalized flow can generate Q12"):
        assert forbidden in body, (
            f"the doc must explicitly disclaim: {forbidden}")
