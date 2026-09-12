import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

BASE = (
    ROOT
    / "verification"
    / "degree12"
    / "modular"
)


def load(name):
    return json.loads(
        (BASE / name).read_text()
    )


def test_degree12_eight_prime_structure_stable():

    x = load("obstruction_scan.json")

    assert (
        x["status"]
        == "stable_eight_prime_modular_obstruction_structure"
    )

    reports = x["reports"]

    assert [r["prime"] for r in reports] == [
        30203,
        30211,
        30223,
        30241,
        30253,
        30259,
        30269,
        30271,
    ]

    reference_leading = reports[0]["leading_pivots"]
    reference_kernel = reports[0]["kernel_pivots"]

    for r in reports:

        assert r["equation_count"] == 36
        assert r["equation_rank"] == 13

        assert r["obstruction_nullity"] == 68
        assert r["leading_rank"] == 68

        assert r["direct_control_rank"] == 10
        assert r["upper_dimension"] == 10
        assert r["dimension_match"] is True

        assert r["leading_pivots"] == reference_leading
        assert r["kernel_pivots"] == reference_kernel


def test_degree12_modular_dimension_sandwich():

    x = load("manifest.json")

    assert x["candidate_weight100_dimension"] == 81
    assert x["reduced_coordinate_dimension"] == 78
    assert x["degree12_coordinate_dimension"] == 72

    assert x["obstruction_dimension"] == 68
    assert x["degree12_leading_rank"] == 68

    assert x["direct_control_rank_at_free_seed"] == 10
    assert x["modular_upper_orbit_dimension"] == 10

    assert (
        x["dimension_sandwich_closes_at_every_tested_prime"]
        is True
    )


def test_degree12_not_falsely_promoted_to_QQ():

    x = load("manifest.json")

    assert (
        "characteristic-zero degree-12 orbit dimension"
        in x["not_yet_claimed"]
    )

    diagnostic = load(
        "rational_reconstruction_diagnostic.json"
    )

    assert diagnostic["failure_count"] > 0
