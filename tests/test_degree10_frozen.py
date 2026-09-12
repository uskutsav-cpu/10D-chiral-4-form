import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "verification" / "degree10"


def load(name):
    return json.loads((BASE / name).read_text())


def test_degree10_rational_model_prime_separation():
    x = load("rational_model.json")

    assert x["fit_primes"] == [
        30139,
        30161,
        30169,
        30181,
        30187,
    ]

    assert x["holdout_prime"] == 30197
    assert x["holdout_prime"] not in x["fit_primes"]
    assert x["coefficients_reconstructed"] == 87


def test_degree10_obstruction_dimension_lock():
    x = load("obstruction_certificate.json")

    assert x["obstruction_nullity"] == 12
    assert x["degree10_leading_rank"] == 12
    assert x["reduced_ambient_dimension"] == 18
    assert x["constraint_upper_dimension"] == 6
    assert x["free_seed_direct_control_rank"] == 6
    assert x["exact_local_dimension_locked"] is True
    assert x["local_orbit_dimension"] == 6

    rows = x["obstructions_integer_coordinates"]

    assert len(rows) == 12
    assert all(len(row) == 19 for row in rows)


def test_degree10_global_constructive_certificate():
    x = load("global_reachability.json")

    assert (
        x["status"]
        == "global_constructive_reachability_in_degree10_truncation"
    )

    assert x["orbit_dimension"] == 6
    assert x["global_parameterization_verified"] is True
    assert x["all_vector_fields_tangent"] is True
    assert x["lower_triangular_control_system_verified"] is True
    assert x["degree10_translation_rank"] == 2

    assert Fraction(x["translation_determinant"]) != 0

    assert set(x["translation_controls"]) == {
        "tr2*tr3",
        "tr5",
    }

    assert set(x["stage_amounts"]) == {
        "tr2",
        "tr3",
        "tr4",
        "tr2*tr2",
        "tr2*tr3",
        "tr5",
    }


def test_degree10_linear_hull_is_not_orbit():
    x = load("reachability_analysis.json")

    assert x["linear_hull"]["rank"] == 8
    assert x["linear_hull"]["stabilized"] is True

    assert x["lie"]["rank"] == 6
    assert x["lie"]["completed_depth"] == 3
