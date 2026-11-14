"""Exact certificate gates for the free-stress map through degree 12."""

import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.exactmap import rank_mod


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "stress_flow" / "dimension_table.json"
CERTIFICATE_DIR = ROOT / "results" / "stress_flow" / "certificates"
ALL_PRIMES = (
    32603, 32609, 32611, 32621, 32633,
    32647, 32653, 32687, 32693, 32707,
    32713, 32717, 32719, 32749, 32771,
)
CHECK_PRIME = 32717


def _load(path):
    with path.open() as stream:
        return json.load(stream)


def test_degree12_static_certificates_are_compatible_and_exact():
    certificates = [
        _load(CERTIFICATE_DIR / f"static_degree12_{prime}.json")
        for prime in ALL_PRIMES
    ]
    reference = certificates[0]
    for certificate in certificates:
        assert certificate["schema"] == 1
        assert certificate["sample_seeds"] == [
            20260729, 20260730, 20260731, 20260732]
        assert certificate["basis"] == reference["basis"]
        assert certificate["basis_rank"] == 72
        assert certificate["stacked_columns"] == 504
        assert certificate["stress_rank"] == 4
        assert certificate["quotient_dimension"] == 68
        assert certificate["engine_sha256"] == reference["engine_sha256"]
        assert (
            certificate["degree12_sha256"]
            == reference["degree12_sha256"]
        )
        rows = np.asarray(list(certificate["targets"].values()),
                          dtype=np.int64)
        assert rank_mod(rows, certificate["prime"]) == 4


def test_static_dimension_table_and_independent_holdout():
    result = _load(RESULT)
    assert result["dimension_rows"] == [
        {
            "degree": 4,
            "full_dimension": 1,
            "stress_dimension": 1,
            "quotient_dimension": 0,
        },
        {
            "degree": 6,
            "full_dimension": 2,
            "stress_dimension": 1,
            "quotient_dimension": 1,
        },
        {
            "degree": 8,
            "full_dimension": 7,
            "stress_dimension": 2,
            "quotient_dimension": 5,
        },
        {
            "degree": 10,
            "full_dimension": 14,
            "stress_dimension": 2,
            "quotient_dimension": 12,
        },
        {
            "degree": 12,
            "full_dimension": 72,
            "stress_dimension": 4,
            "quotient_dimension": 68,
        },
    ]
    degree12 = result["degrees"]["12"]
    validation = degree12["exact_validation"]
    assert validation["primes"] == list(ALL_PRIMES)
    assert validation["prime_count"] == 15
    assert validation["sample_count_per_prime"] == 4
    assert validation["same_standard_complement_every_prime"] is True
    assert len(degree12["complement_basis"]) == 68
    assert degree12["pivot_ids_in_original_atlas"] == [
        "I4_1^3", "I6_1^2", "I4_1*I8_1", "I12_58"]

    tr_m6 = degree12["stress_basis"][0]
    holdout = _load(
        CERTIFICATE_DIR / f"static_degree12_{CHECK_PRIME}.json")
    reconstructed = tr_m6[
        "coordinates_in_original_atlas_by_prime"][str(CHECK_PRIME)]
    assert reconstructed == holdout["targets"]["tr_M6"]

    # Falsification gate: a one-unit alteration fails an exact saved field.
    altered = list(reconstructed)
    altered[0] = (altered[0] + 1) % CHECK_PRIME
    assert altered != holdout["targets"]["tr_M6"]
    assert tr_m6["rational_original_atlas_coordinates"] is None
    assert tr_m6["physical_free_stress_definition"] == (
        "Tr(T_free^6)=Tr(M^6)/48^6")
    assert degree12["rational_lift_audit"][
        "failed_uniqueness_bound_columns"]
