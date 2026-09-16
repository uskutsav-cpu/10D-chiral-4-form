from __future__ import annotations

import inspect

from chiral4form import degree12_exact_targets as targets
from chiral4form import degree12_exact_source_rows as source_rows
from chiral4form.degree12_exact_common import (
    EXPECTED_B11,
    EXPECTED_B13,
    EXPECTED_SOURCE_ROWS,
    step1c_record,
    theorem_record,
)
from chiral4form.degree12_missing_targets_exact import recover_missing_targets


def test_frozen_degree12_theorem_signature():
    theorem = theorem_record()
    assert theorem["qq_theorem"] is True
    assert theorem["obstruction_nullity"] == 68
    assert theorem["degree12_leading_obstruction_rank"] == 68
    assert theorem["local_orbit_dimension"] == 10
    assert len(theorem["primes"]) == 15


def test_corrected_source_row_metadata():
    record = step1c_record()
    basis = {row["basis_id"]: row for row in record["equation"]["basis"]}
    assert EXPECTED_SOURCE_ROWS == (0, 1, 2, 3, 4, 9, 16, 17, 19, 20, 21, 28, 35)
    assert EXPECTED_B11 == (21, "tr2", "A^4")
    assert EXPECTED_B13 == (35, "tr6", "1")
    assert (basis["B11"]["source_row"], basis["B11"]["generator"], basis["B11"]["output_monomial_text"]) == EXPECTED_B11
    assert (basis["B13"]["source_row"], basis["B13"]["generator"], basis["B13"]["output_monomial_text"]) == EXPECTED_B13


def test_legacy_missing_target_entrypoint_routes_to_hardened_code():
    assert recover_missing_targets is targets.recover_exact_targets


def test_b11_uses_orbit_restriction_not_raw_gradient_square_guess():
    from chiral4form import degree12_exact_targets as targets

    reg = targets.Registry.load(targets.REGISTRY_PATH)
    traces = targets._trace_bound_series(reg)
    substitution = targets._degree10_orbit_substitution_qq(reg)

    assert substitution
    assert all(q.n == 6 for q in substitution.values())

    candidates = [
        meta
        for meta in targets._all_source_descriptors()
        if meta.get("basis_id") == "B11"
    ]

    assert len(candidates) == 1
    b11 = candidates[0]

    assert b11.get("generator") == "tr2"
    assert b11.get("output_monomial_text") == "A^4"

    denominator, cleared, proof = targets._target_bound(
        reg,
        b11,
        traces=traces,
        substitution=substitution,
    )

    assert denominator >= 1
    assert cleared > 0
    assert int(proof["cleared_integer_bound"]) == cleared


def test_characteristic_zero_target_proof_has_height_gate():
    import json
    import math
    from pathlib import Path
    from chiral4form import degree12_exact_targets as targets

    reg = targets.Registry.load(targets.REGISTRY_PATH)
    traces = targets._trace_bound_series(reg)
    substitution = targets._degree10_orbit_substitution_qq(reg)

    bounds = []
    for meta in targets._all_source_descriptors():
        _, cleared, proof = targets._target_bound(
            reg,
            meta,
            traces=traces,
            substitution=substitution,
        )
        assert int(proof["cleared_integer_bound"]) == cleared

        if proof.get("structural_zero_leading_target"):
            assert cleared == 0
            zero = proof["zero_certificate"]
            assert zero["all_theorem_primes_leading_zero"] is True
        else:
            assert cleared > 0
            bounds.append(cleared)

    assert bounds
    max_bound = max(bounds)

    cert = Path(targets.CERT_PATH)
    assert cert.exists(), cert
    payload = json.loads(cert.read_text())

    direct_prime_lists = []

    def walk(obj, path=()):
        if isinstance(obj, dict):
            for key, value in obj.items():
                p = path + (str(key),)
                joined = ".".join(p).lower()

                if (
                    isinstance(value, list)
                    and "direct" in joined
                    and "prime" in joined
                    and value
                    and all(isinstance(x, int) for x in value)
                ):
                    direct_prime_lists.append(value)

                walk(value, p)

        elif isinstance(obj, list):
            for i, value in enumerate(obj):
                walk(value, path + (str(i),))

    walk(payload)

    assert direct_prime_lists, (
        "exact target certificate must record the direct physics primes "
        "used for bounded CRT"
    )

    # Use the largest recovered direct-prime set.
    direct_primes = max(direct_prime_lists, key=len)
    modulus = math.prod(int(p) for p in direct_primes)

    # Characteristic-zero uniqueness gate.
    assert modulus > 2 * max_bound


def test_source_rows_consume_QQ_targets_without_integer_cast():
    from chiral4form.degree12_exact_source_rows import reconstruct_source_rows

    source, kernel = reconstruct_source_rows()

    # Successful exact reconstruction is the semantic regression:
    # the QQ source coordinates survive into the exact 81D equations.
    assert int(source["equation_rank"]) == 13
    assert int(kernel["kernel_dimension"]) == 68
    assert int(kernel["leading_rank"]) == 68


def test_b11_height_bound_is_derived_from_registry_and_closes_under_envelope():
    reg = targets.Registry.load(targets.REGISTRY_PATH)
    audit = targets._derived_tau_entry_bound(reg)
    assert audit["derived_max"] <= targets.TAU_TERM_ENTRY_BOUND
    assert audit["gradient_pair_bound"] > 0
    denominator, cleared, proof = targets._b11_characteristic_zero_bound(reg)
    assert denominator >= 1
    assert cleared >= 1
    assert proof["orbit_denominator_clearer"] == str(denominator)
    assert proof["B11_cleared_integer_bound"] == str(cleared)


def test_b11_crt_uses_direct_physics_and_theorem_primes_are_verification_only():
    import json
    from pathlib import Path
    from chiral4form import degree12_exact_targets as targets

    cert = Path(targets.CERT_PATH)
    assert cert.exists(), cert

    payload = json.loads(cert.read_text())

    direct = payload["direct_physical_primes"]
    used_for_reconstruction = payload["theorem_primes_used_in_reconstruction"]
    verification_only = payload["theorem_primes_verification_only"]

    assert isinstance(direct, list)
    assert direct
    assert used_for_reconstruction == []
    assert verification_only

    assert payload["all_36_direct_physical_crosschecks"]
    assert payload["all_36_theorem_prime_crosschecks"]


def test_compatibility_constant_aliases_are_consistent():
    assert targets.GRAPH_VALUE_BOUND == targets.GRAPH12_SAMPLE_BOUND
