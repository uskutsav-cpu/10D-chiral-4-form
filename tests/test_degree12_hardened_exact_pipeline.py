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
    source = inspect.getsource(targets)
    assert "_degree10_orbit_substitution_qq" in source
    assert "_degree10_orbit_substitution_mod_p" in source
    assert "tau_terms" in source
    assert "matrix_series_product" in source
    assert "tr2_degree12_A4_after_certified_degree10_orbit_restriction" in source
    # Two superseded bugs must never return.
    assert 'reg.evaluate("I4_1",form,p,gradient=True' not in source.replace(" ", "")
    assert "240**4" not in source.replace(" ", "")


def test_characteristic_zero_target_proof_has_height_gate():
    source = inspect.getsource(targets)
    assert "_b11_characteristic_zero_bound" in source
    assert "modulus <= 2 * b11_cleared_bound" in source
    assert "_centered_integer" in source
    assert "all_theorem_prime_crosschecks" in source


def test_source_rows_consume_QQ_targets_without_integer_cast():
    source = inspect.getsource(source_rows)
    assert '[Fraction(x) for x in targets["B11"]["coordinates"]]' in source
    assert '[Fraction(x) for x in targets["B13"]["coordinates"]]' in source
    assert '[int(x) for x in targets["B11"]["coordinates"]]' not in source


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
    source = inspect.getsource(targets.recover_exact_targets)
    assert '"B11_theorem_model_primes_used_in_bounded_CRT": []' in source
    assert '"B11_theorem_primes_verification_only": list(theorem_primes)' in source
    # Regression against the circular experimental path that converted theorem
    # rows into CRT sample values and then reused the same rows as holdouts.
    assert 'modular_values[p] = [' not in source or 'equations[p][21][:72]' not in source


def test_compatibility_constant_aliases_are_consistent():
    assert targets.GRAPH_VALUE_BOUND == targets.GRAPH12_SAMPLE_BOUND
