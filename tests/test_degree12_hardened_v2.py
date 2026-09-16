from __future__ import annotations

import inspect

from chiral4form import degree12_exact_targets as targets
from chiral4form import degree12_exact_source_rows as source_rows
from chiral4form.degree12_exact_common import step1c_record


def test_all_36_source_descriptors_are_present_once():
    rows = targets._all_source_descriptors()
    assert len(rows) == 36
    assert [x["source_row"] for x in rows] == list(range(36))
    assert len([x for x in rows if x["kind"] == "basis"]) == 13
    assert len([x for x in rows if x["kind"] == "dependent"]) == 23


def test_source_row_reconstruction_does_not_use_legacy_target_bound_table():
    source = inspect.getsource(source_rows)
    assert "target_bound" not in source
    assert "all_source_rows" in source
    assert "all_23_step1c_dependencies_verified_exactly_over_QQ" in source


def test_exact_degree10_fields_are_six_dimensional():
    fields = source_rows._exact_reduced_degree10_fields()
    assert fields
    for field in fields.values():
        assert len(field) == 6
        assert all(q.n == 6 for q in field)
    assert source_rows.ORBIT_COMPONENTS == (0, 1, 3, 9, 10, 22)


def test_post_orbit_bounds_exist_for_every_source_target():
    reg = targets.Registry.load(targets.REGISTRY_PATH)
    traces = targets._trace_bound_series(reg)
    substitution = targets._degree10_orbit_substitution_qq(reg)
    for meta in targets._all_source_descriptors():
        denominator, cleared, proof = targets._target_bound(
            reg, meta, traces=traces, substitution=substitution
        )
        assert denominator >= 1
        assert (
            cleared >= 1
            or proof.get("structural_zero_leading_target") is True
        )
        assert int(proof["cleared_integer_bound"]) == cleared


def test_generic_source_evaluator_covers_every_required_generator():
    record = step1c_record()["equation"]
    generators = {x["generator"] for x in record["basis"]}
    generators |= {x["generator"] for x in record["dependencies"]}
    factors = {token for g in generators for token in g.split("*")}
    assert factors <= {"tr1", "tr2", "tr3", "tr4", "tr5", "tr6"}


def test_theorem_primes_are_verification_only_in_v2():
    source = inspect.getsource(targets.recover_exact_targets)
    assert '"theorem_primes_used_in_reconstruction": []' in source
    assert '"theorem_primes_verification_only": list(theorem_primes)' in source


def test_old_raw_zero_bound_path_cannot_return():
    source = inspect.getsource(source_rows)
    assert "zero certified bound" not in source
    assert "_basis_target_bounds" not in source
