"""Evidence tracker for obstruction-ideal stabilization.

Finite checks at degree 14, 16, ... can falsify a proposed finite generating
set but cannot by themselves prove all-orders stabilization. This module
enforces that claim boundary.
"""
from __future__ import annotations


def stabilization_evidence(probes):
    """probes: iterable of {'degree': even int, 'new_generator_count': int}."""
    rows = sorted(
        (
            {
                "degree": int(row["degree"]),
                "new_generator_count": int(row["new_generator_count"]),
            }
            for row in probes
        ),
        key=lambda x: x["degree"],
    )
    if any(x["degree"] < 14 or x["degree"] % 2 for x in rows):
        raise ValueError("stabilization probes must be even degrees >=14")
    if any(x["new_generator_count"] < 0 for x in rows):
        raise ValueError("new generator counts must be nonnegative")

    consecutive_zero = 0
    best_zero_run = 0
    for row in rows:
        if row["new_generator_count"] == 0:
            consecutive_zero += 1
            best_zero_run = max(best_zero_run, consecutive_zero)
        else:
            consecutive_zero = 0

    return {
        "schema": 1,
        "status": "finite_stabilization_evidence_only",
        "probes": rows,
        "longest_consecutive_zero_new_generator_run": best_zero_run,
        "finite_generation_proved": False,
        "interpretation": (
            "Zero new-generator quotients at successive finite degrees support "
            "a stabilization conjecture but do not prove an all-orders induction."
        ),
        "required_for_theorem": (
            "a structural induction/Noetherian-module argument showing every "
            "higher-degree obstruction lies in the prolonged lower ideal"
        ),
    }


def finite_generation_gate(*, induction_proof_supplied: bool, proof_reference=None):
    if not induction_proof_supplied:
        return {
            "schema": 1,
            "status": "finite_generation_not_proved",
            "finite_generation_proved": False,
            "remaining_gate": (
                "supply a symbolic induction/module theorem, not merely finite-degree checks"
            ),
        }
    if not proof_reference:
        raise ValueError("a proof reference is required when marking induction supplied")
    return {
        "schema": 1,
        "status": "finite_generation_theorem_external_or_symbolic_proof_supplied",
        "finite_generation_proved": True,
        "proof_reference": str(proof_reference),
    }
