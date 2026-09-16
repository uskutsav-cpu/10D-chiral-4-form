"""Promote the exact finite degree-12 full ideal to an all-orders invariant cylinder."""
from __future__ import annotations

import json
from pathlib import Path

from .formal_filtration import prove_truncation_commutes
from .provenance import atomic_json

D10 = Path("verification/all_orders/degree10_unconditional.json")
TARGETS = Path("verification/all_orders/degree12_exact/exact_source_targets.json")
SOURCE = Path("verification/all_orders/degree12_exact/source_equation_space.json")
KERNEL = Path("verification/all_orders/degree12_exact/kernel_lift.json")
IDEAL = Path("verification/all_orders/degree12_exact/ideal.json")
TANGENCY = Path("verification/all_orders/degree12_exact/tangency.json")
GLOBAL = Path("verification/all_orders/degree12_exact/global_check.json")


def promote(output: str | Path = "verification/all_orders/degree12_exact/all_orders_theorem.json"):
    d10 = json.loads(D10.read_text())
    targets = json.loads(TARGETS.read_text())
    source = json.loads(SOURCE.read_text())
    kernel = json.loads(KERNEL.read_text())
    ideal = json.loads(IDEAL.read_text())
    tangency = json.loads(TANGENCY.read_text())
    global_check = json.loads(GLOBAL.read_text())
    gates = {
        "degree10_unconditional": d10.get("status") == "unconditional_all_orders_degree10_ideal_theorem",
        "all_36_targets_verified_on_all_theorem_primes": targets.get("all_36_theorem_prime_crosschecks") is True,
        "exact_source_rows": source.get("status") == "exact_characteristic_zero_degree12_source_equation_space",
        "exact_kernel_68": kernel.get("exact") is True and kernel.get("kernel_dimension") == 68 and kernel.get("leading_rank") == 68,
        "exact_relative_ideal_68": ideal.get("status") == "exact_characteristic_zero_degree12_triangular_ideal" and ideal.get("constraint_count") == 68,
        "finite_full_ideal_tangency": tangency.get("relative_full_ideal_tangent") is True,
        "finite_global_orbit_equality": global_check.get("global_consistency") is True,
    }
    if not all(gates.values()):
        raise AssertionError("cannot promote degree 12: " + json.dumps(gates, sort_keys=True))
    filtration = prove_truncation_commutes(12)
    if not filtration.all_checks:
        raise AssertionError("field-degree filtration theorem failed at cutoff 12")

    payload = {
        "schema": 3,
        "status": "unconditional_all_orders_degree12_obstruction_ideal_theorem",
        "cutoff": 12,
        "gates": gates,
        "constraint_count_relative_degree12": 68,
        "base_dimension": 6,
        "fiber_dimension": 4,
        "finite_orbit_dimension": 10,
        "finite_full_ideal_definition": "I12 = I10 + preimage(<Q1,...,Q68>)",
        "finite_full_ideal_tangent": True,
        "filtration_verified": True,
        "filtration_certificate": filtration.to_json(),
        "all_orders_cylinder_invariant": True,
        "claim": (
            "For the declared formal analytic derivative-free pure-stress class, pi_<=12^*(I12) "
            "is invariant under every full formal pure-stress flow. At finite cutoff 12, the free-seed "
            "zero set equals the globally reachable 10-dimensional pure-stress orbit."
        ),
        "not_claimed": [
            "degree-14 obstruction stabilization",
            "finite generation or stabilization of the complete all-orders obstruction tower",
            "formal orbit equality beyond the degree-12 cylinder",
            "nonanalytic ModMax sector",
        ],
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(path, payload)
    return payload
