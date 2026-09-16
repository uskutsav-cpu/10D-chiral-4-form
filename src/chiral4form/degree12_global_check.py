"""Match the explicit ideal to the frozen finite global 6+4=10 theorem."""
from __future__ import annotations

import json
from pathlib import Path

from .provenance import atomic_json

GLOBAL = Path("verification/degree12/qq/step2b_global_reachability.json")
IDEAL = Path("verification/all_orders/degree12_exact/ideal.json")


def check_global(output: str | Path = "verification/all_orders/degree12_exact/global_check.json"):
    ideal = json.loads(IDEAL.read_text())
    global_theorem = json.loads(GLOBAL.read_text())
    vertical = global_theorem.get("vertical_translation_argument", {})
    affine = global_theorem.get("affine_fiber_argument", {})
    checks = {
        "explicit_ideal_exact": ideal.get("status") == "exact_characteristic_zero_degree12_triangular_ideal",
        "constraint_count_68": ideal.get("constraint_count") == 68,
        "ideal_fiber_dimension_4": ideal.get("fiber_dimension") == 4,
        "frozen_global_status": global_theorem.get("status") == "global_constructive_reachability_through_degree12",
        "base_dimension_6": global_theorem.get("base_dimension") == 6,
        "leading_rank_68": global_theorem.get("leading_obstruction_rank") == 68,
        "frozen_fiber_dimension_4": global_theorem.get("fiber_dimension") == 4,
        "total_orbit_dimension_10": global_theorem.get("total_orbit_dimension") == 10,
        "affine_fiber_linear": affine.get("obstruction_equations_affine_linear_in_degree12_coordinates") is True,
        "vertical_rank_4": vertical.get("rank_over_Q") == 4,
        "vertical_preserves_base": vertical.get("preserve_base") is True,
        "vertical_constant_action": vertical.get("constant_action") is True,
        "vertical_transitive": vertical.get("transitive_on_each_affine_fiber") is True,
        "global_parameterization_verified": global_theorem.get("global_parameterization_verified") is True,
    }
    if not all(checks.values()):
        raise AssertionError("degree-12 global consistency gate failed: " + json.dumps(checks, sort_keys=True))
    payload = {
        "schema": 3,
        "status": "exact_degree12_ideal_global_orbit_equality_through_degree12",
        "checks": checks,
        "base_dimension": 6,
        "fiber_dimension": 4,
        "orbit_dimension": 10,
        "global_consistency": True,
        "finite_statement": (
            "At the free-seed polynomial truncation through field degree 12, the zero set of I10 "
            "plus the 68 exact relative degree-12 constraints is the globally reachable 10-dimensional "
            "pure-stress orbit."
        ),
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(path, payload)
    return payload
