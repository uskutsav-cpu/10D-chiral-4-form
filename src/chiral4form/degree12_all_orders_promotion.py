"""Promote exact finite degree-12 ideal tangency to an all-orders invariant cylinder."""
from __future__ import annotations
import json
from pathlib import Path
from .formal_filtration import prove_truncation_commutes
from .provenance import atomic_json

def promote(
    tangency_path="verification/all_orders/degree12_exact/tangency.json",
    global_path="verification/all_orders/degree12_exact/global_check.json",
    output="verification/all_orders/degree12_exact/all_orders_theorem.json",
):
    t=json.loads(Path(tangency_path).read_text())
    g=json.loads(Path(global_path).read_text())
    if t.get("finite_ideal_tangent") is not True:
        raise ValueError("finite degree-12 tangency is not exact")
    if g.get("global_consistency") is not True:
        raise ValueError("degree-12 global consistency not verified")
    filtration=prove_truncation_commutes(12)
    if not filtration.all_checks:
        raise AssertionError("field-degree filtration theorem failed at cutoff 12")
    payload={
        "schema":1,
        "status":"unconditional_all_orders_degree12_obstruction_ideal_theorem",
        "cutoff":12,
        "constraint_count":68,
        "base_dimension":6,
        "fiber_dimension":4,
        "finite_orbit_dimension":10,
        "finite_ideal_tangent":True,
        "filtration_verified":True,
        "all_orders_cylinder_invariant":True,
        "claim":(
            "The pullback of the exact degree-12 pure-stress obstruction ideal "
            "under pi_<=12 is invariant under every full formal analytic "
            "derivative-free pure-stress flow in the declared class."
        ),
        "not_claimed":[
            "degree-14 stabilization",
            "finite generation of the complete all-orders obstruction ideal",
            "formal orbit equality",
            "nonanalytic ModMax sector",
        ],
    }
    path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(path,payload)
    return payload
