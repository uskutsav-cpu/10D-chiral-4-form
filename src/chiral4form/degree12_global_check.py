"""Reconcile the explicit ideal with the frozen global 6+4=10 reachability theorem."""
from __future__ import annotations
import json
from pathlib import Path
from .provenance import atomic_json

VERTICAL=Path("verification/degree12/qq/step2a_vertical_rank_theorem.json")
GLOBAL=Path("verification/degree12/qq/step2b_global_reachability.json")

def check_global(
    ideal_path="verification/all_orders/degree12_exact/ideal.json",
    output="verification/all_orders/degree12_exact/global_check.json",
):
    ideal=json.loads(Path(ideal_path).read_text())
    v=json.loads(VERTICAL.read_text())
    g=json.loads(GLOBAL.read_text())
    checks={
        "ideal_fiber_dimension_4":ideal.get("fiber_dimension")==4,
        "vertical_fiber_dimension_4":v.get("fiber_dimension")==4,
        "vertical_rank_over_Q_4":g.get("vertical_translation_argument",{}).get("rank_over_Q")==4,
        "base_dimension_6":g.get("base_dimension")==6,
        "total_orbit_dimension_10":g.get("total_orbit_dimension")==10,
        "global_parameterization_verified":g.get("global_parameterization_verified") is True,
        "affine_fiber_linear":g.get("affine_fiber_argument",{}).get(
            "obstruction_equations_affine_linear_in_degree12_coordinates"
        ) is True,
    }
    if not all(checks.values()):
        raise AssertionError("degree-12 global consistency failed: "+json.dumps(checks))
    payload={
        "schema":1,
        "status":"degree12_explicit_ideal_global_consistency_verified",
        "checks":checks,
        "base_dimension":6,
        "fiber_dimension":4,
        "orbit_dimension":10,
        "global_consistency":True,
    }
    path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(path,payload)
    return payload
