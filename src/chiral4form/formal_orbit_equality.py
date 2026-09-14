"""Fail-closed gate for formal orbit equality."""
from __future__ import annotations
import json
from pathlib import Path

def assess_orbit_equality(stabilization_path="verification/all_orders/stabilization_theorem.json",
                          sufficiency_path="verification/all_orders/orbit_sufficiency_certificate.json"):
    stab=json.loads(Path(stabilization_path).read_text()) if Path(stabilization_path).exists() else {}
    suff=json.loads(Path(sufficiency_path).read_text()) if Path(sufficiency_path).exists() else None
    finite_generation=stab.get("finite_generation_proved") is True
    suff_ok=bool(suff and suff.get("schema")==1 and
                 suff.get("status")=="formal_control_lifting_verified" and
                 suff.get("machine_checkable") is True and
                 suff.get("compatible_finite_jet_surjectivity") is True and
                 suff.get("tangent_distribution_equals_variety_dimension") is True)
    proved=finite_generation and suff_ok
    return {"schema":1,
        "status":"formal_pure_stress_orbit_equality_theorem" if proved else "formal_orbit_equality_not_yet_proved",
        "finite_generation_available":finite_generation,
        "sufficiency_certificate_present":suff is not None,
        "formal_orbit_equality_proved":proved,
        "theorem":"Orbit_stress(0)=Z(I_infinity)_0" if proved else None,
        "remaining_gate":None if proved else
            ("prove compatible constructive control lifting on every finite jet and "
             "formal inverse-limit sufficiency; obstruction containment alone is not equality")}
