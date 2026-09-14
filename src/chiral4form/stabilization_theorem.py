"""Fail-closed theorem gate for finite generation/stabilization."""
from __future__ import annotations
import json
from pathlib import Path

def assess_stabilization(q14_path="verification/all_orders/degree14_quotient.json",
                         q16_path="verification/all_orders/degree16_probe.json",
                         induction_path="verification/all_orders/stabilization_induction_certificate.json"):
    q14=json.loads(Path(q14_path).read_text()) if Path(q14_path).exists() else {}
    q16=json.loads(Path(q16_path).read_text()) if Path(q16_path).exists() else {}
    finite_evidence=(q14.get("stable_new_dimension")==0 and
        q16.get("status") in {"degree16_modular_obstruction_probe","degree16_new_obstruction_quotient_zero"})
    induction=None
    if Path(induction_path).exists():
        induction=json.loads(Path(induction_path).read_text())
    proved=bool(induction and induction.get("schema")==1 and
                induction.get("status")=="symbolic_induction_verified" and
                induction.get("machine_checkable") is True)
    return {"schema":1,
        "status":"all_orders_obstruction_stabilization_theorem" if proved else "stabilization_not_yet_proved",
        "finite_degree_evidence":finite_evidence,
        "degree14_new_dimension":q14.get("stable_new_dimension"),
        "degree16_status":q16.get("status"),
        "induction_certificate_present":induction is not None,
        "finite_generation_proved":proved,
        "remaining_gate":None if proved else
            ("supply a machine-checkable symbolic induction/module certificate proving "
             "every higher obstruction lies in the prolonged finite ideal; finite degree checks alone do not suffice")}
