#!/usr/bin/env python3
import json
from pathlib import Path
from chiral4form.stabilization_theorem import assess_stabilization
if __name__=="__main__":
    r=assess_stabilization()
    p=Path("verification/all_orders/stabilization_theorem.json");p.write_text(json.dumps(r,indent=2)+"\n")
    print("STABILIZATION STATUS:",r["status"])
    print("finite generation proved:",r["finite_generation_proved"])
    if r["remaining_gate"]:print("remaining gate:",r["remaining_gate"])
