#!/usr/bin/env python3
import json
from pathlib import Path
from chiral4form.formal_orbit_equality import assess_orbit_equality
if __name__=="__main__":
    r=assess_orbit_equality()
    p=Path("verification/all_orders/formal_orbit_equality.json");p.write_text(json.dumps(r,indent=2)+"\n")
    print("FORMAL ORBIT EQUALITY STATUS:",r["status"])
    print("proved:",r["formal_orbit_equality_proved"])
    if r["remaining_gate"]:print("remaining gate:",r["remaining_gate"])
