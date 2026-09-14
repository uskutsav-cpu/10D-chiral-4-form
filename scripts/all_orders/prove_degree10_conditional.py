#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
from chiral4form.degree10_all_orders import degree10_conditional_all_orders_certificate

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,
        default=Path("verification/all_orders/degree10_lift.json"))
    a=p.parse_args()
    r=degree10_conditional_all_orders_certificate()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("==============================================")
    print("DEGREE-10 ALL-ORDERS LIFT")
    print("==============================================")
    print("status:",r["status"])
    print("loaded-model cylinder invariant:",r["all_orders_cylinder_invariant_under_loaded_map"])
    print("unconditional physical theorem:",r["unconditional_physical_theorem"])
    print("remaining gate:",r["remaining_gate"])
    print("WROTE:",a.output)
if __name__=="__main__":main()
