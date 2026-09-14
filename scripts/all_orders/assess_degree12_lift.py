#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.degree12_all_orders import degree12_lift_readiness

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,
        default=Path("verification/all_orders/degree12_lift_readiness.json"))
    a=p.parse_args(); r=degree12_lift_readiness()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("==============================================")
    print("DEGREE-12 ALL-ORDERS LIFT GATE")
    print("==============================================")
    print("status:",r["status"])
    print("finite local theorem:",r["finite_characteristic_zero_local_theorem"])
    print("finite global theorem:",r["finite_global_reachability_theorem"])
    print("all-orders lift ready:",r["all_orders_lift_ready"])
    print("WROTE:",a.output)
if __name__=="__main__":main()
