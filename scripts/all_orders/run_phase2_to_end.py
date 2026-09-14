#!/usr/bin/env python3
from __future__ import annotations
import subprocess,sys
from pathlib import Path

def main():
    here=Path(__file__).resolve().parent
    for name in (
        "prove_degree10_conditional.py",
        "assess_degree12_lift.py",
        "prove_abstract_generalized_completion.py",
        "plan_degree14.py",
        "final_all_orders_status.py",
    ):
        print("===== RUN",name,"=====",flush=True)
        subprocess.run([sys.executable,str(here/name)],check=True)
    print("==============================================")
    print("PASS: ALL CODED ALL-ORDERS PHASES EXECUTED")
    print("==============================================")
    print("See verification/all_orders/final_status.json for theorem vs blocker status.")
if __name__=="__main__":main()
