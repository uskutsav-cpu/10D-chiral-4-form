#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path

def load(path):
    return json.loads(Path(path).read_text())

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,
        default=Path("verification/all_orders/final_status.json"))
    a=p.parse_args()

    phase1=load("verification/all_orders/manifest.json")
    d10=load("verification/all_orders/degree10_lift.json")
    d12=load("verification/all_orders/degree12_lift_readiness.json")
    abstract=load("verification/all_orders/abstract_generalized_completion.json")
    d14=load("verification/all_orders/degree14_project_plan.json")

    proven=[
        "formal field-degree filtration theorem",
        "all-orders K6=0 obstruction",
        "all-orders lifted degree-8 obstruction ideal",
        "existence of a finite all-orders polynomial generalized completion under a Hilbert-basis generalized flow class",
    ]
    conditional=[]
    if d10["unconditional_physical_theorem"]:
        proven.append("all-orders lifted degree-10 orbit ideal")
    else:
        conditional.append(
            "degree-10 all-orders ideal lift is exact for the frozen QQ model but conditional on the degree-10 physics/basis map"
        )
    if d12["all_orders_lift_ready"]:
        proven.append("all-orders lifted degree-12 orbit ideal")
    else:
        conditional.append(
            "degree-12 all-orders lift blocked until explicit exact QQ ideal generators and finite ideal tangency are supplied"
        )

    r={
        "schema":1,
        "status":"all_orders_project_status",
        "phase1_complete":phase1["status"]=="all_orders_phase1_complete",
        "proven_results":proven,
        "conditional_results":conditional,
        "degree14_status":d14["status"],
        "major_open_theorem":(
            "find a finitely generated stress-invariant obstruction ideal I_infinity "
            "and prove the formal pure-stress orbit equals the local component Z(I_infinity)_0"
        ),
        "next_hard_math":[
            "unconditional degree-10 physics/basis-map proof",
            "compact exact degree-12 orbit ideal and tangency proof",
            "explicit degree-14 invariant basis/map",
            "degree-14 quotient modulo prolongations",
            "structural induction proving stabilization/finite generation",
            "sufficiency/orbit-equality theorem",
        ],
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("==============================================")
    print("ALL-ORDERS PROJECT STATUS")
    print("==============================================")
    print("proven:")
    for x in proven: print(" -",x)
    print("conditional/open:")
    for x in conditional: print(" -",x)
    print("degree14:",r["degree14_status"])
    print("WROTE:",a.output)
if __name__=="__main__":main()
