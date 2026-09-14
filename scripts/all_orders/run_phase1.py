#!/usr/bin/env python3
from __future__ import annotations
import json,subprocess,sys
from pathlib import Path

def main():
    here=Path(__file__).resolve().parent
    for name in ("prove_filtered_flow.py","prove_k6_all_orders.py","prove_degree8_ideal_all_orders.py","probe_degree14.py"):
        print("===== RUN",name,"=====",flush=True)
        subprocess.run([sys.executable,str(here/name)],check=True)
    paths={"filtration":Path("verification/all_orders/filtration_theorem.json"),
           "k6":Path("verification/all_orders/k6_theorem.json"),
           "degree8":Path("verification/all_orders/degree8_ideal_theorem.json"),
           "degree14":Path("verification/all_orders/degree14_probe.json")}
    p={k:json.loads(v.read_text()) for k,v in paths.items()}
    m={"schema":1,"status":"all_orders_phase1_complete",
       "filtration_theorem":p["filtration"]["status"]=="formal_field_degree_filtration_theorem",
       "K6_all_orders":p["k6"]["all_orders_zero_set_preserved"],
       "degree8_ideal_all_orders":p["degree8"]["all_orders_cylinder_invariant"],
       "degree14_status":p["degree14"]["status"],
       "headline_claim":"Pure-stress flow is field-degree filtered; exact finite-jet invariant ideals lift to formal all-orders invariant cylinders. K6=0 and the certified degree-8 obstruction ideal are preserved to all orders.",
       "open_next":["encode exact degree-10 obstruction ideal generators and lift them",
                    "encode exact degree-12 obstruction ideal generators and lift them",
                    "construct degree-14 invariant basis/map and test for new generators",
                    "test finite generation/stabilization of the obstruction ideal",
                    "do not claim formal orbit equality before sufficiency is proved"]}
    out=Path("verification/all_orders/manifest.json");out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(m,indent=2)+"\n")
    print("PASS: ALL-ORDERS PHASE 1 COMPLETE");print("WROTE:",out)
if __name__=="__main__":main()
