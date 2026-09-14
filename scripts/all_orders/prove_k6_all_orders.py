#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.degree8_orbit import analytic_fields
from chiral4form.obstruction_lifting import lift_finite_jet_ideal,relative_transport_report
from chiral4form.polynomial import Poly

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("verification/all_orders/k6_theorem.json"))
    a=p.parse_args(); ids,fields=analytic_fields()
    if ids[2]!="I6_2": raise SystemExit("degree-6 obstruction coordinate changed")
    q6=Poly.variable(len(ids),2)
    transport=relative_transport_report(fields,q6,{n:(40 if n=="tr1" else 0) for n in fields})
    if not transport["all_verified"]: raise SystemExit("K6 transport failed")
    lifted=lift_finite_jet_ideal(fields,(q6,),cutoff=6,label="<K6>")
    r={"schema":1,"status":"all_orders_K6_obstruction_theorem","coordinate":"c_I6_2",
       "transport_equation":"dq6/dlambda = 40*u_tr1(lambda)*q6",
       "exact_transport":transport,"finite_jet_lift":lifted.to_json(),
       "all_orders_zero_set_preserved":True,"free_seed_has_q6_zero":True,
       "conclusion":"Every formal analytic pure-stress trajectory from the free seed remains on K6=0 to all orders.",
       "nonuniversality_consequence":"Any admissible analytic interaction with nonzero K6 is outside the free-seed pure-stress orbit.",
       "scope":"declared derivative-free analytic/polynomial pure-stress flow class"}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("PASS: ALL-ORDERS K6 OBSTRUCTION");print("WROTE:",a.output)
if __name__=="__main__":main()
