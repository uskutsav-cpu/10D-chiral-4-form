#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.degree8_orbit import analytic_fields,constraints
from chiral4form.obstruction_lifting import lift_finite_jet_ideal,relative_transport_report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("verification/all_orders/degree8_ideal_theorem.json"))
    a=p.parse_args(); ids,fields=analytic_fields(); qs=constraints(); omega=qs[-1]
    tr=relative_transport_report(fields,omega,{n:(60 if n=="tr1" else 0) for n in fields})
    if not tr["all_verified"]: raise SystemExit("Omega8 transport failed")
    lifted=lift_finite_jet_ideal(fields,qs,cutoff=8,label="degree-8 pure-stress obstruction ideal")
    r={"schema":1,"status":"all_orders_degree8_obstruction_ideal_theorem",
       "coordinate_ids":list(ids),
       "finite_generators":["c_I6_2","c_I8_3","c_I8_4","c_I8_5","c_I8_6","Omega8 = c_I8_2 + 16*c_I4_1^3"],
       "omega8_transport":tr,"finite_jet_lift":lifted.to_json(),
       "all_orders_cylinder_invariant":True,
       "conclusion":"The certified degree-8 obstruction variety lifts to an invariant cylinder in the full formal interaction space.",
       "claim_boundary":"This does not prove that the lifted degree-8 ideal defines the entire all-orders orbit."}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("PASS: ALL-ORDERS DEGREE-8 IDEAL");print("WROTE:",a.output)
if __name__=="__main__":main()
