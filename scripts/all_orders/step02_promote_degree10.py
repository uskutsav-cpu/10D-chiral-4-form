#!/usr/bin/env python3
import json
from pathlib import Path
from chiral4form.degree10_all_orders import degree10_conditional_all_orders_certificate
if __name__=="__main__":
    mapcert=Path("verification/all_orders/degree10_map/map_certificate.json")
    if not mapcert.exists():raise SystemExit("run step01 first")
    m=json.loads(mapcert.read_text())
    if m["status"]!="bounded_integer_residual_degree10_physics_map_certificate" or not m["all_integer_residuals_zero"]:
        raise SystemExit("degree-10 map certificate is not unconditional")
    c=degree10_conditional_all_orders_certificate()
    if not c["all_orders_cylinder_invariant_under_loaded_map"]:
        raise SystemExit("degree-10 ideal tangency/lift failed")
    out={"schema":1,"status":"unconditional_all_orders_degree10_ideal_theorem",
         "map_certificate_hash":m["certificate_hash"],"finite_QQ_tangency":c["finite_QQ_tangency"],
         "all_orders_cylinder_invariant":True,
         "conclusion":"The exact degree-10 pure-stress orbit ideal lifts to an all-orders invariant cylinder."}
    p=Path("verification/all_orders/degree10_unconditional.json");p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,indent=2)+"\n")
    print("PASS: UNCONDITIONAL ALL-ORDERS DEGREE-10 IDEAL")
