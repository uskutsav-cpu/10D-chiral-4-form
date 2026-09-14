#!/usr/bin/env python3
import json
from pathlib import Path
if __name__=="__main__":
    t=Path("verification/all_orders/degree12_tangency.json")
    if not t.exists():raise SystemExit("run step04 first")
    r=json.loads(t.read_text())
    unconditional=r.get("conditional") is False and r["status"].startswith("unconditional")
    out={"schema":1,
         "status":"all_orders_degree12_ideal_theorem" if unconditional else "blocked_degree12_all_orders_promotion",
         "all_orders_cylinder_invariant":bool(unconditional),
         "remaining_gate":None if unconditional else r.get("remaining_gate","explicit exact degree-12 ideal/tangency proof")}
    p=Path("verification/all_orders/degree12_promotion.json");p.write_text(json.dumps(out,indent=2)+"\n")
    print("DEGREE-12 PROMOTION STATUS:",out["status"])
