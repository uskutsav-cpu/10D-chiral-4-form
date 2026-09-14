#!/usr/bin/env python3
import json
from pathlib import Path
if __name__=="__main__":
    p=Path("verification/all_orders/degree12_explicit_ideal.json")
    if not p.exists():raise SystemExit("run step03 first")
    r=json.loads(p.read_text())
    if not r["status"].startswith("conditional_exact"):
        out={"schema":1,"status":"blocked_no_explicit_degree12_ideal",
             "reason":r["status"]}
    else:
        out={"schema":1,"status":"conditional_degree12_ideal_tangency_from_exact_lifted_equation_space",
             "constraint_count":r["constraint_count"],"leading_rank":r["leading_rank"],
             "conditional":True,
             "remaining_gate":"characteristic-zero physics exactness of the lifted equation space"}
    q=Path("verification/all_orders/degree12_tangency.json");q.write_text(json.dumps(out,indent=2)+"\n")
    print("DEGREE-12 TANGENCY STATUS:",out["status"])
