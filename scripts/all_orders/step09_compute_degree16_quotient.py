#!/usr/bin/env python3
import json
from pathlib import Path
from chiral4form.higher_obstruction import modular_obstruction_report
from chiral4form.normalization_audit import load_records
from chiral4form.registry import Registry
if __name__=="__main__":
    q14=Path("verification/all_orders/degree14_quotient.json")
    if not q14.exists():raise SystemExit("run step08 first")
    q=json.loads(q14.read_text())
    if q.get("stable_new_dimension")!=0:
        out={"schema":1,"status":"degree16_probe_not_triggered","reason":"degree14 introduced new or unresolved obstructions"}
    else:
        reg=Registry.load("runs/all_orders/degree16_basis/registry_degree16.json")
        records=load_records(sorted(Path("runs/all_orders/degree16_map").glob("fields_prime*.json")))
        reports={str(p):modular_obstruction_report(rec,reg,16,p) for p,rec in sorted(records.items())}
        out={"schema":1,"status":"degree16_modular_obstruction_probe",
             "per_prime":{p:{"equation_rank":r["equation_rank"],"obstruction_nullity":r["obstruction_nullity"],
                             "leading_obstruction_rank":r["leading_obstruction_rank"]} for p,r in reports.items()},
             "note":"full Q16_new quotient requires an explicit prolonged I14; this probe is falsification infrastructure"}
    p=Path("verification/all_orders/degree16_probe.json");p.write_text(json.dumps(out,indent=2)+"\n")
    print("DEGREE-16 STATUS:",out["status"])
