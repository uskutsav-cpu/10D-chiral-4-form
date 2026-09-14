#!/usr/bin/env python3
import json
from pathlib import Path
from chiral4form.higher_obstruction import modular_obstruction_report
from chiral4form.normalization_audit import load_records
from chiral4form.registry import Registry
from chiral4form.obstruction_prolongation import prolongation_vectors,quotient_report

if __name__=="__main__":
    reg=Registry.load("runs/all_orders/degree14_basis/registry_degree14.json")
    records=load_records(sorted(Path("runs/all_orders/degree14_map").glob("fields_prime*.json")))
    ideal_path=Path("verification/all_orders/degree12_explicit_ideal.json")
    ideal=json.loads(ideal_path.read_text()) if ideal_path.exists() else None
    reports={}
    for p,record in sorted(records.items()):
        obs=modular_obstruction_report(record,reg,14,p)
        row={"equation_rank":obs["equation_rank"],"obstruction_nullity":obs["obstruction_nullity"],
             "leading_obstruction_rank":obs["leading_obstruction_rank"]}
        if ideal and ideal["status"].startswith("conditional_exact"):
            ids=tuple(record["coordinate_ids"])
            try:
                P,_=prolongation_vectors(ideal,reg,14,p,ids)
                row["quotient"]=quotient_report(obs["kernel"],P,p)
            except ValueError as exc:
                row["quotient"]={"status":"blocked_missing_full_ring_I12_embedding","reason":str(exc)}
        else:
            row["quotient"]={"status":"blocked_missing_explicit_I12"}
        reports[str(p)]=row
    qdims={r.get("quotient",{}).get("new_obstruction_quotient_dimension") for r in reports.values()}
    qdims.discard(None)
    out={"schema":1,"status":"degree14_modular_new_obstruction_quotient",
         "per_prime":reports,"stable_new_dimension":next(iter(qdims)) if len(qdims)==1 else None,
         "characteristic_zero_theorem":False,
         "interpretation":"finite-prime falsification evidence; requires exact I12 and map proof for theorem status"}
    p=Path("verification/all_orders/degree14_quotient.json");p.write_text(json.dumps(out,indent=2)+"\n")
    print("DEGREE-14 NEW-OBSTRUCTION DIMENSION:",out["stable_new_dimension"])
