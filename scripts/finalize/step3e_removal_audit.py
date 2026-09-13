#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.finite_field import matrix_rank
from chiral4form.fitting import fields_from_json
from chiral4form.normalization_audit import load_records
from chiral4form.registry import Registry
DEGREES=(4,6,8,10,12)

def req(cat): return frozenset(x[2:] for x in cat["factors"] if x.startswith("S:"))
def const(poly):
    z=(0,)*poly.n
    if any(k!=z for k in poly.terms): return None
    return int(poly.terms.get(z,0))
def rows(ids,fields,cat,reg,d,sel):
    idx={x:i for i,x in enumerate(ids)}
    low=[idx[x] for dd in DEGREES if dd<d for x in reg.degree_bases[dd]]
    cols=[idx[x] for x in reg.degree_bases[d]]
    out=[]
    for c in cat:
        if int(c["leading_degree"])!=d or not req(c)<=sel: continue
        f=fields[c["id"]]
        if any(f[j] for j in low): continue
        v=[]; ok=True
        for j in cols:
            q=const(f[j])
            if q is None: ok=False; break
            v.append(q)
        if ok: out.append(v)
    return out
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run",type=Path,default=Path("runs/degree12-generalized-selected-verify"))
    ap.add_argument("--completion",type=Path,default=Path("verification/degree12/generalized_completion_theorem.json"))
    ap.add_argument("--output",type=Path,default=Path("verification/degree12/generalized_removal_minimality.json"))
    args=ap.parse_args()
    th=json.loads(args.completion.read_text()); selected=tuple(th["selected_extras"])
    target={int(k):int(v) for k,v in th["homogeneous_dimensions"].items()}
    reg=Registry.from_json(json.loads((args.run/"registry.json").read_text()))
    rec=load_records(sorted(args.run.glob("fields_prime*.json"))); primes=tuple(sorted(rec))
    result={}
    for removed in selected:
        remaining=frozenset(x for x in selected if x!=removed)
        per={}; all_drop=True
        for p in primes:
            ids,fields=fields_from_json(rec[p]); cat=rec[p]["generator_catalogue"]
            ranks={}
            for d in DEGREES:
                rr=rows(ids,fields,cat,reg,d,remaining)
                ranks[d]=matrix_rank(rr,p,ncols=target[d])
            drops=[d for d in DEGREES if ranks[d]<target[d]]
            if not drops: all_drop=False
            per[str(p)]={"ranks":{str(d):ranks[d] for d in DEGREES},"drop_degrees":drops}
        result[removed]={"all_fresh_primes_lose_completion":all_drop,"per_prime":per}
        print(removed, "NONREDUNDANT" if all_drop else "REMOVABLE")
    nonred=[x for x,v in result.items() if v["all_fresh_primes_lose_completion"]]
    removable=[x for x in selected if x not in nonred]
    payload={"schema":2,"status":"generalized_removal_audit_complete","selected_extras":list(selected),
      "fresh_verification_primes":list(primes),"removal_results":result,
      "fresh_prime_nonredundant_extras":nonred,"removable_under_strong_criterion":removable,
      "all_selected_extras_fresh_prime_nonredundant":len(nonred)==len(selected),
      "claim_scope":"removal audit for strong layerwise global constant-translation completion within the declared primitive graph catalogue",
      "not_a_claim_of":["global minimality over arbitrary scalar combinations","all-orders minimality"]}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(payload,indent=2)+"\n")
    print("PASS: STEP 3E REMOVAL AUDIT")
    print("nonredundant:",len(nonred),"removable:",len(removable))
    print("WROTE:",args.output)
if __name__=="__main__": main()
