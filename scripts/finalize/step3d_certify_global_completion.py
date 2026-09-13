#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.finite_field import determinant,matrix_rank
from chiral4form.fitting import fields_from_json
from chiral4form.normalization_audit import load_records,audit_model_records
from chiral4form.registry import Registry

DEGREES=(4,6,8,10,12)

def const(poly):
    zero=(0,)*poly.n
    if any(k!=zero for k in poly.terms): return None
    return int(poly.terms.get(zero,0))

def layer_rows(ids,fields,catalogue,registry,d):
    idx={x:i for i,x in enumerate(ids)}
    lower=[idx[x] for dd in DEGREES if dd<d for x in registry.degree_bases[dd]]
    cols=[idx[x] for x in registry.degree_bases[d]]
    cat={x["id"]:x for x in catalogue}
    names=[]; rows=[]
    for name in sorted(fields):
        if int(cat[name]["leading_degree"])!=d: continue
        f=fields[name]
        if any(f[j] for j in lower): continue
        row=[]; ok=True
        for j in cols:
            v=const(f[j])
            if v is None: ok=False; break
            row.append(v)
        if ok: names.append(name); rows.append(row)
    return names,rows

def greedy(names,rows,p,n):
    outn=[]; outr=[]; r=0
    for name,row in zip(names,rows):
        nr=matrix_rank([*outr,row],p,ncols=n)
        if nr>r:
            outn.append(name); outr.append(row); r=nr
            if r==n: break
    return outn,outr

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run",type=Path,default=Path("runs/degree12-generalized-selected-verify"))
    ap.add_argument("--discovery",type=Path,default=Path("verification/degree12/generalized_discovery.json"))
    ap.add_argument("--output",type=Path,default=Path("verification/degree12/generalized_completion_theorem.json"))
    args=ap.parse_args()
    discovery=json.loads(args.discovery.read_text())
    reg=Registry.from_json(json.loads((args.run/"registry.json").read_text()))
    rec=load_records(sorted(args.run.glob("fields_prime*.json")))
    audit=audit_model_records(rec); primes=tuple(sorted(rec))
    target={d:len(reg.degree_bases[d]) for d in DEGREES}
    common=None; per={}
    for p in primes:
        ids,fields=fields_from_json(rec[p]); cat=rec[p]["generator_catalogue"]
        data={}
        for d in DEGREES:
            names,rows=layer_rows(ids,fields,cat,reg,d)
            rank=matrix_rank(rows,p,ncols=target[d])
            print(f"{p} d={d}: constant-rank {rank}/{target[d]} from {len(rows)} generators")
            if rank!=target[d]:
                diag={"schema":1,"status":"strong_global_constant_translation_failed","prime":p,"degree":d,
                      "rank":rank,"target":target[d],"generator_count":len(rows)}
                out=args.output.with_name("generalized_completion_failure.json")
                out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(diag,indent=2)+"\n")
                raise SystemExit(f"strong global criterion failed at prime {p}, degree {d}; wrote {out}")
            data[d]=(names,rows)
        if common is None:
            common={}
            for d in DEGREES:
                names,rows=data[d]; bn,br=greedy(names,rows,p,target[d])
                if len(bn)!=target[d]: raise SystemExit(f"basis extraction failed d={d}")
                common[d]=bn
        minors={}
        for d in DEGREES:
            names,rows=data[d]; by=dict(zip(names,rows))
            sq=[by[x] for x in common[d]]
            det=determinant(sq,p)
            if det==0: raise SystemExit(f"common minor vanished p={p} d={d}")
            minors[str(d)]=int(det)
        per[str(p)]={"ranks":{str(d):target[d] for d in DEGREES},"minor_residues":minors}
    result={"schema":2,"status":"characteristic_zero_global_generalized_completion_through_degree12",
      "normalization_audit":audit,"selected_extras":discovery["selected_extras"],
      "selected_count":int(discovery["selected_count"]),"candidate_class":discovery["candidate_class"],
      "fresh_verification_primes":list(primes),"homogeneous_dimensions":{str(d):target[d] for d in DEGREES},
      "common_global_translation_bases":{str(d):common[d] for d in DEGREES},"per_prime_certificates":per,
      "characteristic_zero_rank_argument":"Each layer has a full-size constant-translation minor nonzero modulo good primes; therefore the exact characteristic-zero layer map has full rank.",
      "global_constructive_algorithm":{"degree_order":[4,6,8,10,12],
        "procedure":"At degree d solve the full-rank constant translation system for the desired degree-d displacement. These controls vanish on lower layers; higher-layer spillover is corrected later."},
      "full_polynomial_interaction_space_reachable_through_degree12":True,"global_generalized_completion":True,
      "scope":"finite polynomial interactions through field degree 12; signed piecewise controls; selected primitive graph invariants used as genuine independent S factors in f(tau,S)",
      "not_claimed":["global cardinality minimality over arbitrary scalar combinations","all-orders completion","nonanalytic completion"]}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: STEP 3D GLOBAL GENERALIZED COMPLETION")
    print("homogeneous ranks:",target)
    print("selected extras:",len(discovery["selected_extras"]))
    print("WROTE:",args.output)
if __name__=="__main__": main()
