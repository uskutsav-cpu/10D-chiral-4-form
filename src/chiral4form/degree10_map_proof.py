"""Bounded-residual characteristic-zero proof of the degree-10 physics/basis map."""
from __future__ import annotations
import json
from fractions import Fraction
from math import lcm
from pathlib import Path
import numpy as np
from .checkpoints import CheckpointStore
from .finite_field import matrix_rank
from .fitting import fields_from_json
from .forms import dense,hodge_compact,layout
from .provenance import atomic_json,semantic_hash
from .registry import Registry
from .stress import evaluate_generators
from .tensors import TensorBudget
from .weighted import stress_generators

DEGREE=10
UPPER={4:1,6:2,8:7,10:14}
GRAD_DEN=240
TAU_DEN=GRAD_DEN**2
POINT_SEED=202609141001

def _is_prime(n):
    if n<2:return False
    if n%2==0:return n==2
    d=3
    while d*d<=n:
        if n%d==0:return False
        d+=2
    return True

def proof_primes():
    return tuple(n for n in range(65521,30000,-2) if _is_prime(n) and n not in (2,3,5))

def graph_value_bound(d): return 10**(5*d//2)
def gradient_component_bound(d): return d*10**(5*d//2-5)

def _add(table,key,value): table[key]=table.get(key,0)+int(value)

def coefficient_ids(reg):
    return tuple(i for d in sorted(reg.degree_bases) if d<=DEGREE for i in reg.degree_bases[d])

def tau_scaled_bounds(reg):
    ids=coefficient_ids(reg)
    out={(2,()):TAU_DEN*10**4}
    for name in ids:
        d=reg.items[name].degree
        _add(out,(d,(name,)),TAU_DEN*(d-2)*graph_value_bound(d))
    grad=[name for name in ids if reg.items[name].degree+2<=DEGREE]
    for i,a in enumerate(grad):
        for b in grad[i:]:
            da,db=reg.items[a].degree,reg.items[b].degree
            d=da+db-2
            if d>DEGREE:continue
            symmetry=2 if a!=b else 1
            bound=TAU_DEN*25*symmetry*10**4*gradient_component_bound(da)*gradient_component_bound(db)
            _add(out,(d,tuple(sorted((a,b)))),bound)
    return out

def matrix_bound_product(left,right):
    out={}
    for (da,ma),a in left.items():
        for (db,mb),b in right.items():
            if da+db<=DEGREE:_add(out,(da+db,tuple(sorted(ma+mb))),10*a*b)
    return out

def scalar_bound_product(left,right):
    out={}
    for (da,ma),a in left.items():
        for (db,mb),b in right.items():
            if da+db<=DEGREE:_add(out,(da+db,tuple(sorted(ma+mb))),a*b)
    return out

def trace_scaled_bounds(reg):
    tau=tau_scaled_bounds(reg);ids=coefficient_ids(reg)
    traces={"tr1":{(reg.items[i].degree,(i,)):TAU_DEN*10*(reg.items[i].degree-2)*graph_value_bound(reg.items[i].degree) for i in ids}}
    power={(0,()):1}
    for k in range(1,6):
        power=matrix_bound_product(power,tau)
        if k>=2: traces[f"tr{k}"]={key:10*value for key,value in power.items()}
    return traces

def generator_scaled_bounds(reg):
    traces=trace_scaled_bounds(reg);out={}
    for gen in stress_generators(DEGREE):
        series={(0,()):1};tau_power=0
        for factor in gen.factors:
            k=int(factor[2:]);tau_power+=k
            series=scalar_bound_product(series,traces[factor])
        scale=TAU_DEN**tau_power
        for (d,m),bound in series.items():
            out[(gen.id,d,m)]={"bound":int(bound),"scale":int(scale)}
    return out

def model_coefficients(reg,model):
    ids,fields=fields_from_json(model)
    if tuple(ids)!=coefficient_ids(reg): raise ValueError("model coordinates disagree with registry")
    result={}
    for g,field in fields.items():
        for j,f in enumerate(field):
            d=reg.items[ids[j]].degree
            for ex,q in f.terms.items():
                mono=tuple(sorted(name for name,e in zip(ids,ex) for _ in range(e)))
                result.setdefault((g,d,mono),{})[ids[j]]=Fraction(q)
    return ids,fields,result

def residual_bound_report(reg,model):
    _,_,candidate=model_coefficients(reg,model);actual=generator_scaled_bounds(reg)
    rows={};worst=0
    for key in sorted(set(actual)|set(candidate)):
        g,d,mono=key;a=actual.get(key,{"bound":0,"scale":1});qs=candidate.get(key,{})
        den=1
        for q in qs.values():den=lcm(den,q.denominator)
        D=lcm(int(a["scale"]),den)
        B=(D//int(a["scale"]))*int(a["bound"])+sum(abs(int(q*D))*graph_value_bound(d) for q in qs.values())
        worst=max(worst,B)
        rows["|".join((g,str(d),",".join(mono)))]={
            "generator":g,"degree":d,"monomial":list(mono),"clearing_denominator":D,
            "integer_residual_bound":B,"actual_scale":int(a["scale"]),
            "actual_scaled_bound":int(a["bound"]),
            "candidate_nonzero_coordinates":{k:str(v) for k,v in qs.items()}}
    return {"worst_integer_residual_bound":worst,"rows":rows,"key_count":len(rows)}

def fixed_form(electric,p):
    if len(electric)!=126 or any(type(x)is not int or abs(x)>1 for x in electric):
        raise ValueError("fixed tensors require 126 signed integer electric components")
    raw=np.zeros(252,dtype=np.int64);raw[layout()["electric"]]=electric
    return dense((raw+hodge_compact(raw,p))%p,p)

def candidate_value(key,candidate,basis_values,reg,p):
    _,d,_=key;total=0
    for name,q in candidate.get(key,{}).items():
        j=reg.degree_bases[d].index(name)
        total=(total+(q.numerator*pow(q.denominator,-1,p)%p)*basis_values[str(d)][j])%p
    return total

def choose_points(reg,prime,max_points=40):
    rng=np.random.default_rng(POINT_SEED)
    selected=[];mats={d:[] for d in UPPER}
    budget=TensorBudget(8*1024**3,200_000_000_000,True)
    for _ in range(max_points):
        e=rng.integers(-1,2,126,dtype=np.int64).tolist();form=fixed_form(e,prime);cache={}
        vals={d:reg.basis_values(d,form,prime,cache=cache,budget=budget) for d in UPPER}
        old={d:matrix_rank(mats[d],prime,ncols=UPPER[d]) for d in UPPER}
        new={d:matrix_rank(mats[d]+[vals[d]],prime,ncols=UPPER[d]) for d in UPPER}
        if any(new[d]>old[d] for d in UPPER):
            selected.append(e)
            for d in UPPER:mats[d].append(vals[d])
        if all(matrix_rank(mats[d],prime,ncols=UPPER[d])==UPPER[d] for d in UPPER):
            return selected
    raise RuntimeError("failed to find injective fixed evaluation points")

def run(output,registry_path="data/fixtures/low_degree_registry.json",model_path="verification/degree10/rational_model.json",progress=print):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    reg=Registry.load(registry_path);model=json.loads(Path(model_path).read_text())
    if any(len(reg.degree_bases.get(d,()))!=n for d,n in UPPER.items()):
        raise ValueError("registry dimensions are not 1,2,7,14")
    _,_,candidate=model_coefficients(reg,model);bounds=residual_bound_report(reg,model);B=bounds["worst_integer_residual_bound"]
    denominator_lcm=1
    for row in candidate.values():
        for q in row.values(): denominator_lcm=lcm(denominator_lcm,q.denominator)
    primes=[];M=1
    for p in proof_primes():
        if denominator_lcm%p==0 or TAU_DEN%p==0: continue
        primes.append(p);M*=p
        if len(primes)>=3 and M>2*B:break
    if M<=2*B:raise RuntimeError("insufficient proof-prime product")
    points=choose_points(reg,primes[0])
    budget=TensorBudget(8*1024**3,200_000_000_000,True)
    support=generator_scaled_bounds(reg);keys=sorted(set(support)|set(candidate))
    namespace={"schema":1,"engine":"degree10_bounded_map_v1","registry":reg.fingerprint,
               "model_hash":semantic_hash(model),"points":points,"primes":primes,"bounds":bounds}
    store=CheckpointStore(output/"checkpoints",namespace);samples=[]
    for p in primes:
        for i,electric in enumerate(points):
            def compute():
                form=fixed_form(electric,p);cache={}
                targets,_=evaluate_generators(reg,form,p,DEGREE,(),cache=cache,budget=budget)
                basis={str(d):reg.basis_values(d,form,p,cache=cache,budget=budget) for d in UPPER}
                residuals=[]
                for key in keys:
                    actual=int(targets.get(key,0))%p;expected=candidate_value(key,candidate,basis,reg,p)
                    meta=bounds["rows"]["|".join((key[0],str(key[1]),",".join(key[2])))]
                    residuals.append({"generator":key[0],"degree":key[1],"monomial":list(key[2]),
                        "cleared_residual":meta["clearing_denominator"]*(actual-expected)%p})
                return {"prime":p,"electric":electric,"basis_values":basis,"residuals":residuals}
            s,cached=store.run("map_sample",{"prime":p,"point":i},compute)
            if any(r["cleared_residual"] for r in s["residuals"]):
                raise AssertionError(f"degree-10 map refuted at prime {p}, point {i}")
            samples.append(s);progress(f"degree10 map prime {p}: {i+1}/{len(points)} {'cached' if cached else 'computed'}",flush=True)
    first=[s for s in samples if s["prime"]==primes[0]]
    injective={str(d):matrix_rank([s["basis_values"][str(d)] for s in first],primes[0],ncols=n)==n for d,n in UPPER.items()}
    if not all(injective.values()):raise AssertionError("evaluation set not injective")
    record={"schema":1,"status":"bounded_integer_residual_degree10_physics_map_certificate",
        "external_upper_bounds":{str(k):v for k,v in UPPER.items()},
        "external_source":"Cederwall et al., J. Phys. A 59 (2026) 065203",
        "registry_fingerprint":reg.fingerprint,"model_path":str(model_path),
        "model_original_status":model.get("coefficient_status"),"bounds":bounds,
        "primes":primes,"modulus_product":M,"electric_points":points,
        "evaluation_injective":injective,"sample_count":len(samples),
        "all_integer_residuals_zero":True,
        "claim":"All declared degree-10 polynomial stress-generator basis identities hold over QQ.",
        "limitations":["external invariant-space dimensions are theorem inputs","declared HLS conventions are assumed"]}
    record["certificate_hash"]=semantic_hash(record)
    atomic_json(output/"map_certificate.json",record)
    atomic_json(output/"summary.json",{"status":record["status"],"certificate_hash":record["certificate_hash"],
        "point_count":len(points),"prime_count":len(primes),"modulus_product":M,
        "worst_integer_residual_bound":B,"evaluation_injective":injective})
    return record
