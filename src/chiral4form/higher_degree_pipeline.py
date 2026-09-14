"""Resumable pure-stress tensor->fit pipeline for custom degree-14/16 registries."""
from __future__ import annotations
import json,platform
from pathlib import Path
import numpy as np,sympy
from .checkpoints import CheckpointStore
from .finite_field import _check_prime_field_modulus
from .fitting import fit_degree,vector_fields_from_fits,fields_to_json
from .forms import random_selfdual,compact,hodge,inner
from .graphs import allocation_profile,plan
from .provenance import atomic_json,semantic_hash
from .registry import Registry
from .stress import evaluate_generators
from .tensors import TensorBudget

def _prime(n):
    if n<2:return False
    if n%2==0:return n==2
    d=3
    while d*d<=n:
        if n%d==0:return False
        d+=2
    return True

def default_primes(count=3):
    out=[]
    for n in range(65521,50000,-2):
        if _prime(n) and n not in (2,3,5):
            out.append(n)
            if len(out)==count:return out
    raise RuntimeError("not enough finite-field primes")

def preflight(reg,degree,budget):
    rows=[]
    for item in reg.items.values():
        if item.degree<=degree and item.graph is not None:
            profile=allocation_profile(item.graph,10,gradient=item.degree<=degree-4)
            rows.append({"invariant":item.id,"degree":item.degree,**profile,
                "within_budget":profile["reserved_bytes"]<=budget.max_bytes and
                    profile["multiply_adds"]<=budget.max_multiply_adds})
    return {"degree":degree,"graphs":rows,"all_graphs_within_budget":all(x["within_budget"] for x in rows)}

def compute_sample(reg,p,seed,degree,budget):
    form=random_selfdual(seed,p);cache={}
    if not np.array_equal(hodge(form,p),form) or inner(form,form,p)!=0:
        raise AssertionError("self-duality sample gate failed")
    targets,catalog=evaluate_generators(reg,form,p,degree,(),cache=cache,budget=budget)
    bases={str(d):reg.basis_values(d,form,p,cache=cache,budget=budget)
           for d in sorted(reg.degree_bases) if d<=degree}
    return {"seed":seed,"prime":p,"compact_selfdual_input":compact(form,p).tolist(),
        "basis_values":bases,
        "targets":[{"generator":g,"degree":d,"monomial":list(m),"value":v}
                   for (g,d,m),v in sorted(targets.items())],
        "generator_catalogue":[{"id":g.id,"factors":list(g.factors),
                               "leading_degree":g.leading_degree} for g in catalog]}

def run_higher(registry_path,degree,output,primes=None,seed_start=202609160001,
               fit_margin=2,holdouts=2,max_bytes=8*1024**3,
               max_multiply_adds=200_000_000_000,progress=print):
    if degree not in (14,16):raise ValueError("supported higher cutoffs are 14 and 16")
    reg=Registry.load(registry_path)
    needed=tuple(range(4,degree+1,2))
    if any(d not in reg.degree_bases for d in needed):
        raise ValueError("registry lacks a homogeneous basis through requested degree")
    primes=default_primes(3) if primes is None else list(primes)
    for p in primes:_check_prime_field_modulus(p)
    budget=TensorBudget(max_bytes,max_multiply_adds,True)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    planrec=preflight(reg,degree,budget);atomic_json(output/"execution_plan.json",planrec)
    if not planrec["all_graphs_within_budget"]:
        bad=[x["invariant"] for x in planrec["graphs"] if not x["within_budget"]]
        atomic_json(output/"status.json",{"state":"preflight_refused","over_budget_graphs":bad})
        raise MemoryError("higher-degree graph preflight refused; see execution_plan.json")
    largest=max(len(reg.degree_bases[d]) for d in needed)
    fit_count=largest+fit_margin;count=fit_count+holdouts
    namespace={"schema":1,"registry":reg.fingerprint,"degree":degree,"primes":primes,
               "seed_start":seed_start,"fit_count":fit_count,"count":count}
    store=CheckpointStore(output/"checkpoints",namespace)
    atomic_json(output/"registry.json",reg.to_json())
    atomic_json(output/"configuration.json",namespace)
    status={"state":"running","degree":degree,"primes_finished":[]};atomic_json(output/"status.json",status)
    reports=[]
    for p in primes:
        samples=[]
        for i in range(count):
            seed=seed_start+i
            s,cached=store.run("tensor_sample",{"prime":p,"seed":seed},
                lambda p=p,seed=seed:compute_sample(reg,p,seed,degree,budget))
            samples.append(s);progress(f"degree{degree} prime {p}: sample {i+1}/{count} {'cached' if cached else 'computed'}",flush=True)
        fits=[]
        for d in needed:
            f,_=store.run("basis_fit",{"prime":p,"degree":d,"sample_hash":semantic_hash(samples)},
                lambda d=d:fit_degree(samples,reg.degree_bases[d],d,p,fit_count))
            fits.append(f);atomic_json(output/f"fit_degree{d}_prime{p}.json",f)
        ids,fields=vector_fields_from_fits(fits)
        model=fields_to_json(ids,fields,"finite_field_samples_and_disjoint_holdouts; not symbolic QQ identity")
        model["generator_catalogue"]=samples[0]["generator_catalogue"]
        atomic_json(output/f"fields_prime{p}.json",model)
        reports.append({"prime":p,"fit_count":fit_count,"holdouts":holdouts,
                        "basis_dimensions":{str(d):len(reg.degree_bases[d]) for d in needed}})
        status["primes_finished"].append(p);atomic_json(output/"status.json",status)
    summary={"schema":1,"state":"completed_declared_computations","degree":degree,
        "prime_reports":reports,"registry_fingerprint":reg.fingerprint,
        "warning":"finite-field fitted higher-degree map; characteristic-zero proof is separate"}
    atomic_json(output/"summary.json",summary);status["state"]="completed_declared_computations";atomic_json(output/"status.json",status)
    return summary
