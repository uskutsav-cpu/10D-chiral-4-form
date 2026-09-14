"""Rank-select higher-degree invariant bases from streamed contraction graphs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from .finite_field import _check_prime_field_modulus
from .forms import random_selfdual
from .graphs import allocation_profile,evaluate_graph
from .nauty_multigraphs import stream_multigraphs,graph_label
from .provenance import atomic_json
from .registry import Registry,Invariant
from .tensors import TensorBudget

PUBLISHED_DIMENSIONS={14:247,16:1364}
DEFAULT_PRIME=65521

class IncrementalRowBasis:
    def __init__(self,p,ncols):
        self.p=p;self.ncols=ncols;self.rows=[];self.pivots=[]
    @property
    def rank(self):return len(self.rows)
    def add(self,vector):
        p=self.p;v=[int(x)%p for x in vector]
        for pivot,row in zip(self.pivots,self.rows):
            if v[pivot]:
                f=v[pivot]
                v=[(x-f*y)%p for x,y in zip(v,row)]
        pivot=next((i for i,x in enumerate(v) if x),None)
        if pivot is None:return False
        inv=pow(v[pivot],-1,p);v=[x*inv%p for x in v]
        pos=0
        while pos<len(self.pivots) and self.pivots[pos]<pivot:pos+=1
        self.pivots.insert(pos,pivot);self.rows.insert(pos,v)
        # Maintain reduced pivot columns.
        for i,row in enumerate(self.rows):
            if i!=pos and row[pivot]:
                f=row[pivot];self.rows[i]=[(x-f*y)%p for x,y in zip(row,v)]
        return True

def merge_registry(base,degree,graphs):
    items=list(base.items.values())
    ids=[]
    for i,M in enumerate(graphs,1):
        name=f"I{degree}_{i}";ids.append(name);items.append(Invariant(name,degree,M,()))
    bases=dict(base.degree_bases);bases[degree]=ids
    meta=dict(base.metadata)
    meta[f"degree{degree}_basis"]={
        "status":"evaluation-independent contraction basis conditional on published dimension",
        "published_dimension":len(ids),
        "enumerator":"nauty genbg incidence-graph encoding",
    }
    return Registry(items,bases,meta)

def build_higher_basis(base_registry_path,degree,output,prime=DEFAULT_PRIME,
                       max_candidates=200000,max_bytes=8*1024**3,
                       max_multiply_adds=200_000_000_000,seed_start=202609150001,
                       progress=print):
    if degree not in PUBLISHED_DIMENSIONS:raise ValueError("supported new degrees are 14 and 16")
    _check_prime_field_modulus(prime)
    target=PUBLISHED_DIMENSIONS[degree]
    base=Registry.load(base_registry_path)
    if degree in base.degree_bases and len(base.degree_bases[degree])==target:
        return base
    sample_count=target
    forms=[random_selfdual(seed_start+i,prime) for i in range(sample_count)]
    budget=TensorBudget(max_bytes,max_multiply_adds,True)
    basis=IncrementalRowBasis(prime,sample_count);selected=[];scanned=0;over_budget=0
    checkpoint=Path(output);checkpoint.mkdir(parents=True,exist_ok=True)
    for M in stream_multigraphs(degree):
        scanned+=1
        if scanned>max_candidates:break
        profile=allocation_profile(M,10,gradient=False)
        if profile["reserved_bytes"]>max_bytes or profile["multiply_adds"]>max_multiply_adds:
            over_budget+=1;continue
        values=[evaluate_graph(M,F,prime,gradient=False,budget=budget) for F in forms]
        if basis.add(values):
            selected.append(M)
            progress(f"degree {degree}: rank {basis.rank}/{target}, scanned={scanned}, over_budget={over_budget}",flush=True)
            atomic_json(checkpoint/"progress.json",{
                "degree":degree,"prime":prime,"target":target,"rank":basis.rank,
                "scanned":scanned,"over_budget":over_budget,
                "selected_labels":[graph_label(x) for x in selected]})
            if basis.rank==target:break
    if basis.rank!=target:
        atomic_json(checkpoint/"failure.json",{
            "status":"basis_not_completed","degree":degree,"target":target,"rank":basis.rank,
            "scanned":scanned,"over_budget":over_budget,
            "suggestion":"increase max_candidates/resources or use a faster/spinor invariant evaluator"})
        raise RuntimeError(f"degree-{degree} basis incomplete: rank {basis.rank}/{target}")
    reg=merge_registry(base,degree,selected)
    atomic_json(checkpoint/f"registry_degree{degree}.json",reg.to_json())
    atomic_json(checkpoint/"summary.json",{
        "status":"higher_degree_basis_completed","degree":degree,"dimension":target,
        "prime":prime,"evaluation_rank":basis.rank,"scanned":scanned,"over_budget":over_budget,
        "published_dimension_input":target,"registry_fingerprint":reg.fingerprint})
    return reg
