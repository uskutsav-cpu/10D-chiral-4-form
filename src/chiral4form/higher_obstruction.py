"""Homogeneous higher-jet obstruction spaces and prolongation quotients."""
from __future__ import annotations
import json
from pathlib import Path
from .finite_field import matrix_rank,nullspace
from .fitting import fields_from_json
from .polynomial import Poly,derivation
from .registry import Registry

def coupling_weight(field_degree): return 10*(field_degree-2)

def enumerate_weight_monomials(weights,target,allowed_indices):
    n=len(weights);out=[]
    def rec(start,remain,ex):
        if remain==0:
            out.append(tuple(ex));return
        for pos in range(start,len(allowed_indices)):
            i=allowed_indices[pos];w=weights[i]
            if w>remain:continue
            ex[i]+=1;rec(pos,remain-w,ex);ex[i]-=1
    rec(0,target,[0]*n)
    return out

def candidate_basis(ids,reg,degree,p):
    n=len(ids);idx={name:i for i,name in enumerate(ids)}
    weights=[coupling_weight(reg.items[name].degree) for name in ids]
    target=coupling_weight(degree)
    target_ids=list(reg.degree_bases[degree])
    Z=[Poly.variable(n,idx[name],p) for name in target_ids]
    lower_indices=[i for i,name in enumerate(ids) if reg.items[name].degree<degree]
    exps=enumerate_weight_monomials(weights,target,lower_indices)
    lower=[Poly(n,{e:1},p) for e in exps]
    labels=target_ids+[f"lower:{','.join(f'{ids[i]}^{e}' for i,e in enumerate(x) if e)}" for x in exps]
    return Z+lower,labels,len(Z)

def obstruction_equations(record,reg,degree,p):
    ids,fields=fields_from_json(record)
    if tuple(ids)!=tuple(i for d in sorted(reg.degree_bases) if d<=degree for i in reg.degree_bases[d]):
        raise ValueError("field model IDs disagree with registry")
    basis,labels,leading_count=candidate_basis(ids,reg,degree,p)
    target=coupling_weight(degree);rows=[]
    for gname in sorted(fields):
        field=fields[gname]
        derivs=[derivation(field,q)-(target*q if gname=="tr1" else Poly(q.n,p=p)) for q in basis]
        mons=sorted(set().union(*(f.terms.keys() for f in derivs)))
        for m in mons:
            row=[int(f.terms.get(m,0))%p for f in derivs]
            if any(row):rows.append(row)
    return ids,labels,leading_count,rows

def modular_obstruction_report(record,reg,degree,p):
    ids,labels,leading_count,E=obstruction_equations(record,reg,degree,p)
    K=nullspace(E,p,ncols=len(labels))
    return {"prime":p,"degree":degree,"ansatz_dimension":len(labels),
        "equation_rank":matrix_rank(E,p,ncols=len(labels)),
        "obstruction_nullity":len(K),
        "leading_obstruction_rank":matrix_rank([r[:leading_count] for r in K],p,ncols=leading_count),
        "labels":labels,"kernel":K}
