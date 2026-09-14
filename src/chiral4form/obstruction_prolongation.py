"""Linear graded prolongation of an explicit lower obstruction ideal."""
from __future__ import annotations
import json
from pathlib import Path
from .finite_field import matrix_rank,matvec
from .polynomial import Poly
from .higher_obstruction import coupling_weight,enumerate_weight_monomials,candidate_basis

def _poly_vector(poly,basis):
    lookup={}
    for j,q in enumerate(basis):
        if len(q.terms)!=1:continue
        (m,c),=q.terms.items()
        if c==1:lookup[m]=j
    row=[0]*len(basis)
    for m,c in poly.terms.items():
        if m not in lookup:raise ValueError("polynomial term is outside target homogeneous ansatz")
        row[lookup[m]]=int(c)%poly.p
    return row

def prolongation_vectors(ideal_record,reg,degree,p,ids):
    lower_ids=tuple(ideal_record.get("coordinate_ids",()))
    if not lower_ids:
        raise ValueError("ideal record must declare coordinate_ids")
    if tuple(ids[:len(lower_ids)]) != lower_ids:
        raise ValueError(
            "lower ideal coordinates are not a prefix of the higher-jet coordinate ring; "
            "an explicit embedding/restriction map is required"
        )
    constraints=[Poly.from_json(x) for x in ideal_record["constraints"]]
    n=len(ids);weights=[coupling_weight(reg.items[name].degree) for name in ids]
    basis,labels,_=candidate_basis(ids,reg,degree,p)
    rows=[]
    for q0 in constraints:
        # Embed lower ring into the current ring.
        q=Poly(n,{tuple(m)+(0,)*(n-q0.n):c for m,c in q0.terms.items()},p)
        term_weights={sum(e*w for e,w in zip(m,weights)) for m in q.terms}
        if len(term_weights)!=1:raise ValueError("lower ideal generator is not homogeneous")
        w0=term_weights.pop();remaining=coupling_weight(degree)-w0
        if remaining<0:continue
        allowed=[i for i,name in enumerate(ids) if reg.items[name].degree<degree]
        for e in enumerate_weight_monomials(weights,remaining,allowed):
            multiplier=Poly(n,{e:1},p)
            rows.append(_poly_vector(q*multiplier,basis))
    return rows,labels

def quotient_report(kernel,prolongations,p):
    if not kernel:raise ValueError("nonempty obstruction kernel required")
    n=len(kernel[0])
    rk=matrix_rank(kernel,p,ncols=n);rp=matrix_rank(prolongations,p,ncols=n)
    ru=matrix_rank([*prolongations,*kernel],p,ncols=n)
    # If every prolongation is a genuine obstruction, span(P) must lie in span(K).
    inclusion=(ru==rk)
    newdim=None if not inclusion else rk-rp
    return {"obstruction_rank":rk,"prolongation_rank":rp,"union_rank":ru,
            "prolongations_inside_obstruction_space":inclusion,
            "new_obstruction_quotient_dimension":newdim}
