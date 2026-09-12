"""HLS equation (2.33), in mixed components with tau = 48 T.

Uses constrained anti-self-dual derivatives, never derivatives of a naive
unconstrained F^2 action. The trace sign is derived directly from (2.33):
Tr(tau)|V_d = 10(d-2) V_d. The sign printed in (2.36) is not silently used.
"""
from __future__ import annotations
import numpy as np
from .exact import rational
from .forms import mixed_bilinear
from .tensors import DEFAULT_BUDGET,matrix_product
from .weighted import stress_generators


def residue(q,p):
    q=rational(q)
    if q.denominator%p==0:
        raise ValueError('prime divides a normalization denominator')
    return q.numerator*pow(q.denominator,-1,p)%p


def _add(table,key,value,p):
    table[key]=(table.get(key,0)+value)%p


def coefficient_ids(reg,max_degree):
    return tuple(i for d in sorted(reg.degree_bases) if d<=max_degree for i in reg.degree_bases[d])


def tau_terms(reg,form,p,max_degree,*,cache=None,budget=DEFAULT_BUDGET):
    """Exact homogeneous matrix coefficients indexed by (field degree, c-monomial)."""
    cache={} if cache is None else cache
    ids=coefficient_ids(reg,max_degree)
    terms={(2,()):mixed_bilinear(form,form,p,budget=budget)}
    identity=np.eye(10,dtype=np.int64)
    gradients={}
    for name in ids:
        d=reg.items[name].degree
        v=reg.evaluate(name,form,p,cache=cache,budget=budget)
        _add(terms,(d,(name,)),(d-2)*v*identity,p)
        # No derivative of degree d can enter a pair if d+4-2 exceeds cutoff.
        if d+2<=max_degree:
            _,gradients[name]=reg.evaluate(name,form,p,gradient=True,cache=cache,budget=budget)
    names=tuple(gradients)
    for pos,a in enumerate(names):
        for b in names[pos:]:
            d=reg.items[a].degree+reg.items[b].degree-2
            if d>max_degree:
                continue
            B=mixed_bilinear(gradients[a],gradients[b],p,budget=budget)
            if a!=b:
                B=(B+mixed_bilinear(gradients[b],gradients[a],p,budget=budget))%p
            _add(terms,(d,tuple(sorted((a,b)))),-25*B,p)
    return {k:v for k,v in terms.items() if np.any(v)}


def matrix_series_product(left,right,p,max_degree,budget=DEFAULT_BUDGET,max_terms=50000):
    out={}
    for (da,ma),a in left.items():
        for (db,mb),b in right.items():
            if da+db<=max_degree:
                key=(da+db,tuple(sorted(ma+mb)))
                _add(out,key,matrix_product(a,b,p,budget),p)
                if len(out)>max_terms:
                    raise RuntimeError('formal tensor-series term budget exceeded')
    return {k:v for k,v in out.items() if np.any(v)}


def scalar_series_product(left,right,p,max_degree,max_terms=50000):
    out={}
    for (da,ma),a in left.items():
        for (db,mb),b in right.items():
            if da+db<=max_degree:
                key=(da+db,tuple(sorted(ma+mb)))
                out[key]=(out.get(key,0)+int(a)*int(b))%p
                if len(out)>max_terms:
                    raise RuntimeError('formal scalar-series term budget exceeded')
    return {k:v for k,v in out.items() if v}


def trace_series(reg,form,p,max_degree,*,cache=None,budget=DEFAULT_BUDGET):
    cache={} if cache is None else cache
    traces={}
    traces['tr1']={(reg.items[i].degree,(i,)):10*(reg.items[i].degree-2)*reg.evaluate(i,form,p,cache=cache,budget=budget)%p
                   for i in coefficient_ids(reg,max_degree)}
    # In Tr(tau^k), k>=2, any individual tensor factor can have degree at
    # most max_degree-2. Tr(tau) above uses its independent exact Euler identity.
    tau=tau_terms(reg,form,p,max(2,max_degree-2),cache=cache,budget=budget)
    power={(0,()):np.eye(10,dtype=np.int64)}
    for k in range(1,min(10,max_degree//2)+1):
        power=matrix_series_product(power,tau,p,max_degree,budget)
        if k>=2:
            traces[f'tr{k}']={key:int(np.trace(A))%p for key,A in power.items() if int(np.trace(A))%p}
    return traces


def evaluate_generators(reg,form,p,max_degree,extras=(),*,cache=None,budget=DEFAULT_BUDGET):
    """Evaluate all f(tau,S) monomials, not merely changes to the initial seed."""
    cache={} if cache is None else cache
    scalar=trace_series(reg,form,p,max_degree,cache=cache,budget=budget)
    extra_degrees={}
    for name in extras:
        if name not in reg.items:
            raise ValueError(f'unknown extra invariant {name}')
        token='S:'+name;d=reg.items[name].degree
        extra_degrees[token]=d
        scalar[token]={(d,()):reg.evaluate(name,form,p,cache=cache,budget=budget)}
    catalog=stress_generators(max_degree,extra_degrees)
    out={}
    for gen in catalog:
        series={(0,()):1}
        for factor in gen.factors:
            series=scalar_series_product(series,scalar.get(factor,{}),p,max_degree)
        for (d,monomial),value in series.items():
            out[(gen.id,d,monomial)]=int(value)
    return out,catalog


def evaluate_tau_at_coefficients(terms,coefficients,p):
    out=np.zeros((10,10),dtype=np.int64)
    for (_,mono),A in terms.items():
        c=1
        for name in mono:
            c=c*residue(coefficients.get(name,0),p)%p
        out=(out+c*A)%p
    return out
