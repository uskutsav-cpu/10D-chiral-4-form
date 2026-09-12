"""Distinct, explicitly scoped notions of flow reachability.

A coefficient hull is a LINEAR over-approximation, not a nonlinear orbit.
Finite-depth Lie brackets give a lower bound on a distribution rank, not a
complete classification. Jet coefficients are realizable for a fixed
polynomial vector field, but their linear span is not an orbit dimension.
"""
from __future__ import annotations
from itertools import combinations
from .polynomial import Poly,bracket,derivation,field_at
from .finite_field import rref as rref_p
from .exact import rref_q


def row_basis(rows,n,p):
    rr,piv=rref_p(rows,p,ncols=n) if p is not None else rref_q(rows,n)
    return rr[:len(piv)]


def linear_invariant_hull(fields,seed=None,max_iterations=100):
    if type(max_iterations) is not int or max_iterations < 1:
        raise ValueError('positive iteration budget required')
    if not fields:
        raise ValueError('at least one vector field is required')
    n=len(fields[0]);p=fields[0][0].p
    W=[] if seed is None else row_basis([seed],n,p)
    history=[]
    for step in range(max_iterations):
        k=len(W)
        z=[Poly.variable(k,i,p) for i in range(k)]
        substitutions=[sum((W[j][i]*z[j] for j in range(k)),Poly(k,p=p)) for i in range(n)]
        rows=list(W)
        for field in fields:
            restricted=[f.substitute(substitutions,n_out=k) for f in field]
            monomials=sorted(set().union(*(f.terms for f in restricted)))
            rows += [[f.terms.get(m,0) for f in restricted] for m in monomials]
        new=row_basis(rows,n,p)
        history.append({'iteration':step,'rank':len(new)})
        if new==W:
            return {'interpretation':'minimal_linear_flow_invariant_hull_not_orbit',
                    'basis':[[str(x) for x in row] for row in W],'rank':len(W),
                    'history':history,'stabilized':True,'domain':'QQ' if p is None else f'GF({p})'}
        W=new
    return {'interpretation':'partial_linear_hull_not_yet_invariant','rank':len(W),
            'history':history,'stabilized':False,'basis':[[str(x) for x in row] for row in W]}


def lie_rank(fields,point=None,max_depth=3,max_fields=10000):
    if type(max_depth) is not int or max_depth < 1 or max_fields < 1:
        raise ValueError('positive Lie depth and field budget required')
    if not fields:
        raise ValueError('at least one vector field is required')
    n=len(fields[0]);p=fields[0][0].p
    point=[0]*n if point is None else point
    frontier=[(str(i),tuple(f)) for i,f in enumerate(fields)]
    rows=[];witnesses=[];history=[];count=0
    seen=set()
    for depth in range(1,max_depth+1):
        next_frontier=[]
        for label,f in frontier:
            key=tuple(tuple(sorted(x.terms.items())) for x in f)
            if key in seen:
                continue
            seen.add(key);count+=1
            if count>max_fields:
                return {'interpretation':'finite_depth_lie_rank_lower_bound','rank':len(row_basis(rows,n,p)),
                    'completed_depth':depth-1,'terminated':'resource_limit','history':history,'witnesses':witnesses}
            value=field_at(f,point)
            previous=len(row_basis(rows,n,p));candidate=rows+[value]
            if len(row_basis(candidate,n,p))>previous:
                rows=candidate;witnesses.append({'word':label,'value':[str(x) for x in value]})
            # CRITICAL: do not discard a field just because its value is zero.
            if depth<max_depth and any(f):
                for i,g in enumerate(fields):
                    b=bracket(g,f)
                    if any(b):
                        next_frontier.append((f'[{i},{label}]',b))
        history.append({'depth':depth,'rank':len(row_basis(rows,n,p)),'fields_visited':count})
        frontier=next_frontier
        if not frontier:
            break
    return {'interpretation':'finite_depth_lie_rank_lower_bound','rank':len(row_basis(rows,n,p)),
            'completed_depth':depth,'terminated':'depth_or_bracket_exhaustion',
            'history':history,'witnesses':witnesses,'not_a_global_orbit_classification':True}


def trajectory_jets(field,order,point=None):
    """Formal autonomous solution coefficients c(t), via repeated derivations.

    Exactly solves the supplied finite polynomial ODE to the requested time
    order. No inference about an untruncated field theory is made.
    """
    from math import factorial
    n=len(field);p=field[0].p
    if order<0 or (p is not None and order>=p):
        raise ValueError('invalid jet order, or factorial is not invertible')
    point=[0]*n if point is None else point
    derivatives=[Poly.variable(n,i,p) for i in range(n)]
    out=[]
    for j in range(order+1):
        f=factorial(j)
        row=field_at(derivatives,point)
        coeff=[(x*pow(f,-1,p)%p if p else x/f) for x in row]
        out.append([str(x) for x in coeff])
        derivatives=[derivation(field,x) for x in derivatives]
    return {'interpretation':'autonomous_formal_time_jet','coefficients':out,'order':order}


def completion_search(base_fields,candidate_fields,*,target_rank,metric='linear_hull',max_subsets=10000,max_depth=3):
    """Exhaustive cardinality search ONLY inside the supplied finite catalogue.

    For physics f(T,S), candidate_fields must contain all mixed generator
    monomials for each selected S set; use the pipeline's per-subset builder
    for that case. This helper handles a finite catalogue of vector fields.
    """
    if metric not in ('linear_hull','lie_rank'):
        raise ValueError('unknown completion metric')
    names=sorted(candidate_fields);tested=0;records=[]
    for size in range(len(names)+1):
        for subset in combinations(names,size):
            if tested>=max_subsets:
                return {'status':'budget_exhausted','tested':tested,'global_minimality':False,'records':records}
            fields=list(base_fields)+[candidate_fields[i] for i in subset]
            report=linear_invariant_hull(fields) if metric=='linear_hull' else lie_rank(fields,max_depth=max_depth)
            tested+=1;records.append({'subset':list(subset),'rank':report['rank']})
            if report['rank']>=target_rank and (metric!='linear_hull' or report['stabilized']):
                return {'status':'minimum_in_supplied_finite_catalogue','metric':metric,'subset':list(subset),
                        'cardinality':size,'tested':tested,'global_minimality':False,'records':records}
    return {'status':'not_completed_in_catalogue','tested':tested,'global_minimality':False,'records':records}
