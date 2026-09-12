"""Exact polynomial ideal tangency in the declared finite coefficient model."""
from __future__ import annotations
import sympy as sp
from .polynomial import Poly,derivation


def sympy_expr(f,variables):
    out=sp.S.Zero
    for m,c in f.terms.items():
        c=sp.Rational(str(c))
        term=c
        for x,e in zip(variables,m):
            term*=x**e
        out+=term
    return out


def invariant_ideal_report(fields,constraints):
    if not fields:
        raise ValueError('nonempty vector-field catalogue required')
    first=next(iter(fields.values()))[0];n,p=first.n,first.p
    variables=sp.symbols(f'c0:{n}')
    if any(g.n!=n or g.p!=p for g in constraints):
        raise ValueError('constraint ring mismatch')
    if not constraints:
        return {'status':'whole_space','all_tangent':True,'remainders':{}}
    exprs=[sympy_expr(g,variables) for g in constraints]
    options={'domain':sp.QQ} if p is None else {'modulus':p}
    G=sp.groebner(exprs,*variables,order='grevlex',**options)
    rows={}
    for name,field in fields.items():
        rows[name]=[str(G.reduce(sympy_expr(derivation(field,g),variables))[1]) for g in constraints]
    tangent=all(x=='0' for row in rows.values() for x in row)
    return {'status':'exact_finite_model_ideal_tangency','domain':'QQ' if p is None else f'GF({p})',
        'all_tangent':tangent,'groebner_basis':[str(x.as_expr()) for x in G.polys],
        'remainders':rows,'not_an_all_orders_field_theory_theorem':True,
        'logical_scope':'X(I) subset I is sufficient for preservation under a locally unique smooth real ODE when coefficients are over QQ'}
