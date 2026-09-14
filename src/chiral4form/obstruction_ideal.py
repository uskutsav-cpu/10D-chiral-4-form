"""Ideal membership and the degree-14 stabilization falsification scaffold."""
from __future__ import annotations
import sympy as sp
from .ideals import sympy_expr
from .polynomial import Poly

def ideal_membership_report(generators,candidates)->dict:
    if not generators: raise ValueError("nonempty generators required")
    if not candidates:
        return {"status":"no_candidates","new_generator_count":0,"remainders":[]}
    first=generators[0]; n,p=first.n,first.p
    if any(x.n!=n or x.p!=p for x in [*generators,*candidates]):
        raise ValueError("ring mismatch")
    vars=sp.symbols(f"c0:{n}")
    opts={"domain":sp.QQ} if p is None else {"modulus":p}
    G=sp.groebner([sympy_expr(x,vars) for x in generators],*vars,order="grevlex",**opts)
    rows=[]; new=0
    for f in candidates:
        rem=sp.expand(G.reduce(sympy_expr(f,vars))[1])
        is_new=rem!=0; new+=int(is_new)
        rows.append({"candidate":str(sympy_expr(f,vars)),
                     "remainder":str(rem),
                     "is_new_mod_lower_ideal":is_new})
    return {"status":"exact_ideal_membership_report",
            "domain":"QQ" if p is None else f"GF({p})",
            "groebner_basis":[str(x.as_expr()) for x in G.polys],
            "new_generator_count":new,"remainders":rows}

def degree14_probe_readiness(registry)->dict:
    has_basis=14 in registry.degree_bases and bool(registry.degree_bases[14])
    has_items=any(item.degree==14 for item in registry.items.values())
    if not has_basis:
        return {"schema":1,"status":"blocked_missing_degree14_invariant_basis",
                "degree14_registry_items_present":has_items,
                "degree14_basis_present":False,
                "required_next_input":("explicit degree-14 invariant basis or justified "
                                       "evaluation representation plus degree-14 stress map"),
                "scientific_question":("Do degree-14 obstruction candidates reduce to zero "
                                       "modulo prolongation of the lifted lower ideal?"),
                "do_not_claim":["degree-14 stabilization",
                                "finite generation of the all-orders obstruction ideal"]}
    return {"schema":1,"status":"degree14_probe_ready",
            "degree14_registry_items_present":has_items,
            "degree14_basis_present":True,
            "degree14_basis_dimension":len(registry.degree_bases[14])}
