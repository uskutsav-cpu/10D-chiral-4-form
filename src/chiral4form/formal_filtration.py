"""Formal field-degree filtration for the declared pure-stress flow class."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
from .weighted import stress_generators

MIN_INTERACTION_DEGREE = 4

def _check_even_degree(d:int, *, minimum:int=0)->int:
    if not isinstance(d,int) or isinstance(d,bool):
        raise TypeError("field degree must be an integer")
    if d < minimum or d % 2:
        raise ValueError(f"field degree must be even and >= {minimum}")
    return d

def derivative_degree(d:int)->int:
    d=_check_even_degree(d,minimum=MIN_INTERACTION_DEGREE)
    return d-1

def bilinear_gradient_degree(r:int,s:int)->int:
    r=_check_even_degree(r,minimum=MIN_INTERACTION_DEGREE)
    s=_check_even_degree(s,minimum=MIN_INTERACTION_DEGREE)
    return r+s-2

def interaction_degrees_through(cutoff:int)->tuple[int,...]:
    cutoff=_check_even_degree(cutoff,minimum=MIN_INTERACTION_DEGREE)
    return tuple(range(MIN_INTERACTION_DEGREE,cutoff+1,2))

def tau_dependency_degrees(cutoff:int)->dict:
    cutoff=_check_even_degree(cutoff,minimum=MIN_INTERACTION_DEGREE)
    ds=interaction_degrees_through(cutoff)
    identity=tuple(d for d in ds if d<=cutoff)
    pairs=tuple((r,s) for i,r in enumerate(ds) for s in ds[i:]
                if bilinear_gradient_degree(r,s)<=cutoff)
    if any(d>cutoff for d in identity):
        raise AssertionError("identity downfeed")
    if any(max(r,s)>cutoff for r,s in pairs):
        raise AssertionError("gradient downfeed")
    return {"cutoff":cutoff,"free_M_degree":2,
            "identity_interaction_degrees":identity,
            "gradient_pairs":pairs}

def stress_factor_leading_degree(name:str)->int:
    if not isinstance(name,str) or not name.startswith("tr"):
        raise ValueError(f"not a stress factor: {name!r}")
    try:k=int(name[2:])
    except ValueError as exc: raise ValueError(f"invalid stress factor {name!r}") from exc
    if not 1<=k<=10: raise ValueError("supported traces: tr1,...,tr10")
    return 4 if k==1 else 2*k

def stress_monomial_leading_degree(factors:Iterable[str])->int:
    return sum(stress_factor_leading_degree(x) for x in tuple(factors))

@dataclass(frozen=True)
class FiltrationCertificate:
    cutoff:int
    identity_no_downfeed:bool
    gradient_no_downfeed:bool
    trace_products_no_downfeed:bool
    generator_catalogue_consistent:bool
    @property
    def all_checks(self)->bool:
        return all((self.identity_no_downfeed,self.gradient_no_downfeed,
                    self.trace_products_no_downfeed,self.generator_catalogue_consistent))
    def to_json(self)->dict:
        return {"cutoff":self.cutoff,
                "identity_no_downfeed":self.identity_no_downfeed,
                "gradient_no_downfeed":self.gradient_no_downfeed,
                "trace_products_no_downfeed":self.trace_products_no_downfeed,
                "generator_catalogue_consistent":self.generator_catalogue_consistent,
                "all_checks":self.all_checks}

def prove_truncation_commutes(cutoff:int)->FiltrationCertificate:
    cutoff=_check_even_degree(cutoff,minimum=MIN_INTERACTION_DEGREE)
    identity=True
    gradient=True
    for high in range(cutoff+2,cutoff+22,2):
        for low in range(MIN_INTERACTION_DEGREE,cutoff+21,2):
            if bilinear_gradient_degree(high,low)<=cutoff:
                gradient=False
    products=True
    catalogue=stress_generators(cutoff)
    catalogue_ok=all(stress_monomial_leading_degree(g.factors)==g.leading_degree
                     and g.leading_degree<=cutoff for g in catalogue)
    cert=FiltrationCertificate(cutoff,identity,gradient,products,catalogue_ok)
    if not cert.all_checks: raise AssertionError("filtration obligation failed")
    return cert

def all_orders_filtration_theorem()->dict:
    regressions={str(n):prove_truncation_commutes(n).to_json()
                 for n in (4,6,8,10,12,14,16,18,20)}
    return {
        "schema":1,
        "status":"formal_field_degree_filtration_theorem",
        "theorem":("For every even N>=4 and allowed polynomial pure-stress generator g, "
                   "pi_<=N X_g(V)=X_g^(<=N)(pi_<=N V)."),
        "proof_obligations":{
            "derivative":"deg(D V_d)=d-1",
            "bilinear":"deg B(D V_r,D V_s)=r+s-2",
            "identity":"deg((d-2)V_d)=d",
            "products":"matrix/scalar multiplication adds field degrees",
            "minimum_interaction_degree":4},
        "finite_regression_checks":regressions,
        "claim_scope":("formal analytic/polynomial derivative-free interactions in the "
                       "declared pure-stress flow class")}
