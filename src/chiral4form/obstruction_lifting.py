"""Lift exact finite-jet invariant ideals to formal all-orders cylinders."""
from __future__ import annotations
from dataclasses import dataclass
from .formal_filtration import prove_truncation_commutes
from .ideals import invariant_ideal_report
from .polynomial import Poly,derivation

def extend_poly(poly:Poly,new_n:int)->Poly:
    if new_n<poly.n: raise ValueError("cannot shrink polynomial ring")
    if new_n==poly.n: return poly
    extra=(0,)*(new_n-poly.n)
    return Poly(new_n,{tuple(m)+extra:c for m,c in poly.terms.items()},poly.p)

@dataclass(frozen=True)
class LiftedIdealCertificate:
    cutoff:int
    label:str
    finite_tangent:bool
    filtration_verified:bool
    all_orders_cylinder_invariant:bool
    finite_report:dict
    def to_json(self)->dict:
        return {"cutoff":self.cutoff,"label":self.label,
                "finite_tangent":self.finite_tangent,
                "filtration_verified":self.filtration_verified,
                "all_orders_cylinder_invariant":self.all_orders_cylinder_invariant,
                "finite_report":self.finite_report,
                "lifting_argument":("pi_<=N X=X_N pi_<=N and X_N(J_N) subset J_N imply "
                                    "X(pi^*J_N) subset pi^*J_N. Time-dependent signed "
                                    "controls preserve tangency by linearity.")}

def lift_finite_jet_ideal(fields,constraints,*,cutoff:int,label:str)->LiftedIdealCertificate:
    if not fields: raise ValueError("nonempty fields required")
    if not constraints: raise ValueError("nonempty ideal required")
    first=next(iter(fields.values()))[0]
    if first.p is not None: raise ValueError("exact QQ fields required")
    finite=invariant_ideal_report(fields,constraints)
    tangent=bool(finite["all_tangent"])
    filtration=prove_truncation_commutes(cutoff)
    cert=LiftedIdealCertificate(cutoff,label,tangent,filtration.all_checks,
        tangent and filtration.all_checks,finite)
    if not cert.all_orders_cylinder_invariant:
        raise AssertionError(f"cannot lift non-invariant ideal {label}")
    return cert

def relative_transport_report(fields,obstruction:Poly,multipliers)->dict:
    rows={}; ok_all=True
    for name,field in fields.items():
        if name not in multipliers: raise ValueError(f"missing multiplier {name}")
        h=multipliers[name]
        if not isinstance(h,Poly):
            h=Poly.constant(obstruction.n,h,obstruction.p)
        lhs=derivation(field,obstruction); rhs=h*obstruction
        ok=lhs==rhs; ok_all &= ok
        rows[name]={"verified":ok}
    return {"status":"exact_relative_transport","all_verified":ok_all,
            "per_generator":rows,"zero_set_preserved":ok_all}
