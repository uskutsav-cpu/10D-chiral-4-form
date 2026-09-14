"""Attempt an explicit degree-12 orbit ideal from the modular equation-space lift."""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path
from .degree12_certificate import certify_records,LOWER_LABELS
from .exact import nullspace_q,rref_q
from .normalization_audit import load_records
from .polynomial import Poly
from .provenance import atomic_json

def candidate_basis_q():
    n=78
    A,B,C,D,U,V=[Poly.variable(n,i) for i in range(6)]
    Z=[Poly.variable(n,6+i) for i in range(72)]
    lower=[A**5,A**3*B,A**2*C,A**2*D,A*B**2,A*U,A*V,B*C,B*D]
    return Z+lower

def build_degree12_ideal(models_dir="runs/degree12-reduced-extended",
                         output="verification/all_orders/degree12_explicit_ideal.json"):
    paths=sorted(Path(models_dir).glob("fields_prime*.json"))
    if len(paths)<2:raise ValueError("need at least two reduced degree-12 models")
    records=load_records(paths)
    result=certify_records(records,holdout_primes=(),max_bad_primes=0)
    lift=result["equation_space_lift"]
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    if not lift.get("qq_verified"):
        blocker={"schema":1,"status":"blocked_degree12_equation_space_not_lifted",
            "model_count":len(paths),"primes":sorted(records),
            "equation_lift_status":lift.get("status"),
            "unresolved_count":len(lift.get("unresolved",[])),
            "next_action":("add independent reduced degree-12 prime models, then rerun. "
                           "Do not replace this with a modular-only ideal claim.")}
        atomic_json(out,blocker);return blocker
    eq=[[Fraction(x) for x in row] for row in lift["rational_rows"]]
    kernel=nullspace_q(eq,81)
    if len(kernel)!=68:raise AssertionError("expected 68 degree-12 obstructions")
    _,piv=rref_q([row[:72] for row in kernel],ncols=72)
    if len(piv)!=68:raise AssertionError("degree-12 leading obstruction rank changed")
    basis=candidate_basis_q()
    constraints=[]
    for row in kernel:
        q=Poly(78)
        for c,f in zip(row,basis):
            q=q+c*f
        constraints.append(q)
    payload={"schema":1,"status":"conditional_exact_degree12_orbit_ideal_from_lifted_equation_space",
        "coordinate_ids":list(next(iter(records.values()))["coordinate_ids"]),
        "candidate_labels":list(next(iter(result["modular"]["reports"]),{})) if False else
            list(next(iter(records.values()))["coordinate_ids"][6:])+list(LOWER_LABELS),
        "constraint_count":68,"leading_rank":68,
        "constraints":[q.to_json() for q in constraints],
        "equation_space_lift_status":lift["status"],
        "conditional_on_exactness_of_lifted_equation_space":True,
        "unconditional_physics_tangency":False,
        "remaining_gate":("certify the lifted equation-space coefficients as the exact "
                          "characteristic-zero physics-derived equation space")}
    atomic_json(out,payload);return payload
