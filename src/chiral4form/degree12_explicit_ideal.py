"""Build the explicit triangular degree-12 obstruction ideal from the lifted kernel."""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path
from .polynomial import Poly
from .provenance import atomic_json

LOWER_LABELS=("A^5","A^3*B","A^2*C","A^2*D","A*B^2","A*U","A*V","B*C","B*D")

def ansatz_polynomials():
    n=78
    A,B,C,D,U,V=[Poly.variable(n,i) for i in range(6)]
    Z=[Poly.variable(n,6+i) for i in range(72)]
    lower=[A**5,A**3*B,A**2*C,A**2*D,A*B**2,A*U,A*V,B*C,B*D]
    labels=[f"Z{i+1}" for i in range(72)]+list(LOWER_LABELS)
    return Z+lower,labels

def build_ideal(kernel_path="verification/all_orders/degree12_exact/kernel_lift.json",
                output="verification/all_orders/degree12_exact/ideal.json"):
    rec=json.loads(Path(kernel_path).read_text())
    lift=rec["lift"]
    if not lift.get("qq_verified"):
        raise ValueError("degree-12 kernel is not an exact verified QQ lift")
    rows=[[Fraction(x) for x in row] for row in lift["rational_rows"]]
    if len(rows)!=68 or any(len(row)!=81 for row in rows):
        raise ValueError("lifted kernel matrix is not 68x81")
    basis,labels=ansatz_polynomials()
    constraints=[]
    for row in rows:
        q=Poly(78)
        for c,b in zip(row,basis):
            if c:q=q+c*b
        constraints.append(q)

    pivots=list(rec["pivots"])
    free=list(rec["free_degree12_columns"])
    if len(pivots)!=68 or len(free)!=4 or any(p>=72 for p in pivots):
        raise ValueError("invalid triangular kernel metadata")

    substitutions=[Poly.variable(78,i) for i in range(78)]
    triangular=[]
    for r,pivot in enumerate(pivots):
        row=rows[r]
        if row[pivot]!=1:
            raise ValueError("kernel row is not normalized at pivot")
        if any(row[p] for p in pivots if p!=pivot):
            raise ValueError("kernel is not reduced in pivot columns")
        rhs=Poly(78)
        for j,c in enumerate(row):
            if j!=pivot and c:
                rhs=rhs-c*basis[j]
        substitutions[6+pivot]=rhs
        triangular.append({
            "pivot_column":pivot,
            "pivot_label":labels[pivot],
            "rhs":rhs.to_json(),
        })

    if any(q.substitute(substitutions).terms for q in constraints):
        raise AssertionError("triangular substitutions do not annihilate all constraints")

    payload={
        "schema":1,
        "status":"exact_degree12_triangular_ideal",
        "constraint_count":68,
        "leading_rank":68,
        "degree12_coordinate_count":72,
        "fiber_dimension":4,
        "pivots":pivots,
        "free_degree12_columns":free,
        "ansatz_labels":labels,
        "constraints":[q.to_json() for q in constraints],
        "triangular_equations":triangular,
        "all_generators_reduce_to_zero":True,
        "scope":"degree-12 obstruction ideal on the certified degree-10 orbit quotient",
    }
    path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(path,payload)
    return payload
