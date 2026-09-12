"""Exact covariance of coefficient-space vector fields.

Given old couplings c=P z, the transformed field is P^{-1} X(P z).
This P is a COUPLING transformation. If invariant columns change as J=B I,
then c=B^T a. Do not confuse an invariant basis map with a coupling map.
"""
from fractions import Fraction
from .polynomial import Poly
from .exact import rref_q
from .finite_field import solve_many


def inverse(matrix,p=None):
    n=len(matrix)
    if n<1 or any(len(row)!=n for row in matrix):
        raise ValueError('nonempty square transformation required')
    identity=[[int(i==j) for j in range(n)] for i in range(n)]
    if p is not None:
        return solve_many(matrix,identity,p)
    augmented=[list(row)+eye for row,eye in zip(matrix,identity)]
    rr,piv=rref_q(augmented)
    if piv[:n]!=list(range(n)):
        raise ValueError('singular transformation')
    return [row[n:] for row in rr]


def transform_fields(fields,new_to_old,field_degrees=None):
    if not fields:
        raise ValueError('nonempty field catalogue required')
    first=next(iter(fields.values()));n=len(first);p=first[0].p
    if len(new_to_old)!=n:
        raise ValueError('transformation dimension mismatch')
    inv=inverse(new_to_old,p)
    if field_degrees is not None:
        if len(field_degrees)!=n:
            raise ValueError('field-degree list has wrong length')
        if any(new_to_old[i][j] and field_degrees[i]!=field_degrees[j] for i in range(n) for j in range(n)):
            raise ValueError('transformation mixes different homogeneous field degrees')
    z=[Poly.variable(n,i,p) for i in range(n)]
    substitutions=[sum((new_to_old[i][j]*z[j] for j in range(n)),Poly(n,p=p)) for i in range(n)]
    transformed={}
    for name,field in fields.items():
        if len(field)!=n or any(f.n!=n or f.p!=p for f in field):
            raise ValueError('inconsistent vector-field ring')
        composed=[f.substitute(substitutions) for f in field]
        transformed[name]=tuple(sum((inv[i][j]*composed[j] for j in range(n)),Poly(n,p=p)) for i in range(n))
    return transformed,substitutions
