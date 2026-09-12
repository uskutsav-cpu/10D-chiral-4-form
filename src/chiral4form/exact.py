"""Characteristic-zero linear algebra and bounded modular reconstruction."""
from __future__ import annotations
from fractions import Fraction
from math import gcd, isqrt
from .finite_field import _check_prime_field_modulus

def rational(x) -> Fraction:
    if isinstance(x,(float,bool)):
        raise TypeError("use integers or exact rational strings, not floating point")
    return Fraction(x)


def rref_q(rows, ncols=None):
    n = len(rows[0]) if rows else (ncols or 0)
    if n < 0 or (ncols is not None and ncols != n):
        raise ValueError('column count mismatch')
    if any(len(row)!=n for row in rows):
        raise ValueError("ragged matrix")
    a = [[rational(x) for x in row] for row in rows]
    piv,r = [],0
    for c in range(n):
        k = next((i for i in range(r,len(a)) if a[i][c]),None)
        if k is None:
            continue
        a[r],a[k]=a[k],a[r]
        f = a[r][c]
        a[r]=[x/f for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                f=a[i][c]
                a[i]=[x-f*y for x,y in zip(a[i],a[r])]
        piv.append(c)
        r+=1
        if r==len(a):
            break
    return a,piv


def nullspace_q(rows,ncols):
    rr,piv=rref_q(rows,ncols)
    out=[]
    for f in range(ncols):
        if f in piv:
            continue
        v=[Fraction(0)]*ncols
        v[f]=Fraction(1)
        for i,c in enumerate(piv):
            v[c]=-rr[i][f]
        out.append(v)
    return out


def crt(residues,primes):
    if len(residues)!=len(primes) or not primes:
        raise ValueError("nonempty aligned CRT inputs required")
    x,m=0,1
    for a,p in zip(residues,primes):
        _check_prime_field_modulus(p)
        if gcd(m,p)!=1:
            raise ValueError("CRT primes must be distinct")
        x+=((a-x)*pow(m,-1,p)%p)*m
        m*=p
    return x%m,m


def rational_reconstruct(a,m,numerator_bound,denominator_bound):
    """Unique bounded fraction congruent to a modulo m, if it exists.

    A user-supplied height bound constrains the *candidate*. It is NOT proof
    that unknown coefficients obey that bound; verify identities separately.
    """
    N,D=numerator_bound,denominator_bound
    if N<0 or D<1 or 2*N*D>=m:
        raise ValueError("reconstruction requires 2*N*D < modulus")
    a%=m
    if a==0:
        return Fraction(0)
    r0,r1,t0,t1=m,a,0,1
    while abs(r1)>N:
        if r1==0:
            raise ValueError("no bounded rational reconstruction")
        q=r0//r1
        r0,r1=r1,r0-q*r1
        t0,t1=t1,t0-q*t1
    if t1==0:
        raise ValueError("zero reconstructed denominator")
    f=Fraction(r1,t1)
    if abs(f.numerator)>N or f.denominator>D or gcd(f.denominator,m)!=1:
        raise ValueError("candidate exceeds bounds or has bad denominator")
    if (f.numerator-a*f.denominator)%m:
        raise ValueError("reconstruction congruence failed")
    return f


def reconstruct_with_holdout(residues,primes,holdout_residue,holdout_prime,bound=None):
    if holdout_prime in primes:
        raise ValueError("holdout prime must not be a fit prime")
    _check_prime_field_modulus(holdout_prime)
    a,m=crt(residues,primes)
    B=isqrt((m-1)//2) if bound is None else bound
    f=rational_reconstruct(a,m,B,B)
    if f.denominator%holdout_prime==0 or f.numerator*pow(f.denominator,-1,holdout_prime)%holdout_prime!=holdout_residue%holdout_prime:
        raise ValueError("independent prime rejects the reconstructed candidate")
    return {"value":str(f),"status":"bounded_candidate_with_modular_holdout",
            "not_a_symbolic_identity_proof":True,"modulus":str(m),"bound":B}
