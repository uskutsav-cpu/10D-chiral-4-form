from fractions import Fraction
from dataclasses import replace
from copy import deepcopy
import numpy as np
import pytest
from chiral4form.finite_field import (
    _check_prime_field_modulus, rref, matrix_rank, nullspace, matvec,
    matmul, solve_many, determinant, independent_rows, span_contains)
from chiral4form.certificates import build_quotient_certificate
from chiral4form.exact import rational, rref_q, nullspace_q, crt, rational_reconstruct, reconstruct_with_holdout

@pytest.mark.parametrize("p", [3,5,101,30011,32749,65521,2147483647])
def test_prime_validation(p):
    _check_prime_field_modulus(p)

@pytest.mark.parametrize("p", [0,1,2,4,9,15,341,561,1105,3215031751,-101,2**64+13])
def test_reject_bad_moduli(p):
    with pytest.raises((ValueError,TypeError)):
        _check_prime_field_modulus(p)

@pytest.mark.parametrize("p", [True,101.0,"101"])
def test_reject_noninteger_moduli(p):
    with pytest.raises((ValueError,TypeError)):
        _check_prime_field_modulus(p)

@pytest.mark.parametrize("rows", [[[1,2],[3]], [[1.0]], [[True]]])
def test_strict_matrix_inputs(rows):
    with pytest.raises((TypeError,ValueError)):
        rref(rows,101)

@pytest.mark.parametrize("p", [101,103,30011])
def test_rank_nullity_and_many_rhs(p):
    A=[[1,2,3],[0,1,4],[0,0,1],[1,3,7]]
    X=[[2,3],[5,7],[11,13]]
    B=matmul(A,X,p)
    assert solve_many(A,B,p)==X
    assert matrix_rank(A,p)==3
    assert determinant(A[:3],p)==1
    assert len(independent_rows(A,p))==3
    assert nullspace(A,p)==[]
    assert span_contains(A,[7,8,9],p)

@pytest.mark.parametrize("A,B", [([[1,0]],[[1]]), ([[1],[1]],[[1],[2]])])
def test_nonunique_or_inconsistent_system_rejected(A,B):
    with pytest.raises(ValueError):
        solve_many(A,B,101)

@pytest.mark.parametrize("n", [0,1,5])
def test_empty_matrix_preserves_ambient_dimension(n):
    N=nullspace([],101,ncols=n)
    assert len(N)==n
    assert matrix_rank([],101,ncols=n)==0
    assert determinant([],101)==1

@pytest.mark.parametrize("mutation", ["annihilators","reachable_rank","minor_determinant","quotient_dim"])
def test_forged_certificates_rejected(mutation):
    c=build_quotient_certificate([[1,0,1,0],[0,1,0,1]],4,101)
    assert c.verify([[1,0,1,0],[0,1,0,1]])
    changes={"annihilators": [c.annihilators[0],c.annihilators[0]],
             "reachable_rank":1,"minor_determinant":0,"quotient_dim":1}
    forged=replace(c,**{mutation:changes[mutation]})
    assert not forged.verify([[1,0,1,0],[0,1,0,1]])

def test_zero_row_certificate():
    c=build_quotient_certificate([],3,101)
    assert c.verify([]) and c.quotient_dim==3

def test_numpy_integer_matrix():
    assert matrix_rank(np.eye(3,dtype=np.int64).tolist(),101)==3

@pytest.mark.parametrize("value", [0,Fraction(32,3),Fraction(-17,13),400,-96])
def test_reconstruction_and_disjoint_holdout(value):
    value=Fraction(value);ps=[30011,30013]
    residues=[value.numerator*pow(value.denominator,-1,p)%p for p in ps]
    a,m=crt(residues,ps)
    assert rational_reconstruct(a,m,500,500)==value
    hp=30029;hr=value.numerator*pow(value.denominator,-1,hp)%hp
    result=reconstruct_with_holdout(residues,ps,hr,hp,500)
    assert Fraction(result['value'])==value
    assert result['not_a_symbolic_identity_proof']
    with pytest.raises(ValueError):
        reconstruct_with_holdout(residues,ps,(hr+1)%hp,hp,500)

def test_crt_rejects_repeated_primes_and_missing_height():
    with pytest.raises(ValueError):crt([1,2],[101,101])
    with pytest.raises(ValueError):rational_reconstruct(12,101,10,10)
    with pytest.raises(ValueError):reconstruct_with_holdout([1],[101],1,101)

@pytest.mark.parametrize("bad", [0.5,True])
def test_exact_rationals_reject_floats(bad):
    with pytest.raises(TypeError):rational(bad)

def test_rational_linear_algebra():
    a,piv=rref_q([[1,2,3],[2,4,6]])
    assert piv==[0] and a[0]==[1,2,3]
    assert len(nullspace_q([[1,2,3]],3))==2
