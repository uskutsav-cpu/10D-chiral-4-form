from fractions import Fraction
import numpy as np
import pytest
from chiral4form.forms import random_selfdual,hodge,inner,boost,compact,dense,mixed_bilinear,layout
from chiral4form.tensors import TensorBudget,binary_contract,matrix_product,as_field_array
from chiral4form.graphs import parse_graph,evaluate_graph,validate_graph
from chiral4form.stress import residue,evaluate_generators,trace_series
from chiral4form.weighted import stress_generators

@pytest.mark.parametrize("p", [101,30011])
def test_lorentzian_hodge(p):
    f=random_selfdual(1964,p)
    assert f.shape==(10,)*5
    assert np.array_equal(hodge(f,p),f)
    assert np.array_equal(hodge(hodge(f,p),p),f)
    assert inner(f,f,p)==0
    assert np.array_equal(dense(compact(f,p),p),f)
    assert len(layout()['electric'])==126
    assert np.array_equal(f.swapaxes(0,1),(-f)%p)

@pytest.mark.parametrize("p", [101,30011])
def test_real_boost_preserves_scalar(registry,p):
    f=random_selfdual(9871,p);g=boost(f,2,p)
    assert np.array_equal(hodge(g,p),g)
    assert registry.evaluate('I4_1',f,p)==registry.evaluate('I4_1',g,p)

@pytest.mark.parametrize("name", ['I4_1','I6_1'])
def test_constrained_derivative_and_trace_normalization(registry,name):
    p=101;f=random_selfdual(1776,p)
    value,gradient=registry.evaluate(name,f,p,gradient=True)
    assert np.array_equal(hodge(gradient,p),(-gradient)%p)
    assert inner(f,gradient,p)==registry.items[name].degree*value%p
    M=mixed_bilinear(f,f,p);M2=matrix_product(M,M,p)
    trace=int(np.trace(M2 if name=='I4_1' else matrix_product(M2,M,p)))%p
    factor=2 if name=='I4_1' else residue(Fraction(32,3),p)
    assert trace==factor*value%p

@pytest.mark.parametrize("shape", [(3,5,4),(49,60,51)])
def test_blas_path_matches_integer_reference(shape):
    m,k,n=shape;p=30011;rng=np.random.default_rng(5)
    a=rng.integers(0,p,(m,k),dtype=np.int64);b=rng.integers(0,p,(k,n),dtype=np.int64)
    exact=binary_contract(a,(0,1),b,(1,2),(0,2),p,TensorBudget(exact_blas=False))
    fast=binary_contract(a,(0,1),b,(1,2),(0,2),p,TensorBudget(exact_blas=True))
    assert np.array_equal(exact,fast)
    assert int(fast[0,0])==sum(int(x)*int(y) for x,y in zip(a[0],b[:,0]))%p

def test_tensor_budget_and_labels_fail_early():
    a=np.ones((10,10),dtype=np.int64)
    with pytest.raises(MemoryError):matrix_product(a,a,101,TensorBudget(max_bytes=10))
    with pytest.raises(MemoryError):matrix_product(a,a,101,TensorBudget(max_multiply_adds=1))
    with pytest.raises(ValueError):binary_contract(a,(0,0),a,(0,1),(1,),101)
    with pytest.raises(ValueError):binary_contract(a,(0,1),a,(1,2),(0,0),101)
    with pytest.raises(TypeError):as_field_array(np.ones(3),101)
    with pytest.raises(ValueError):as_field_array(np.ones(3,dtype=int),1000003)

@pytest.mark.parametrize("s", ['nope','n4[02^8]','n4[00^5]','n4[02^3,02^3]'])
def test_bad_graph_rejected(s):
    with pytest.raises(ValueError):parse_graph(s)

def test_high_order_graph_delimiters():
    # Disconnected valid 12-vertex graph used only to check parser precision.
    label='n12['+','.join(f'{i}-{i+1}^3,{i}-{i+2}^2,{i+1}-{i+3}^2,{i+2}-{i+3}^3' for i in (0,4,8))+']'
    g=parse_graph(label)
    assert len(g)==12 and all(sum(row)==5 for row in g)

@pytest.mark.parametrize("degree,count", [(6,3),(8,7),(12,18)])
def test_complete_trace_generator_catalogue(degree,count):
    g=stress_generators(degree)
    assert len(g)==count
    assert len({x.id for x in g})==count
    assert all(x.leading_degree<=degree for x in g)

def test_genuine_extra_generator_not_seed(registry):
    f=random_selfdual(155,p:=101)
    targets,catalog= evaluate_generators(registry,f,p,6,('I6_2',))
    assert targets[('S:I6_2',6,())]==registry.evaluate('I6_2',f,p)
    assert any(g.id=='S:I6_2' for g in catalog)
