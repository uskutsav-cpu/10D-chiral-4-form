"""Independent Lorentzian Hodge/constrained-derivative implementation.

Conventions: eta=diag(-,+,...,+), epsilon_(0...d-1)=+1. Dense tensors store
lower components. A self-dual odd middle form is parameterized by components
containing the time index. Gradients are anti-self-dual LOWER tensors D with
sum Lambda^I D_I = homogeneous_degree * invariant.
"""
from __future__ import annotations
from functools import lru_cache
from itertools import combinations,permutations
from math import factorial,prod
import numpy as np
from .tensors import as_field_array,binary_contract,DEFAULT_BUDGET


def permutation_sign(seq):
    return -1 if sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))%2 else 1

@lru_cache(maxsize=8)
def layout(d=10,k=5):
    if d!=2*k or k%2!=1:
        raise ValueError('this real Lorentzian middle-form backend requires d=2k and odd k')
    basis=tuple(combinations(range(d),k));index={b:i for i,b in enumerate(basis)}
    perm=tuple(permutations(range(k)))
    ps=np.array([permutation_sign(q) for q in perm],dtype=np.int64)
    flat=np.empty((len(basis),len(perm)),dtype=np.int64)
    complements=[];star_sign=[];metric=[]
    for i,I in enumerate(basis):
        J=tuple(x for x in range(d) if x not in I)
        complements.append(index[J])
        star_sign.append(permutation_sign(I+J)*(-1 if 0 in J else 1))
        metric.append(-1 if 0 in I else 1)
        for j,q in enumerate(perm):
            flat[i,j]=np.ravel_multi_index(tuple(I[n] for n in q),(d,)*k)
    return {'basis':basis,'flat':flat,'perm_sign':ps,'complements':np.array(complements),
            'star_sign':np.array(star_sign),'metric':np.array(metric),
            'electric':np.array([i for i,I in enumerate(basis) if 0 in I])}


def compact(tensor,p,d=10,k=5):
    tensor=as_field_array(tensor,p)
    if tensor.shape!=(d,)*k:
        raise ValueError('wrong dense form shape')
    return tensor.reshape(-1)[layout(d,k)['flat'][:,0]].copy()


def dense(values,p,d=10,k=5):
    values=as_field_array(values,p)
    L=layout(d,k)
    if values.shape!=(len(L['basis']),):
        raise ValueError('wrong compact form dimension')
    result=np.zeros(d**k,dtype=np.int64)
    result[L['flat']]=(values[:,None]*L['perm_sign'][None,:])%p
    return result.reshape((d,)*k)


def hodge_compact(values,p,d=10,k=5):
    values=as_field_array(values,p);L=layout(d,k)
    if values.shape!=(len(L['basis']),):
        raise ValueError('wrong compact form dimension')
    return L['star_sign']*values[L['complements']]%p


def hodge(tensor,p,d=10,k=5):
    return dense(hodge_compact(compact(tensor,p,d,k),p,d,k),p,d,k)


def random_selfdual(seed,p,d=10,k=5):
    rng=np.random.default_rng(seed);L=layout(d,k)
    raw=np.zeros(len(L['basis']),dtype=np.int64)
    raw[L['electric']]=rng.integers(0,p,size=len(L['electric']),dtype=np.int64)
    # This gives unit electric coordinates, not half-electric coordinates.
    return dense((raw+hodge_compact(raw,p,d,k))%p,p,d,k)


def raise_all(tensor,p,d=10,k=5):
    values=compact(tensor,p,d,k)*layout(d,k)['metric']%p
    return dense(values,p,d,k)


def constrained_gradient(ambient_lower_input_gradient,p,d=10,k=5):
    ambient=as_field_array(ambient_lower_input_gradient,p)
    if ambient.shape!=(d,)*k:
        raise ValueError('wrong ambient derivative shape')
    L=layout(d,k)
    # Antisymmetrize, then lower all indices. Canonical slots only suffice.
    values=(ambient.reshape(-1)[L['flat']]*L['perm_sign'][None,:]).sum(axis=1)%p
    values=values*pow(factorial(k),-1,p)%p
    values=values*L['metric']%p
    values=(values-hodge_compact(values,p,d,k))*pow(2,-1,p)%p
    return dense(values,p,d,k)


def inner(a,b,p,d=10,k=5):
    L=layout(d,k)
    # Python int accumulation avoids overflow regardless of form dimension.
    return factorial(k)*sum(int(x)*int(y)*int(s) for x,y,s in
        zip(compact(a,p,d,k),compact(b,p,d,k),L['metric']))%p


def mixed_bilinear(a,b,p,d=10,k=5,budget=DEFAULT_BUDGET):
    """A_mu,rho... B^nu,rho... as a mixed d-by-d matrix."""
    a=as_field_array(a,p);b=raise_all(b,p,d,k)
    shared=tuple(range(1,k))
    return binary_contract(a,(0,)+shared,b,(k,)+shared,(0,k),p,budget)


def boost(tensor,t,p,d=10,k=5,axis=1,budget=DEFAULT_BUDGET):
    """Rational proper Lorentz boost of covariant slots; singular t rejected."""
    if not 1<=axis<d or (1-t*t)%p==0:
        raise ValueError('invalid boost parameter or spatial axis')
    inv=pow((1-t*t)%p,-1,p)
    c=(1+t*t)*inv%p;s=2*t*inv%p
    B=np.eye(d,dtype=np.int64);B[0,0]=B[axis,axis]=c;B[0,axis]=B[axis,0]=s
    result=as_field_array(tensor,p)
    for slot in range(k):
        labels=list(range(k));new=k+1
        out=labels.copy();out[slot]=new
        result=binary_contract(result,labels,B,(new,slot),out,p,budget)
    return result
