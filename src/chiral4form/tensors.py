"""Memory-bounded, overflow-checked binary tensor contractions over F_p.

Every contraction reduces modulo p before another multiplication. Float64 BLAS
is used only with a rigorous exact-integer dot-product bound; otherwise int64.
Arrays must have unique labels; hyperedges are not supported here.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import prod
import numpy as np
from .finite_field import _check_prime_field_modulus

@dataclass(frozen=True)
class TensorBudget:
    max_bytes: int = 512 * 1024**2
    max_multiply_adds: int = 5_000_000_000
    exact_blas: bool = True
    def __post_init__(self):
        if type(self.max_bytes) is not int or type(self.max_multiply_adds) is not int or type(self.exact_blas) is not bool:
            raise TypeError('integer budgets and a boolean exact_blas flag required')
        if self.max_bytes<1 or self.max_multiply_adds<1:
            raise ValueError('tensor resource budgets must be positive')

DEFAULT_BUDGET=TensorBudget()


def as_field_array(value,p):
    _check_prime_field_modulus(p)
    if p>65521:
        raise ValueError('int64 tensor backend is restricted to primes <= 65521')
    arr=np.asarray(value)
    if arr.dtype.kind not in 'iu':
        raise TypeError('integer tensor arrays required')
    # Unsigned values may not fit signed int64 before modular reduction.
    if arr.dtype.kind=='u':
        arr=arr % np.uint64(p)
    return np.asarray(arr,dtype=np.int64)%p


def binary_contract(a,la,b,lb,out,p,budget=DEFAULT_BUDGET):
    """Einstein contraction of two tensors, summing every common label.

    This restricted primitive suffices for degree-regular contraction graphs
    and reverse differentiation. Restrictions are checked, not assumed.
    """
    a=as_field_array(a,p);b=as_field_array(b,p)
    la,lb,out=tuple(la),tuple(lb),tuple(out)
    if a.ndim!=len(la) or b.ndim!=len(lb) or len(set(la))!=len(la) or len(set(lb))!=len(lb):
        raise ValueError('tensor rank/label mismatch or repeated label')
    common=[x for x in la if x in lb]
    left=[x for x in la if x not in lb]
    right=[x for x in lb if x not in la]
    if len(out)!=len(left)+len(right) or set(out)!=set(left+right):
        raise ValueError('output must be the uncontracted labels exactly once')
    for x in common:
        if a.shape[la.index(x)]!=b.shape[lb.index(x)]:
            raise ValueError('contracted dimensions differ')
    ashape=[a.shape[la.index(x)] for x in left]
    bshape=[b.shape[lb.index(x)] for x in right]
    k=prod(a.shape[la.index(x)] for x in common)
    m,n=prod(ashape),prod(bshape)
    if k*(p-1)**2>np.iinfo(np.int64).max:
        raise OverflowError('integer dot-product bound exceeds int64')
    # Reserve transposed copies, output, and modular-reduction temporaries.
    requested=8*(4*a.size+4*b.size+4*m*n)
    if requested>budget.max_bytes:
        raise MemoryError(f'contraction requires <= {requested} bytes; budget {budget.max_bytes}')
    if m*n*k>budget.max_multiply_adds:
        raise MemoryError(f'contraction work {m*n*k} exceeds limit {budget.max_multiply_adds}')
    aa=np.transpose(a,[la.index(x) for x in left+common]).reshape(m,k)
    bb=np.transpose(b,[lb.index(x) for x in common+right]).reshape(k,n)
    if budget.exact_blas and k*(p-1)**2 < 2**52 and m*n*k > 100_000:
        # Integer factors/products and every partial nonnegative dot sum are
        # exactly representable in float64. The extra factor-two margin is
        # conservative. Do not reuse this path for arbitrary real operands.
        value=(aa.astype(np.float64)@bb.astype(np.float64)).astype(np.int64)%p
    else:
        value=(aa@bb)%p
    value=value.reshape(tuple(ashape+bshape))
    order=left+right
    return np.transpose(value,[order.index(x) for x in out])


def matrix_product(a,b,p,budget=DEFAULT_BUDGET):
    return binary_contract(a,(0,1),b,(1,2),(0,2),p,budget)
