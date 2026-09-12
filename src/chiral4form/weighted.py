"""Finite weighted generator catalogues with explicit truncation semantics."""
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class Generator:
    factors: tuple[str,...]
    leading_degree: int
    @property
    def id(self):
        return '*'.join(self.factors)


def weighted_monomials(names,weights,max_degree,*,include_constant=False,max_count=100000):
    if len(names)!=len(weights) or len(set(names))!=len(names) or any(w<=0 for w in weights):
        raise ValueError('unique names and strictly positive aligned weights required')
    if max_degree<0:
        raise ValueError('nonnegative cutoff required')
    out=[]
    def rec(start,degree,factors):
        if factors or include_constant:
            out.append(Generator(tuple(factors),degree))
            if len(out)>max_count:
                raise RuntimeError('weighted generator catalogue exceeds resource limit')
        for i in range(start,len(names)):
            if degree+weights[i]<=max_degree:
                rec(i,degree+weights[i],factors+[names[i]])
    rec(0,0,[])
    return tuple(sorted(out,key=lambda x:(x.leading_degree,len(x.factors),x.factors)))


def stress_generators(max_degree,extras=None):
    # Free tau is traceless. Tr(tau) begins at degree FOUR, not two.
    names=[f'tr{k}' for k in range(1,11)]
    weights=[4]+[2*k for k in range(2,11)]
    for name,degree in sorted((extras or {}).items()):
        if name in names or degree<4 or degree%2:
            raise ValueError('extra generators must be distinct homogeneous invariants of even degree >=4')
        names.append(name);weights.append(degree)
    return weighted_monomials(names,weights,max_degree)
