"""Independent Python integer checks of the orbit-compressed character tables.

Input weights are derived by enumerating the 252 exterior basis vectors.
The denominator is multiplied from20positive roots, not the Weyl permutation
sum used by the C++ engine. No graph data or expected Hilbert counts is read.
"""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, prod
from pathlib import Path
import json
import random
import time

ZERO=(0,)*5

def dominant(w):
    a=sorted(map(abs,w),reverse=True)
    if sum(x<0 for x in w)%2:a[-1]=-a[-1]
    return tuple(a)

def exterior_weights():
    vector=[tuple(sign*(i==j) for i in range(5)) for j in range(5) for sign in(-1,1)]
    full=Counter(tuple(sum(vector[j][i] for j in indices) for i in range(5)) for indices in combinations(range(10),5))
    chiral={}
    for w,m in full.items():
        if all(w):
            if prod(w)==1:chiral[w]=m
        else:
            assert m%2==0
            chiral[w]=m//2
    assert sum(full.values())==252 and sum(chiral.values())==126
    return chiral

def root_denominator():
    roots=[]
    for i,j in combinations(range(5),2):
        for sign in(-1,1):
            root=[0]*5;root[i]=1;root[j]=sign;roots.append(tuple(root))
    polynomial={ZERO:1}
    for root in roots:
        updated=Counter(polynomial)
        for w,m in polynomial.items():
            shifted=tuple(x-y for x,y in zip(w,root))
            updated[shifted]-=m
        polynomial={w:m for w,m in updated.items() if m}
    assert len(polynomial)==1920 and set(polynomial.values())=={-1,1}
    rho=tuple(sum(r[i] for r in roots)//2 for i in range(5))
    assert rho==(4,3,2,1,0)
    dimension=prod(Fraction(sum((rho[i]+1)*r[i] for i in range(5)),sum(rho[i]*r[i] for i in range(5))) for r in roots)
    assert dimension==126
    folded=Counter()
    for w,m in polynomial.items():folded[dominant(tuple(-x for x in w))]+=m
    return {w:m for w,m in folded.items() if m}

def orbit_size(w):
    a=tuple(map(abs,w));z=a.count(0)
    return factorial(5)*2**(5-z if z else 4)//prod(factorial(n) for n in Counter(a).values())

def candidates(n):
    for a in combinations_with_replacement(range(n+1),5):
        if (sum(a)-n)%2:continue
        w=tuple(reversed(a));yield w
        if a[0]:yield (*w[:4],-w[4])

def recurrence(n,w,characters,weights):
    result=0
    for k in range(1,n+1):
        previous=characters[n-k];limit=n-k
        for nu,m in weights.items():
            shifted=tuple(w[i]-k*nu[i] for i in range(5))
            if max(map(abs,shifted))>limit:continue
            result+=m*previous.get(dominant(shifted),0)
    assert result%n==0
    return result//n

def main():
    ap=argparse.ArgumentParser();ap.add_argument('prefix',type=Path);ap.add_argument('--full-through',type=int,default=12);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    chars=[{} for _ in range(23)]
    for line in Path(str(args.prefix)+'-weights.tsv').read_text().splitlines():
        n,*row=map(int,line.split());w=tuple(row[:5]);m=row[5]
        assert w==dominant(w) and m>0 and w not in chars[n]
        assert max(map(abs,w))<=n and (sum(w)-n)%2==0
        chars[n][w]=m
    weights=exterior_weights();denominator=root_denominator()
    assert chars[0]=={ZERO:1}
    assert chars[1]=={w:m for w,m in weights.items() if w==dominant(w)}
    results=[];rng=random.Random(202609079711)
    summary=json.loads(Path(str(args.prefix)+'.json').read_text())
    for n,character in enumerate(chars):
        start=time.monotonic()
        assert sum(m*orbit_size(w) for w,m in character.items())==comb(125+n,n)
        inv=sum(sign*character.get(w,0) for w,sign in denominator.items())
        assert inv==int(summary['degrees'][n]['invariants'])
        checked=0
        if n:
            targets=list(candidates(n))
            if n>args.full_through:targets=rng.sample(targets,min(128,len(targets)))
            for w in targets:
                assert recurrence(n,w,chars,weights)==character.get(w,0),(n,w)
                checked+=1
        results.append({'degree':n,'invariants':inv,'recurrence_entries_checked':checked,'recurrence_check_exhaustive':n<=args.full_through,'seconds':time.monotonic()-start})
        print(json.dumps(results[-1]),flush=True)
    # Every intermediate fits signed128: a conservative bound uses all1920
    # denominator terms and all126*n Newton contributions, before cancellation.
    bound=1920*126*22*comb(147,22)
    assert bound<2**127
    args.output.write_text(json.dumps({'passed':True,'weight_derivation':'exterior basis enumeration and chiral splitting','denominator_derivation':'product over20positive D5roots','signed128_conservative_bound':str(bound),'degrees':results},indent=2)+'\n')
if __name__=='__main__':main()
