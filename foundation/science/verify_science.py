#!/usr/bin/env python3
"""Check frozen evidence binding and the new logical closure; no tensor rerun."""
import hashlib,json
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parent

def rank(rows,prime=None):
    a=[[int(x)%prime if prime else Fraction(x) for x in row] for row in rows]
    n=0
    for c in range(len(a[0])):
        k=next((k for k in range(n,len(a)) if a[k][c]),None)
        if k is None:continue
        a[n],a[k]=a[k],a[n]
        v=pow(a[n][c],-1,prime) if prime else 1/a[n][c]
        a[n]=[(x*v)%prime if prime else x*v for x in a[n]]
        for k in range(n+1,len(a)):
            v=a[k][c]
            a[k]=[(x-v*y)%prime if prime else x-v*y for x,y in zip(a[k],a[n])]
        n+=1
        if n==len(a):break
    return n

def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    load=lambda p:json.loads((ROOT/p).read_text())
    for p,h in load('manifest.json').items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    h={r['degree']:int(r['invariants']) for r in load('character/independent-d5.json')['degrees']}
    checks=load('character/python-character-check.json');assert checks['passed']
    for d in checks['degrees']:
        assert d['invariants']==h[d['degree']]
        if d['degree']<=12:assert d['recurrence_check_exhaustive']
    assert [h[n] for n in (4,6,8,10,12)]==[1,2,7,14,72]
    source=load('evidence/source-certificate.json')
    assert source['ranks']['atlas8']==h[8] and source['ranks']['atlas10']==h[10]
    for degree in ('degree8','degree10'):
        assert source['modulus']>2*max(source['bounds'][degree]['integer_numerator_absolute_bounds'])
    m=source['maps']['degree10'];a=m['matrix_12x14']
    products=[[int(i==j) for i in range(14)] for j in (12,13)]
    assert rank(a)==12 and rank(products)==2 and rank(a+products)==13
    assert rank([r[:12] for r in a])==11
    for w in m['product_intersection_witnesses']:
        row=[sum(Fraction(c)*Fraction(r[j]) for c,r in zip(w['source_combination'],a)) for j in range(14)]
        assert row==[Fraction(0)]*12+list(map(Fraction,w['product_coefficients']))
    for f in load('evidence/degree12-certificate.json')['fields']:
        assert rank(f['stacked_gradient'],f['prime'])==h[12]==72
        assert int(f['integer_minor_determinant'])%f['prime']==f['determinant_mod_prime']!=0
    cert=load('rank81/certificate.json');assert cert['passed']
    for p,key in [('rank81/input.json','input_sha256'),('rank81/verify.py','implementation_sha256')]:
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==cert[key]
    for f in cert['points']:
        p=f['prime'];assert rank(f['jacobian'],p)==81 and rank(f['orbit_rows'],p)==45
        assert int(f['integer_minor_determinant'])%p==f['determinant_mod_p']!=0
        assert int(f['orbit_integer_minor_determinant'])%p==f['orbit_determinant_mod_p']!=0
    print(json.dumps({'passed':True,'dimensions':h,'degree12_basis':72,'source10_ranks':[12,2,13,1],'rank81':True,'scope':'hash binding and saved arithmetic; full fresh tensor computation is rank81/verify.py'},indent=2))
if __name__=='__main__':main()
