"""Small exact sparse multivariate polynomials over Q or a validated F_p.

Expressions are loaded as coefficient/exponent JSON, never with eval().
"""
from __future__ import annotations
from fractions import Fraction
from types import MappingProxyType
from .exact import rational
from .finite_field import _check_prime_field_modulus

MAX_TERMS=100000

class Poly:
    def __init__(self,n,terms=None,p=None):
        if not isinstance(n,int) or n<0:
            raise ValueError('nonnegative variable count required')
        if p is not None:
            _check_prime_field_modulus(p)
        self.n,self.p=n,p
        out={}
        for monomial,c in (terms or {}).items():
            m=tuple(monomial)
            if len(m)!=n or any(not isinstance(e,int) or isinstance(e,bool) or e<0 for e in m):
                raise ValueError('invalid polynomial exponent tuple')
            q=rational(c)
            if p is not None:
                if q.denominator%p==0:
                    raise ValueError('prime divides coefficient denominator')
                q=q.numerator*pow(q.denominator,-1,p)%p
            if q:
                out[m]=q
        if len(out)>MAX_TERMS:
            raise RuntimeError('polynomial term resource limit exceeded')
        self.terms=MappingProxyType(out)

    @classmethod
    def constant(cls,n,c,p=None):
        return cls(n,{(0,)*n:c},p)
    @classmethod
    def variable(cls,n,i,p=None):
        if not 0<=i<n:
            raise ValueError('variable index out of range')
        m=[0]*n;m[i]=1
        return cls(n,{tuple(m):1},p)
    def _other(self,x):
        if not isinstance(x,Poly):
            return Poly.constant(self.n,x,self.p)
        if (self.n,self.p)!=(x.n,x.p):
            raise ValueError('polynomial rings do not match')
        return x
    def __add__(self,other):
        other=self._other(other);out=dict(self.terms)
        for m,c in other.terms.items():
            out[m]=out.get(m,0)+c
        return Poly(self.n,out,self.p)
    __radd__=__add__
    def __neg__(self):
        return Poly(self.n,{m:-c for m,c in self.terms.items()},self.p)
    def __sub__(self,other):
        return self+-self._other(other)
    def __rsub__(self,other):
        return self._other(other)+-self
    def __mul__(self,other):
        other=self._other(other);out={}
        for a,c in self.terms.items():
            for b,d in other.terms.items():
                m=tuple(x+y for x,y in zip(a,b))
                out[m]=out.get(m,0)+c*d
                if self.p:
                    out[m]%=self.p
                if len(out)>MAX_TERMS:
                    raise RuntimeError('polynomial multiplication term budget exceeded')
        return Poly(self.n,out,self.p)
    __rmul__=__mul__
    def __pow__(self,e):
        if not isinstance(e,int) or e<0:
            raise ValueError('polynomial power must be a nonnegative integer')
        result=Poly.constant(self.n,1,self.p);base=self
        while e:
            if e&1:
                result=result*base
            e//=2
            if e:
                base=base*base
        return result
    def __bool__(self):
        return bool(self.terms)
    def __eq__(self,other):
        try:
            other=self._other(other)
        except (ValueError,TypeError):
            return False
        return self.terms==other.terms
    def derivative(self,i):
        if not 0<=i<self.n:
            raise ValueError('derivative variable index out of range')
        out={}
        for m,c in self.terms.items():
            if m[i]:
                a=list(m);a[i]-=1
                out[tuple(a)]=c*m[i]
        return Poly(self.n,out,self.p)
    def evaluate(self,point):
        if len(point)!=self.n:
            raise ValueError('point has wrong dimension')
        point=[self._other(x).terms.get((0,)*self.n,0) if not isinstance(x,Poly) else x for x in point]
        if any(isinstance(x,Poly) for x in point):
            raise TypeError('use substitute for polynomial points')
        value=0
        for m,c in self.terms.items():
            term=c
            for x,e in zip(point,m):
                term=term*(pow(int(x),e,self.p) if self.p else x**e)
                if self.p:
                    term%=self.p
            value+=term
        return value%self.p if self.p else Fraction(value)
    def substitute(self,values,*,n_out=None):
        if len(values)!=self.n:
            raise ValueError('substitution length mismatch')
        n_out=values[0].n if values else (0 if n_out is None else n_out)
        if any(not isinstance(v,Poly) or v.n!=n_out or v.p!=self.p for v in values):
            raise ValueError('substitution rings mismatch')
        result=Poly(n_out,p=self.p)
        for m,c in self.terms.items():
            term=Poly.constant(n_out,c,self.p)
            for v,e in zip(values,m):
                if e:
                    term=term*(v**e)
            result=result+term
        return result
    def truncate(self,degree):
        return Poly(self.n,{m:c for m,c in self.terms.items() if sum(m)<=degree},self.p)
    def to_json(self):
        return {'nvars':self.n,'prime':self.p,'terms':[{'powers':list(m),'coefficient':str(c)} for m,c in sorted(self.terms.items())]}
    @classmethod
    def from_json(cls,record):
        terms={}
        for row in record['terms']:
            m=tuple(row['powers'])
            if m in terms:
                raise ValueError('duplicate polynomial term')
            terms[m]=row['coefficient']
        return cls(record['nvars'],terms,record.get('prime'))


def derivation(field,polynomial):
    if len(field)!=polynomial.n:
        raise ValueError('vector field dimension mismatch')
    return sum((a*polynomial.derivative(i) for i,a in enumerate(field)),Poly(polynomial.n,p=polynomial.p))


def bracket(x,y):
    """[X,Y] acting on coordinate functions, including fields zero at the seed."""
    if len(x)!=len(y) or not x:
        raise ValueError('nonempty vector fields of equal dimension required')
    return tuple(derivation(x,f)-derivation(y,g) for f,g in zip(y,x))


def field_at(field,point):
    return [f.evaluate(point) for f in field]
