"""Quadratic localization a+b*r with r^2=I4, on the stated I4!=0 patch.

Base symbols are formal algebraic variables. This is NOT a claim that the
81 functionally independent invariants generate a freely polynomial ring.
No arbitrary input strings are evaluated and no floating square roots occur.
"""
from __future__ import annotations
from dataclasses import dataclass
import sympy as sp


def _expr(value):
    if isinstance(value,(str,float,bool)):
        raise TypeError('supply exact numeric or SymPy expression objects')
    return sp.sympify(value)

@dataclass(frozen=True)
class QuadraticElement:
    a: object
    b: object
    square: object
    def __post_init__(self):
        object.__setattr__(self,'a',sp.cancel(_expr(self.a)))
        object.__setattr__(self,'b',sp.cancel(_expr(self.b)))
        object.__setattr__(self,'square',_expr(self.square))
    def _other(self,x):
        if not isinstance(x,QuadraticElement):
            return QuadraticElement(x,0,self.square)
        if x.square!=self.square:
            raise ValueError('different quadratic extensions')
        return x
    def __add__(self,other):
        other=self._other(other)
        return QuadraticElement(self.a+other.a,self.b+other.b,self.square)
    __radd__=__add__
    def __neg__(self):
        return QuadraticElement(-self.a,-self.b,self.square)
    def __sub__(self,other):
        return self+-self._other(other)
    def __rsub__(self,other):
        return self._other(other)+-self
    def __mul__(self,other):
        other=self._other(other)
        return QuadraticElement(self.a*other.a+self.b*other.b*self.square,
                                self.a*other.b+self.b*other.a,self.square)
    __rmul__=__mul__
    def inverse(self):
        norm=sp.cancel(self.a**2-self.square*self.b**2)
        if norm==0:
            raise ZeroDivisionError('zero norm; cannot invert on this algebraic patch')
        return QuadraticElement(self.a/norm,-self.b/norm,self.square)
    def derivative(self,variable):
        return QuadraticElement(sp.diff(self.a,variable),
            sp.diff(self.b,variable)+self.b*sp.diff(self.square,variable)/(2*self.square),self.square)
    def is_zero(self):
        return self.a==0 and self.b==0
    def euler(self,weights):
        return sum((weight*var*self.derivative(var) for var,weight in weights.items()),
                   QuadraticElement(0,0,self.square))


def modmax_symbolic_report():
    I4,I8,I12,b=sp.symbols('I4 I8 I12 b')
    root=QuadraticElement(0,1,I4);V=b*root
    conformal=(V.euler({I4:4})-2*V).is_zero()
    a=(1-sp.Rational(24,7)*b*b)/48;c=b*b/(12*I4)
    stress_square=sp.expand(a*a*I4+2*a*c*I8+c*c*I12)
    correction=b*b*(1-sp.Rational(24,7)*b*b)*I8/(288*I4)+b**4*I12/(144*I4**2)
    target=I4*(1-sp.Rational(24,7)*b*b)**2/(4*24**2)+correction
    exact=sp.cancel(stress_square-target)==0
    return {'status':'published_identity_reproduction_not_new_reachability_result','patch':'I4 != 0; r^2=I4',
        'conformal_homogeneity_verified':conformal,'stress_square_identity_verified':exact,
        'stress_square':str(stress_square),'correction':str(correction),
        'source':'arXiv:2509.14351v2 equations (3.12)-(3.16)',
        'alternative_pure_stress_flow':'UNRESOLVED; reproduction does not exclude another generator',
        'real_branch':'not selected; algebraic calculation only'}
