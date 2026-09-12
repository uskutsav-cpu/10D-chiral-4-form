"""Conditional analytic sextic obstruction, separately from sampled fits.

Hypotheses: classical HLS (2.33), analytic scalar generators at the origin,
no derivatives, V starts at field degree 4, tau=48T, fixed homogeneous
coordinate basis. Negative powers and ModMax square roots are excluded.
"""
from __future__ import annotations
from fractions import Fraction
from .polynomial import Poly,derivation
from .weighted import stress_generators


def analytic_fields():
    # Coordinates c4, c6J, c6K multiply the committed graph invariants.
    n=3;c4,c6J,c6K=[Poly.variable(n,i) for i in range(n)];z=Poly(n)
    return ('I4_1','I6_1','I6_2'),{
        'tr1':(20*c4,40*c6J,40*c6K),
        'tr2':(Poly.constant(n,2),z,z),
        'tr3':(z,Poly.constant(n,Fraction(32,3)),z)}


def sextic_report():
    ids,fields=analytic_fields();q=Fraction(125,3)*Poly.variable(3,2)
    expected={'tr1':40*q,'tr2':Poly(3),'tr3':Poly(3)}
    checks={name:derivation(field,q)==expected[name] for name,field in fields.items()}
    if not all(checks.values()):
        raise AssertionError('sextic analytic algebra check failed')
    catalog=[g.id for g in stress_generators(6)]
    if set(catalog)!=set(fields):
        raise AssertionError('degree-six generator catalogue is incomplete')
    return {'status':'analytic_derivation_under_explicit_physics_hypotheses',
        'coordinate_ids':list(ids),'generator_catalogue':catalog,'checks':checks,
        'equation':'dq6/dlambda = 40 * u_tr1(lambda) * q6',
        'q6':'(125/3) * c_I6_2; this is a coupling-space coordinate, not K6(F)=0',
        'basis_relation':'K6 = -I6_1/1125 + 3*I6_2/125',
        'proof_steps':[
            'All polynomial stress monomials eligible through degree six are tr1,tr2,tr3.',
            'Equation (2.33) and anti-self-duality imply Tr(tau)=sum_d 10(d-2)V_d.',
            'The degree-six term of Tr(tau^2) is proportional to Tr(M) and vanishes.',
            'The degree-six part of Tr(tau^3) is Tr(M^3)=(32/3)I6_1.',
            'Consequently only tr1 acts on q6, and acts homogeneously by 40q6.'
        ],'higher_field_degrees':'cannot feed degree six under positive-degree analytic polynomial grading',
        'does_not_prove':['full orbit reachability','minimal generalized completion','nonanalytic ModMax classification'],
        'external_inputs':['HLS arXiv:2509.14351v2 eq. (2.33)',
                           'predecessor sextic normalization; tensor smoke tests cross-check Tr(M^3)'],
        'mentor_review_required':['conventions','flow-law class','intrinsic K6 identification']}


def verify_sampled_sextic(fields,ids,prime):
    if tuple(ids)!=('I4_1','I6_1','I6_2'):
        raise ValueError('sextic cross-check requires the pinned three-coordinate order')
    _,analytic=analytic_fields()
    for name,field in analytic.items():
        expected=tuple(Poly(3,dict(f.terms),prime) for f in field)
        if name not in fields or fields[name]!=expected:
            raise ValueError(f'sampled field {name} does not match the independent analytic sextic model')
    return True
