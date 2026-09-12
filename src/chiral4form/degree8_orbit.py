"""An exact degree-eight REDUCED MODEL orbit and a nonlinear obstruction.

The physics-to-polynomial map is cross-checked against independent tensor
samples. The algebraic proof below is exact over QQ conditional on that map;
it is not an unconditional all-orders result or a novelty assertion.
"""
from fractions import Fraction
from .exact import rational
from .polynomial import Poly,derivation
from .reachability import linear_invariant_hull,lie_rank
from .ideals import invariant_ideal_report

IDS=('I4_1','I6_1','I6_2','I8_1','I8_2','I8_3','I8_4','I8_5','I8_6','I4_1^2')


def analytic_fields():
    n=10;c=[Poly.variable(n,i) for i in range(n)];z=Poly(n)
    def vector(entries):
        v=[z]*n
        for i,x in entries.items():
            v[i]=x if isinstance(x,Poly) else Poly.constant(n,x)
        return tuple(v)
    fields={
        'tr1':vector({i:(20 if i==0 else 40 if i<3 else 60)*x for i,x in enumerate(c)}),
        'tr1*tr1':vector({9:400*c[0]**2}),
        'tr1*tr2':vector({9:40*c[0]}),
        'tr2':vector({0:2,3:16*c[0]**2,4:-96*c[0]**2,9:-8*c[0]**2}),
        'tr2*tr2':vector({9:4}),
        'tr3':vector({1:Fraction(32,3),9:12*c[0]}),
        'tr4':vector({3:1})}
    return IDS,fields


def constraints():
    c=[Poly.variable(10,i) for i in range(10)]
    return (c[2],c[5],c[6],c[7],c[8],c[4]+16*c[0]**3)


def construct_endpoint(a,b,c,d):
    """Four controlled segments realize every point of the constraint graph.

    a,b,c,d are exact target c4,c6J,c8_1,c8_product coordinates. Signed
    generator amounts are allowed. No positivity/causality condition is imposed.
    """
    a,b,c,d=map(rational,(a,b,c,d))
    x=[Fraction(0)]*10
    # Integrate tr2 from c4=0 to c4=a exactly.
    x[0]=a;x[3]=Fraction(8,3)*a**3;x[4]=-16*a**3;x[9]=-Fraction(4,3)*a**3
    stage1=a/2
    # tr3 transports only c6J and c8_product in this truncation.
    stage2=Fraction(3,32)*b;x[1]=b;x[9]+=12*a*stage2
    stage3=c-x[3];x[3]+=stage3
    stage4=(d-x[9])/4;x[9]+=4*stage4
    if any(q.evaluate(x) for q in constraints()):
        raise AssertionError('constructed endpoint violates the orbit ideal')
    return {'endpoint':[str(q) for q in x],'stages':[
        {'generator':'tr2','amount':str(stage1)},
        {'generator':'tr3','amount':str(stage2)},
        {'generator':'tr4','amount':str(stage3)},
        {'generator':'tr2*tr2','amount':str(stage4)}]}


def report():
    ids,fields=analytic_fields();qs=constraints();q8=qs[-1]
    transport={name:derivation(f,q8)==(60*q8 if name=='tr1' else Poly(10)) for name,f in fields.items()}
    if not all(transport.values()):
        raise AssertionError('nonlinear octic transport identity failed')
    tangent=invariant_ideal_report(fields,qs)
    hull=linear_invariant_hull(list(fields.values()))
    lie=lie_rank(list(fields.values()),max_depth=3)
    return {'status':'exact_QQ_reduced_model_result_conditional_on_physics_basis_map',
        'coordinate_ids':list(ids),'nonlinear_obstruction':'c_I8_2 + 16*c_I4_1^3',
        'transport_equation':'dOmega8/dlambda = 60*u_tr1(lambda)*Omega8',
        'transport_checks':transport,'ideal_tangency':tangent,
        'linear_hull_rank':hull['rank'],'lie_rank_through_depth_3':lie['rank'],
        'reduced_model_orbit_dimension':4,
        'nonzero_control_minor':'256/3 on (c4,c6J,c8_1,c8_product)',
        'orbit_constraints':['c_I6_2=0','c_I8_3=0','c_I8_4=0','c_I8_5=0','c_I8_6=0','c_I8_2=-16*c_I4_1^3'],
        'sufficiency_witness':construct_endpoint(2,Fraction(3,7),5,Fraction(-11,13)),
        'proof_scope':'finite degree-eight polynomial model, signed independent scalar controls, free seed',
        'interpretation_warning':'five-dimensional linear hull is not a five-dimensional orbit',
        'not_established':['priority/novelty','degree-ten or degree-twelve orbit classification',
                           'conformal/nonanalytic reachability','causality or unitary physical admissibility'],
        'external_inputs':['HLS (2.33),(3.3) and the stated graph-basis identities',
                           'independent tensor fitting checks the map on finite samples, not all fields']}


def verify_sampled(fields,ids,p):
    expected_ids,analytic=analytic_fields()
    if tuple(ids)!=expected_ids:
        raise ValueError('degree-eight coordinate ordering differs')
    if set(fields)!=set(analytic):
        raise ValueError('degree-eight stress generator catalogue differs')
    for name,f in analytic.items():
        expected=tuple(Poly(10,dict(x.terms),p) for x in f)
        if fields[name]!=expected:
            raise ValueError(f'fresh tensor fit contradicts the analytic degree-eight model: {name}')
    return True
