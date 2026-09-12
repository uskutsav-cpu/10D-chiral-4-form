from fractions import Fraction
import pytest
from chiral4form.polynomial import Poly,bracket,derivation
from chiral4form.reachability import linear_invariant_hull,lie_rank,trajectory_jets,completion_search
from chiral4form.ideals import invariant_ideal_report
from chiral4form.sextic import sextic_report
from chiral4form.degree8_orbit import analytic_fields,constraints,construct_endpoint,report
from chiral4form.localization import QuadraticElement,modmax_symbolic_report
import sympy as sp

@pytest.mark.parametrize("p", [None,101])
def test_polynomial_exact_arithmetic_and_serialization(p):
    x,y=Poly.variable(2,0,p),Poly.variable(2,1,p)
    f=(x+y)**3
    assert f.derivative(0)==3*(x+y)**2
    assert Poly.from_json(f.to_json())==f
    assert f.evaluate([2,3])==(125 if p is None else 24)
    assert f.substitute([y,x])==f
    assert (x+y)**0==1
    assert (x+y).truncate(0)==0
    assert not (x-x)

@pytest.mark.parametrize("record", [
    {'nvars':1,'terms':[{'powers':[-1],'coefficient':'2'}]},
    {'nvars':1,'terms':[{'powers':[1],'coefficient':'2'},{'powers':[1],'coefficient':'3'}]},
    {'nvars':1,'prime':3,'terms':[{'powers':[1],'coefficient':'1/3'}]}])
def test_reject_malformed_polynomial(record):
    with pytest.raises((ValueError,TypeError)):Poly.from_json(record)

@pytest.mark.parametrize("p", [None,101])
def test_linear_hull_is_not_orbit(p):
    x=Poly.variable(2,0,p);one=Poly.constant(2,1,p)
    field=(one,x)
    assert linear_invariant_hull([field])['rank']==2
    assert lie_rank([field],max_depth=4)['rank']==1
    jets=trajectory_jets(field,3)
    assert jets['coefficients'][1]==['1','0']
    assert jets['coefficients'][2][1]==('1/2' if p is None else '51')

@pytest.mark.parametrize("p", [None,101])
def test_zero_at_seed_field_can_generate_lie_direction(p):
    x=Poly.variable(2,0,p);z=Poly(2,p=p);one=Poly.constant(2,1,p)
    X=(one,z);Y=(z,x)
    assert lie_rank([X,Y],max_depth=1)['rank']==1
    assert lie_rank([X,Y],max_depth=2)['rank']==2
    assert bracket(X,Y)==(z,one)

@pytest.mark.parametrize("depth", [0,-1])
def test_invalid_lie_depth_rejected(depth):
    with pytest.raises(ValueError):lie_rank([(Poly.constant(1,1),)],max_depth=depth)

def test_non_tangent_ideal_rejected():
    x=Poly.variable(1,0);one=Poly.constant(1,1)
    assert not invariant_ideal_report({'translation':(one,)},[x])['all_tangent']
    assert invariant_ideal_report({'scaling':(x,)},[x])['all_tangent']

@pytest.mark.parametrize("a,b,c,d", [(0,0,0,0),(2,3,5,7),(-2,Fraction(3,7),11,Fraction(-2,9))])
def test_degree8_constructive_orbit(a,b,c,d):
    result=construct_endpoint(a,b,c,d)
    endpoint=list(map(Fraction,result['endpoint']))
    assert [endpoint[i] for i in (0,1,3,9)]==list(map(Fraction,(a,b,c,d)))
    assert all(q.evaluate(endpoint)==0 for q in constraints())

def test_exact_degree8_transport_and_ranks():
    ids,fields=analytic_fields();omega=constraints()[-1]
    for name,field in fields.items():
        assert derivation(field,omega)==(60*omega if name=='tr1' else Poly(10))
    r=report()
    assert r['reduced_model_orbit_dimension']==4
    assert r['linear_hull_rank']==5
    assert r['ideal_tangency']['all_tangent']

def test_sextic_conditional_report():
    r=sextic_report()
    assert '40' in str(r) and 'hypotheses' in str(r).lower()

def test_quadratic_extension_and_modmax():
    I=sp.Symbol('I4',nonzero=True)
    r=QuadraticElement(0,1,I)
    assert (r*r).a==I and (r*r).b==0
    assert (r*r.inverse()).a==1
    assert sp.simplify(r.derivative(I).b-1/(2*I))==0
    result=modmax_symbolic_report()
    assert 'UNRESOLVED' in str(result)

def test_finite_catalogue_completion_scope():
    z=Poly(2);one=Poly.constant(2,1)
    r=completion_search([(one,z)],{'y':(z,one)},target_rank=2,metric='lie_rank')
    assert r['cardinality']==1 and not r['global_minimality']

def test_exact_catalogue_completion_requires_nonlinear_octic_direction():
    from chiral4form.degree8_completion import report
    r=report()
    assert r['necessary_catalogue_cardinality']==6
    assert r['local_orbit_dimension']==10
    assert r['five_extra_model_lie_rank']==9
    assert Fraction(r['full_rank_control_minor'])!=0
    assert all(x['ideal_tangency']['all_tangent'] for x in r['removal_obstructions'])
