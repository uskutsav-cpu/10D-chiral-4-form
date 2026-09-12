from fractions import Fraction
import pytest
from chiral4form.polynomial import Poly
from chiral4form.basis_change import inverse,transform_fields
from chiral4form.reachability import lie_rank,linear_invariant_hull
from chiral4form.degree8_orbit import analytic_fields,constraints
from chiral4form.ideals import invariant_ideal_report

@pytest.mark.parametrize('p',[None,101,30011])
def test_field_covariance_preserves_orbit_not_just_hull(p):
    x=Poly.variable(2,0,p);one=Poly.constant(2,1,p)
    fields={'parabola':(one,x)}
    converted,_=transform_fields(fields,[[1,2],[0,1]])
    assert lie_rank(list(converted.values()))['rank']==1
    assert linear_invariant_hull(list(converted.values()))['rank']==2

@pytest.mark.parametrize('p',[None,101])
def test_singular_change_rejected(p):
    with pytest.raises(ValueError):inverse([[1,1],[1,1]],p)

def test_octic_obstruction_ideal_survives_homogeneous_change():
    ids,fields=analytic_fields();P=[[int(i==j) for j in range(10)] for i in range(10)]
    P[3][4]=2;P[6][9]=-3;P[1][2]=7
    transformed,subs=transform_fields(fields,P,[4,6,6]+[8]*7)
    qs=[q.substitute(subs) for q in constraints()]
    assert invariant_ideal_report(transformed,qs)['all_tangent']
    assert lie_rank(list(transformed.values()))['rank']==4

def test_grade_mixing_rejected():
    x=Poly.variable(2,0);y=Poly.variable(2,1)
    with pytest.raises(ValueError):transform_fields({'X':(x,y)},[[1,1],[0,1]],[4,6])
