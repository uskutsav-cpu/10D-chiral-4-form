"""Exact completion inside a declared finite homogeneous generator catalogue.

The six extras below have degree 6 or 8. Consequently no mixed product with
any stress scalar can contribute through field degree eight. Adding them to
f(tau,S) therefore gives genuine constant coefficient-space vector fields.
This is NOT a minimum over all imaginable F5 structures or inhomogeneous S.
"""
from fractions import Fraction
from itertools import combinations
import sympy as sp
from .polynomial import Poly,derivation
from .degree8_orbit import IDS,analytic_fields,constraints
from .ideals import invariant_ideal_report
from .reachability import lie_rank

CANDIDATES=('I6_2','I8_2','I8_3','I8_4','I8_5','I8_6')


def fields_with_extras(extras=CANDIDATES):
    if len(set(extras))!=len(extras) or set(extras)-set(CANDIDATES):
        raise ValueError('use a distinct subset of the declared degree-eight catalogue')
    ids,fields=analytic_fields()
    for name in extras:
        fields['S:'+name]=tuple(Poly.constant(10,1 if i==IDS.index(name) else 0) for i in range(10))
    return ids,fields


def report():
    ids,fields=fields_with_extras()
    # This determinant is a polynomial identity at arbitrary coupling values.
    c=sp.symbols('c0:10')
    from .ideals import sympy_expr
    chosen=['tr2','tr3','tr4','tr2*tr2']+['S:'+name for name in CANDIDATES]
    matrix=sp.Matrix([[sympy_expr(f,c) for f in fields[name]] for name in chosen])
    det=sp.factor(matrix.det())
    if det==0:
        raise AssertionError('generalized controls fail the full-rank test')
    qs=constraints()
    by_name={'I6_2':qs[0],'I8_2':qs[-1],
             'I8_3':qs[1],'I8_4':qs[2],'I8_5':qs[3],'I8_6':qs[4]}
    removal=[]
    for omitted in CANDIDATES:
        _,reduced=fields_with_extras(tuple(x for x in CANDIDATES if x!=omitted))
        q=by_name[omitted]
        tangency=invariant_ideal_report(reduced,[q])
        if not tangency['all_tangent']:
            raise AssertionError('claimed catalogue necessity is false')
        removal.append({'omitted':omitted,'surviving_obstruction':q.to_json(),
                        'ideal_tangency':tangency,'nonzero_obstruction':bool(q)})
    # With the five previously proposed extras, the nonlinear Omega8 survives.
    _,old=fields_with_extras(tuple(x for x in CANDIDATES if x!='I8_2'))
    old_rank=lie_rank(list(old.values()),max_depth=3)['rank']
    return {'status':'exact_QQ_completion_in_declared_six_element_catalogue_conditional_on_degree8_map',
        'scope':'degree-eight finite coefficient model; homogeneous extras and signed controls; free seed',
        'candidate_catalogue':list(CANDIDATES),'necessary_catalogue_cardinality':6,
        'full_rank_control_minor':str(det),'local_orbit_dimension':10,
        'removal_obstructions':removal,'five_extra_model_lie_rank':old_rank,
        'five_extra_model_missing_constraint':'c_I8_2 + 16*c_I4_1^3 = 0',
        'global_minimum_over_all_possible_F5_structures':'NOT_CLAIMED',
        'why_no_mixed_generators_missing':'min extra degree 6 plus min stress degree 4 exceeds cutoff 8',
        'not_established':['degree-ten/twelve completion','minimality outside this catalogue','nonanalytic or physically admissible completion']}
