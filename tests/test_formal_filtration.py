from chiral4form.formal_filtration import (
    bilinear_gradient_degree,derivative_degree,prove_truncation_commutes,
    stress_factor_leading_degree,stress_monomial_leading_degree)
from chiral4form.weighted import stress_generators

def test_basic_degree_rules():
    assert derivative_degree(4)==3
    assert derivative_degree(12)==11
    assert bilinear_gradient_degree(4,4)==6
    assert bilinear_gradient_degree(6,8)==12

def test_no_higher_degree_downfeed():
    for cutoff in range(4,22,2):
        assert prove_truncation_commutes(cutoff).all_checks

def test_catalogue_weights_match():
    for cutoff in range(4,22,2):
        for g in stress_generators(cutoff):
            assert stress_monomial_leading_degree(g.factors)==g.leading_degree

def test_trace_weights():
    assert stress_factor_leading_degree("tr1")==4
    assert stress_factor_leading_degree("tr2")==4
    assert stress_factor_leading_degree("tr3")==6
    assert stress_factor_leading_degree("tr6")==12
