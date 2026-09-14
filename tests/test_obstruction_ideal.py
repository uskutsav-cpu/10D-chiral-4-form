from chiral4form.obstruction_ideal import ideal_membership_report
from chiral4form.polynomial import Poly

def test_membership_detects_new_generator():
    n=2;x=Poly.variable(n,0);y=Poly.variable(n,1)
    r=ideal_membership_report((x,),(x*y,y))
    assert r["new_generator_count"]==1
    assert r["remainders"][0]["is_new_mod_lower_ideal"] is False
    assert r["remainders"][1]["is_new_mod_lower_ideal"] is True
