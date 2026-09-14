from chiral4form.degree8_orbit import analytic_fields,constraints
from chiral4form.obstruction_lifting import lift_finite_jet_ideal,relative_transport_report
from chiral4form.polynomial import Poly

def test_k6_transport_and_lift():
    ids,fields=analytic_fields()
    q6=Poly.variable(len(ids),2)
    r=relative_transport_report(fields,q6,{n:(40 if n=="tr1" else 0) for n in fields})
    assert r["all_verified"]
    c=lift_finite_jet_ideal(fields,(q6,),cutoff=6,label="<K6>")
    assert c.all_orders_cylinder_invariant

def test_degree8_full_ideal_lifts():
    _,fields=analytic_fields()
    c=lift_finite_jet_ideal(fields,constraints(),cutoff=8,label="degree8")
    assert c.all_orders_cylinder_invariant

def test_omega8_transport():
    _,fields=analytic_fields()
    omega=constraints()[-1]
    r=relative_transport_report(fields,omega,{n:(60 if n=="tr1" else 0) for n in fields})
    assert r["all_verified"]
