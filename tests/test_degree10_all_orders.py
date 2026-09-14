from chiral4form.degree10_all_orders import (
    degree10_constraints,
    degree10_orbit_parameterization,
)

def test_degree10_constraint_count_and_graph():
    qs=degree10_constraints()
    assert len(qs)==18
    names,values=degree10_orbit_parameterization()
    assert names==("A","B","C","D","U","V")
    assert len(values)==24
    for q in qs:
        assert not q.substitute(values,n_out=6)
