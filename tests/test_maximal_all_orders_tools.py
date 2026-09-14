from chiral4form.nauty_multigraphs import collapse_incidence_graph
from chiral4form.higher_degree_basis import IncrementalRowBasis
from chiral4form.higher_obstruction import coupling_weight,enumerate_weight_monomials
from chiral4form.stabilization_theorem import assess_stabilization
from chiral4form.formal_orbit_equality import assess_orbit_equality

def test_incidence_collapse_toy_multigraph():
    # Two left vertices joined by five parallel edges.
    n=7;A=[[0]*n for _ in range(n)]
    for r in range(2,7):
        for i in (0,1):A[i][r]=A[r][i]=1
    M=collapse_incidence_graph(A,2)
    assert M==((0,5),(5,0))

def test_incremental_row_basis():
    B=IncrementalRowBasis(101,3)
    assert B.add([1,0,0])
    assert not B.add([2,0,0])
    assert B.add([0,1,0])
    assert B.rank==2

def test_weight_monomials_small():
    weights=[20,40]
    rows=enumerate_weight_monomials(weights,60,[0,1])
    assert (3,0) in rows
    assert (1,1) in rows

def test_stabilization_gate_never_promotes_finite_evidence(tmp_path):
    q14=tmp_path/"q14.json";q16=tmp_path/"q16.json"
    q14.write_text('{"stable_new_dimension":0}')
    q16.write_text('{"status":"degree16_modular_obstruction_probe"}')
    r=assess_stabilization(q14,q16,tmp_path/"missing.json")
    assert r["finite_degree_evidence"] is True
    assert r["finite_generation_proved"] is False

def test_orbit_equality_requires_sufficiency(tmp_path):
    s=tmp_path/"s.json";s.write_text('{"finite_generation_proved":true}')
    r=assess_orbit_equality(s,tmp_path/"missing.json")
    assert r["formal_orbit_equality_proved"] is False
