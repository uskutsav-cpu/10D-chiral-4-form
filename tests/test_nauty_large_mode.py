from chiral4form.nauty_multigraphs import (
    _construct_known_incidence,
    collapse_incidence_graph,
    genbg_command,
    local_incidence_smoke_test,
)

def test_constructed_degree14_incidence_is_valid():
    A=_construct_known_incidence(14)
    assert len(A)==49
    M=collapse_incidence_graph(A,14)
    assert len(M)==14
    assert all(sum(row)==5 for row in M)

def test_degree16_incidence_size():
    A=_construct_known_incidence(16)
    assert len(A)==56

def test_large_commands_prefer_genbgL_when_available():
    cmd=genbg_command(14)
    assert cmd[0].endswith("genbgL")
    assert cmd[-3:]==["14","35","70:70"]

def test_local_smoke():
    assert local_incidence_smoke_test(14)["passed"]
