from chiral4form.predecessor_import import GENERATOR_MAP

def test_predecessor_generator_map_complete_and_unique():
    assert len(GENERATOR_MAP) == 18
    assert len(set(GENERATOR_MAP.values())) == 18
    assert GENERATOR_MAP["tr_tau"] == "tr1"
    assert GENERATOR_MAP["tr_tau^2"] == "tr1*tr1"
    assert GENERATOR_MAP["tr_tau2*tr_tau3"] == "tr2*tr3"
    assert GENERATOR_MAP["tr_tau2^3"] == "tr2*tr2*tr2"
