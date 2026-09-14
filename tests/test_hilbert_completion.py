from chiral4form.hilbert_completion import abstract_hilbert_completion_theorem
from chiral4form.stabilization import stabilization_evidence,finite_generation_gate

def test_abstract_completion_is_existence_only():
    r=abstract_hilbert_completion_theorem()
    assert r["existence_only"] is True
    assert "finite" in r["conclusion"].lower()

def test_finite_probes_do_not_overclaim():
    r=stabilization_evidence([
        {"degree":14,"new_generator_count":0},
        {"degree":16,"new_generator_count":0},
    ])
    assert r["longest_consecutive_zero_new_generator_run"]==2
    assert r["finite_generation_proved"] is False
    gate=finite_generation_gate(induction_proof_supplied=False)
    assert gate["finite_generation_proved"] is False
