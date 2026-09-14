from chiral4form.registry import Registry
from chiral4form.registry_resolution import registry_has_cutoff

def test_low_fixture_is_intentionally_not_degree10():
    reg=Registry.load("data/fixtures/low_degree_registry.json")
    assert registry_has_cutoff(reg,8)
    assert not registry_has_cutoff(reg,10)
