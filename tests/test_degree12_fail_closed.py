from pathlib import Path


def test_qq_theorem_is_not_promoted_by_design():
    path = Path(__file__).resolve().parents[1] / "src" / "chiral4form" / "degree12_certificate.py"
    source = path.read_text()
    assert '"qq_theorem": False' in source
    assert '"direct_characteristic_zero_physics_certificate": False' in source
