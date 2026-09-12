from pathlib import Path

from chiral4form import load_baseline

ROOT = Path(__file__).resolve().parents[1]


def test_imported_baseline_is_internally_consistent():
    baseline = load_baseline(ROOT / "data" / "baseline" / "paper2_baseline.json")
    assert [r.degree for r in baseline.records] == [4, 6, 8, 10, 12]
    assert [r.dynamic_reachable_dim for r in baseline.records] == [1, 1, 3, 11, 67]
    assert [r.quotient_dim for r in baseline.records] == [0, 1, 4, 3, 5]
    assert "imported" in baseline.status
