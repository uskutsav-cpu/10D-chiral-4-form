import pytest

from chiral4form.normalization_audit import audit_model_records


def model(p, ids=("A", "B")):
    z = {"nvars": len(ids), "prime": p, "terms": []}
    return {
        "coordinate_ids": list(ids),
        "restriction": "same",
        "fields": {"tr1": [dict(z) for _ in ids]},
    }


def test_normalization_audit_passes_fixed_structure():
    out = audit_model_records({101: model(101), 103: model(103)})
    assert out["status"] == "normalization_audit_passed"


def test_normalization_audit_rejects_coordinate_drift():
    with pytest.raises(ValueError):
        audit_model_records({101: model(101), 103: model(103, ("B", "A"))})
