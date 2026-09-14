"""Exactness gates for promoting fitted finite models to physics theorems."""
from __future__ import annotations

import json
from pathlib import Path


UNCONDITIONAL_STATUSES = {
    "symbolic_physics_basis_identity_proved",
    "exact_characteristic_zero_physics_map",
    "bounded_integer_residual_certificate",
}


def model_exactness_report(path):
    path = Path(path)
    record = json.loads(path.read_text())
    status = record.get("coefficient_status", "unknown")
    exact = status in UNCONDITIONAL_STATUSES

    return {
        "schema": 1,
        "model": str(path),
        "coefficient_status": status,
        "unconditional": exact,
        "acceptable_unconditional_statuses": sorted(UNCONDITIONAL_STATUSES),
        "remaining_gate": (
            None
            if exact
            else (
                "Construct an unconditional characteristic-zero physics/basis-map "
                "certificate. Modular holdouts and bounded rational reconstruction "
                "alone do not prove the unknown physical coefficient height."
            )
        ),
    }
