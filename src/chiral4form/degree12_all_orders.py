"""Degree-12 all-orders promotion gate.

The current repository proves a characteristic-zero local/global finite theorem
through degree 12, but it does not store an explicit 68-generator exact QQ
orbit ideal with a finite-model ideal-tangency certificate.

This module prevents accidental overclaiming: the finite theorem is accepted,
but the all-orders ideal lift stays blocked until explicit polynomial ideal
generators and exact finite tangency are provided.
"""
from __future__ import annotations

import json
from pathlib import Path


def _load(path):
    return json.loads(Path(path).read_text())


def degree12_lift_readiness(
    local_path="verification/degree12/qq/degree12_local_theorem_certificate.json",
    global_path="verification/degree12/global_reachability_manifest.json",
):
    local = _load(local_path)
    global_ = _load(global_path)

    finite_local = (
        local.get("status")
        == "characteristic_zero_degree12_local_orbit_theorem_certificate"
        and local.get("qq_theorem") is True
    )
    finite_global = (
        global_.get("status")
        == "frozen_degree12_global_constructive_reachability"
        and global_["theorem_summary"].get(
            "degree12_global_constructive_reachability"
        ) is True
    )

    explicit_ideal = local.get("explicit_orbit_ideal_QQ")
    exact_tangency = local.get("exact_orbit_ideal_tangency")

    ready = bool(
        finite_local
        and finite_global
        and explicit_ideal
        and exact_tangency
    )

    return {
        "schema": 1,
        "status": (
            "degree12_all_orders_lift_ready"
            if ready
            else "blocked_missing_explicit_degree12_ideal_tangency"
        ),
        "finite_characteristic_zero_local_theorem": finite_local,
        "finite_global_reachability_theorem": finite_global,
        "explicit_exact_QQ_orbit_ideal_present": bool(explicit_ideal),
        "exact_finite_ideal_tangency_present": bool(exact_tangency),
        "all_orders_lift_ready": ready,
        "required_next_input": (
            None
            if ready
            else (
                "A compact exact QQ generating set for the degree-12 orbit ideal "
                "and an exact finite-model proof X_g(I_12) subset I_12. "
                "The rank/nullity certificate alone is not an ideal-tangency proof."
            )
        ),
        "do_not_use_as_ideal_generators": [
            "the 13 source-equation basis rows from the PDE classification",
            "RREF rows from modular reconstruction",
            "raw annihilator vectors without orbit-ideal interpretation",
        ],
    }
