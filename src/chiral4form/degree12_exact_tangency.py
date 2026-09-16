"""Exact tangency of the degree-12 obstruction ideal under all 36 sources."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from .provenance import atomic_json

SOURCE = Path("verification/all_orders/degree12_exact/source_equation_space.json")
KERNEL = Path("verification/all_orders/degree12_exact/kernel_lift.json")
IDEAL = Path("verification/all_orders/degree12_exact/ideal.json")


def _dot(left, right) -> Fraction:
    if len(left) != len(right):
        raise ValueError("tangency row lengths differ")
    return sum((Fraction(a) * Fraction(b) for a, b in zip(left, right)), Fraction(0))


def certify_tangency(output: str | Path = "verification/all_orders/degree12_exact/tangency.json"):
    source = json.loads(SOURCE.read_text())
    kernel = json.loads(KERNEL.read_text())
    ideal = json.loads(IDEAL.read_text())
    if source.get("status") != "exact_characteristic_zero_degree12_source_equation_space":
        raise ValueError("exact characteristic-zero degree-12 source equations are missing")
    if source.get("equation_count") != 36 or source.get("equation_rank") != 13:
        raise ValueError("source equation certificate does not contain rank-13 all-36 rows")
    if source.get("all_36_rows_verified_mod_every_theorem_prime") is not True:
        raise ValueError("all-36 source-row theorem-prime gate is false")
    if kernel.get("exact") is not True or kernel.get("kernel_dimension") != 68:
        raise ValueError("exact 68-dimensional degree-12 kernel is missing")
    if ideal.get("status") != "exact_characteristic_zero_degree12_triangular_ideal":
        raise ValueError("exact 68-generator degree-12 ideal is missing")
    if ideal.get("all_generators_reduce_to_zero") is not True:
        raise ValueError("degree-12 triangular normal form did not verify")

    equations = [[Fraction(x) for x in row["row"]] for row in source["all_source_rows"]]
    constraints = [[Fraction(x) for x in row] for row in kernel["lift"]["rational_rows"]]
    if len(equations) != 36 or any(len(row) != 81 for row in equations):
        raise ValueError("all exact source equations are not 36x81")
    if len(constraints) != 68 or any(len(row) != 81 for row in constraints):
        raise ValueError("exact obstruction matrix is not 68x81")

    failures = []
    for i, constraint in enumerate(constraints):
        for j, equation in enumerate(equations):
            if _dot(equation, constraint):
                failures.append((i, j))
    if failures:
        raise AssertionError(f"exact degree-12 tangency residuals are nonzero: {failures[:10]}")

    payload = {
        "schema": 4,
        "status": "exact_characteristic_zero_degree12_ideal_tangency",
        "constraint_count": 68,
        "source_equation_count": 36,
        "independent_source_equation_count": 13,
        "all_36_exact_source_residuals_zero": True,
        "finite_quotient_ideal_tangent": True,
        "relative_full_ideal_tangent": True,
        "argument": (
            "All 36 reduced pure-stress source equations were reconstructed directly over QQ. "
            "The 68 constraints are the exact QQ nullspace of that rank-13 equation space, so "
            "every one of the 36 physical sources annihilates every relative constraint. "
            "Unconditional invariance of I10 therefore gives tangency of the full preimage I12."
        ),
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(path, payload)
    return payload
