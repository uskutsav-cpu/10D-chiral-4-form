"""Build the explicit triangular 68-generator degree-12 quotient ideal."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from .degree12_certificate import LOWER_LABELS
from .polynomial import Poly
from .provenance import atomic_json


def ansatz_polynomials(labels=None):
    # Preserve the earlier helper API for downstream tests/tools while adding
    # strict validation when the physical coordinate labels are supplied.
    if labels is None:
        labels = [f"Z{i+1}" for i in range(72)] + list(LOWER_LABELS)
    labels = list(labels)
    if len(labels) != 81:
        raise ValueError("degree-12 obstruction ansatz must have exactly 81 labels")
    if tuple(labels[-9:]) != tuple(LOWER_LABELS):
        raise ValueError("degree-12 lower-monomial label ordering changed")
    n = 78
    A, B, C, D, U, V = [Poly.variable(n, i) for i in range(6)]
    Z = [Poly.variable(n, 6 + i) for i in range(72)]
    lower = [A**5, A**3 * B, A**2 * C, A**2 * D, A * B**2, A * U, A * V, B * C, B * D]
    return Z + lower, labels


def build_ideal(
    kernel_path: str | Path = "verification/all_orders/degree12_exact/kernel_lift.json",
    output: str | Path = "verification/all_orders/degree12_exact/ideal.json",
):
    record = json.loads(Path(kernel_path).read_text())
    if record.get("status") != "exact_degree12_kernel_lift" or record.get("exact") is not True:
        raise ValueError("degree-12 exact QQ kernel certificate is missing")
    if record.get("exact_proof_method") not in {
        "direct_characteristic_zero_source_equation_reconstruction",
        "nullspace_of_all_36_direct_characteristic_zero_source_rows",
    }:
        raise ValueError("degree-12 kernel was not produced by an accepted direct characteristic-zero proof path")
    if record.get("kernel_dimension") != 68 or record.get("leading_rank") != 68:
        raise ValueError("degree-12 exact kernel signature is not 68/68")
    lift = record.get("lift", {})
    if lift.get("qq_verified") is not True or lift.get("all_theorem_primes_verified") is not True:
        raise ValueError("degree-12 kernel lacks exact/holdout verification")

    rows = [[Fraction(x) for x in row] for row in lift["rational_rows"]]
    if len(rows) != 68 or any(len(row) != 81 for row in rows):
        raise ValueError("exact degree-12 kernel matrix is not 68x81")
    basis, labels = ansatz_polynomials(record["ansatz_labels"])

    pivots = [int(x) for x in record["pivots"]]
    free = [int(x) for x in record["free_degree12_columns"]]
    if len(pivots) != 68 or len(set(pivots)) != 68 or any(not 0 <= p < 72 for p in pivots):
        raise ValueError("kernel does not have 68 distinct degree-12 pivot columns")
    if free != [c for c in range(72) if c not in pivots] or len(free) != 4:
        raise ValueError("kernel free-coordinate metadata is inconsistent")

    constraints = []
    substitutions = [Poly.variable(78, i) for i in range(78)]
    triangular = []
    for row_index, (row, pivot) in enumerate(zip(rows, pivots)):
        if row[pivot] != 1:
            raise ValueError(f"kernel RREF row {row_index} is not normalized at pivot {pivot}")
        if any(row[p] for p in pivots if p != pivot):
            raise ValueError(f"kernel RREF row {row_index} is not reduced in pivot columns")
        q = Poly(78)
        for coefficient, monomial in zip(row, basis):
            if coefficient:
                q = q + coefficient * monomial
        constraints.append(q)

        rhs = Poly(78)
        for column, coefficient in enumerate(row):
            if column != pivot and coefficient:
                rhs = rhs - coefficient * basis[column]
        substitutions[6 + pivot] = rhs
        triangular.append({
            "pivot_column": pivot,
            "pivot_label": labels[pivot],
            "rhs": rhs.to_json(),
        })

    residuals = [q.substitute(substitutions) for q in constraints]
    if any(residual.terms for residual in residuals):
        raise AssertionError("triangular normal form does not annihilate all 68 exact generators")

    payload = {
        "schema": 3,
        "status": "exact_characteristic_zero_degree12_triangular_ideal",
        "constraint_count": 68,
        "leading_rank": 68,
        "degree12_coordinate_count": 72,
        "fiber_dimension": 4,
        "pivots": pivots,
        "free_degree12_columns": free,
        "free_degree12_labels": [labels[i] for i in free],
        "ansatz_labels": labels,
        "constraints": [q.to_json() for q in constraints],
        "triangular_equations": triangular,
        "all_generators_reduce_to_zero": True,
        "scope": "relative obstruction ideal on the exact degree-10 pure-stress orbit quotient",
        "full_ideal_definition": "I12 := I10 + preimage(<Q1,...,Q68>)",
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(path, payload)
    return payload
