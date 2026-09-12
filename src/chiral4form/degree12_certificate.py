"""Degree-12 modular/QQ certification orchestration.

The module is fail-closed: stable modular ranks are reported independently
from any characteristic-zero lift, and no ``qq_theorem`` flag is emitted unless
both equation/control row spaces lift and verify exactly.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Sequence

from .exact import nullspace_q, rref_q
from .finite_field import matrix_rank, nullspace, rref
from .fitting import fields_from_json
from .normalization_audit import audit_model_records, audit_registry_records, load_records
from .polynomial import Poly, derivation, field_at
from .subspace_lift import canonical_row_spaces, lift_stable_row_space

LOWER_LABELS = (
    "A^5", "A^3*B", "A^2*C", "A^2*D", "A*B^2",
    "A*U", "A*V", "B*C", "B*D",
)


def _weight100_basis(ids: Sequence[str], p: int):
    if len(ids) != 78 or tuple(ids[:6]) != ("A", "B", "C", "D", "U", "V"):
        raise ValueError("expected reduced 78-coordinate model ordered A,B,C,D,U,V,+72 degree-12 coordinates")
    n = 78
    A, B, C, D, U, V = [Poly.variable(n, i, p) for i in range(6)]
    Z = [Poly.variable(n, 6 + i, p) for i in range(72)]
    lower = [
        A**5,
        A**3 * B,
        A**2 * C,
        A**2 * D,
        A * B**2,
        A * U,
        A * V,
        B * C,
        B * D,
    ]
    return Z + lower, list(ids[6:]) + list(LOWER_LABELS)


def obstruction_equations(record: Mapping, p: int):
    ids, fields = fields_from_json(record)
    basis, labels = _weight100_basis(ids, p)
    equations: list[list[int]] = []
    equation_counts: dict[str, int] = {}
    for gname in sorted(fields):
        field = fields[gname]
        derivs = (
            [derivation(field, q) - 100 * q for q in basis]
            if gname == "tr1"
            else [derivation(field, q) for q in basis]
        )
        monomials = sorted(set().union(*(f.terms.keys() for f in derivs)))
        before = len(equations)
        for monomial in monomials:
            row = [int(f.terms.get(monomial, 0)) % p for f in derivs]
            if any(row):
                equations.append(row)
        equation_counts[gname] = len(equations) - before
    return tuple(ids), fields, labels, equations, equation_counts


def modular_certificate(records_by_prime: Mapping[int, Mapping]) -> dict:
    normalization = audit_model_records(records_by_prime)
    reports = []
    equation_rows: dict[int, list[list[int]]] = {}
    control_rows: dict[int, list[list[int]]] = {}
    eq_piv_ref = None
    control_piv_ref = None
    labels_ref = None
    for p, record in sorted(records_by_prime.items()):
        ids, fields, labels, equations, counts = obstruction_equations(record, int(p))
        if labels_ref is None:
            labels_ref = labels
        elif labels != labels_ref:
            raise ValueError(f"obstruction ansatz labels differ at prime {p}")
        eq_rr, eq_piv = rref(equations, int(p), ncols=81)
        eq_rank = len(eq_piv)
        if eq_rank != 13:
            raise ValueError(f"expected degree-12 equation rank 13 at prime {p}, got {eq_rank}")
        kernel = nullspace(equations, int(p), ncols=81)
        nullity = len(kernel)
        leading_rank = matrix_rank([row[:72] for row in kernel], int(p), ncols=72)
        seed = [0] * 78
        controls = [field_at(fields[name], seed) for name in sorted(fields)]
        c_rr, c_piv = rref(controls, int(p), ncols=78)
        control_rank = len(c_piv)
        upper = 78 - leading_rank
        if nullity != 68 or leading_rank != 68 or control_rank != 10 or upper != 10:
            raise ValueError(
                f"degree-12 modular signature changed at prime {p}: "
                f"nullity={nullity}, leading={leading_rank}, control={control_rank}, upper={upper}"
            )
        if eq_piv_ref is None:
            eq_piv_ref = tuple(eq_piv); control_piv_ref = tuple(c_piv)
        elif tuple(eq_piv) != eq_piv_ref or tuple(c_piv) != control_piv_ref:
            raise ValueError(f"pivot pattern changed at prime {p}")
        equation_rows[int(p)] = equations
        control_rows[int(p)] = controls
        reports.append({
            "prime": int(p),
            "equation_count": len(equations),
            "equation_rank": eq_rank,
            "obstruction_nullity": nullity,
            "degree12_leading_rank": leading_rank,
            "control_rank": control_rank,
            "upper_orbit_dimension": upper,
            "dimension_match": upper == control_rank,
            "equations_by_generator": counts,
        })
    return {
        "status": "stable_modular_degree12_certificate",
        "normalization": normalization,
        "labels": labels_ref,
        "reports": reports,
        "equation_rows_by_prime": equation_rows,
        "control_rows_by_prime": control_rows,
        "equation_pivots": list(eq_piv_ref or ()),
        "control_pivots": list(control_piv_ref or ()),
        "modular_orbit_dimension": 10,
        "qq_theorem": False,
    }


def certify_records(
    records_by_prime: Mapping[int, Mapping],
    *,
    holdout_primes: Sequence[int] = (),
    max_bad_primes: int = 0,
    registries: Sequence[Mapping] = (),
) -> dict:
    modular = modular_certificate(records_by_prime)
    equation_spaces, eq_pivots = canonical_row_spaces(
        modular["equation_rows_by_prime"], ncols=81, expected_rank=13
    )
    control_spaces, control_pivots = canonical_row_spaces(
        modular["control_rows_by_prime"], ncols=78, expected_rank=10
    )
    eq_lift = lift_stable_row_space(
        equation_spaces,
        eq_pivots,
        ncols=81,
        holdout_primes=holdout_primes,
        max_bad_primes=max_bad_primes,
    )
    control_lift = lift_stable_row_space(
        control_spaces,
        control_pivots,
        ncols=78,
        holdout_primes=holdout_primes,
        max_bad_primes=max_bad_primes,
    )
    result = {
        "schema": 1,
        "status": "degree12_certification_complete" if eq_lift.get("qq_verified") and control_lift.get("qq_verified") else "degree12_modular_certified_qq_pending",
        "registry_audit": audit_registry_records(registries),
        "modular": {k: v for k, v in modular.items() if not k.endswith("_rows_by_prime")},
        "equation_space_lift": eq_lift,
        "control_space_lift": control_lift,
        "qq_algebraic_candidate": False,
        "qq_theorem": False,
        "direct_characteristic_zero_physics_certificate": False,
    }
    if not (eq_lift.get("qq_verified") and control_lift.get("qq_verified")):
        return result
    eq_q = [[x for x in row] for row in eq_lift["rational_rows"]]
    eq_rr, eq_piv_q = rref_q(eq_q, ncols=81)
    if len(eq_piv_q) != 13:
        raise ValueError("lifted QQ equation rank is not 13")
    kernel_q = nullspace_q(eq_q, ncols=81)
    if len(kernel_q) != 68:
        raise ValueError("lifted QQ obstruction nullity is not 68")
    _, leading_piv_q = rref_q([row[:72] for row in kernel_q], ncols=72)
    control_q = [[x for x in row] for row in control_lift["rational_rows"]]
    _, control_piv_q = rref_q(control_q, ncols=78)
    if len(leading_piv_q) != 68 or len(control_piv_q) != 10:
        raise ValueError("lifted QQ dimension sandwich does not close")
    result.update({
        "qq_algebraic_candidate": True,
        "qq_equation_rank": 13,
        "qq_obstruction_nullity": 68,
        "qq_degree12_leading_rank": 68,
        "qq_control_rank": 10,
        "qq_orbit_dimension_candidate": 10,
        # This remains false by design: a modularly lifted algebraic candidate
        # is not automatically a direct characteristic-zero tensor proof.
        "qq_theorem": False,
    })
    return result


def certify_directory(
    models_dir: Path,
    *,
    holdout_prime: int | None = None,
    max_bad_primes: int = 0,
) -> dict:
    paths = sorted(Path(models_dir).glob("fields_prime*.json"))
    if len(paths) < 2:
        raise ValueError("at least two fields_prime*.json models required")
    records = load_records(paths)
    holdout = () if holdout_prime is None else (int(holdout_prime),)
    return certify_records(records, holdout_primes=holdout, max_bad_primes=max_bad_primes)
