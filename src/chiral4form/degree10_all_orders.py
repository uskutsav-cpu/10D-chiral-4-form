"""Degree-10 obstruction ideal and conditional all-orders lift.

The finite polynomial vector field stored in verification/degree10/rational_model.json
is exact over QQ as a polynomial object once loaded, but the repository explicitly
labels its physics/basis-map status as a bounded rational candidate with modular
validation, not as an unconditional symbolic physics-map proof.

Accordingly this module:
1. reconstructs the exact finite degree-10 orbit ideal from the frozen obstruction
   certificate;
2. verifies exact QQ tangency of the loaded finite model to that ideal;
3. uses the already-proved filtration theorem to obtain a CONDITIONAL all-orders
   cylinder theorem;
4. refuses to relabel that result unconditional until an exact degree-10 physics-map
   certificate exists.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

from .degree8_orbit import constraints as degree8_constraints
from .exact import rational
from .fitting import fields_from_json
from .obstruction_lifting import extend_poly, lift_finite_jet_ideal
from .polynomial import Poly, derivation


DEGREE10_FREE_POSITIONS = (0, 12)  # I10_1 and I4_1*I6_1 in the 14 d10 slots.


def _load_json(path) -> dict:
    return json.loads(Path(path).read_text())


def load_degree10_artifacts(
    model_path="verification/degree10/rational_model.json",
    obstruction_path="verification/degree10/obstruction_certificate.json",
):
    model = _load_json(model_path)
    obstruction = _load_json(obstruction_path)
    ids, fields = fields_from_json(model)
    if tuple(ids) != tuple(obstruction["coordinate_ids"]):
        raise ValueError("degree-10 model/certificate coordinate ordering differs")
    if len(ids) != 24:
        raise ValueError("expected 24 cumulative coordinates through degree 10")
    if obstruction["obstruction_nullity"] != 12:
        raise ValueError("expected 12 independent degree-10 obstruction relations")
    return model, obstruction, ids, fields


def degree10_candidate_basis_polynomials(ids, obstruction):
    n = len(ids)
    c = [Poly.variable(n, i) for i in range(n)]
    index = {name: i for i, name in enumerate(ids)}

    lower_aliases = {
        "a^4": c[index["I4_1"]] ** 4,
        "a^2*b": c[index["I4_1"]] ** 2 * c[index["I6_1"]],
        "a*c": c[index["I4_1"]] * c[index["I8_1"]],
        "a*d": c[index["I4_1"]] * c[index["I4_1^2"]],
        "b^2": c[index["I6_1"]] ** 2,
    }

    out = []
    for label in obstruction["candidate_basis"]:
        if label in index:
            out.append(c[index[label]])
        elif label in lower_aliases:
            out.append(lower_aliases[label])
        else:
            raise ValueError(f"unknown degree-10 obstruction basis label {label!r}")
    if len(out) != 19:
        raise ValueError("expected 19-function obstruction ansatz")
    return tuple(out)


def degree10_constraints(
    model_path="verification/degree10/rational_model.json",
    obstruction_path="verification/degree10/obstruction_certificate.json",
):
    model, obstruction, ids, fields = load_degree10_artifacts(
        model_path, obstruction_path
    )
    n = len(ids)
    lower = tuple(extend_poly(q, n) for q in degree8_constraints())
    basis = degree10_candidate_basis_polynomials(ids, obstruction)

    q10 = []
    for row in obstruction["obstructions_integer_coordinates"]:
        if len(row) != len(basis):
            raise ValueError("obstruction row length mismatch")
        q = Poly(n)
        for coefficient, function in zip(row, basis):
            q = q + rational(coefficient) * function
        q10.append(q)

    if len(q10) != 12:
        raise ValueError("expected twelve degree-10 constraints")
    return tuple(lower + tuple(q10))


def degree10_orbit_parameterization(
    obstruction_path="verification/degree10/obstruction_certificate.json",
):
    """Return the exact six-parameter graph solving all 18 constraints.

    The six free variables are A,B,C,D,U,V:
      A=c_I4_1
      B=c_I6_1
      C=c_I8_1
      D=c_I4_1^2
      U=c_I10_1
      V=c_I4_1*I6_1
    """
    obstruction = _load_json(obstruction_path)
    rows = [
        [Fraction(x) for x in row]
        for row in obstruction["obstructions_integer_coordinates"]
    ]
    leading = [row[:14] for row in rows]
    tails = [row[14:] for row in rows]

    dep_positions = [
        j for j in range(14)
        if j not in DEGREE10_FREE_POSITIONS
    ]
    if len(dep_positions) != 12:
        raise AssertionError("expected twelve dependent degree-10 positions")

    M = sp.Matrix(
        [[sp.Rational(leading[i][j].numerator, leading[i][j].denominator)
          for j in dep_positions]
         for i in range(12)]
    )
    if M.det() == 0:
        raise ValueError("degree-10 obstruction solve matrix is singular")
    Minv = M.inv()

    nfree = 6
    A, B, C, D, U, V = [Poly.variable(nfree, i) for i in range(nfree)]
    lower_tail = [A**4, A**2 * B, A * C, A * D, B**2]

    rhs = []
    for i in range(12):
        q = (
            -leading[i][DEGREE10_FREE_POSITIONS[0]] * U
            -leading[i][DEGREE10_FREE_POSITIONS[1]] * V
        )
        for coefficient, monomial in zip(tails[i], lower_tail):
            q = q - coefficient * monomial
        rhs.append(q)

    dependent = []
    for j in range(12):
        q = Poly(nfree)
        for i in range(12):
            coeff = Fraction(str(Minv[j, i]))
            q = q + coeff * rhs[i]
        dependent.append(q)

    z = [None] * 14
    z[DEGREE10_FREE_POSITIONS[0]] = U
    z[DEGREE10_FREE_POSITIONS[1]] = V
    for pos, q in zip(dep_positions, dependent):
        z[pos] = q
    if any(x is None for x in z):
        raise AssertionError("incomplete degree-10 graph parameterization")

    zero = Poly(nfree)
    lower_values = [
        A,            # I4_1
        B,            # I6_1
        zero,         # I6_2
        C,            # I8_1
        -16 * A**3,   # I8_2
        zero, zero, zero, zero,  # I8_3..I8_6
        D,            # I4_1^2
    ]
    values = tuple(lower_values + z)
    if len(values) != 24:
        raise AssertionError("parameterization dimension mismatch")
    return ("A", "B", "C", "D", "U", "V"), values


def degree10_tangency_report(
    model_path="verification/degree10/rational_model.json",
    obstruction_path="verification/degree10/obstruction_certificate.json",
):
    model, obstruction, ids, fields = load_degree10_artifacts(
        model_path, obstruction_path
    )
    constraints = degree10_constraints(model_path, obstruction_path)
    free_names, substitution = degree10_orbit_parameterization(obstruction_path)

    constraint_zero = []
    for i, q in enumerate(constraints):
        reduced = q.substitute(substitution, n_out=6)
        constraint_zero.append(not bool(reduced))
        if reduced:
            raise AssertionError(f"constraint {i} does not vanish on graph")

    tangency = {}
    all_tangent = True
    for name, field in fields.items():
        rows = []
        for i, q in enumerate(constraints):
            residual = derivation(field, q).substitute(substitution, n_out=6)
            ok = not bool(residual)
            rows.append(ok)
            all_tangent &= ok
            if not ok:
                raise AssertionError(
                    f"degree-10 tangency failed: generator={name}, constraint={i}"
                )
        tangency[name] = rows

    map_status = model.get("coefficient_status", "unknown")
    exact_map = map_status in {
        "symbolic_physics_basis_identity_proved",
        "exact_characteristic_zero_physics_map",
        "bounded_integer_residual_certificate",
    }

    return {
        "schema": 1,
        "status": "exact_QQ_degree10_ideal_tangency_under_loaded_map",
        "coordinate_ids": list(ids),
        "free_parameters": list(free_names),
        "constraint_count": len(constraints),
        "all_constraints_vanish_on_graph": all(constraint_zero),
        "all_loaded_vector_fields_tangent": all_tangent,
        "per_generator_tangency": tangency,
        "physics_map_status": map_status,
        "physics_map_unconditional": exact_map,
        "conditional_all_orders_lift_available": all_tangent,
        "unconditional_all_orders_lift_available": all_tangent and exact_map,
        "claim_boundary": (
            "Exact algebra over QQ is proved for the loaded frozen degree-10 "
            "polynomial model. Promotion to an unconditional physical all-orders "
            "theorem requires an unconditional degree-10 physics/basis-map proof."
        ),
    }


def degree10_conditional_all_orders_certificate(
    model_path="verification/degree10/rational_model.json",
    obstruction_path="verification/degree10/obstruction_certificate.json",
):
    report = degree10_tangency_report(model_path, obstruction_path)
    model, obstruction, ids, fields = load_degree10_artifacts(
        model_path, obstruction_path
    )
    constraints = degree10_constraints(model_path, obstruction_path)

    lifted = lift_finite_jet_ideal(
        fields,
        constraints,
        cutoff=10,
        label="degree-10 pure-stress orbit ideal under frozen QQ model",
    )

    return {
        "schema": 1,
        "status": (
            "all_orders_degree10_ideal_theorem"
            if report["physics_map_unconditional"]
            else "conditional_all_orders_degree10_ideal_theorem"
        ),
        "finite_QQ_tangency": report,
        "finite_jet_lift": lifted.to_json(),
        "all_orders_cylinder_invariant_under_loaded_map": True,
        "unconditional_physical_theorem": report["physics_map_unconditional"],
        "remaining_gate": (
            None
            if report["physics_map_unconditional"]
            else "prove the degree-10 physics/basis map by a bounded-residual or symbolic identity certificate"
        ),
    }
