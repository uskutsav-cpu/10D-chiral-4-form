"""Exact all-36 characteristic-zero degree-12 source equation reconstruction.

The leading 72 entries of every source row come from the direct physical
post-orbit target certificate.  The nine lower-monomial entries are derived
exactly over QQ from the unconditional degree-10 quotient fields.  Every one
of the 36 resulting 81-vectors is then reduced modulo every frozen theorem
prime and required to match the independently constructed reduced models.

This removes reliance on the superseded raw-target Step1C height shortcut.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from .degree12_exact_common import load_theorem_equations, step1c_record
from .exact import nullspace_q, rref_q
from .fitting import fields_from_json
from .modular_reconstruction import reduce_fraction
from .polynomial import Poly, derivation
from .provenance import atomic_json

DEG10_MODEL = Path("verification/degree10/rational_model.json")
DEG10_OBS = Path("verification/degree10/obstruction_certificate.json")
DEG10_UNCOND = Path("verification/all_orders/degree10_unconditional.json")
TARGET_CERT = Path("verification/all_orders/degree12_exact/exact_source_targets.json")
VAR_INDEX = {"A": 0, "B": 1, "C": 2, "D": 3, "U": 4, "V": 5}
ORBIT_COMPONENTS = (0, 1, 3, 9, 10, 22)


def _inverse_q(matrix):
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        raise ValueError("square rational matrix required")
    augmented = [
        [Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    rr, pivots = rref_q(augmented, ncols=2 * n)
    if pivots[:n] != list(range(n)):
        raise ValueError("degree-10 orbit leading matrix is singular over QQ")
    return [row[n:] for row in rr[:n]]


def _exact_degree10_orbit_substitution(ids: tuple[str, ...]) -> dict[str, Poly]:
    obs = json.loads(DEG10_OBS.read_text())
    if tuple(obs.get("coordinate_ids", ())) != ids:
        raise ValueError("degree-10 rational model and obstruction coordinate orders differ")
    rows = [[Fraction(int(x)) for x in row] for row in obs["obstructions_integer_coordinates"]]
    if len(rows) != 12 or any(len(row) != 19 for row in rows):
        raise ValueError("malformed degree-10 obstruction certificate")
    leading = [row[:14] for row in rows]
    tails = [row[14:] for row in rows]
    n = 6
    A, B, C, D, U, V = [Poly.variable(n, i) for i in range(n)]
    zero = Poly(n)
    lower = [A**4, A**2 * B, A * C, A * D, B**2]
    free = {0: U, 12: V}
    dependent_positions = [j for j in range(14) if j not in free]
    inverse = _inverse_q([[row[j] for j in dependent_positions] for row in leading])

    rhs = []
    for i in range(12):
        value = zero - leading[i][0] * U - leading[i][12] * V
        for coefficient, monomial in zip(tails[i], lower):
            value = value - coefficient * monomial
        rhs.append(value)
    dependent = []
    for j in range(12):
        value = zero
        for i in range(12):
            value = value + inverse[j][i] * rhs[i]
        dependent.append(value)
    z10 = [None] * 14
    z10[0], z10[12] = U, V
    for pos, value in zip(dependent_positions, dependent):
        z10[pos] = value
    if any(value is None for value in z10):
        raise AssertionError("exact degree-10 orbit substitution is incomplete")

    substitution: dict[str, Poly] = {
        "I4_1": A,
        "I6_1": B,
        "I6_2": zero,
        "I8_1": C,
        "I8_2": -16 * A**3,
        "I8_3": zero,
        "I8_4": zero,
        "I8_5": zero,
        "I8_6": zero,
        "I4_1^2": D,
    }
    for name, value in zip(ids[10:24], z10):
        substitution[name] = value
    if set(substitution) != set(ids):
        raise AssertionError("exact degree-10 substitution does not cover all 24 coordinates")
    return substitution


def _exact_reduced_degree10_fields() -> dict[str, tuple[Poly, ...]]:
    uncond = json.loads(DEG10_UNCOND.read_text())
    if uncond.get("status") != "unconditional_all_orders_degree10_ideal_theorem":
        raise ValueError("unconditional degree-10 theorem certificate is missing")
    model = json.loads(DEG10_MODEL.read_text())
    ids, fields = fields_from_json(model)
    ids = tuple(ids)
    if len(ids) != 24:
        raise ValueError("frozen degree-10 rational model does not have 24 coordinates")
    sample = next(iter(next(iter(fields.values()))))
    if sample.p is not None:
        raise ValueError("frozen degree-10 rational model is not over QQ")
    substitution = _exact_degree10_orbit_substitution(ids)
    values = [substitution[name] for name in ids]
    reduced = {}
    for name, field in fields.items():
        full = tuple(component.substitute(values, n_out=6) for component in field)
        if len(full) != 24:
            raise AssertionError(f"{name}: substituted degree-10 field has {len(full)} components, expected 24")
        reduced[name] = tuple(full[i] for i in ORBIT_COMPONENTS)
        if len(reduced[name]) != 6 or any(q.n != 6 for q in reduced[name]):
            raise AssertionError(f"{name}: reduced degree-10 quotient field is not 6D")
    return reduced


def _monomial_exp(meta: dict) -> tuple[int, ...]:
    exponent = [0] * 6
    for factor in meta.get("output_monomial", ()):
        variable = str(factor["variable"])
        if variable not in VAR_INDEX:
            raise ValueError(f"unknown reduced output variable {variable}")
        exponent[VAR_INDEX[variable]] += int(factor["power"])
    return tuple(exponent)


def _all_descriptors(step1c: dict) -> tuple[dict, ...]:
    rows = []
    for meta in step1c["equation"]["basis"]:
        rows.append({**meta, "source_row": int(meta["source_row"]), "kind": "basis"})
    for meta in step1c["equation"]["dependencies"]:
        rows.append({
            "source_row": int(meta["source_row"]),
            "generator": meta["generator"],
            "output_monomial": meta.get("output_monomial", []),
            "output_monomial_text": meta["output_monomial_text"],
            "kind": "dependent",
        })
    rows.sort(key=lambda x: x["source_row"])
    if [x["source_row"] for x in rows] != list(range(36)):
        raise ValueError("Step1C metadata does not label exactly source rows 0..35")
    return tuple(rows)


def _tail_row(meta: dict, fields: dict[str, tuple[Poly, ...]]) -> list[Fraction]:
    n = 6
    A, B, C, D, U, V = [Poly.variable(n, i) for i in range(n)]
    lower = [A**5, A**3 * B, A**2 * C, A**2 * D, A * B**2, A * U, A * V, B * C, B * D]
    generator = str(meta["generator"])
    exponent = _monomial_exp(meta)
    # Generators whose leading field degree is 12 preserve the degree<=10 base,
    # so their action on every lower monomial is identically zero.
    if generator not in fields:
        return [Fraction(0)] * 9
    row = []
    for q in lower:
        d = derivation(fields[generator], q)
        if generator == "tr1":
            d = d - 100 * q
        row.append(Fraction(d.terms.get(exponent, 0)))
    return row


def _exact_linear_combination(rows: dict[int, list[Fraction]], dependency: dict) -> list[Fraction]:
    out = [Fraction(0)] * 81
    for term in dependency["terms"]:
        coefficient = Fraction(term["coefficient"])
        basis_id = str(term["basis_id"])
        source_row = next(
            int(meta["source_row"])
            for meta in step1c_record()["equation"]["basis"]
            if str(meta["basis_id"]) == basis_id
        )
        out = [a + coefficient * b for a, b in zip(out, rows[source_row])]
    return out


def reconstruct_source_rows(
    output: str | Path = "verification/all_orders/degree12_exact/source_equation_space.json",
    kernel_output: str | Path = "verification/all_orders/degree12_exact/kernel_lift.json",
):
    step1c = step1c_record()
    targets = json.loads(TARGET_CERT.read_text())
    if targets.get("status") != "exact_characteristic_zero_all_degree12_source_targets":
        raise ValueError("exact all-36 characteristic-zero source-target certificate is missing")
    if targets.get("all_36_direct_physical_crosschecks") is not True or targets.get("all_36_theorem_prime_crosschecks") is not True:
        raise ValueError("all-36 source-target certificate did not pass its physical/theorem-prime gates")

    target_rows = {int(x["source_row"]): x for x in targets.get("rows", ())}
    if sorted(target_rows) != list(range(36)):
        raise ValueError("source-target certificate does not contain exact rows 0..35")

    model_dir, primes, equations, labels = load_theorem_equations()
    fields = _exact_reduced_degree10_fields()
    descriptors = _all_descriptors(step1c)
    full_rows: dict[int, list[Fraction]] = {}
    for meta in descriptors:
        source_row = int(meta["source_row"])
        leading = [Fraction(x) for x in target_rows[source_row]["coordinates"]]
        if len(leading) != 72:
            raise ValueError(f"source target row {source_row} is not 72-dimensional")
        row = leading + _tail_row(meta, fields)
        if len(row) != 81:
            raise AssertionError(f"exact source row {source_row} has length {len(row)}, expected 81")
        for p in primes:
            got = [reduce_fraction(q, p) for q in row]
            want = [int(x) % p for x in equations[p][source_row]]
            if got != want:
                raise AssertionError(f"exact source row {source_row} disagrees with theorem model at prime {p}")
        full_rows[source_row] = row

    ordered_all = [full_rows[i] for i in range(36)]
    all_rref, all_pivots = rref_q(ordered_all, ncols=81)
    if len(all_pivots) != 13:
        raise AssertionError(f"exact all-source row rank is {len(all_pivots)}, expected 13")

    basis_meta = step1c["equation"]["basis"]
    basis_rows = [full_rows[int(meta["source_row"])] for meta in basis_meta]
    _, basis_pivots = rref_q(basis_rows, ncols=81)
    if len(basis_pivots) != 13:
        raise AssertionError(f"corrected 13-row basis has rank {len(basis_pivots)}, expected 13")

    # Now the legacy dependency formulas are checked as exact QQ identities
    # against directly reconstructed characteristic-zero rows; no height argument
    # is needed for them anymore.
    dependency_rows_verified = []
    basis_id_to_row = {str(meta["basis_id"]): int(meta["source_row"]) for meta in basis_meta}
    for dependency in step1c["equation"]["dependencies"]:
        source_row = int(dependency["source_row"])
        combination = [Fraction(0)] * 81
        for term in dependency["terms"]:
            coefficient = Fraction(term["coefficient"])
            source = full_rows[basis_id_to_row[str(term["basis_id"])]]
            combination = [a + coefficient * b for a, b in zip(combination, source)]
        if combination != full_rows[source_row]:
            raise AssertionError(f"Step1C dependency row {source_row} is not an exact QQ identity")
        dependency_rows_verified.append(source_row)

    kernel = nullspace_q(ordered_all, ncols=81)
    if len(kernel) != 68:
        raise AssertionError(f"exact obstruction nullity is {len(kernel)}, expected 68")
    _, leading_pivots = rref_q([row[:72] for row in kernel], ncols=72)
    if len(leading_pivots) != 68:
        raise AssertionError(f"exact leading obstruction rank is {len(leading_pivots)}, expected 68")
    kernel_rref, kernel_pivots = rref_q(kernel, ncols=81)
    if len(kernel_pivots) != 68 or any(column >= 72 for column in kernel_pivots):
        raise AssertionError("exact kernel does not have 68 pivots in the degree-12 coordinate block")
    free = [column for column in range(72) if column not in kernel_pivots]
    if len(free) != 4:
        raise AssertionError(f"exact degree-12 fiber has {len(free)} free coordinates, expected 4")

    source_payload = {
        "schema": 4,
        "status": "exact_characteristic_zero_degree12_source_equation_space",
        "proof_method": "all_36_direct_physical_characteristic_zero_target_rows_plus_exact_QQ_degree10_tails",
        "model_directory": str(model_dir),
        "equation_count": 36,
        "equation_rank": 13,
        "independent_source_rows": [int(meta["source_row"]) for meta in basis_meta],
        "dependency_rows_verified_exactly_over_QQ": dependency_rows_verified,
        "theorem_prime_reduction_checks": list(primes),
        "all_source_rows": [
            {
                **meta,
                "row": [str(q) for q in full_rows[int(meta["source_row"])]],
                "leading_target_certificate": str(TARGET_CERT),
                "tail_method": "direct_QQ_derivation_from_unconditional_degree10_quotient_fields",
            }
            for meta in descriptors
        ],
        "source_basis": [
            {
                **meta,
                "row": [str(q) for q in full_rows[int(meta["source_row"])]],
            }
            for meta in basis_meta
        ],
        "all_36_rows_verified_mod_every_theorem_prime": True,
        "all_23_step1c_dependencies_verified_exactly_over_QQ": True,
    }
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(out, source_payload)

    kernel_payload = {
        "schema": 4,
        "status": "exact_degree12_kernel_lift",
        "exact": True,
        "exact_proof_method": "nullspace_of_all_36_direct_characteristic_zero_source_rows",
        "source_equation_certificate": str(out),
        "ansatz_labels": list(labels),
        "ansatz_dimension": 81,
        "kernel_dimension": 68,
        "leading_rank": 68,
        "pivots": list(kernel_pivots),
        "free_degree12_columns": free,
        "free_degree12_count": 4,
        "lift": {
            "status": "exact_QQ_nullspace_of_all_36_exact_source_rows",
            "rank": 68,
            "pivots": list(kernel_pivots),
            "rational_rows": [[str(q) for q in row] for row in kernel_rref[:68]],
            "qq_candidate": True,
            "qq_verified": True,
            "all_theorem_primes_verified": True,
        },
    }
    kernel_path = Path(kernel_output)
    kernel_path.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(kernel_path, kernel_payload)
    return source_payload, kernel_payload
